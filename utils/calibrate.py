from __future__ import annotations

# from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast, get_args
from pydantic import ConfigDict
from pydantic.dataclasses import dataclass as pydantic_dataclass

import csv
import math
import re
import shutil
import json
import traceback
import os
import signal
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from postprocess import PostprocessConfig, run as _postprocess
from joblib import Parallel, delayed

from calc_utils import (
    mean_abs_error_aligned,
    mean_sq_error_aligned,
    mean_rel_err_aligned,
    convert_load,
)
from config import (
    COMPLGAS_BOUNDS,
    GasElastic,
    GasPlastic,
    SolidElastic,
    SolidPlastic,
    WS,
    WC,
    MICROSTRUCTURE_PATH,
    NTHREADS,
    CALIBRATIONS_DIR,
    CALIBRATION_INDEX_FILE,
    CALIBRATIONS_NEXT_ID_FILE,
    SRJ_OFFSET,
    DEFAULT_rb,
)
from data_utils import (
    SrjData,
    CreepData,
    SimResults,
    MicrostructureInfo,
)
from run import run_solver as run_solver_impl
from run import RunConfig
from write_fft import write_fft


CalibMode = Literal["srj", "creep", "all", "multi-load"]
ElasticMode = Literal["opt", "fixed"]
GasElasticMode = Literal[None, "opt", "fixed"]
CalibCase = Literal[0, 1, 2, 3, 4, 5, 6]
GasMode = Literal[0, 1, 2]
Metric = Literal["MAE", "MSE", "MRE"]
SrjAlign = Literal["time", "strain"]
Load = Literal[475, 525, 575]
MAX_LOSS = 1e2


def _write_json(path: Path, payload: dict[str, object]) -> None:

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _ensure_trailing_newline(path: Path) -> None:
    """Ensure a text file ends with a newline (so appends don't corrupt CSVs)."""

    path = Path(path)
    try:
        if not path.exists():
            return
        if path.stat().st_size == 0:
            return
        with path.open("rb+") as f:
            f.seek(-1, os.SEEK_END)
            last = f.read(1)
            if last != b"\n":
                f.seek(0, os.SEEK_END)
                f.write(b"\n")
    except OSError:
        pass


def _next_eval_file(results_dir: Path) -> Path:
    return Path(results_dir) / "_next_eval.txt"


def _max_existing_eval_dir_id(results_dir: Path) -> int:
    results_dir = Path(results_dir)
    pat = re.compile(r"^eval_(\d{4})$")
    mx = 0
    try:
        for p in results_dir.iterdir():
            if not p.is_dir():
                continue
            m = pat.match(p.name)
            if m is None:
                continue
            mx = max(mx, int(m.group(1)))
    except OSError:
        pass
    return mx


def _completed_eval_rows(results_dir: Path) -> int:
    results_path = Path(results_dir) / "calibration_results.csv"
    if not results_path.exists():
        return 0
    df = pd.read_csv(results_path)
    return int(len(df))


def _read_or_init_next_eval_id(results_dir: Path) -> int:
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    path = _next_eval_file(results_dir)

    file_next = 1
    if path.exists():
        try:
            file_next = int(path.read_text().strip())
        except Exception:
            file_next = 1

    repaired_next = max(
        1,
        file_next,
        _completed_eval_rows(results_dir) + 1,
        _max_existing_eval_dir_id(results_dir) + 1,
    )

    if (not path.exists()) or repaired_next != file_next:
        path.write_text(str(repaired_next))

    return repaired_next


def _allocate_eval_root(results_dir: Path) -> tuple[int, Path]:
    """Allocate the next sequential eval dir from `_next_eval.txt`."""
    results_dir = Path(results_dir)
    eval_dir_id = _read_or_init_next_eval_id(results_dir)
    eval_root = results_dir / f"eval_{eval_dir_id:04d}"

    if eval_root.exists():
        raise RuntimeError(
            f"Expected fresh eval dir {eval_root}, but it already exists. "
            "This indicates a broken _next_eval.txt state."
        )

    eval_root.mkdir(parents=True, exist_ok=False)
    _next_eval_file(results_dir).write_text(str(eval_dir_id + 1))
    return eval_dir_id, eval_root


def _log_eval_failure(
    *,
    eval_root: Path,
    where: str,
    exc: BaseException,
    extra: dict[str, object] | None = None,
) -> None:
    """Persist failure context for an evaluation.

    Calibration failures are costly to reproduce. We keep the eval directory on
    failure and write a short status JSON + a full traceback.
    """

    eval_root = Path(eval_root)
    eval_root.mkdir(parents=True, exist_ok=True)

    status: dict[str, object] = {
        "where": where,
        "exception_type": type(exc).__name__,
        "exception_message": str(exc),
    }
    if extra:
        status.update(extra)

    try:
        _write_json(eval_root / "eval_status.json", status)
    except Exception:
        pass

    try:
        (eval_root / "traceback.txt").write_text(
            "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        )
    except Exception:
        pass


def _decode_os_system_status(status: int) -> tuple[int, int | None]:
    """Decode the status returned by `os.system`.

    Returns
    -------
    (exit_code, term_signal)

    Notes
    -----
    - When a program exits normally, term_signal is None and exit_code is the
      process exit code.
    - When a program is terminated by a signal, term_signal is that signal
      number, and exit_code is conventionally 128 + signal.
    """

    try:
        if os.WIFSIGNALED(status):
            sig = int(os.WTERMSIG(status))
            return 128 + sig, sig
        if os.WIFEXITED(status):
            return int(os.WEXITSTATUS(status)), None
    except Exception:
        # Fall back to legacy behavior: treat as a raw exit code.
        pass
    # Best-effort fallback.
    if status > 255:
        return int(status >> 8), None
    return int(status), None


def _is_interrupt_exit(exit_code: int, term_signal: int | None) -> bool:
    """Return True if a solver status looks like user interrupt (Ctrl+C)."""

    if term_signal == int(signal.SIGINT):
        return True
    # Common conventions:
    # - bash returns 130 for SIGINT
    # - some wrappers return 2
    return int(exit_code) in {2, 130}


def allocate_next_calibration_id() -> int:
    """Allocate the next sequential numeric calibration id.

    Convention
    ----------
    - results are written under: results/calibrations/<ID>/
    - `_next_id.txt` stores the *next* ID to allocate.

    This is intentionally simple for a single-user workstation workflow.
    """

    CALIBRATIONS_DIR.mkdir(parents=True, exist_ok=True)

    # Initialize from directory scan if the file doesn't exist.
    if not CALIBRATIONS_NEXT_ID_FILE.exists():
        max_existing = 0
        for p in CALIBRATIONS_DIR.iterdir():
            if p.is_dir() and p.name.isdigit():
                max_existing = max(max_existing, int(p.name))
        CALIBRATIONS_NEXT_ID_FILE.write_text(str(max_existing + 1))

    nxt = int(CALIBRATIONS_NEXT_ID_FILE.read_text().strip())
    CALIBRATIONS_NEXT_ID_FILE.write_text(str(nxt + 1))
    return nxt


def _results_dir(cfg: CalibrateConfig) -> Path:
    return CALIBRATIONS_DIR / str(cfg.calib_id)


def _fftw_save_path(cfg: CalibrateConfig) -> Path:

    if cfg.fftw_save is None:

        fftw_save = _results_dir(cfg) / "fftw_wisdom_omp.save"
    else:

        fftw_save = Path(cfg.fftw_save)
        if fftw_save.exists() and fftw_save.is_dir():
            fftw_save = fftw_save / "fftw_wisdom_omp.save"

    return fftw_save


def calibration_dir(calib_id: int) -> Path:

    return CALIBRATIONS_DIR / str(int(calib_id))


def _results_file(cfg: CalibrateConfig) -> Path:
    return _results_dir(cfg) / "calibration_results.csv"


def _evpfft_dir(run_dir: Path) -> Path:
    """Return the directory where LS-EVPFFT writes str_str.out/vtk/etc."""

    return Path(run_dir) / "evpfft_outputs"


@pydantic_dataclass(config=ConfigDict(frozen=True, extra="forbid", strict=True))
class CalibrateConfig(RunConfig):
    calib_id: int
    mode: CalibMode = "srj"
    case: CalibCase = 0
    start_calib_id: int | None = None
    gas_elastic: GasElasticMode = None
    solid_elastic: ElasticMode = "fixed"
    w_s: float = WS  # srj
    w_c: float = WC  # creep
    w_iters: float = 0
    n_calls: int = 15
    random_state: int = 42
    keep_vtk: bool = False
    metric: Metric = "MSE"
    srj_align: SrjAlign = "strain"
    fixed: tuple | None = None
    micro_info: MicrostructureInfo | None = None
    ncalls: int = 15
    random_state: int = 42


@pydantic_dataclass(config=ConfigDict(frozen=True, extra="forbid", strict=True))
class Candidate:
    """Concrete parameter set used to write input files for one evaluation."""

    # Solid plasticity
    tau0xf: float
    tau0xb: float
    tau1x: float
    thet0: float
    thet1: float
    nrsx: float
    gamd0x: float

    # Solid cubic elastic constants (MPa)
    c11: float
    c12: float
    c44: float

    # Gas isotropic elasticity (used for igas=2)
    young: float
    nu: float

    # Gas kinetics (igas=2-relevant)
    gas_tau0xf: float
    gas_tau0xb: float
    gas_nrsx: float
    gas_gamd0x: float

    # complgas (only meaningful when igas==1)
    complgas: float


def _default_candidate(sim=None) -> Candidate:
    solp = SolidPlastic()
    sole = SolidElastic()
    gasp = GasPlastic()
    gase = GasElastic()
    return Candidate(
        tau0xf=solp.tau0xf,
        tau0xb=DEFAULT_rb,
        tau1x=solp.tau1x,
        thet0=solp.thet0,
        thet1=solp.thet1,
        nrsx=solp.nrsx,
        gamd0x=solp.gamd0x,
        c11=sole.c11,
        c12=sole.c12,
        c44=sole.c44,
        young=gase.young,
        nu=gase.nu,
        gas_tau0xf=gasp.tau0xf,
        gas_tau0xb=gasp.tau0xb,
        gas_nrsx=gasp.nrsx,
        gas_gamd0x=gasp.gamd0x,
        complgas=0.01,
    )


def _load_resume_history(
    cfg: CalibrateConfig,
) -> tuple[list[list[float]], list[float]]:
    """Load (x0, y0) from an existing calibration_results.csv for resume.

    The optimizer's state is represented by the set of already-evaluated points.
    We treat the per-evaluation CSV as the source of truth so long calibrations
    can be resumed across processes/machines without needing pickled objects.
    """

    results_path = _results_file(cfg)
    if not results_path.exists():
        raise FileNotFoundError(
            f"--resume_id requested but results file not found: {results_path}"
        )

    df = pd.read_csv(results_path)
    if len(df) == 0:
        return [], []

    layout = _x_layout(cfg)
    missing = [c for c in (layout + ["loss"]) if c not in df.columns]
    if missing:
        raise ValueError(
            f"Cannot resume from {results_path}: missing columns {missing} "
            f"(expected layout={layout} + ['base','loss'])."
        )
    if "tau0xb" in layout and "tau0xf" in layout:
        df["tau0xb"] = df["tau0xb"] / df["tau0xf"]
    x0 = df[layout].astype(float).values.tolist()
    y0 = df["loss"].astype(float).tolist()
    return x0, y0


def _load_calibration_index_row(calib_id: int) -> dict[str, object]:
    """Load a single row from the global calibration index.

    Used for post-hoc calibration cases where some parameters are frozen to a
    previous calibration id.
    """

    path = CALIBRATION_INDEX_FILE
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    mask = df["calib_id"].astype(int) == int(calib_id)
    if not mask.any():
        raise KeyError(f"calib_id={calib_id} not found in {path}")
    return cast(dict[str, object], df.loc[mask].iloc[-1].to_dict())


def _parameterize_run_cfg(
    base: CalibrateConfig,
    *,
    run_dir: Path,
    sim_case: str,
    cand: Candidate | None = None,
) -> CalibrateConfig:
    """Create a fully-parameterized CalibrateConfig for one run directory.

    This function is the only place where run_dir-dependent paths are assigned.
    The returned cfg is safe to pass to:
      - write_solid_files / write_gas_files
      - write_fft
      - run_solver
    """

    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    evpfft_dir = _evpfft_dir(run_dir)
    evpfft_dir.mkdir(parents=True, exist_ok=True)
    fft = run_dir / "fft.in"
    out_txt = run_dir / "output.txt"
    cuel = run_dir / "cuel.sx"
    cupl2 = run_dir / "cupl2.sx"
    cuel_gas = run_dir / "cuel_gas.sx"
    cupl2_gas = run_dir / "cupl2_gas.sx"

    if int(base.igas) in {0, 1}:
        cuel_gas = cuel
        cupl2_gas = cupl2

    fftw_save = _fftw_save_path(base)
    micro_info = MicrostructureInfo.load(base.microstructure)

    cfg = CalibrateConfig(
        **{
            **base.__dict__,
            "sim_case": sim_case,
            "run_dir": run_dir,
            "evpfft_dir": evpfft_dir,
            "fft": fft,
            "out_txt": out_txt,
            "cuel": cuel,
            "cupl2": cupl2,
            "cuel_gas": cuel_gas,
            "cupl2_gas": cupl2_gas,
            "fftw_save": fftw_save,
            "nphase": micro_info.nphase,
            "nx": micro_info.nx,
            "ny": micro_info.ny,
            "nz": micro_info.nz,
        }
    )

    if cand is not None:
        cfg = CalibrateConfig(
            **{
                **cfg.__dict__,
                "complgas": cand.complgas,
            }
        )
    return cfg


def run_solver(cfg: CalibrateConfig) -> int:
    """Run the solver inside a calibration evaluation directory.

    This wraps `utils.run.run_solver` so calibrations don't depend on the
    run.py CLI semantics.
    """

    return int(run_solver_impl(cfg, v=False, tee=True))


def _delete_vtk(out_dir: Path) -> None:
    """Delete VTK files in an output directory to keep calibration runs small."""

    if not out_dir.exists():
        return
    for p in out_dir.glob("*.vtk"):
        try:
            p.unlink()
        except OSError:
            pass


def _x_layout(cfg: CalibrateConfig) -> list[str]:
    """Return the list of optimization variables for this calibration case.

    Case mapping (integer):
      0: solid only, igas=0
      1: solid only, igas=1 (complgas fixed)
      2: solid + complgas, igas=1
      3: solid only, igas=2
      4: posthoc complgas only, igas=1 (solid frozen)
      5: posthoc gas kinetics only, igas=2 (solid frozen)
      6: joint igas=1 (solid + complgas)
      7: joint igas=2 (solid + gas kinetics)
    """

    names: list[str] = []

    sol_names: list[str] = [
        "tau0xf",
        "tau0xb",
        "tau1x",
        "thet0",
        "thet1",
        "nrsx",
        "gamd0x",
    ]
    gas_names: list[str] = [
        "gas_tau0xf",
        "gas_tau0xb",
        "gas_gamd0x",
        "gas_nrsx",
    ]

    if cfg.case in {0, 1, 2, 3, 6, 7}:
        names += sol_names

    if cfg.case in {2, 4, 6}:
        names += ["complgas"]

    if cfg.case in {5, 7}:
        names += gas_names
        if cfg.gas_elastic == "opt":
            names += ["young", "nu"]

    if cfg.solid_elastic == "opt" and cfg.case not in {4, 5}:
        names += ["c11", "c12", "c44"]

    fixed = cfg.fixed
    return [name for name in names if name not in fixed]


def _initial_guess(cfg: CalibrateConfig) -> list[float]:

    solp = SolidPlastic()
    sole = SolidElastic()
    gasp = GasPlastic()
    gase = GasElastic()

    defaults: dict[str, float] = {
        "tau0xf": solp.tau0xf,
        "tau0xb": solp.tau0xb,
        "tau1x": solp.tau1x,
        "thet0": solp.thet0,
        "thet1": solp.thet1,
        "gamd0x": solp.gamd0x,
        "nrsx": solp.nrsx,
        "c11": sole.c11,
        "c12": sole.c12,
        "c44": sole.c44,
        "complgas": cfg.complgas,
        "gas_tau0xf": gasp.tau0xf,
        "gas_gamd0x": gasp.gamd0x,
        "gas_nrsx": gasp.nrsx,
        "young": gase.young,
        "nu": gase.nu,
    }

    # If post-hoc, use the baseline's best parameters for fixed values.
    if cfg.start_calib_id is not None:
        row = _load_calibration_index_row(cfg.start_calib_id)
        defaults.update(
            {
                "tau0xf": row.get("tau0xf", defaults["tau0xf"]),
                "tau0xb": row.get("tau0xb", defaults["tau0xb"]),
                "tau1x": row.get("tau1x", defaults["tau1x"]),
                "thet0": row.get("thet0", defaults["thet0"]),
                "thet1": row.get("thet1", defaults["thet1"]),
                "gamd0x": row.get("gamd0x", defaults["gamd0x"]),
                "nrsx": row.get("nrsx", defaults["nrsx"]),
                "c11": row.get("c11", defaults["c11"]),
                "c12": row.get("c12", defaults["c12"]),
                "c44": row.get("c44", defaults["c44"]),
                "complgas": row.get("complgas", defaults["complgas"]),
                "gas_tau0xf": row.get("gas_tau0xf", defaults["gas_tau0xf"]),
                "gas_tau0xb": row.get("gas_tau0xb", defaults["gas_tau0xb"]),
                "gas_gamd0x": row.get("gas_gamd0x", defaults["gas_gamd0x"]),
                "gas_nrsx": row.get("gas_nrsx", defaults["gas_nrsx"]),
                "young": row.get("young", defaults["young"]),
                "nu": row.get("nu", defaults["nu"]),
            }
        )

    defaults["tau0xb"] = defaults["tau0xb"] / defaults["tau0xf"]
    return [defaults[name] for name in _x_layout(cfg)]


def _space(cfg: CalibrateConfig):
    # Local import so non-calibration runs don't require skopt.
    from skopt.space import Real

    all_bounds: dict[str, Real] = {}
    all_bounds.update(SolidPlastic.bounds())
    all_bounds.update(SolidElastic.bounds())
    all_bounds.update(GasPlastic.bounds())
    all_bounds.update(GasElastic.bounds())
    all_bounds["complgas"] = COMPLGAS_BOUNDS

    layout = _x_layout(cfg)
    missing = [n for n in layout if n not in all_bounds]
    if missing:
        raise ValueError(f"Missing bounds for: {missing}")

    return [all_bounds[name] for name in layout]


def candidate_from_x(
    cfg: CalibrateConfig, x: list[float], from_results: bool = False
) -> Candidate:
    layout = _x_layout(cfg)
    vals = {name: x[i] for i, name in enumerate(layout)}

    # Start from defaults then overwrite the calibrating subset.
    base = _default_candidate()

    # If post-hoc, overwrite fixed solid from baseline.
    if cfg.start_calib_id is not None and cfg.case in {4, 5}:
        # Posthoc cases freeze the solid to the baseline.
        row = _load_calibration_index_row(cfg.start_calib_id)
        base = Candidate(
            **{
                **base.__dict__,
                "tau0xf": row.get("tau0xf", base.tau0xf),
                "tau0xb": row.get("tau0xb", base.tau0xb),
                "tau1x": row.get("tau1x", base.tau1x),
                "thet0": row.get("thet0", base.thet0),
                "thet1": row.get("thet1", base.thet1),
                "nrsx": row.get("nrsx", base.nrsx),
                "gamd0x": row.get("gamd0x", base.gamd0x),
                "c11": row.get("c11", base.c11),
                "c12": row.get("c12", base.c12),
                "c44": row.get("c44", base.c44),
            }
        )

    # Solid calibration vars

    tau0xf = vals.get("tau0xf", base.tau0xf)
    tau0xb = vals.get("tau0xb", base.tau0xb)
    if not from_results:
        tau0xb *= tau0xf
    tau1x = vals.get("tau1x", base.tau1x)
    thet0 = vals.get("thet0", base.thet0)
    thet1 = vals.get("thet1", base.thet1)
    gamd0x = vals.get("gamd0x", base.gamd0x)
    nrsx = vals.get("nrsx", base.nrsx)

    # Solid elastic
    c11 = vals.get("c11", base.c11)
    c12 = vals.get("c12", base.c12)
    c44 = vals.get("c44", base.c44)
    # complgas
    complgas = vals.get("complgas", base.complgas)
    # Gas kinetics vars
    gas_tau0xf = vals.get("gas_tau0xf", base.gas_tau0xf)
    gas_tau0xb = vals.get("gas_tau0xb", base.gas_tau0xb)
    gas_gamd0x = vals.get("gas_gamd0x", base.gas_nrsx)
    gas_nrsx = vals.get("gas_nrsx", base.gas_gamd0x)

    young = vals.get("young", base.young)
    nu = vals.get("nu", base.nu)

    return Candidate(
        tau0xf=tau0xf,
        tau0xb=tau0xb,
        tau1x=tau1x,
        thet0=thet0,
        thet1=thet1,
        nrsx=nrsx,
        gamd0x=gamd0x,
        c11=c11,
        c12=c12,
        c44=c44,
        young=young,
        nu=nu,
        gas_tau0xf=gas_tau0xf,
        gas_tau0xb=gas_tau0xb,
        gas_nrsx=gas_nrsx,
        gas_gamd0x=gas_gamd0x,
        complgas=complgas,
    )


def _write_inputs_for_candidate(
    cfg: CalibrateConfig, cand: Candidate, *, work_dir: Path
) -> tuple[Path, Path, Path, Path]:
    # work_dir is a run directory root (eval_XXXX[/srj|creep]).
    run_cfg = _parameterize_run_cfg(
        cfg, run_dir=Path(work_dir), sim_case=str(cfg.sim_case)
    )

    from run import write_gas_files, write_solid_files

    solid_elastic = SolidElastic(c11=cand.c11, c12=cand.c12, c44=cand.c44)
    solid_plastic = SolidPlastic(
        tau0xf=cand.tau0xf,
        tau0xb=cand.tau0xb,
        tau1x=cand.tau1x,
        thet0=cand.thet0,
        thet1=cand.thet1,
        nrsx=cand.nrsx,
        gamd0x=cand.gamd0x,
    )
    write_solid_files(run_cfg, solid_elastic=solid_elastic, solid_plastic=solid_plastic)
    # Gas files are only physically used for igas=2; for igas=1 they are ignored.
    # For igas=0, point gas files to the solid phase files to match run.py behavior.
    if int(run_cfg.igas) == 2:
        gas_elastic = GasElastic(young=cand.young, nu=cand.nu)
        gas_plastic = GasPlastic(
            tau0xf=cand.gas_tau0xf,
            tau0xb=cand.gas_tau0xb,
            nrsx=cand.gas_nrsx,
            gamd0x=cand.gas_gamd0x,
        )
        write_gas_files(run_cfg, gas_elastic=gas_elastic, gas_plastic=gas_plastic)

    return run_cfg.cuel, run_cfg.cupl2, run_cfg.cuel_gas, run_cfg.cupl2_gas


def _srj_objective(cfg: CalibrateConfig, out, srj_exp: SrjData) -> float:

    exp_time, exp_mean_stress, exp_mean_strain = (
        srj_exp.time.copy(),
        srj_exp.mean_stress.copy(),
        srj_exp.mean_strain.copy(),
    )
    sim_time, sim_stress, sim_strain = (
        out.sim_time.copy(),
        out.sav33.copy(),
        out.eav33.copy(),
    )

    exp_time -= SRJ_OFFSET
    sim_time -= SRJ_OFFSET

    exp_mask = exp_time > 0
    sim_mask = sim_time > 0

    exp_time = exp_time[exp_mask]
    exp_mean_stress = exp_mean_stress[exp_mask]
    exp_mean_strain = exp_mean_strain[exp_mask]

    sim_time = sim_time[sim_mask]
    sim_stress = sim_stress[sim_mask]
    sim_strain = sim_strain[sim_mask]

    metric = cfg.metric
    if metric == "MAE":
        err_func = mean_abs_error_aligned
    elif metric == "MSE":
        err_func = mean_sq_error_aligned
    elif metric == "MRE":
        err_func = mean_rel_err_aligned
    else:
        raise ValueError(f"Unknown metric: {metric}")

    sim_y = sim_stress
    exp_y = exp_mean_stress

    if cfg.srj_align == "strain":
        exp_x = exp_mean_strain
        sim_x = sim_strain
    elif cfg.srj_align == "time":
        exp_x = exp_time
        sim_x = sim_time
    else:
        raise ValueError(f"Unknown srj_align: {cfg.srj_align}")

    # Restrict to overlapping x-range to avoid endpoint artifacts.
    if len(sim_x) == 0 or len(exp_x) == 0:
        stress_err = float("inf")
    else:

        x_lo = max(np.min(sim_x), np.min(exp_x))
        x_hi = min(np.max(sim_x), np.max(exp_x))
        exp_ok = (exp_x >= x_lo) & (exp_x <= x_hi)
        sim_ok = (sim_x >= x_lo) & (sim_x <= x_hi)

        stress_err = err_func(
            sim_x=sim_x[sim_ok],
            sim_y=sim_y[sim_ok],
            exp_x=exp_x[exp_ok],
            exp_y=exp_y[exp_ok],
        )

    return stress_err


def _creep_objective(
    cfg: CalibrateConfig,
    out,
    crp_exp: CreepData,
) -> float:

    load_vals = crp_exp.get_load_vals(cfg.load)
    exp_strain = load_vals["mean_plastic_strain"].copy()
    exp_time = load_vals["time"][1:].copy()

    sim_strain = out.epav33.copy()
    sim_time = out.sim_time.copy()

    metric = cfg.metric

    if metric == "MAE":
        err_func = mean_abs_error_aligned
    elif metric == "MSE":
        err_func = mean_sq_error_aligned
    elif metric == "MRE":
        err_func = mean_rel_err_aligned
    else:
        raise ValueError(f"Unknown metric: {metric}")

    if len(sim_time) == 0 or len(exp_time) == 0:
        strain_err = float("inf")
    else:

        x_lo = max(np.min(sim_time), np.min(exp_time))
        x_hi = min(np.max(sim_time), np.max(exp_time))
        exp_ok = (exp_time >= x_lo) & (exp_time <= x_hi)
        sim_ok = (sim_time >= x_lo) & (sim_time <= x_hi)

        strain_err = err_func(
            sim_x=sim_time[sim_ok],
            sim_y=sim_strain[sim_ok],
            exp_x=exp_time[exp_ok],
            exp_y=exp_strain[exp_ok],
        )

    return strain_err


def _results_header(cfg: CalibrateConfig) -> str:
    """CSV header for per-evaluation results.

    Note: this is intentionally case-driven so post-hoc runs log only the
    parameters that are being evaluated/optimized.
    """

    cols: list[str] = []
    layout = _x_layout(cfg)

    # Log evaluation-space variables first (these match bounds/layout).
    cols += list(layout)

    cols += ["loss"]
    return ",".join(cols)


def _ensure_results_header(cfg: CalibrateConfig) -> None:
    """Ensure calibration_results.csv exists and matches the expected schema."""

    results_path = _results_file(cfg)
    results_path.parent.mkdir(parents=True, exist_ok=True)
    expected = _results_header(cfg)

    if not results_path.exists() or results_path.stat().st_size == 0:
        results_path.write_text(expected + "\n")
        return

    with results_path.open("r", errors="ignore") as f:
        header = (f.readline() or "").strip()
    if header != expected:
        raise ValueError(
            f"Unexpected header in {results_path}.\n"
            f"Expected:\n  {expected}\n"
            f"Got:\n  {header}\n"
            "Start a new calibration id (or rewrite the header) to proceed."
        )


def _append_eval(
    cfg: CalibrateConfig,
    *,
    cand: Candidate,
    loss: float,
) -> None:
    results_path = _results_file(cfg)
    _ensure_results_header(cfg)
    _ensure_trailing_newline(results_path)

    layout = _x_layout(cfg)
    eval_vars: dict[str, float] = {
        # solid calibration vars
        "tau0xf": cand.tau0xf,
        "tau0xb": cand.tau0xb,
        "tau1x": cand.tau1x,
        "thet0": cand.thet0,
        "thet1": cand.thet1,
        "gamd0x": cand.gamd0x,
        "nrsx": cand.nrsx,
        "c11": cand.c11,
        "c12": cand.c12,
        "c44": cand.c44,
        # gas calibration vars
        "complgas": cand.complgas,
        "gas_tau0xf": cand.gas_tau0xf,
        "gas_tau0xb": cand.gas_tau0xb,
        "gas_gamd0x": cand.gas_gamd0x,
        "gas_log10_m": cand.gas_nrsx,
        "young": cand.young,
        "nu": cand.nu,
    }

    fields = [eval_vars[n] for n in layout]
    fields += [loss]

    with results_path.open("a") as f:
        f.write(",".join(f"{v:.6g}" for v in fields) + "\n")


def prepare_one(
    cfg: CalibrateConfig,
    eval_root: Path,
    cand: Candidate,
    sim_case: str,
    *,
    subdir: str | None = None,
) -> CalibrateConfig:
    run_dir = eval_root if subdir is None else (eval_root / subdir)
    run_dir.mkdir(parents=True, exist_ok=True)

    cuel, cupl2, cuel_gas, cupl2_gas = _write_inputs_for_candidate(
        cfg, cand, work_dir=run_dir
    )

    run_cfg = _parameterize_run_cfg(cfg, run_dir=run_dir, sim_case=sim_case, cand=cand)
    run_cfg = CalibrateConfig(
        **{
            **run_cfg.__dict__,
            "complgas": cand.complgas,
            "cuel": cuel,
            "cupl2": cupl2,
            "cuel_gas": cuel_gas,
            "cupl2_gas": cupl2_gas,
        }
    )
    write_fft(run_cfg)
    return run_cfg, run_dir


srj_exp = srj_exp = SrjData.load()
crp_exp = CreepData.load(skip_roughness=True)


def run_one(
    cfg: CalibrateConfig,
    eval_root: Path,
    cand: Candidate,
    sim_case: str,
    *,
    keep_vtk: bool = False,
    skip_vtk: bool = True,
    subdir: str | None = None,
    postprocess: bool = False,
) -> float:

    global srj_exp
    global crp_exp

    run_cfg, run_dir = prepare_one(cfg, eval_root, cand, sim_case, subdir=subdir)
    rc_status = run_solver(run_cfg)
    exit_code, term_sig = _decode_os_system_status(int(rc_status))

    if _is_interrupt_exit(exit_code, term_sig):
        raise KeyboardInterrupt()

    try:
        assert exit_code == 0
    except:

        (eval_root / "solver_failed.txt").write_text(
            f"Solver failed. exit_code={exit_code} term_signal={term_sig}\n"
        )
        raise

    out_dir = _evpfft_dir(run_dir)
    out = SimResults.load(out_dir, skip_vtk=True)
    iters = out.iter
    avg_iters = np.mean(iters)

    if sim_case.split("_")[0] == "srj":

        base = _srj_objective(cfg, out, srj_exp)
    elif sim_case.split("_")[0] == "creep":

        base = _creep_objective(cfg, out, crp_exp)
    else:
        raise ValueError(f"Unknown sim_case {sim_case}")

    loss = base + cfg.w_iters * avg_iters

    if not keep_vtk:
        _delete_vtk(out_dir)

    if postprocess:

        postprocess_cfg = PostprocessConfig(
            out_dir=run_dir,
            sim_case=sim_case,
            igas=int(cfg.igas),
            microstructure=cfg.microstructure,
            load=int(cfg.sol_load),
            skip_vtk=skip_vtk,
        )
        _postprocess(postprocess_cfg)

    return loss


def run_srj_creep(
    cfg: CalibrateConfig,
    eval_root: Path,
    cand: Candidate,
    keep_vtk=False,
    skip_vtk=True,
    postprocess=False,
) -> float:
    srj_loss, creep_loss = Parallel(
        n_jobs=2,
        backend="loky",
        prefer="processes",
    )(
        [
            delayed(run_one)(
                cfg,
                eval_root,
                cand,
                "srj_calibration",
                keep_vtk=keep_vtk,
                subdir="srj",
                skip_vtk=skip_vtk,
                postprocess=postprocess,
            ),
            delayed(run_one)(
                cfg,
                eval_root,
                cand,
                "creep_calibration",
                keep_vtk=keep_vtk,
                subdir="creep",
                skip_vtk=skip_vtk,
                postprocess=postprocess,
            ),
        ]
    )

    loss = cfg.w_s * srj_loss + cfg.w_c * creep_loss

    return loss


def run_multi_load(
    cfg: CalibrateConfig,
    eval_root: Path,
    cand: Candidate,
    keep_vtk=False,
    skip_vtk=True,
    postprocess=False,
) -> float:

    loads = get_args(Load)
    nloads = len(loads)
    jobs = []

    ws = 1 / nloads

    for load in loads:

        load_cfg = CalibrateConfig(
            **{
                **cfg.__dict__,
                "load": convert_load(load, cfg.micro_info),
                "sol_load": load,
            }
        )
        job = delayed(run_one)(
            load_cfg,
            eval_root,
            cand,
            "creep_calibration",
            keep_vtk=keep_vtk,
            subdir=f"creep{load}",
            skip_vtk=skip_vtk,
            postprocess=postprocess,
        )

        jobs.append(job)

    losses = Parallel(
        n_jobs=len(loads),
        backend="loky",
        prefer="processes",
    )(jobs)

    loss = np.sum(losses) * ws

    return loss


def run_mode(
    cfg: CalibrateConfig,
    eval_root: Path,
    cand: Candidate,
    keep_vtk=False,
    skip_vtk=True,
    postprocess=False,
) -> float:

    if cfg.mode == "srj":
        loss = run_one(
            cfg,
            eval_root,
            cand,
            "srj_calibration",
            keep_vtk=keep_vtk,
            skip_vtk=skip_vtk,
            postprocess=postprocess,
        )
    elif cfg.mode == "creep":
        loss = run_one(
            cfg,
            eval_root,
            cand,
            "creep_calibration",
            keep_vtk=keep_vtk,
            skip_vtk=skip_vtk,
            postprocess=postprocess,
        )
    elif cfg.mode == "srj_creep":
        loss = run_srj_creep(
            cfg,
            eval_root,
            cand,
            keep_vtk=keep_vtk,
            skip_vtk=skip_vtk,
            postprocess=postprocess,
        )
    elif cfg.mode == "multi-load":
        loss = run_multi_load(
            cfg,
            eval_root,
            cand,
            keep_vtk=keep_vtk,
            skip_vtk=skip_vtk,
            postprocess=postprocess,
        )

    else:
        raise ValueError(f"Unknown calibration mode: {cfg.mode}")

    return loss


def calibration_core(cfg: CalibrateConfig, x: list[float]) -> float:
    cand = candidate_from_x(cfg, x)
    eval_dir_id, eval_root = _allocate_eval_root(_results_dir(cfg))
    eval_idx = eval_dir_id - 1

    keep_eval_root = False

    try:
        loss = run_mode(cfg, eval_root, cand)

    except Exception as e:
        keep_eval_root = True
        _log_eval_failure(
            eval_root=eval_root,
            where="calibration_core",
            exc=e,
            extra={"eval_idx": eval_idx, "eval_dir_id": eval_dir_id, "mode": cfg.mode},
        )

        loss = MAX_LOSS

    finally:
        if not keep_eval_root:
            try:
                shutil.rmtree(eval_root)
            except OSError:
                pass

    if math.isnan(loss) or math.isinf(loss):
        loss = MAX_LOSS

    _append_eval(
        cfg,
        cand=cand,
        loss=loss,
    )

    return loss


INDEX_COLUMNS: tuple[str, ...] = (
    "calib_id",
    "created_at",
    "notes",
    "start_calib_id",
    "mode",
    "srj_align",
    "case",
    "igas",
    "complgas",
    "load",
    "gas_elastic",
    "solid_elastic",
    "nthreads",
    "modal",
    "mix",
    "microstructure",
    "fftw_save",
    "w_s",
    "w_c",
    "w_iters",
    "ncalls",
    "random_state",
    # solid plastic (both reduced variables and derived physical)
    "tau0xf",
    "tau0xb",
    "tau1x",
    "thet0",
    "thet1",
    "gamd0x",
    "nrsx",
    # solid elastic
    "c11",
    "c12",
    "c44",
    # gas phase
    "complgas",
    "gas_tau0xf",
    "gas_tau0xb",
    "gas_tau0xb",
    "gas_gamd0x",
    "gas_nrsx",
    "young",
    "nu",
    "loss",
)


def _case_default_igas(case: int) -> int:
    if case == 0:
        return 0
    if case in {1, 2, 4, 6}:
        return 1
    return 2


def _append_index_row(row: dict[str, object]) -> None:
    path = CALIBRATION_INDEX_FILE
    path.parent.mkdir(parents=True, exist_ok=True)

    # If the index doesn't exist, create it with the current column set.
    desired = list(INDEX_COLUMNS)

    try:
        target_id = int(row.get("calib_id", ""))
    except Exception as e:
        raise ValueError(
            f"Index row missing/invalid calib_id: {row.get('calib_id')}"
        ) from e

    if not path.exists() or path.stat().st_size == 0:
        new_header = desired
        new_rows: list[dict[str, object]] = [{k: row.get(k, "") for k in new_header}]
    else:
        with path.open("r", newline="") as f:
            reader = csv.DictReader(f)
            header = list(reader.fieldnames or [])
            old_rows = list(reader)

        extras = [c for c in header if c not in desired]
        new_header = desired + extras

        # Remove any existing rows for this calib_id (resume should replace).
        kept_rows: list[dict[str, object]] = []
        insert_at: int | None = None
        for r in old_rows:
            try:
                rid = int(r.get("calib_id", ""))
            except Exception:
                # Preserve malformed/legacy rows rather than dropping.
                rid = None
            if rid == target_id:
                if insert_at is None:
                    insert_at = len(kept_rows)
                continue
            kept_rows.append(r)

        if insert_at is None:
            insert_at = len(kept_rows)

        kept_rows.insert(insert_at, {k: row.get(k, "") for k in new_header})
        new_rows = kept_rows

    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=new_header)
        w.writeheader()
        for r in new_rows:
            w.writerow({k: r.get(k, "") for k in new_header})
    tmp.replace(path)
    return None


def _row_to_candidate(cfg: CalibrateConfig, row: dict[str, object]) -> Candidate:
    """Convert a calibration_results row (x-layout variables) into a Candidate."""

    # Build an x vector in layout order.
    x = []
    for name in _x_layout(cfg):
        x.append(float(row[name]))
    return candidate_from_x(cfg, x, from_results=True)


def _write_best_row_csv(dst_dir: Path, row: dict[str, object]) -> None:
    dst_dir.mkdir(parents=True, exist_ok=True)
    p = dst_dir / "best_row.csv"
    keys = row.keys()
    with p.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        w.writerow([row.get(k, "") for k in keys])


def finalize(cfg: CalibrateConfig) -> None:

    results_path = _results_file(cfg)
    if not results_path.exists():
        raise RuntimeError(f"Expected calibration results at {results_path}")

    df = pd.read_csv(results_path)
    if len(df) == 0:
        raise RuntimeError(f"Empty calibration results at {results_path}")

    best_row = df.iloc[int(df["loss"].astype(float).argmin())].to_dict()
    cand = _row_to_candidate(cfg, best_row)
    _write_best_row_csv(calibration_dir(cfg.calib_id), best_row)

    idx_row: dict[str, object] = {k: "" for k in INDEX_COLUMNS}
    idx_row.update(
        {
            "calib_id": int(cfg.calib_id),
            "created_at": datetime.now(timezone.utc).isoformat(),
            "start_calib_id": cfg.start_calib_id,
            "mode": cfg.mode,
            "srj_align": cfg.srj_align,
            "case": cfg.case,
            "igas": cfg.igas,
            "complgas": cfg.complgas,
            "load": cfg.load,
            "gas_elastic": cfg.gas_elastic,
            "solid_elastic": cfg.solid_elastic,
            "nthreads": cfg.nthreads,
            "modal": cfg.modal,
            "mix": cfg.mix,
            "microstructure": cfg.microstructure,
            "fftw_save": cfg.fftw_save,
            "w_s": cfg.w_s,
            "w_c": cfg.w_c,
            "w_iters": cfg.w_iters,
            "ncalls": cfg.ncalls,
            "random_state": cfg.random_state,
            "notes": cfg.notes,
            # solid
            "tau0xf": cand.tau0xf,
            "tau0xb": cand.tau0xb,
            "tau1x": cand.tau1x,
            "thet0": cand.thet0,
            "thet1": cand.thet1,
            "nrsx": cand.nrsx,
            "gamd0x": cand.gamd0x,
            "c11": cand.c11,
            "c12": cand.c12,
            "c44": cand.c44,
            # gas
            "complgas": cand.complgas,
            "gas_tau0xf": cand.gas_tau0xf,
            "gas_tau0xb": cand.gas_tau0xb,
            "gas_nrsx": cand.gas_nrsx,
            "gas_gamd0x": cand.gas_gamd0x,
            "young": cand.young,
            "nu": cand.nu,
            # metrics
            "loss": best_row.get("loss", np.inf),
        }
    )

    _append_index_row(idx_row)

    final_root = _results_dir(cfg) / "final"
    run_mode(cfg, final_root, cand, keep_vtk=True, skip_vtk=False, postprocess=True)


def main(argv: list[str] | None = None) -> int:
    import argparse
    import textwrap

    p = argparse.ArgumentParser(
        prog="calibrate.py",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    p.add_argument("--mode", choices=get_args(CalibMode), default="srj")
    p.add_argument(
        "--case",
        type=int,
        default=0,
        choices=get_args(CalibCase),
        help=textwrap.dedent(
            """\
            Integer calibration case:
                0: solid only                      igas=0    (optional --gas_elastic opt)
                1: solid only, complgas fixed      igas=1    (optional --gas_elastic opt)
                2: solid only, gas fixed           igas=2    (optional --gas_elastic opt)
                3: complgas only                   igas=1    (solid params from --paramset)
                4: gas only                        igas=2    (solid frozen from --paramset) (optional --gas_elastic opt)
                5: joint solid & complgas          igas=1    (optional --solid_elastic opt)
                6: joint solid & gas               igas=2    (optional --solid_elastic opt) (optional --gas_elastic opt)
            
            Note: All igas=2 simulations use isotropic elastic parameters and kinetic parameters (no hardening)\n\n\n
            """
        ),
    )

    p.add_argument(
        "--paramset",
        type=int,
        default=None,
        help="optional starting calibration id",
    )
    p.add_argument("--complgas", type=float, default=None)
    p.add_argument("--load", type=int, choices=get_args(Load), default=None)
    p.add_argument("--nthreads", type=int, default=NTHREADS)
    p.add_argument(
        "--notes",
        type=str,
        default="",
        help="free-form notes recorded in results/calibrations/calibration_index.csv",
    )
    p.add_argument("--micro", type=Path, default=MICROSTRUCTURE_PATH)
    p.add_argument("--fftw_save", type=Path, default=None)
    p.add_argument("--gas_elastic", choices=get_args(GasElasticMode), default=None)
    p.add_argument(
        "--solid_elastic",
        choices=get_args(ElasticMode),
        default="fixed",
    )
    p.add_argument("--w_iters", type=float, default=0)
    p.add_argument("--w_s", type=float, default=WS)
    p.add_argument("--w_c", type=float, default=WC)
    p.add_argument(
        "--ncalls",
        type=int,
        default=15,
        help=(
            "Number of solver evaluations to run. When resuming (via --resume_id), "
            "this is the *additional* number of evaluations to run in this invocation."
        ),
    )
    p.add_argument("--random_state", type=int, default=42)
    p.add_argument(
        "--resume_id",
        type=int,
        default=None,
        help=(
            "Resume an existing calibration id by reusing its calibration_results.csv. "
            "When resuming, --n_calls is interpreted as the additional number of evaluations to run."
        ),
    )
    p.add_argument(
        "--keep_vtk",
        action="store_true",
        help="keep VTKs for the best run (default: delete)",
    )
    p.add_argument("--metric", choices=get_args(Metric), default="MSE")
    p.add_argument(
        "--srj_align",
        choices=get_args(SrjAlign),
        default="time",
        help="SRJ objective alignment axis: 'strain'or 'time'.",
    )
    p.add_argument(
        "--mix",
        type=int,
        default=1,
        help="apply mixing to velgrad update",
    )
    p.add_argument(
        "--modal",
        type=int,
        default=0,
        help="use modified augmented lagrangian",
    )

    FIXABLE_PARAMS = tuple(Candidate.__annotations__.keys())

    p.add_argument(
        "--fixed",
        nargs="+",
        metavar="PARAM",
        choices=FIXABLE_PARAMS,
        default=[],
        help=(
            "Material parameters to hold fixed during calibration. "
            f"Allowed values: {', '.join(FIXABLE_PARAMS)}"
        ),
    )

    args = p.parse_args(argv)

    if args.case in {3, 4} and args.paramset is None:
        print(
            "WARNING: posthoc case requested but no --paramset provided; using repo defaults as baseline"
        )

    igas = _case_default_igas(int(args.case))

    if args.resume_id is None:
        work_id = allocate_next_calibration_id()
        work_dir = calibration_dir(work_id)
        work_dir.mkdir(parents=True, exist_ok=True)
        print(f"Starting calibration {work_id} in {work_dir}")
    else:
        work_id = int(args.resume_id)
        work_dir = calibration_dir(work_id)
        if not work_dir.exists():
            raise FileNotFoundError(
                f"--resume_id={work_id} requested but directory not found: {work_dir}"
            )
        print(f"Resuming calibration {work_id} in {work_dir}")

    sol_load = args.load
    micro_info = MicrostructureInfo.load(args.micro)
    if args.load is not None and args.case > 0:
        args.load = convert_load(args.load, micro_info)

    cfg = CalibrateConfig(
        calib_id=work_id,
        out_dir=work_dir,
        mode=args.mode,
        case=args.case,
        start_calib_id=args.paramset,
        igas=igas,
        complgas=args.complgas,
        nthreads=args.nthreads,
        notes=str(args.notes) if str(args.notes).strip() else None,
        microstructure=args.micro,
        fftw_save=args.fftw_save,
        load=args.load,
        gas_elastic=args.gas_elastic,
        solid_elastic=args.solid_elastic,
        w_s=args.w_s,
        w_c=args.w_c,
        w_iters=args.w_iters,
        srj_align=args.srj_align,
        modal=args.modal,
        mix=args.mix,
        sol_load=sol_load,
        fixed=tuple(args.fixed),
        micro_info=micro_info,
        ncalls=args.ncalls,
        metric=args.metric,
    )

    dims = _space(cfg)
    x0 = _initial_guess(cfg)
    resume_x0: list[list[float]] | None = None
    resume_y0: list[float] | None = None
    gp_n_calls = int(args.ncalls)
    if args.resume_id is not None:
        prev_x0, prev_y0 = _load_resume_history(cfg)
        completed = len(prev_x0)
        additional = int(args.ncalls)
        if additional < 0:
            raise ValueError("--n_calls must be >= 0 when resuming.")

        gp_n_calls = additional
        _read_or_init_next_eval_id(_results_dir(cfg))

        print(
            f"Resuming calibration {work_id}: {completed} completed evals; "
            f"running {additional} more (gp_minimize total n_calls={gp_n_calls})."
        )
        if completed > 0:
            resume_x0 = prev_x0
            resume_y0 = prev_y0

    from skopt import gp_minimize

    try:

        kwargs = {
            "func": lambda x: calibration_core(cfg, list(x)),
            "dimensions": dims,
            "n_calls": int(gp_n_calls),
            "random_state": int(args.random_state),
        }
        if resume_x0 is not None and resume_y0 is not None:
            kwargs.update({"x0": resume_x0, "y0": resume_y0})
        else:
            kwargs.update({"x0": x0})

        gp_minimize(**kwargs)
        finalize(cfg)

        return 0

    except KeyboardInterrupt:

        try:

            finalize(cfg)

        except Exception:
            (calibration_dir(work_id) / "finalize_traceback.txt").write_text(
                traceback.format_exc()
            )

            raise

        raise
    except Exception:

        try:

            finalize(cfg)

        except Exception:

            (calibration_dir(work_id) / "finalize_traceback.txt").write_text(
                traceback.format_exc()
            )

        raise


if __name__ == "__main__":
    raise SystemExit(main())


# def _isotropic_young_from_cubic(c11: float, c12: float, c44: float) -> float:
#     """Voigt-Reuss-Hill isotropic Young's modulus (MPa) for cubic stiffnesses."""

#     c11 = c11
#     c12 = c12
#     c44 = c44
#     bulk_k = (c11 + 2.0 * c12) / 3.0
#     shear_g_v = (c11 - c12 + 3.0 * c44) / 5.0
#     denom = 4.0 * c44 + 3.0 * (c11 - c12)
#     shear_g_r = shear_g_v if denom == 0.0 else (5.0 * (c11 - c12) * c44) / denom
#     shear_g = 0.5 * (shear_g_v + shear_g_r)
#     denom_e = 3.0 * bulk_k + shear_g
#     if denom_e == 0.0:
#         return np.inf
#     return 9.0 * bulk_k * shear_g / denom_e


# @functools.lru_cache(maxsize=1)
# def _default_solid_young() -> float:
#     sole = SolidElastic()
#     return _isotropic_young_from_cubic(sole.c11, sole.c12, sole.c44)

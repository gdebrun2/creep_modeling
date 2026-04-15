from __future__ import annotations

from pathlib import Path
from typing import Any
import pandas as pd


from io_utils import write_cuel, write_cupl2, allocate_run_dir, register_run, os_system
from config import (
    GasElastic,
    GasPlastic,
    SolidElastic,
    SolidPlastic,
    SOLVER_EXE,
    NTHREADS,
    DEFAULT_PARAMSET,
    MICROSTRUCTURE_PATH,
    DEFAULT_SIM,
    CALIBRATION_INDEX_FILE,
)
from data_utils import MicrostructureInfo
from calc_utils import convert_load
from write_fft import FftConfig, GasMode, SimCase, write_fft
from typing import get_args
from pydantic import ConfigDict
from pydantic.dataclasses import dataclass as pydantic_dataclass
from postprocess import PostprocessConfig, run as _postprocess


@pydantic_dataclass(config=ConfigDict(frozen=True, extra="forbid", strict=True))
class RunConfig(FftConfig):

    calib_id: int | None = DEFAULT_PARAMSET
    nthreads: int = NTHREADS
    notes: str | None = None
    run_dir: Path | None = None
    out_dir: Path | None = None
    evpfft_dir: Path | None = None
    out_txt: Path | None = None
    fftw_save: Path | None = None
    modal: int | None = None
    mix: int | None = None
    sol_load: int | None = None


def load_calib_index(calib_id: int | None) -> dict[str, Any]:
    """Load a single row from the global calibration index.

    Returns an empty dict if the index file does not exist or if the requested id
    is not present.
    """
    if not calib_id:
        print("Continuing without calibrated parameters...")
        return {}
    path = CALIBRATION_INDEX_FILE
    if not path.exists():
        print(
            f"No calibration index found at {path}. \nContinuing with default parameters..."
        )
        return {}

    df = pd.read_csv(path)
    mask = df["calib_id"].astype(int) == int(calib_id)
    if not mask.any():
        print(
            f"No calibration id={calib_id} found in {path}. \nContinuing with default parameters..."
        )
        return {}

    row = df.loc[mask].iloc[-1]

    return row.to_dict()


def _solid_plastic_from_row(
    row: dict[str, Any], defaults: SolidPlastic = SolidPlastic()
) -> SolidPlastic:
    """Extract solid plastic params from a calibration row."""

    return SolidPlastic(
        tau0xf=float(row.get("tau0xf", defaults.tau0xf)),
        tau0xb=float(row.get("tau0xb", defaults.tau0xb)),
        tau1x=float(row.get("tau1x", defaults.tau1x)),
        thet0=float(row.get("thet0", defaults.thet0)),
        thet1=float(row.get("thet1", defaults.thet1)),
        nrsx=float(row.get("nrsx", defaults.nrsx)),
        gamd0x=float(row.get("gamd0x", defaults.gamd0x)),
    )


def _solid_elastic_from_row(
    row: dict[str, Any], defaults: SolidElastic = SolidElastic()
) -> SolidElastic:
    """Extract solid elastic constants if present; otherwise return defaults."""

    # Accept both legacy and explicit names.
    c11 = row.get("c11", defaults.c11)
    c12 = row.get("c12", defaults.c12)
    c44 = row.get("c44", defaults.c44)
    return SolidElastic(c11=c11, c12=c12, c44=c44)


def _gas_elastic_from_row(
    row: dict[str, Any], defaults: GasElastic = GasElastic()
) -> GasElastic:

    return GasElastic(
        young=float(row.get("young", defaults.young)),
        nu=float(row.get("nu", defaults.nu)),
    )


def _gas_plastic_from_row(
    row: dict[str, Any], defaults: GasPlastic = GasPlastic()
) -> GasPlastic:

    return GasPlastic(
        tau0xf=float(row.get("gas_tau0xf", defaults.tau0xf)),
        tau0xb=float(row.get("gas_tau0xb", defaults.tau0xb)),
        tau1x=float(row.get("gas_tau1x", defaults.tau1x)),
        thet0=float(row.get("gas_thet0", defaults.thet0)),
        thet1=float(row.get("gas_thet1", defaults.thet1)),
        nrsx=float(row.get("gas_nrsx", defaults.nrsx)),
        gamd0x=float(row.get("gas_gamd0x", defaults.gamd0x)),
    )


def write_solid_files(
    cfg: RunConfig, *, solid_elastic: SolidElastic, solid_plastic: SolidPlastic
) -> None:
    """Write solid `cuel.sx` and `cupl2.sx` into run_dir."""

    write_cuel(cuel_path=cfg.cuel, elastic=solid_elastic, iso=0)
    write_cupl2(cupl2_path=cfg.cupl2, plastic=solid_plastic)
    return None


def write_gas_files(
    cfg: RunConfig, *, gas_elastic: GasElastic, gas_plastic: GasPlastic
) -> None:
    """Write gas/damper phase files into run_dir (`cuel_gas.sx`, `cupl2_gas.sx`)."""

    write_cuel(cuel_path=cfg.cuel_gas, elastic=gas_elastic, iso=1)
    write_cupl2(cupl2_path=cfg.cupl2_gas, plastic=gas_plastic)
    return None


def write_inputs(cfg: RunConfig, skip_fft: bool = False) -> None:
    """Write all required inputs into the per-run directory and generate fft.in."""

    row = load_calib_index(cfg.calib_id)
    # solid material params
    sim = cfg.sim_case.split("_")[0]
    solid_plastic = _solid_plastic_from_row(row, defaults=SolidPlastic())
    solid_elastic = _solid_elastic_from_row(row)
    write_solid_files(cfg, solid_elastic=solid_elastic, solid_plastic=solid_plastic)

    if int(cfg.igas) == 2:
        gas_elastic = _gas_elastic_from_row(row)
        gas_plastic = _gas_plastic_from_row(row)
        write_gas_files(cfg, gas_elastic=gas_elastic, gas_plastic=gas_plastic)

    # fft.in
    if not skip_fft:
        write_fft(cfg)
    return None


def run_solver(cfg: RunConfig, *, v: bool = False, tee: bool = True) -> int:
    """Run the solver with fft.in in run_dir, teeing output.txt.

    The solver is invoked with:
    - --inputfile <run_dir>/fft.in
    - --outdir <run_dir>/evpfft_outputs
    so all generated output files (including FFTW wisdom) land inside the
    per-run folder.
    """

    if cfg.evpfft_dir is None or cfg.out_txt is None:
        raise ValueError("cfg.evpfft_dir and cfg.out_txt must be set before run_solver")

    prefix = "set -o pipefail; " if tee else ""
    cmd = prefix + (
        f"{SOLVER_EXE} --nthreads {cfg.nthreads} "
        f"--in {cfg.fft} "
        f"--out {cfg.evpfft_dir} "
        f"--fftw_save {cfg.fftw_save} "
        f"--mix {cfg.mix} "
        f"--modal {cfg.modal} "
    )
    if tee:
        cmd += f"| tee {cfg.out_txt}; "

    if v:
        cmd += " say -v Daniel done "
    return os_system(cmd)


def run_solver_nlive(
    cfg: RunConfig,
    *,
    v: bool = False,
    tee: bool = True,
    n_minutes: float = 2.0,
) -> int:
    """Run the solver and trigger postprocessing every `n_minutes`.

    The solver itself is launched through `run_solver(...)` so its stdout/stderr
    behavior remains unchanged. Periodic plot refreshes are launched as separate
    subprocesses while the solver runs.
    """

    if cfg.out_dir is None:
        raise ValueError("cfg.out_dir must be set before run_solver_nlive")

    import subprocess
    import sys
    import threading
    import time
    from contextlib import suppress

    interval_s = max(float(n_minutes) * 60.0, 30.0)
    poll_s = min(5.0, interval_s)

    load = int(cfg.sol_load) if cfg.sol_load is not None else int(cfg.load or 0)
    postprocess_cmd = [
        sys.executable,
        "/Users/gtdebru/creep_modeling/utils/postprocess.py",
        "--out_dir",
        str(cfg.out_dir),
        "--sim_case",
        str(cfg.sim_case),
        "--igas",
        str(int(cfg.igas)),
        "--micro",
        str(cfg.microstructure),
        "--load",
        str(load),
        "--skip_vtk",
        "1",
    ]

    result: dict[str, Any] = {"rc": 1, "exc": None}
    done = threading.Event()

    def _solver_target() -> None:
        try:
            result["rc"] = int(run_solver(cfg, v=v, tee=tee))
        except BaseException as exc:  # preserve KeyboardInterrupt/SystemExit
            result["exc"] = exc
        finally:
            done.set()

    worker = threading.Thread(target=_solver_target, daemon=True)
    worker.start()

    next_refresh = time.monotonic() + interval_s
    try:
        while not done.wait(timeout=poll_s):
            now = time.monotonic()
            if now >= next_refresh:
                with suppress(Exception):
                    subprocess.run(
                        postprocess_cmd,
                        check=False,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                    )
                next_refresh = now + interval_s
    except KeyboardInterrupt:
        raise
    finally:
        worker.join(timeout=1.0)

    if result["exc"] is not None:
        raise result["exc"]
    return int(result["rc"])


def main(argv: list[str] | None = None) -> int:
    import argparse
    import textwrap

    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--calibrate", action="store_true")
    pre_args, remaining = pre.parse_known_args(argv)
    if pre_args.calibrate:
        from calibrate import main as calibrate_main

        # Delegate: run.py is the single entrypoint, but calibrate.py is callable too.
        return int(calibrate_main(remaining))

    parser = argparse.ArgumentParser(
        prog="run.py",
        description="Run LS-EVPFFT simulations (writes inputs to a run directory and launches solver).",
        formatter_class=argparse.RawTextHelpFormatter,  # <-- preserves \n in help
    )

    parser.add_argument(
        "--sim_case",
        default=DEFAULT_SIM,
        choices=get_args(SimCase),
        type=str,
        help=textwrap.dedent(
            f"""\
            Type of simulation to run.

            One of:
              - creep_strain : strain-controlled creep
              - creep_stress, creep : stress-controlled creep
              - tension      : tension simulation
              - test         : short test case
              - srj          : strain-rate jump calibration
              - srj_old      : old srj calibration
              - creep_strain_calibration : strain controlled creep calibration
              - creep_calibration : stress controlled creep calibration

            Default: {DEFAULT_SIM}
            """
        ),
    )

    parser.add_argument(
        "--calib_id",
        default=None,
        type=int,
        help="Calibration id to pull parameters from",
    )
    parser.add_argument(
        "--load", type=int, default=0, help="max load in MPa (creep cases)"
    )
    parser.add_argument(
        "--complgas",
        type=float,
        default=None,
        help="gas compliance scalar",
    )
    parser.add_argument(
        "--igas",
        default=0,
        type=int,
        choices=get_args(GasMode),
    )
    parser.add_argument(
        "--nthreads", type=int, default=NTHREADS, help="solver thread count"
    )

    parser.add_argument(
        "--micro",
        type=Path,
        default=MICROSTRUCTURE_PATH,
        help="path to microstructure file",
    )

    parser.add_argument(
        "--mix",
        type=int,
        default=1,
        help="apply mixing to velgrad update",
    )
    parser.add_argument(
        "--modal",
        type=int,
        default=0,
        help="use modified augmented lagrangian",
    )

    parser.add_argument(
        "--notes",
        type=str,
        default="",
        help="free-form notes recorded in results/run_index.csv",
    )

    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="write inputs + fft.in but do not run the solver",
    )

    parser.add_argument(
        "--nlive",
        action="store_true",
        help="regenerate diagnostics plots during the run (every 5 minutes)",
    )
    parser.add_argument(
        "--out_dir",
        type=Path,
        default=None,
        help="path of the output directory (default: creep_modeling/results/SIM/INT/)",
    )

    parser.add_argument(
        "--evpfft_dir",
        type=Path,
        default=None,
        help="Path to the evpfft_outputs (default: out_dir/evpfft_outputs",
    )

    parser.add_argument(
        "--fftw_save",
        type=Path,
        default=None,
        help="path to FFTW wisdom save file (default: <out_dir>/evpfft_outputs/fftw_wisdom_omp.save)",
    )

    parser.add_argument(
        "--fft",
        type=Path,
        default=None,
        help="path to an fft.in instead of generating one",
    )

    parser.add_argument(
        "--out_txt",
        type=Path,
        default=None,
        help="path to output.txt log (default: <out_dir>/output.txt",
    )

    parser.add_argument(
        "--cuel",
        type=Path,
        default=None,
        help="path to cuel.sx instead of pulling one from calib_id",
    )

    parser.add_argument(
        "--cupl2",
        type=Path,
        default=None,
        help="path to cupl2.sx instead of pulling one from calib_id",
    )

    parser.add_argument(
        "--cuel_gas",
        type=Path,
        default=None,
        help="path to cuel_gas.sx instead of pulling one from calib_id",
    )

    parser.add_argument(
        "--cupl2_gas",
        type=Path,
        default=None,
        help="path to cupl2_gas.sx instead of pulling one from calib_id",
    )

    parser.add_argument(
        "--skip_vtk",
        type=bool,
        default=True,
        help="Skip adding vtk height data",
    )

    args = parser.parse_args(remaining)

    # TODO if calib_id is none and material files are set, use those
    # TODO: if calib_id is not none and material files are set,
    # log a warning and use the specified material files
    # TODO: if only some material files are set and calib_id is none,
    # use defaults for those that are not set
    # TODO: if some material files are set and calib_id is not none,
    # log a warning and use the set material files and the calib_id parameters for others

    run_id, rd = allocate_run_dir(args.sim_case, out_dir=args.out_dir)

    # write_cuel, write_cupl2, write_cupl2_gas, write_cuel_gas = [True, True, True, True]
    skip_fft = False
    if args.evpfft_dir is None:
        args.evpfft_dir = rd / "evpfft_outputs"
    if args.out_txt is None:
        args.out_txt = rd / "output.txt"
    if args.fft is None:
        args.fft = rd / "fft.in"
    else:
        skip_fft = True
    if args.cuel is None:
        args.cuel = rd / "cuel.sx"
    if args.cupl2 is None:
        args.cupl2 = rd / "cupl2.sx"
    if args.cuel_gas is None:
        args.cuel_gas = rd / "cuel_gas.sx"
    if args.cupl2_gas is None:
        args.cupl2_gas = rd / "cupl2_gas.sx"
    if args.igas == 0 or args.igas == 1:
        args.cuel_gas = args.cuel
        args.cupl2_gas = args.cupl2
    if args.micro is None:
        args.micro = MICROSTRUCTURE_PATH
    if args.fftw_save is None:
        # args.fftw_save = Path(args.evpfft_dir) / "fftw_wisdom_omp.save"
        args.fftw_save = Path("/Users/gtdebru/creep_modeling/results/fftw_wisdom_omp.save")

    micro_info = MicrostructureInfo.load(args.micro)
    sol_load = args.load
    if args.load is not None and args.igas > 0:
        args.load = convert_load(args.load, micro_info)
    if args.complgas is None:
        args.complgas = 0.01

    # Register run in the global index.
    register_run(
        run_id=run_id,
        sim_case=args.sim_case,
        calib_id=args.calib_id,
        igas=args.igas,
        nthreads=args.nthreads,
        notes=str(args.notes) if str(args.notes).strip() else None,
        microstructure=args.micro,
        fft=args.fft,
        run_dir=rd,
        load=args.load,
        complgas=args.complgas,
        mix=args.mix,
        modal=args.modal,
    )

    cfg = RunConfig(
        run_dir=rd,
        out_dir=rd,
        sim_case=args.sim_case,
        calib_id=args.calib_id,
        igas=args.igas,
        nthreads=args.nthreads,
        notes=str(args.notes) if str(args.notes).strip() else None,
        load=args.load,
        complgas=args.complgas,
        microstructure=args.micro,
        mix=args.mix,
        modal=args.modal,
        fft=args.fft,
        fftw_save=args.fftw_save,
        cuel=args.cuel,
        cupl2=args.cupl2,
        cuel_gas=args.cuel_gas,
        cupl2_gas=args.cupl2_gas,
        evpfft_dir=args.evpfft_dir,
        out_txt=args.out_txt,
        nx=micro_info.nx,
        ny=micro_info.ny,
        nz=micro_info.nz,
        nphase=micro_info.nphase,
        sol_load=sol_load,
    )

    write_inputs(cfg, skip_fft)
    print(f"run_dir={rd}")
    if args.dry_run:
        return 0

    rc: int = 1
    interrupted = False
    try:
        if args.nlive:
            rc = int(run_solver_nlive(cfg, v=True))
        else:
            rc = int(run_solver(cfg, v=True))
    except KeyboardInterrupt:
        interrupted = True
    finally:

        load = int(cfg.sol_load)

        postprocess_cfg = PostprocessConfig(
            out_dir=cfg.out_dir,
            sim_case=cfg.sim_case,
            microstructure=cfg.microstructure,
            igas=cfg.igas,
            load=load,
            skip_vtk=args.skip_vtk,
        )

        _postprocess(postprocess_cfg)

    if interrupted:
        raise KeyboardInterrupt
    return int(rc)


if __name__ == "__main__":

    raise SystemExit(main())

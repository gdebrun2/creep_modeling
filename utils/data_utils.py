from __future__ import annotations

from pathlib import Path
from typing import Iterable, Literal

import numpy as np
import pandas as pd
from pydantic import ConfigDict
from pydantic.dataclasses import dataclass as pydantic_dataclass
from calc_utils import (
    to_true_stress,
    to_true_strain,
    total_strain_to_plastic,
    calc_displacement,
    calc_max_slip,
    calc_height,
    calc_mean_slip,
    calculate_sa,
    calculate_sz,
    subtract_initial_plane,
    fit_plane,
    calc_exp_roughness,
)
from glob import glob
from io_utils import append_fields, write_vtk_series
from config import (
    CREEP_CSV_PATH,
    SRJ_CSV_PATH,
    TENSILE_CSV_PATH,
    MICROSTRUCTURE_PATH,
    PHASE_GAS,
    PHASE_SOLID,
    CREEP_PROFILOMETRY_PATH,
)

SRJ_SECTION_COLS: tuple[str, ...] = (
    "Section1Strain",
    "Section2Strain",
    "Section3Strain",
    "Section4Strain",
    "Section5Strain",
)

SIM_RESULTS_COLS: tuple[str, ...] = (
    "iter",
    "macrostep",
    "ipr",
    "sim_time",
    "wall_time",
    "tdot",
    "tdot_ref",
    "eav11",
    "eav22",
    "eav33",
    "eav23",
    "eav13",
    "eav12",
    "sav11",
    "sav22",
    "sav33",
    "sav23",
    "sav13",
    "sav12",
    "epav11",
    "epav22",
    "epav33",
    "epav23",
    "epav13",
    "epav12",
    "eelav11",
    "eelav22",
    "eelav33",
    "eelav23",
    "eelav13",
    "eelav12",
    "vgrad11",
    "vgrad22",
    "vgrad33",
    "vgrad23",
    "vgrad13",
    "vgrad12",
    "edotp11",
    "edotp22",
    "edotp33",
    "edotp23",
    "edotp13",
    "edotp12",
    "evm",
    "evmp",
    "dvm",
    "dvmp",
    "svm",
    "gamma",
    "g",
    "geff",
)

SIM_RESULTS_INT_COLS = frozenset({"iter", "macrostep", "ipr"})

DIAGNOSTICS_COLS: tuple[str, ...] = (
    "iter",
    "macrostep",
    "ipr",
    "sim_time",
    "wall_time",
    "tdot",
    "vm_max_inc",
    "tdotref",
    "erre",
    "errs",
    "errsbc",
    "s33_target",
    "s33",
    "s33_sol",
    "detF_min_sol",
    "detF_min_gas",
    "detF_min_sol_rel",
    "detF_min_gas_rel",
    "volinc_max_sol",
    "volinc_max_gas",
    "volinc_max_sol_rel",
    "volinc_max_gas_rel",
    "stretchinc_max_sol",
    "stretchinc_max_gas",
    "stretchinc_max_sol_rel",
    "stretchinc_max_gas_rel",
    "rotang_max_sol",
    "rotang_max_gas",
    "dL_step_rms_sol",
    "dL_step_rms_gas",
    "dL_step_rms",
    "dL_step_max_sol",
    "dL_step_max_gas",
    "dL_step_max",
    "dvelgradavg_mag",
    "cond_c066mod_max_sol",
    "cond_c066mod_max_gas",
    "cond_cgas_max",
    "imin_sol_i",
    "imin_sol_j",
    "imin_sol_k",
    "imin_gas_i",
    "imin_gas_j",
    "imin_gas_k",
)

DIAGNOSTICS_INT_COLS = frozenset(
    {
        "iter",
        "macrostep",
        "ipr",
        "imin_sol_i",
        "imin_sol_j",
        "imin_sol_k",
        "imin_gas_i",
        "imin_gas_j",
        "imin_gas_k",
    }
)


def _load_csv_arrays(
    csv_path: Path | str,
    required_cols: tuple[str, ...],
    int_cols: frozenset[str],
) -> dict[str, np.ndarray]:
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()

    missing = tuple(col for col in required_cols if col not in df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {csv_path}: {missing}")

    arrays: dict[str, np.ndarray] = {}
    for col in required_cols:
        dtype = int if col in int_cols else float
        arrays[col] = df[col].to_numpy(dtype=dtype)
    return arrays


@pydantic_dataclass(
    config=ConfigDict(
        frozen=True, extra="forbid", strict=True, arbitrary_types_allowed=True
    )
)
class SimResults:
    output_path: Path | str
    csv_path: Path | str
    iter: np.ndarray
    macrostep: np.ndarray
    ipr: np.ndarray
    sim_time: np.ndarray
    wall_time: np.ndarray
    tdot: np.ndarray
    tdot_ref: np.ndarray
    eav11: np.ndarray
    eav22: np.ndarray
    eav33: np.ndarray
    eav23: np.ndarray
    eav13: np.ndarray
    eav12: np.ndarray
    sav11: np.ndarray
    sav22: np.ndarray
    sav33: np.ndarray
    sav23: np.ndarray
    sav13: np.ndarray
    sav12: np.ndarray
    epav11: np.ndarray
    epav22: np.ndarray
    epav33: np.ndarray
    epav23: np.ndarray
    epav13: np.ndarray
    epav12: np.ndarray
    eelav11: np.ndarray
    eelav22: np.ndarray
    eelav33: np.ndarray
    eelav23: np.ndarray
    eelav13: np.ndarray
    eelav12: np.ndarray
    vgrad11: np.ndarray
    vgrad22: np.ndarray
    vgrad33: np.ndarray
    vgrad23: np.ndarray
    vgrad13: np.ndarray
    vgrad12: np.ndarray
    edotp11: np.ndarray
    edotp22: np.ndarray
    edotp33: np.ndarray
    edotp23: np.ndarray
    edotp13: np.ndarray
    edotp12: np.ndarray
    evm: np.ndarray
    evmp: np.ndarray
    dvm: np.ndarray
    dvmp: np.ndarray
    svm: np.ndarray
    vtk_macrostep: np.ndarray | None
    vtk_time: np.ndarray | None
    mean_slip: np.ndarray | None
    max_slip: np.ndarray | None
    sa: np.ndarray | None
    sz: np.ndarray | None
    info: MicrostructureInfo | None
    gamma: np.ndarray | None
    g: np.ndarray
    geff: np.ndarray

    @classmethod
    def load(
        cls,
        output_path: Path | str,  # evpfft dir
        micro_path: Path | str = MICROSTRUCTURE_PATH,
        skip_vtk: bool = False,
    ) -> SimResults:
        """Load a solver `results.csv` file."""
        output_path = Path(output_path)
        csv_path = output_path / "results.csv"

        arrays = _load_csv_arrays(csv_path, SIM_RESULTS_COLS, SIM_RESULTS_INT_COLS)

        # discard diverged steps
        n_steps = arrays["sim_time"].size
        step_ok = (
            np.isfinite(arrays["eav33"])
            & np.isfinite(arrays["sav33"])
            & np.isfinite(arrays["epav33"])
        )
        bad = np.where(~step_ok)[0]

        if skip_vtk:
            info = None
            vtk_macrostep = None
            vtk_time = None
            sa = None
            sz = None
            max_slip = None
            mean_slip = None
        else:
            info = MicrostructureInfo.load(micro_path)
            vtk_macrostep, vtk_time, sa, sz, max_slip, mean_slip = vtk_sweep(
                output_path,
                arrays["sim_time"],
                info,
            )

        if len(bad) > 0:
            n_steps = bad[0]
            for key, arr in list(arrays.items()):
                arrays[key] = arr[:n_steps]

            if vtk_macrostep is not None:
                vtk_macrostep = vtk_macrostep[:n_steps]
                vtk_time = vtk_time[:n_steps]
                sa = sa[:n_steps]
                sz = sz[:n_steps]
                max_slip = max_slip[:n_steps]
                mean_slip = mean_slip[:n_steps]

        return cls(
            output_path=output_path,
            csv_path=csv_path,
            info=info,
            vtk_macrostep=vtk_macrostep,
            vtk_time=vtk_time,
            sa=sa,
            sz=sz,
            max_slip=max_slip,
            mean_slip=mean_slip,
            **arrays,
        )


@pydantic_dataclass(
    config=ConfigDict(
        frozen=True, extra="forbid", strict=True, arbitrary_types_allowed=True
    )
)
class Diagnostics:
    output_dir: Path | str
    csv_path: Path | str
    iter: np.ndarray
    macrostep: np.ndarray
    ipr: np.ndarray
    sim_time: np.ndarray
    wall_time: np.ndarray
    tdot: np.ndarray
    vm_max_inc: np.ndarray
    tdotref: np.ndarray
    erre: np.ndarray
    errs: np.ndarray
    errsbc: np.ndarray
    s33_target: np.ndarray
    s33: np.ndarray
    s33_sol: np.ndarray
    detF_min_sol: np.ndarray
    detF_min_gas: np.ndarray
    detF_min_sol_rel: np.ndarray
    detF_min_gas_rel: np.ndarray
    volinc_max_sol: np.ndarray
    volinc_max_gas: np.ndarray
    volinc_max_sol_rel: np.ndarray
    volinc_max_gas_rel: np.ndarray
    stretchinc_max_sol: np.ndarray
    stretchinc_max_gas: np.ndarray
    stretchinc_max_sol_rel: np.ndarray
    stretchinc_max_gas_rel: np.ndarray
    rotang_max_sol: np.ndarray
    rotang_max_gas: np.ndarray
    dL_step_rms_sol: np.ndarray
    dL_step_rms_gas: np.ndarray
    dL_step_rms: np.ndarray
    dL_step_max_sol: np.ndarray
    dL_step_max_gas: np.ndarray
    dL_step_max: np.ndarray
    dvelgradavg_mag: np.ndarray
    cond_c066mod_max_sol: np.ndarray
    cond_c066mod_max_gas: np.ndarray
    cond_cgas_max: np.ndarray
    imin_sol_i: np.ndarray
    imin_sol_j: np.ndarray
    imin_sol_k: np.ndarray
    imin_gas_i: np.ndarray
    imin_gas_j: np.ndarray
    imin_gas_k: np.ndarray

    @classmethod
    def load(cls, output_path: Path | str) -> Diagnostics:
        """Load a solver `diagnostics.csv` file."""
        output_path = Path(output_path)
        csv_path = output_path / "diagnostics.csv"
        arrays = _load_csv_arrays(csv_path, DIAGNOSTICS_COLS, DIAGNOSTICS_INT_COLS)
        return cls(output_dir=output_path, csv_path=csv_path, **arrays)


@pydantic_dataclass(
    config=ConfigDict(
        frozen=True, extra="forbid", strict=True, arbitrary_types_allowed=True
    )
)
class TensileData:
    csv_path: Path | str
    time: np.ndarray
    strain1: np.ndarray
    strain2: np.ndarray
    strain3: np.ndarray
    stress1: np.ndarray
    stress2: np.ndarray
    stress3: np.ndarray
    mean_stress: np.ndarray
    mean_strain: np.ndarray
    std_stress: np.ndarray
    std_strain: np.ndarray
    force: np.ndarray

    @classmethod
    def load(
        cls,
        csv_path: Path = TENSILE_CSV_PATH,
    ) -> TensileData:
        """Load a tensile sample and compute engineering stress for high/mid/low areas.

        Ports the exact logic in `stress_strain_calc.py` but with repo-relative paths.

        Expected CSV column convention (legacy):
        - col 1: force (N)
        - col 6: strain high
        - col 5: strain mid
        - col 4: strain low

        Cross-sectional areas (m^2) from legacy script:
        - high: 0.00216 * 0.00255
        - mid:  0.00243 * 0.00252
        - low:  0.00263 * 0.00249

        Stress is returned in MPa.
        """

        df = pd.read_csv(csv_path, dtype=float)
        time = df["Time"].values
        force = df["force"].values
        strain1 = df["Strain0"].values
        strain2 = df["Strain1"].values
        strain3 = df["Strain2"].values

        area1 = 0.00216 * 0.00255
        area2 = 0.00243 * 0.00252
        area3 = 0.00263 * 0.00249

        # N/m^2 -> Pa -> MPa
        stress1 = (force / area1) / 1e6
        stress2 = (force / area2) / 1e6
        stress3 = (force / area3) / 1e6

        stress1 = to_true_stress(stress1, strain1)
        stress2 = to_true_stress(stress2, strain2)
        stress3 = to_true_stress(stress3, strain3)

        strain1 = to_true_strain(strain1)
        strain2 = to_true_strain(strain2)
        strain3 = to_true_strain(strain3)

        stacked_stress = np.vstack([stress1, stress2, stress3])
        stacked_strain = np.vstack([strain1, strain2, strain3])
        mean_stress = np.mean(stacked_stress, axis=0)
        mean_strain = np.mean(stacked_strain, axis=0)
        std_stress = np.std(stacked_stress, axis=0)
        std_strain = np.std(stacked_strain, axis=0)

        return cls(
            time=time,
            strain1=strain1,
            strain2=strain2,
            strain3=strain3,
            stress1=stress1,
            stress2=stress2,
            stress3=stress3,
            mean_strain=mean_strain,
            mean_stress=mean_stress,
            std_stress=std_stress,
            std_strain=std_strain,
            force=force,
        )


@pydantic_dataclass(
    config=ConfigDict(
        frozen=True, extra="forbid", strict=True, arbitrary_types_allowed=True
    )
)
class CreepData:

    time: np.ndarray
    csv_path: Path | str
    mean_strain_475: np.ndarray
    mean_strain_525: np.ndarray
    mean_strain_575: np.ndarray
    std_strain_475: np.ndarray
    std_strain_525: np.ndarray
    std_strain_575: np.ndarray
    mean_plastic_strain_475: np.ndarray
    mean_plastic_strain_525: np.ndarray
    mean_plastic_strain_575: np.ndarray
    std_plastic_strain_475: np.ndarray
    std_plastic_strain_525: np.ndarray
    std_plastic_strain_575: np.ndarray
    stress_475: np.ndarray
    stress_525: np.ndarray
    stress_575: np.ndarray
    roughness_time: np.ndarray
    mean_sz_475_10x: np.ndarray
    mean_sz_525_10x: np.ndarray
    mean_sz_575_10x: np.ndarray
    std_sz_475_10x: np.ndarray
    std_sz_525_10x: np.ndarray
    std_sz_575_10x: np.ndarray
    mean_sa_475_10x: np.ndarray
    mean_sa_525_10x: np.ndarray
    mean_sa_575_10x: np.ndarray
    std_sa_475_10x: np.ndarray
    std_sa_525_10x: np.ndarray
    std_sa_575_10x: np.ndarray
    mean_sz_475_50x: np.ndarray
    mean_sz_525_50x: np.ndarray
    mean_sz_575_50x: np.ndarray
    std_sz_475_50x: np.ndarray
    std_sz_525_50x: np.ndarray
    std_sz_575_50x: np.ndarray
    mean_sa_475_50x: np.ndarray
    mean_sa_525_50x: np.ndarray
    mean_sa_575_50x: np.ndarray
    std_sa_475_50x: np.ndarray
    std_sa_525_50x: np.ndarray
    std_sa_575_50x: np.ndarray

    @classmethod
    def load(
        cls,
        csv_path: Path | str = CREEP_CSV_PATH,
        profilometry_path: Path | str = CREEP_PROFILOMETRY_PATH,
        skip_roughness: bool = False,
    ) -> CreepData:
        """Compute true-strain mean±std time histories for the creep experiment.

        The creep CSV contains multiple samples/columns for each load level.
        - columns ending in "a" correspond to the 95% condition (475 MPa)
        - columns ending in "b" correspond to the 105% condition (525 MPa)
        - columns ending in "c" correspond to the 115% condition (575 MPa)

        Missing entries (NaN) are ignored and the mean is computed over available
        samples at each time.
        """

        df = pd.read_csv(csv_path, dtype=float)
        df.columns = df.columns.str.strip()
        cols = df.columns.values
        time = df["Time"].values * 3600

        strain_475: list[np.ndarray] = []
        strain_525: list[np.ndarray] = []
        strain_575: list[np.ndarray] = []

        for col in cols:
            if col == "Time":
                continue
            if col.endswith("a"):
                strain_475.append(df[col].values)
            elif col.endswith("b"):
                strain_525.append(df[col].values)
            elif col.endswith("c"):
                strain_575.append(df[col].values)

        stacked_strain_475 = np.vstack(strain_475).T  # (Nt, Na)
        stacked_strain_525 = np.vstack(strain_525).T  # (Nt, Nb)
        stacked_strain_575 = np.vstack(strain_575).T  # (Nt, Nc)

        # it's the same row for all loads
        bad_row = np.all(np.isnan(stacked_strain_475), axis=1)
        stacked_strain_475 = stacked_strain_475[~bad_row, :]
        stacked_strain_525 = stacked_strain_525[~bad_row, :]
        stacked_strain_575 = stacked_strain_575[~bad_row, :]
        time = time[~bad_row]

        mean_strain_475 = np.nanmean(stacked_strain_475, axis=1)
        mean_strain_525 = np.nanmean(stacked_strain_525, axis=1)
        mean_strain_575 = np.nanmean(stacked_strain_575, axis=1)
        std_strain_475 = np.nanstd(stacked_strain_475, axis=1)
        std_strain_525 = np.nanstd(stacked_strain_525, axis=1)
        std_strain_575 = np.nanstd(stacked_strain_575, axis=1)

        stress_475 = np.ones_like(mean_strain_475) * 475
        stress_525 = np.ones_like(mean_strain_525) * 525
        stress_575 = np.ones_like(mean_strain_575) * 575
        stacked_plastic_strain_475 = total_strain_to_plastic(
            stacked_strain_475, ref_index=1
        )
        stacked_plastic_strain_525 = total_strain_to_plastic(
            stacked_strain_525, ref_index=1
        )
        stacked_plastic_strain_575 = total_strain_to_plastic(
            stacked_strain_575, ref_index=1
        )

        mean_plastic_strain_475 = np.nanmean(stacked_plastic_strain_475, axis=1)
        mean_plastic_strain_525 = np.nanmean(stacked_plastic_strain_525, axis=1)
        mean_plastic_strain_575 = np.nanmean(stacked_plastic_strain_575, axis=1)
        std_plastic_strain_475 = np.nanstd(stacked_plastic_strain_475, axis=1)
        std_plastic_strain_525 = np.nanstd(stacked_plastic_strain_525, axis=1)
        std_plastic_strain_575 = np.nanstd(stacked_plastic_strain_575, axis=1)

        if not skip_roughness:
            roughness = load_exp_roughness(path=profilometry_path)
        else:
            roughness = {
                "roughness_time": np.zeros(1),
                "mean_sz_475_10x": np.zeros(1),
                "mean_sz_525_10x": np.zeros(1),
                "mean_sz_575_10x": np.zeros(1),
                "std_sz_475_10x": np.zeros(1),
                "std_sz_525_10x": np.zeros(1),
                "std_sz_575_10x": np.zeros(1),
                "mean_sa_475_10x": np.zeros(1),
                "mean_sa_525_10x": np.zeros(1),
                "mean_sa_575_10x": np.zeros(1),
                "std_sa_475_10x": np.zeros(1),
                "std_sa_525_10x": np.zeros(1),
                "std_sa_575_10x": np.zeros(1),
                "mean_sz_475_50x": np.zeros(1),
                "mean_sz_525_50x": np.zeros(1),
                "mean_sz_575_50x": np.zeros(1),
                "std_sz_475_50x": np.zeros(1),
                "std_sz_525_50x": np.zeros(1),
                "std_sz_575_50x": np.zeros(1),
                "mean_sa_475_50x": np.zeros(1),
                "mean_sa_525_50x": np.zeros(1),
                "mean_sa_575_50x": np.zeros(1),
                "std_sa_475_50x": np.zeros(1),
                "std_sa_525_50x": np.zeros(1),
                "std_sa_575_50x": np.zeros(1),
            }

        return cls(
            time=time,
            csv_path=csv_path,
            mean_strain_475=mean_strain_475,
            mean_strain_525=mean_strain_525,
            mean_strain_575=mean_strain_575,
            std_strain_475=std_strain_475,
            std_strain_525=std_strain_525,
            std_strain_575=std_strain_575,
            stress_475=stress_475,
            stress_525=stress_525,
            stress_575=stress_575,
            mean_plastic_strain_475=mean_plastic_strain_475,
            mean_plastic_strain_525=mean_plastic_strain_525,
            mean_plastic_strain_575=mean_plastic_strain_575,
            std_plastic_strain_475=std_plastic_strain_475,
            std_plastic_strain_525=std_plastic_strain_525,
            std_plastic_strain_575=std_plastic_strain_575,
            **roughness,
        )

    def get_load_vals(self, load: int) -> dict[str, np.ndarray]:

        if load == 475:

            return {
                "time": self.time,
                "roughness_time": self.roughness_time,
                "mean_strain": self.mean_strain_475,
                "std_strain": self.std_strain_475,
                "mean_plastic_strain": self.mean_plastic_strain_475,
                "std_plastic_strain": self.std_plastic_strain_475,
                "stress": self.stress_475,
                "mean_sa_10x": self.mean_sa_475_10x,
                "std_sa_10x": self.std_sa_475_10x,
                "mean_sa_50x": self.mean_sa_475_50x,
                "std_sa_50x": self.std_sa_475_50x,
                "mean_sz_10x": self.mean_sz_475_10x,
                "std_sz_10x": self.std_sz_475_10x,
                "mean_sz_50x": self.mean_sz_475_50x,
                "std_sz_50x": self.std_sz_475_50x,
            }
        elif load == 525:

            return {
                "time": self.time,
                "roughness_time": self.roughness_time,
                "mean_strain": self.mean_strain_525,
                "std_strain": self.std_strain_525,
                "mean_plastic_strain": self.mean_plastic_strain_525,
                "std_plastic_strain": self.std_plastic_strain_525,
                "stress": self.stress_525,
                "mean_sa_10x": self.mean_sa_525_10x,
                "std_sa_10x": self.std_sa_525_10x,
                "mean_sa_50x": self.mean_sa_525_50x,
                "std_sa_50x": self.std_sa_525_50x,
                "mean_sz_10x": self.mean_sz_525_10x,
                "std_sz_10x": self.std_sz_525_10x,
                "mean_sz_50x": self.mean_sz_525_50x,
                "std_sz_50x": self.std_sz_525_50x,
            }

        elif load == 575:

            return {
                "time": self.time,
                "roughness_time": self.roughness_time,
                "mean_strain": self.mean_strain_575,
                "std_strain": self.std_strain_575,
                "mean_plastic_strain": self.mean_plastic_strain_575,
                "std_plastic_strain": self.std_plastic_strain_575,
                "stress": self.stress_575,
                "mean_sa_10x": self.mean_sa_575_10x,
                "std_sa_10x": self.std_sa_575_10x,
                "mean_sa_50x": self.mean_sa_575_50x,
                "std_sa_50x": self.std_sa_575_50x,
                "mean_sz_10x": self.mean_sz_575_10x,
                "std_sz_10x": self.std_sz_575_10x,
                "mean_sz_50x": self.mean_sz_575_50x,
                "std_sz_50x": self.std_sz_575_50x,
            }

        else:
            raise ValueError("load must be one of 475, 525 or 575")


@pydantic_dataclass(
    config=ConfigDict(
        frozen=True, extra="forbid", strict=True, arbitrary_types_allowed=True
    )
)
class SrjData:
    csv_path: Path | str
    time: np.ndarray
    mean_strain: np.ndarray
    mean_stress: np.ndarray
    std_strain: np.ndarray
    std_stress: np.ndarray
    section_stress: dict[str, np.ndarray]
    section_strain: dict[str, np.ndarray]

    @classmethod
    def load(
        cls,
        csv_path: Path | str = SRJ_CSV_PATH,
        section_strain_cols: Iterable[str] = SRJ_SECTION_COLS,
    ) -> SrjData:

        # the SRJ sample broke after t=643
        df = pd.read_csv(csv_path, dtype=float).iloc[:643]
        df.columns = df.columns.str.strip()

        time = df["Time"].values
        stress = df["Stress(MPa)"].values

        section_true_stress: dict[str, np.ndarray] = {}
        section_true_strain: dict[str, np.ndarray] = {}

        for col in section_strain_cols:
            strain = df[col].values
            section_true_stress[col] = to_true_stress(stress, strain)
            section_true_strain[col] = to_true_strain(strain)

        stacked_strain = np.vstack(
            [section_true_strain[c] for c in section_strain_cols]
        )
        stacked_stress = np.vstack(
            [section_true_stress[c] for c in section_strain_cols]
        )

        mean_true_strain = np.mean(stacked_strain, axis=0)
        mean_true_stress = np.mean(stacked_stress, axis=0)
        std_true_strain = np.std(stacked_strain, axis=0)
        std_true_stress = np.std(stacked_stress, axis=0)

        return cls(
            csv_path=csv_path,
            time=time,
            mean_strain=mean_true_strain,
            mean_stress=mean_true_stress,
            std_strain=std_true_strain,
            std_stress=std_true_stress,
            section_stress=section_true_stress,
            section_strain=section_true_strain,
        )


@pydantic_dataclass(
    config=ConfigDict(
        frozen=True, extra="forbid", strict=True, arbitrary_types_allowed=True
    )
)
class MicrostructureInfo:
    """Basic microstructure metadata inferred from `microstr_*.dat`."""

    path: str | Path
    nx: int
    ny: int
    nz: int
    nphase: int
    # 1D arrays length = nx*ny*nz
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    grain_id: np.ndarray
    ph_id: np.ndarray
    solid_frac: float
    gas_frac: float
    gas_phase: int
    sol_phase: int

    @classmethod
    def load(cls, path: Path | str = MICROSTRUCTURE_PATH) -> MicrostructureInfo:
        """Read microstructure file and infer nx/2/3 and phase ids.

        Expected column meanings per your provided snippet:
        - col0-2: Euler angles (ignored here)
        - col3: x (1-based integer)
        - col4: y (1-based integer)
        - col5: z (1-based integer)
        - col6: grain_id (int)
        - col7: ph_id (int) where 1=solid and 2=gas
        """

        arr = np.loadtxt(path)
        if arr.ndim != 2 or arr.shape[1] < 8:
            raise ValueError(f"Unexpected microstructure format: shape={arr.shape}")

        x = arr[:, 3].astype(int)
        y = arr[:, 4].astype(int)
        z = arr[:, 5].astype(int)
        grain_id = arr[:, 6].astype(int)
        ph_id = arr[:, 7].astype(int)

        nx = int(x.max())
        ny = int(y.max())
        nz = int(z.max())
        nphase = np.unique(ph_id).size

        solid_frac = float(np.mean(ph_id == PHASE_SOLID))
        gas_frac = float(np.mean(ph_id == PHASE_GAS))

        return cls(
            path=path,
            nx=nx,
            ny=ny,
            nz=nz,
            nphase=nphase,
            x=x,
            y=y,
            z=z,
            grain_id=grain_id,
            ph_id=ph_id,
            solid_frac=solid_frac,
            gas_frac=gas_frac,
            gas_phase=PHASE_GAS,
            sol_phase=PHASE_SOLID,
        )


def srj_time_offset(
    *,
    time: np.ndarray,
    stress: np.ndarray,
    stress_lo: float = 200.0,
    stress_hi: float = 400.0,
) -> float:
    """Estimate an SRJ experimental time offset by fitting the elastic ramp.

    Motivation
    ----------
    The SRJ stress-time curve can exhibit an initial nonlinear "toe" (e.g.,
    machine compliance / seating). This can delay the recorded curve in time
    relative to the imposed strain-rate schedule.

    Method
    ------
    Take experimental stress values within a specified window (default:
    [200, 300] MPa), fit a line:

        stress(t) = m*t + b

    and define the offset as the extrapolated time where that fitted line would
    reach zero stress:

        t0 = -b/m

    The experiment can then be plotted/scored on a shifted axis:

        t_shift = t - t0

    so the elastic ramp starts near t_shift = 0 and the initial toe appears at
    negative time.

    Returns
    -------
    t0_s:
        Estimated offset time in seconds. Returns 0.0 if the fit is ill-posed
        (insufficient points, nonpositive slope, or non-finite results).
    """

    t = time.copy()
    s = stress.copy()
    if t.ndim != 1 or s.ndim != 1 or len(t) != len(s) or len(t) == 0:
        return 0.0

    ok = np.isfinite(t) & np.isfinite(s)
    t = t[ok]
    s = s[ok]
    if len(t) < 2:
        return 0.0

    in_win = (s >= float(stress_lo)) & (s <= float(stress_hi))
    idx = np.where(in_win)[0]
    if len(idx) < 2:
        return 0.0

    # Prefer the earliest contiguous segment in the window (initial loading).
    splits = np.where(np.diff(idx) > 1)[0]
    segments = np.split(idx, splits + 1)
    seg = None
    for cand in segments:
        if len(cand) >= 2:
            seg = cand
            break
    if seg is None:
        return 0.0

    t_seg = t[seg]
    s_seg = s[seg]
    if len(t_seg) < 2:
        return 0.0

    try:
        m, b = np.polyfit(t_seg, s_seg, deg=1)
    except Exception:
        return 0.0

    if not np.isfinite(m) or not np.isfinite(b) or float(m) <= 0.0:
        return 0.0

    t0 = float(-b / m)
    if not np.isfinite(t0):
        return 0.0

    # Sanity: offset should be nonnegative and should precede the window start.
    t0 = max(0.0, t0)
    t_min = float(np.min(t_seg))
    if t0 > t_min:
        t0 = t_min
    return t0


def read_vtk_points_slip(
    vtk_path: str, info: MicrostructureInfo
) -> tuple[np.ndarray, np.ndarray]:

    nx = info.nx
    ny = info.ny
    nz = info.nz
    Nnodes = (nz + 1) * (ny + 1) * (nx + 1)
    Ncells = nz * ny * nx

    with open(vtk_path, "r") as f:
        it = iter(f)
        for line in it:
            if line.strip().startswith("POINTS"):
                # (Nnodes, 3) (z,y,x)
                points = np.loadtxt(f, max_rows=Nnodes, dtype=np.float64).reshape(
                    nz + 1, ny + 1, nx + 1, 3
                )
            elif line.strip().startswith("gamdot_acum"):
                slip = np.loadtxt(f, max_rows=Ncells, dtype=np.float64).reshape(
                    nz, ny, nx
                )

                return points, slip


def vtk_sweep(
    output_path: str | Path,
    sim_time: np.ndarray,
    info: MicrostructureInfo,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    output_path = Path(output_path)
    files = sorted(glob(rf"{output_path}/*.vtk"))
    sas = np.zeros((len(files), 4))
    szs = np.zeros((len(files), 4))
    max_slips = np.zeros((len(files), 4))
    mean_slips = np.zeros((len(files), 4))
    vtk_macrostep = np.zeros(len(files), dtype=int)

    for i, file in enumerate(files):
        macrostep = int(file.split("_")[-1].split(".")[0])
        vtk_macrostep[i] = macrostep
        points, slip = read_vtk_points_slip(file, info)

        if i == 0:
            points0 = points.copy()

        disp = calc_displacement(points, points0, info, mode="node")
        height = calc_height(points, info, mode="node")
        sa = calculate_sa(height, info, mode="node")
        sz = calculate_sz(height, info, mode="node")
        max_slip = calc_max_slip(slip, info, mode="cell")
        mean_slip = calc_mean_slip(slip, info, mode="cell")

        sas[i] = np.array(sa)
        szs[i] = np.array(sz)
        max_slips[i] = np.array(max_slip)
        mean_slips[i] = np.array(mean_slip)

        append_fields(
            file,
            point_fields={
                "node_disp": disp.flatten(),
                "node_height": height.flatten(),
            },
        )

    vtk_time = sim_time[vtk_macrostep]

    write_vtk_series(output_path, vtk_time)
    return vtk_macrostep, vtk_time, sas.T, szs.T, max_slips.T, mean_slips.T


def read_creep_profilometry(
    *, path: str | Path = CREEP_PROFILOMETRY_PATH, resolution: Literal[10, 50]
) -> dict[str, np.ndarray]:

    files = sorted(
        glob(f"{str(path)}/*{resolution}x_Height.csv"),
        key=lambda x: (
            float(Path(x).name.split("_")[0].strip("hrs")),
            Path(x).name.split("_")[1],
        ),
    )

    samples_to_idx = {"32": 0, "34": 1, "35": 2, "36": 3, "38": 4, "40": 5}

    times = np.array(
        sorted(
            {float(Path(f).name.split("_")[0].strip("hrs")) * 3600.0 for f in files}
        ),
        dtype=float,
    )
    time_to_idx = {t: i for i, t in enumerate(times)}

    ntimes = len(times)
    nsamples = len(samples_to_idx)
    ny, nx = 768, 1024

    height_475 = np.full((ntimes, nsamples, ny, nx), np.nan, dtype=float)
    height_525 = np.full((ntimes, nsamples, ny, nx), np.nan, dtype=float)
    height_575 = np.full((ntimes, nsamples, ny, nx), np.nan, dtype=float)

    for file in files:
        name = Path(file).name
        time_str, sample = name.split("_")[:2]

        t = float(time_str.strip("hrs")) * 3600.0
        t_idx = time_to_idx[t]

        sample_num = sample[:-1]  # "32a" -> "32"
        suffix = sample[-1]  # "32a" -> "a"
        sample_idx = samples_to_idx[sample_num]

        height = pd.read_csv(file, skiprows=19, header=None, dtype=float).to_numpy()
        height -= fit_plane(height)

        if height.shape != (ny, nx):
            raise ValueError(f"{name} has shape {height.shape}, expected {(ny, nx)}")

        if suffix == "a":
            height_475[t_idx, sample_idx] = height
        elif suffix == "b":
            height_525[t_idx, sample_idx] = height
        elif suffix == "c":
            height_575[t_idx, sample_idx] = height
        else:
            raise ValueError(f"Unexpected suffix in filename: {name}")

    # height_475 = subtract_initial_plane(height_475)
    # height_525 = subtract_initial_plane(height_525)
    # height_575 = subtract_initial_plane(height_575)

    heights = {
        "height_475": height_475,
        "height_525": height_525,
        "height_575": height_575,
    }

    return times, heights


def load_exp_roughness(
    *,
    path: str | Path = CREEP_PROFILOMETRY_PATH,
) -> dict[str, np.ndarray]:
    times, heights_10x = read_creep_profilometry(path=path, resolution=10)
    roughness_10x = calc_exp_roughness(heights_10x, resolution=10)
    times, heights_50x = read_creep_profilometry(path=path, resolution=50)
    roughness_50x = calc_exp_roughness(heights_50x, resolution=50)

    roughness = {"roughness_time": times, **roughness_10x, **roughness_50x}

    return roughness

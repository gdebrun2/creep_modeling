from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal
import os
import datetime as _dt
import numpy as np
from pathlib import Path
from skopt.space import Real

REPO_ROOT = Path(__file__).resolve().parents[1]
# If this is running from a tmp/ mirror, prefer the parent repo root.
if not (REPO_ROOT / "lsevpfft").exists() and (REPO_ROOT.parent / "lsevpfft").exists():
    REPO_ROOT = REPO_ROOT.parent

LSEVPFFT_DIR = REPO_ROOT / "lsevpfft"
LSEVPFFT_SRC_DIR = LSEVPFFT_DIR / "src"
MICROSTRUCTURE_DIR = REPO_ROOT / "microstructures"
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = REPO_ROOT / "results"
RESULTS_NEXT_ID_FILE = RESULTS_DIR / "_next_id.txt"
RUN_INDEX_FILE = RESULTS_DIR / "run_index.csv"
CALIBRATIONS_DIR = RESULTS_DIR / "calibrations"
CALIBRATIONS_NEXT_ID_FILE = CALIBRATIONS_DIR / "_next_id.txt"
CALIBRATION_INDEX_FILE = CALIBRATIONS_DIR / "calibration_index.csv"

# Data Files
SRJ_CSV_PATH = DATA_DIR / "CPValidation_316L" / "StrainRateJump.csv"
CREEP_CSV_PATH = DATA_DIR / "creep_int_95_105_115.csv"
TENSILE_CSV_PATH = DATA_DIR / "33_Tensile_Strain" / "11.csv"
CREEP_PROFILOMETRY_PATH = DATA_DIR / "creep_profilometry"

creep_stress_RESULTS_DIR = RESULTS_DIR / "creep_stress"
creep_strain_RESULTS_DIR = RESULTS_DIR / "creep_strain"
TENSION_RESULTS_DIR = RESULTS_DIR / "tension"
SRJ_RESULTS_DIR = RESULTS_DIR / "srj"
TEST_RESULTS_DIR = RESULTS_DIR / "test"

SIMCASE = Literal[
    "creep",
    "tension",
    "test",
    "srj_calibration",
    "creep_calibration",
]
GASMODE = Literal[0, 1, 2]

MICROSTRUCTURE_PATH = MICROSTRUCTURE_DIR / "micro1.dat"
SOLVER_EXE = LSEVPFFT_DIR / "LS-EVPFFT"

SRJ_OFFSET = 5.25  # s
SRJ_OFFSET_STRESS = 191.90175809309824
NTHREADS = os.cpu_count()


# Default calibration parameter set id used by run.py CLI
DEFAULT_PARAMSET = 1
DEFAULT_SIM = "creep"

PHASE_SOLID = 1
PHASE_GAS = 2


######################### DEFAULT ######################
DEFAULT_C11 = 206000  # C11, axial yield stress / axial stiffness
DEFAULT_C12 = 133000  # C12, lateral yield stress / lateral coupling
DEFAULT_C44 = 119000  # C44, yield shear stress / shear stiffness

DEFAULT_rb = 0.85

DEFAULT_tau0xf = 160.0
DEFAULT_tau0xb = DEFAULT_rb * DEFAULT_tau0xf
DEFAULT_tau1x = 2.0
DEFAULT_thet0 = 1000.0
DEFAULT_thet1 = 200
DEFAULT_nrsx = 30
DEFAULT_gamd0x = 1e-4


DEFAULT_hselfx = 1
DEFAULT_hlatex = 1.4  # Ma et al 2022, peirce et al 1983

WS = 1.0  # srj multiobjective weight
WC = 1.0  # creep multiobjective weight


@dataclass(frozen=True)
class SolidElastic:
    c11: float = DEFAULT_C11
    c12: float = DEFAULT_C12
    c44: float = DEFAULT_C44

    @classmethod
    def bounds(cls) -> dict[str, Real]:
        return {
            "c11": Real(0.95 * cls.c11, 1.05 * cls.c11, name="c11", prior="uniform"),
            "c12": Real(0.95 * cls.c12, 1.05 * cls.c12, name="c12", prior="uniform"),
            "c44": Real(0.95 * cls.c44, 1.05 * cls.c44, name="c44", prior="uniform"),
        }


@dataclass(frozen=True)
class SolidPlastic:
    tau0xf: float = DEFAULT_tau0xf
    tau0xb: float = DEFAULT_tau0xb
    tau1x: float = DEFAULT_tau1x
    thet0: float = DEFAULT_thet0
    thet1: float = DEFAULT_thet1
    nrsx: float = DEFAULT_nrsx
    gamd0x: float = DEFAULT_gamd0x
    hselfx: float = DEFAULT_hselfx
    hlatex: float = DEFAULT_hlatex

    @classmethod
    def bounds(cls) -> dict[str, Real]:

        return {
            "tau0xf": Real(145.0, 175.0, name="tau0xf", prior="uniform"),
            "tau0xb": Real(0.85, 1.0, name="tau0xb", prior="uniform"),
            "tau1x": Real(1.0, 30.0, name="tau1x", prior="uniform"),
            "thet0": Real(300.0, 4000.0, name="thet0", prior="uniform"),
            "thet1": Real(150.0, 300.0, name="thet1", prior="uniform"),
            "gamd0x": Real(1e-6, 1e-3, name="gamd0x", prior="log-uniform"),
            "nrsx": Real(15, 40, name="nrsx", prior="uniform"),
        }


# complgas is a scalar compliance factor (fraction of solid compliance).
COMPLGAS_BOUNDS = Real(0.0005, 0.05, name="complgas", prior="uniform")

######################### DEFAULT GAS ######################
GAS_DEFAULT_tau0xf = 10.0
GAS_DEFAULT_tau0xb = GAS_DEFAULT_tau0xf
GAS_DEFAULT_nrsx = 5
GAS_DEFAULT_gamd0x = 2e-5
GAS_DEFAULT_young = 2000.0
GAS_DEFAULT_nu = 0.49
# written/read for igas = 2 but never used. Do not calibrate
GAS_DEFAULT_tau1x = 1.0e4
GAS_DEFAULT_thet0 = 1.0
GAS_DEFAULT_thet1 = 1.0
GAS_DEFAULT_hselfx = 1.0
GAS_DEFAULT_hlatex = 1.0


@dataclass(frozen=True)
class GasElastic:
    young: float = GAS_DEFAULT_young
    nu: float = GAS_DEFAULT_nu

    @classmethod
    def bounds(cls) -> dict[str, Real]:
        return {
            "young": Real(10.0, 2.0e4, name="young", prior="uniform"),
            "nu": Real(0.0, 0.49, name="nu", prior="uniform"),
        }


@dataclass(frozen=True)
class GasPlastic:
    tau0xf: float = GAS_DEFAULT_tau0xf
    tau0xb: float = GAS_DEFAULT_tau0xb
    tau1x: float = GAS_DEFAULT_tau1x
    thet0: float = GAS_DEFAULT_thet0
    thet1: float = GAS_DEFAULT_thet1
    nrsx: float = GAS_DEFAULT_nrsx
    gamd0x: float = GAS_DEFAULT_gamd0x
    hselfx: float = GAS_DEFAULT_hselfx
    hlatex: float = GAS_DEFAULT_hlatex

    @classmethod
    def bounds(cls) -> dict[str, Real]:

        return {
            "gas_tau0xf": Real(50.0, 2000.0, name="gas_tau0xf", prior="uniform"),
            "gas_tau0xf": Real(
                0.85 * GAS_DEFAULT_tau0xf,
                1.0 * GAS_DEFAULT_tau0xf,
                name="gas_tau0xb",
                prior="uniform",
            ),
            "gas_gamd0x": Real(1e-10, 1e-1, name="gas_gamd0x", prior="log-uniform"),
            "gas_nrsx": Real(10, 200, name="gas_nrsx", prior="uniform"),
        }


######################### LEGACY DEFAULT ######################
# DEFAULT_C11 = 97153  # C11, axial yield stress / axial stiffness
# DEFAULT_C12 = 91292  # C12, lateral yield stress / lateral coupling
# DEFAULT_C44 = 52658  # C44, yield shear stress / shear stiffness
# DEFAULT_tau0xf = 254
# DEFAULT_rb = 0.83
# DEFAULT_tau0xb = DEFAULT_rb * DEFAULT_tau0xf
# DEFAULT_tau0xb = 211
# DEFAULT_tau1x = 21
# DEFAULT_thet0 = 2167
# DEFAULT_thet1 = 102
# DEFAULT_log10_gamd0x = np.log10(1)  # gamd0x = 1e-4 1/s
# DEFAULT_log10_m = np.log10(0.014)  # m=0.05 -> nrsx=20
# DEFAULT_nrsx = 70
# DEFAULT_gamd0x = 1.0
# DEFAULT_hselfx = 1.0
# DEFAULT_hlatex = 1.0  # Ma et al 2022, peirce et al 1983

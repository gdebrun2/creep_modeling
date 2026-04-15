from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import ConfigDict
from pydantic.dataclasses import dataclass as pydantic_dataclass


from config import (
    MICROSTRUCTURE_PATH,
    DEFAULT_SIM,
)

from io_utils import write_text

SimCase = Literal[
    "creep_strain",
    "creep_stress",
    "tension",
    "test",
    "creep_calibration",
    "srj",
    "creep",
    "srj_old",
    "srj_offset",
]
GasMode = Literal[0, 1, 2]


@pydantic_dataclass(config=ConfigDict(frozen=True, extra="forbid", strict=True))
class FftConfig:

    sim_case: SimCase = DEFAULT_SIM
    igas: GasMode = 0
    load: float | None = None
    complgas: float | None = 0
    microstructure: Path = Path(MICROSTRUCTURE_PATH)
    cuel: Path = Path("cuel.sx")
    cupl2: Path = Path("cupl2.sx")
    cuel_gas: Path = Path("cuel_gas.sx")
    cupl2_gas: Path = Path("cupl2_gas.sx")
    fft: Path = Path("fft.in")
    nphase: int | None = None
    nx: int | None = None
    ny: int | None = None
    nz: int | None = None


def write_fft(cfg: FftConfig) -> None:
    """Write `fft.in` into `fft`."""

    if cfg.sim_case == "creep":
        from sim_cases.creep_stress import build_fft
    elif cfg.sim_case == "creep_calibration":
        from sim_cases.creep_calibration import build_fft
    elif cfg.sim_case == "srj":
        from sim_cases.srj_calibration_new import build_fft
    elif cfg.sim_case == "srj_old":
        from sim_cases.srj_calibration import build_fft
    elif cfg.sim_case == "srj_offset":
        from sim_cases.srj_offset import build_fft
    elif cfg.sim_case == "tension":
        from sim_cases.tension import build_fft
    elif cfg.sim_case == "test":
        from sim_cases.test import build_fft

    else:
        raise ValueError(f"Unsupported sim_case: {cfg.sim_case}")

    text = build_fft(cfg)
    write_text(cfg.fft, text)

    return None

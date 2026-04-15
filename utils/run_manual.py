from __future__ import annotations

import os
from pathlib import Path

from config import (
    MICROSTRUCTURE_PATH,
    SOLVER_EXE,
)
from config import GasElastic, GasPlastic, SolidElastic, SolidPlastic
from io_utils import os_system, write_cuel, write_cupl2, write_text
from write_fft import FftConfig
from data_utils import MicrostructureInfo
from postprocess import run as run_postprocess, PostprocessConfig
import numpy as np
import pandas as pd
from run import _solid_elastic_from_row, _solid_plastic_from_row
import sys

if __name__ == "__main__":

    microstructure = Path("/Users/gtdebru/creep_modeling/microstructures/micro1.dat")
    fftw_save = Path("/Users/gtdebru/creep_modeling/microstructures/micro1_fftw.save")
    micro_info = MicrostructureInfo.load(microstructure)

    nx = micro_info.nx
    ny = micro_info.ny
    nz = micro_info.nz
    nphase = micro_info.nphase

    sim_case = "creep_calibration"
    load = 575
    igas = 0
    modal = 0
    mix = 1

    calib_id = 1
    max_loss = 0.003

    df = pd.read_csv(
        f"/Users/gtdebru/creep_modeling/results/calibrations/{calib_id}/calibration_results.csv"
    )
    rows = df[df["loss"] <= max_loss]
    row_nums = rows.index

    for i, row in rows.iterrows():
        out_dir = Path(
            f"~/creep_modeling/results/calibrations/{calib_id}/results_sweep/row{i}/{load}"
        ).expanduser()

        out_dir.mkdir(parents=True, exist_ok=True)
        complgas = 0.01  # only used for igas=1 (dummy gas)
        fft = out_dir / "fft.in"
        cuel = out_dir / "cuel.sx"
        cupl2 = out_dir / "cupl2.sx"
        evpfft_outputs = out_dir / "evpfft_outputs"
        log = out_dir / "output.txt"
        nthreads = os.cpu_count()

        # Only used for igas=2 (damper) or if you explicitly want to write gas files.
        # cuel_gas = out_dir / "cuel_gas.sx"
        # cupl2_gas = out_dir / "cupl2_gas.sx"

        cfg = FftConfig(
            sim_case=sim_case,
            igas=igas,
            nphase=nphase,
            nx=nx,
            ny=ny,
            nz=nz,
            load=load,
            complgas=complgas,
            fft=fft,
            cuel=cuel,
            cupl2=cupl2,
            # cuel_gas=cuel_gas,
            # cupl2_gas=cupl2_gas,
            microstructure=microstructure,
        )

        # solid_elastic = SolidElastic()
        # solid_plastic = SolidPlastic()

        solid_elastic = _solid_elastic_from_row(row)
        solid_plastic = _solid_plastic_from_row(row)
        write_cuel(cuel, elastic=solid_elastic, iso=0)
        write_cupl2(cupl2, plastic=solid_plastic)

        # gas_elastic = GasElastic()
        # gas_plastic = GasPlastic()
        # write_cuel(cuel_gas, elastic=gas_elastic, iso=1)
        # write_cupl2(cupl2_gas, plastic=gas_plastic)

        if cfg.sim_case == "creep":
            from sim_cases.creep import build_fft
        elif cfg.sim_case == "creep_calibration":
            from sim_cases.creep_calibration import build_fft
        elif cfg.sim_case == "srj":
            from sim_cases.srj import build_fft
        elif cfg.sim_case == "srj_calibration":
            from sim_cases.srj_calibration import build_fft
        elif cfg.sim_case == "tension":
            from sim_cases.tension import build_fft
        elif cfg.sim_case == "test":
            from sim_cases.test import build_fft
        text = build_fft(cfg)
        write_text(fft, text)

        cmd = (
            f"{SOLVER_EXE} --nthreads {nthreads} "
            f"--in {fft} "
            f"--out {evpfft_outputs} "  # not the same as the sim output dir
            f"--fftw_save {fftw_save} "
            f"--mix {mix} "
            f"--modal {modal} "
            f"| tee {log}; "
        )

        os_system(cmd)

        cfg = PostprocessConfig(
            out_dir=Path(out_dir),
            sim_case=str(sim_case),
            igas=int(igas),
            microstructure=Path(microstructure),
            load=load,
            skip_vtk=True,
        )

        written = run_postprocess(cfg)
        for pth in written:
            print(str(pth))

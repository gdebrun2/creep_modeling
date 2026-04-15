from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from data_utils import SimResults, Diagnostics
from config import MICROSTRUCTURE_PATH


@dataclass(frozen=True)
class PostprocessConfig:
    out_dir: Path
    sim_case: str
    microstructure: str | Path
    igas: int = 0
    load: int | None = 0
    skip_vtk: bool = False


def _exp_kind_from_sim_case(sim_case: str) -> str | None:

    if sim_case.split("_")[0] == "creep":
        return "creep"
    if sim_case.split("_")[0] == "srj":
        return "srj"
    if sim_case == "tension":
        return "tension"
    return None


def run(cfg: PostprocessConfig) -> list[Path]:
    """Generate the standardized plot suite for an existing run directory."""

    from plotting.mechanical import save_mechanical_response_plots
    from plotting.diagnostics import save_diagnostics_plots

    out_dir = Path(cfg.out_dir)
    evpfft_outputs = out_dir / "evpfft_outputs"
    plots_root = out_dir / "plots"
    plots_root.mkdir(parents=True, exist_ok=True)

    exp_kind = _exp_kind_from_sim_case(cfg.sim_case)

    sim = None
    diag = None
    try:
        sim = SimResults.load(
            evpfft_outputs, micro_path=cfg.microstructure, skip_vtk=cfg.skip_vtk
        )
    except Exception as e:
        (plots_root / "mech").mkdir(parents=True, exist_ok=True)
        (plots_root / "mech" / "mechanical_error.txt").write_text(str(e))

    try:
        diag = Diagnostics.load(evpfft_outputs)
    except Exception as e:

        (plots_root / "diagnostics").mkdir(parents=True, exist_ok=True)
        (plots_root / "diagnostics" / "diagnostics_error.txt").write_text(str(e))

    written: list[Path] = []
    if diag is not None:
        written += save_diagnostics_plots(
            diag=diag,
            plot_dir=plots_root / "diagnostics",
            load=cfg.load,
            igas=cfg.igas,
            exp_kind=exp_kind,
        )

    if sim is not None:
        written += save_mechanical_response_plots(
            sim=sim,
            plot_dir=plots_root / "mech",
            exp_kind=exp_kind,
            load=cfg.load,
            igas=cfg.igas,
            skip_vtk=cfg.skip_vtk,
        )

    return written


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in {"true", "1", "yes", "y"}:
            return True
        if v.lower() in {"false", "0", "no", "n"}:
            return False
        raise argparse.ArgumentTypeError("Expected true/false")

    p = argparse.ArgumentParser(
        prog="python -m utils.postprocess",
        description="Post-process an existing LS-EVPFFT run directory (no solver run).",
    )
    p.add_argument(
        "--out_dir",
        type=Path,
        required=True,
        help="run directory containing output.txt and evpfft_outputs/",
    )
    p.add_argument(
        "--sim_case",
        type=str,
        required=True,
        help="sim_case name used for selecting experimental overlays and time units",
    )
    p.add_argument(
        "--igas",
        type=int,
        default=0,
        help="igas used for the run (roughness/slip only when 1 or 2)",
    )

    p.add_argument(
        "--micro",
        type=Path,
        default=MICROSTRUCTURE_PATH,
        help="path to microstructure",
    )

    p.add_argument(
        "--load",
        type=int,
        default=0,
        help="load for creep cases",
    )

    p.add_argument(
        "--skip_vtk",
        type=str2bool,
        default=False,
        help="Skip adding vtk height data",
    )

    args = p.parse_args(list(argv) if argv is not None else None)

    cfg = PostprocessConfig(
        out_dir=Path(args.out_dir),
        sim_case=str(args.sim_case),
        igas=int(args.igas),
        microstructure=Path(args.micro),
        load=args.load,
        skip_vtk=args.skip_vtk,
    )

    written = run(cfg)
    for pth in written:
        print(str(pth))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

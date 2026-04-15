from __future__ import annotations

from pathlib import Path

import numpy as np

from data_utils import (
    SimResults,
    CreepData,
    SrjData,
    TensileData,
)
from plotting.common import new_fig_ax, save_fig, style_axes, material_params_suptitle
from config import SRJ_OFFSET, SRJ_OFFSET_STRESS

EXP_COLOR = "tab:blue"  # blue
# EXP_PLASTIC_COLOR = "darkred"
SIM_COLOR = "tab:orange"  # orange
# SIM_PLASTIC_COLOR = "tab:green"  # green


def plot_mean_with_std(
    ax,
    x: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    label: str,
    color: str = EXP_COLOR,
) -> None:
    ax.plot(x, mean, linewidth=2.5, color=color, label=label)
    ax.fill_between(x, mean - std, mean + std, color=color, alpha=0.25)


def write_strain_vs_time(
    sim: SimResults,
    exp: TensileData | SrjData | CreepData | dict[str, np.ndarray],
    plot_dir: str | Path,
    exp_kind: str,
    igas: int,
    load: int,
):

    fig, ax = new_fig_ax()
    title = "Strain vs Time"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"
    if exp_kind == "creep":
        sim_strain = sim.epav33
        sim_time = sim.sim_time.copy() / 3600
        exp_strain = exp["mean_plastic_strain"]
        std_strain = exp["std_plastic_strain"]
        exp_time = exp["time"][1:].copy() / 3600
        ylabel = "Plastic Strain (%)"
        xlabel = "Time (hrs)"

    else:
        sim_strain = sim.eav33
        sim_time = sim.sim_time.copy()
        exp_strain = exp.mean_strain
        std_strain = exp.std_strain
        exp_time = exp.time.copy()
        
        if exp_kind == "srj":
            exp_time -= SRJ_OFFSET
            sim_time -= SRJ_OFFSET

        ylabel = "Strain (%)"
        xlabel = "Time (s)"

    ax.plot(
        sim_time,
        sim_strain * 100,
        label="Sim.",
        color=SIM_COLOR,
        zorder=10,
        linestyle="--",
    )
    plot_mean_with_std(
        ax, exp_time, exp_strain * 100, std_strain * 100, label="Exp.", color=EXP_COLOR
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, Path(plot_dir) / "strain_vs_time.png")


def write_stress_vs_time(
    sim: SimResults,
    exp: TensileData | SrjData | CreepData | dict[str, np.ndarray],
    plot_dir: str | Path,
    exp_kind: str,
    igas: int,
    load: int,
):

    fig, ax = new_fig_ax()
    sim_stress = sim.sav33
    ylabel = "Stress (MPa)"

    title = "Stress vs Time"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    if exp_kind == "creep":
        sim_time = sim.sim_time / 3600
        exp_stress = exp["stress"]
        exp_time = exp["time"] / 3600
        xlabel = "Time (hrs)"

        ax.plot(exp_time, exp_stress, label="Exp.", color=EXP_COLOR)

    else:

        sim_time = sim.sim_time.copy()
        exp_stress = exp.mean_stress
        std_stress = exp.std_stress
        exp_time = exp.time.copy()

        if exp_kind == "srj":
            exp_time -= SRJ_OFFSET
            sim_time -= SRJ_OFFSET

        plot_mean_with_std(
            ax, exp_time, exp_stress, std_stress, label="Exp.", color=EXP_COLOR
        )

        xlabel = "Time (s)"

    ax.plot(sim_time, sim_stress, label="Sim.", color=SIM_COLOR)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, Path(plot_dir) / "stress_vs_time.png")


def write_stress_vs_strain(
    sim: SimResults,
    exp: TensileData | SrjData | CreepData | dict[str, np.ndarray],
    plot_dir: str | Path,
    exp_kind: str,
    igas: int,
    load: int,
):

    fig, ax = new_fig_ax()
    title = "Stress vs Strain"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    sim_strain = sim.epav33.copy()
    sim_stress = sim.sav33.copy()

    if exp_kind == "creep":

        exp_strain = exp["mean_plastic_strain"].copy()
        exp_stress = exp["stress"].copy()
        xlabel = "Plastic Strain (%)"
        ylabel = "Stress (MPa)"

    else:
        exp_strain = exp.mean_strain.copy()
        exp_stress = exp.mean_stress.copy()
        xlabel = "Strain (%)"
        ylabel = "Stress (MPa)"

    if "srj" in exp_kind:
        sim_mask = np.where(sim_stress > SRJ_OFFSET_STRESS)
        exp_mask = np.where(exp_stress > SRJ_OFFSET_STRESS)

        exp_strain = exp_strain[exp_mask]
        exp_stress = exp_stress[exp_mask]
        sim_strain = sim_strain[sim_mask]
        sim_stress = sim_stress[sim_mask]

    exp_argsort = np.argsort(exp_strain)
    sim_argsort = np.argsort(sim_strain)

    ax.plot(
        sim_strain[sim_argsort] * 100,
        sim_stress[sim_argsort],
        label="Sim.",
        color=SIM_COLOR,
        zorder=10,
        linestyle="--",
    )
    ax.plot(
        exp_strain[exp_argsort] * 100,
        exp_stress[exp_argsort],
        label="Exp.",
        color=EXP_COLOR,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, Path(plot_dir) / "stress_vs_strain.png")


def write_sa_vs_time_norm(
    sim: SimResults,
    exp: TensileData | SrjData | CreepData | dict[str, np.ndarray],
    plot_dir: str | Path,
    exp_kind: str,
    igas: int,
    load: int,
):

    fig, ax = new_fig_ax()
    ylabel = "Sa (max-norm)"
    xlabel = "Time (s)"

    title = "Sa vs Time"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    sim_time = sim.vtk_time.copy()
    if exp_kind == "creep":
        sim_time /= 3600

    sa_xmin, sa_xmax, sa_ymin, sa_ymax = sim.sa
    ax.plot(sim_time, sa_xmin / sa_xmin.max(), label="-X", color="tab:red")
    ax.plot(sim_time, sa_xmax / sa_xmax.max(), label="+X", color="tab:orange")
    ax.plot(sim_time, sa_ymin / sa_ymin.max(), label="-Y", color="tab:green")
    ax.plot(sim_time, sa_ymax / sa_ymax.max(), label="+Y", color="tab:blue")

    if exp_kind == "creep":

        exp_time = exp["roughness_time"] / 3600
        exp_mean_sa_10x = exp["mean_sa_10x"]
        exp_std_sa_10x = exp["std_sa_10x"]
        exp_mean_sa_50x = exp["mean_sa_50x"]
        exp_std_sa_50x = exp["std_sa_50x"]
        xlabel = "Time (hrs)"

        ax.plot(
            exp_time,
            exp_mean_sa_10x / exp_mean_sa_10x.max(),
            label="10x",
            color="black",
        )
        # ax.plot(
        #     exp_time,
        #     exp_mean_sa_50x / exp_mean_sa_50x.max(),
        #     label="50x",
        #     color="black",
        # )

        # plot_mean_with_std(
        #     ax, exp_time, exp_mean_sa_10x, exp_std_sa_10x, label="10x", color="black"
        # )

        # plot_mean_with_std(
        #     ax, exp_time, exp_mean_sa_50x, exp_std_sa_50x, label="50x", color="yellow"
        # )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, Path(plot_dir) / "sa_vs_time_norm.png")


def write_sz_vs_time_norm(
    sim: SimResults,
    exp: TensileData | SrjData | CreepData | dict[str, np.ndarray],
    plot_dir: str | Path,
    exp_kind: str,
    igas: int,
    load: int,
):

    fig, ax = new_fig_ax()
    ylabel = "Sz (max-norm)"
    xlabel = "Time (s)"
    title = "Sz vs Time"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    time = sim.vtk_time.copy()
    if exp_kind == "creep":
        time /= 3600
    sz_xmin, sz_xmax, sz_ymin, sz_ymax = sim.sa

    ax.plot(time, sz_xmin / sz_xmin.max(), label="-X", color="tab:red")
    ax.plot(time, sz_xmax / sz_xmax.max(), label="+X", color="tab:orange")
    ax.plot(time, sz_ymin / sz_ymin.max(), label="-Y", color="tab:green")
    ax.plot(time, sz_ymax / sz_ymax.max(), label="+Y", color="tab:blue")

    if exp_kind == "creep":
        xlabel = "Time (hrs)"
        exp_time = exp["roughness_time"] / 3600
        exp_mean_sz_10x = exp["mean_sz_10x"]
        exp_std_sz_10x = exp["std_sz_10x"]
        exp_mean_sz_50x = exp["mean_sz_50x"]
        exp_std_sz_50x = exp["std_sz_50x"]

        ax.plot(
            exp_time,
            exp_mean_sz_10x / exp_mean_sz_10x.max(),
            label="10x",
            color="black",
        )
        ax.plot(
            exp_time,
            exp_mean_sz_50x / exp_mean_sz_50x.max(),
            label="10x",
            color="black",
        )

        # plot_mean_with_std(
        #     ax, exp_time, exp_mean_sz_10x, exp_std_sz_10x, label="10x", color="black"
        # )

        # plot_mean_with_std(
        #     ax, exp_time, exp_mean_sz_50x, exp_std_sz_50x, label="50x", color="yellow"
        # )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, Path(plot_dir) / "sz_vs_time_norm.png")


def write_sa_vs_time(
    sim: SimResults,
    exp: TensileData | SrjData | CreepData | dict[str, np.ndarray],
    plot_dir: str | Path,
    exp_kind: str,
    igas: int,
    load: int,
):

    fig, ax = new_fig_ax()
    ylabel = "Sa"
    xlabel = "Time (s)"

    title = "Sa vs Time"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    sim_time = sim.vtk_time.copy()
    if exp_kind == "creep":
        sim_time /= 3600

    sa_xmin, sa_xmax, sa_ymin, sa_ymax = sim.sa

    if exp_kind == "creep":

        exp_time = exp["roughness_time"] / 3600
        exp_mean_sa_10x = exp["mean_sa_10x"] / 10
        exp_std_sa_10x = exp["std_sa_10x"] / 10
        exp_mean_sa_50x = exp["mean_sa_50x"]
        exp_std_sa_50x = exp["std_sa_50x"]
        xlabel = "Time (hrs)"

        plot_mean_with_std(
            ax, exp_time, exp_mean_sa_10x, exp_std_sa_10x, label="10x", color="black"
        )

        # plot_mean_with_std(
        #     ax, exp_time, exp_mean_sa_50x, exp_std_sa_50x, label="50x", color="yellow"
        # )
    else:
        exp_mean_sa_10x = [0]

    ax.plot(
        sim_time,
        exp_mean_sa_10x[0] + sa_xmin * 4,
        label="-X",
        color="tab:red",
        zorder=10,
    )
    ax.plot(
        sim_time,
        exp_mean_sa_10x[0] + sa_xmax * 4,
        label="+X",
        color="tab:orange",
        zorder=11,
    )
    ax.plot(
        sim_time,
        exp_mean_sa_10x[0] + sa_ymin * 4,
        label="-Y",
        color="tab:green",
        zorder=12,
    )
    ax.plot(
        sim_time,
        exp_mean_sa_10x[0] + sa_ymax * 4,
        label="+Y",
        color="tab:blue",
        zorder=13,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, Path(plot_dir) / "sa_vs_time.png")


def write_sz_vs_time(
    sim: SimResults,
    exp: TensileData | SrjData | CreepData | dict[str, np.ndarray],
    plot_dir: str | Path,
    exp_kind: str,
    igas: int,
    load: int,
):

    fig, ax = new_fig_ax()
    ylabel = "Sz"
    xlabel = "Time (s)"
    title = "Sz vs Time"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    time = sim.vtk_time.copy()
    if exp_kind == "creep":
        time /= 3600
    sz_xmin, sz_xmax, sz_ymin, sz_ymax = sim.sa

    ax.plot(time, sz_xmin, label="-X", color="tab:red")
    ax.plot(time, sz_xmax, label="+X", color="tab:orange")
    ax.plot(time, sz_ymin, label="-Y", color="tab:green")
    ax.plot(time, sz_ymax, label="+Y", color="tab:blue")

    if exp_kind == "creep":
        xlabel = "Time (hrs)"
        exp_time = exp["roughness_time"] / 3600
        exp_mean_sz_10x = exp["mean_sz_10x"]
        exp_std_sz_10x = exp["std_sz_10x"]
        exp_mean_sz_50x = exp["mean_sz_50x"]
        exp_std_sz_50x = exp["std_sz_50x"]

        plot_mean_with_std(
            ax, exp_time, exp_mean_sz_10x, exp_std_sz_10x, label="10x", color="black"
        )

        plot_mean_with_std(
            ax, exp_time, exp_mean_sz_50x, exp_std_sz_50x, label="50x", color="yellow"
        )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, Path(plot_dir) / "sz_vs_time.png")


def write_max_slip_vs_time(
    sim: SimResults, plot_dir: str | Path, igas: int, load: int, exp_kind: str
):

    fig, ax = new_fig_ax()
    ylabel = "Slip"
    xlabel = "Time (hrs)"
    title = "Max Slip vs Time"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    sim_time = sim.vtk_time.copy()
    if exp_kind == "creep":
        sim_time /= 3600
    max_slip_xmin, max_slip_xmax, max_slip_ymin, max_slip_ymax = sim.max_slip

    ax.plot(sim_time, max_slip_xmin, label="-X", color="tab:red")
    ax.plot(sim_time, max_slip_xmax, label="+X", color="tab:orange")
    ax.plot(sim_time, max_slip_ymin, label="-Y", color="tab:green")
    ax.plot(sim_time, max_slip_ymax, label="+Y", color="tab:blue")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)

    return save_fig(fig, Path(plot_dir) / "max_slip_vs_time.png")


def write_mean_slip_vs_time(
    sim: SimResults,
    plot_dir: str | Path,
    igas: int,
    load: int,
    exp_kind: str,
):

    fig, ax = new_fig_ax()
    ylabel = "Slip"
    xlabel = "Time (hrs)"
    title = "Mean Slip vs Time"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    sim_time = sim.vtk_time.copy()
    if exp_kind == "creep":
        sim_time /= 3600
    mean_slip_xmin, mean_slip_xmax, mean_slip_ymin, mean_slip_ymax = sim.mean_slip

    ax.plot(sim_time, mean_slip_xmin, label="-X", color="tab:red")
    ax.plot(sim_time, mean_slip_xmax, label="+X", color="tab:orange")
    ax.plot(sim_time, mean_slip_ymin, label="-Y", color="tab:green")
    ax.plot(sim_time, mean_slip_ymax, label="+Y", color="tab:blue")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)

    return save_fig(fig, Path(plot_dir) / "mean_slip_vs_time.png")


def write_g_vs_time(
    sim: SimResults,
    plot_dir: str | Path,
    exp_kind: str,
    igas: int,
    load: int,
):

    fig, ax = new_fig_ax()
    ylabel = r"$g$"

    title = "g/CRSS vs Time"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    if exp_kind == "creep":

        sim_time = sim.sim_time / 3600
        xlabel = "Time (hrs)"

    else:

        sim_time = sim.sim_time.copy()
        if exp_kind == "srj":
            sim_time -= SRJ_OFFSET

        xlabel = "Time (s)"

    ax.plot(sim_time, sim.geff, label="Sim.", color=SIM_COLOR)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, Path(plot_dir) / "geff_vs_time.png")


def write_slip_vs_time(
    sim: SimResults,
    plot_dir: str | Path,
    exp_kind: str,
    igas: int,
    load: int,
):

    fig, ax = new_fig_ax()
    ylabel = r"$\Gamma$"

    title = "Accum. Slip vs Time"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    if exp_kind == "creep":

        sim_time = sim.sim_time / 3600
        xlabel = "Time (hrs)"

    else:

        sim_time = sim.sim_time.copy()

        if exp_kind == "srj":
            sim_time -= SRJ_OFFSET

        xlabel = "Time (s)"

    ax.plot(sim_time, sim.gamma, label="Sim.", color=SIM_COLOR)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, Path(plot_dir) / "slip_vs_time.png")


# def write_max_slip_vs_strain(
#     sim: SimResults, plot_dir: str | Path, igas: int, load: int
# ):

#     fig, ax = new_fig_ax()
#     ylabel = "Slip"
#     xlabel = "Plastic Strain"
#     title = "Max Slilp vs Strain"
#     material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
#     if load:
#         title += f" {load} MPa"
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)

#     slip = sim.epav33
#     max_slip_xmin, max_slip_xmax, max_slip_ymin, max_slip_ymax = sim.max_slip

#     ax.plot(slip, max_slip_xmin, label="-X", color="tab:red")
#     ax.plot(slip, max_slip_xmax, label="+X", color="tab:orange")
#     ax.plot(slip, max_slip_ymin, label="-Y", color="tab:green")
#     ax.plot(slip, max_slip_ymax, label="+Y", color="tab:blue")
#     ax.set_title(title)
#     ax.legend()
#     style_axes(ax)

#     return save_fig(fig, Path(plot_dir) / "max_slip_vs_strain.png")


# def write_mean_slip_vs_strain(
#     sim: SimResults, plot_dir: str | Path, igas: int, load: int
# ):

#     fig, ax = new_fig_ax()
#     ylabel = "Slip"
#     xlabel = "Plastic Strain"
#     title = "Mean Slip vs Strain"
#     material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
#     if load:
#         title += f" {load} MPa"
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)

#     strain = sim.epav33
#     mean_slip_xmin, mean_slip_xmax, mean_slip_ymin, mean_slip_ymax = sim.mean_slip

#     ax.plot(strain, mean_slip_xmin, label="-X", color="tab:red")
#     ax.plot(strain, mean_slip_xmax, label="+X", color="tab:orange")
#     ax.plot(strain, mean_slip_ymin, label="-Y", color="tab:green")
#     ax.plot(strain, mean_slip_ymax, label="+Y", color="tab:blue")
#     ax.set_title(title)
#     ax.legend()
#     style_axes(ax)

#     return save_fig(fig, Path(plot_dir) / "mean_slip_vs_strain.png")


def save_mechanical_response_plots(
    *,
    sim: SimResults,
    plot_dir: str | Path,
    exp_kind: str | None = None,
    load: int | None = 0,
    igas: int | None = 0,
    skip_vtk: bool = False,
) -> list[Path]:

    plot_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    if exp_kind == "creep":
        exp = CreepData.load()
        exp = exp.get_load_vals(load)
    if exp_kind == "srj":
        exp = SrjData.load()

    if exp_kind == "tension":
        exp = TensileData.load()

    written.append(
        write_stress_vs_time(
            sim=sim,
            exp=exp,
            plot_dir=plot_dir,
            exp_kind=exp_kind,
            igas=igas,
            load=load,
        )
    )
    written.append(
        write_strain_vs_time(
            sim=sim,
            exp=exp,
            plot_dir=plot_dir,
            exp_kind=exp_kind,
            igas=igas,
            load=load,
        )
    )

    written.append(
        write_stress_vs_strain(
            sim=sim,
            exp=exp,
            plot_dir=plot_dir,
            exp_kind=exp_kind,
            igas=igas,
            load=load,
        )
    )

    written.append(
        write_slip_vs_time(
            sim, plot_dir, igas=igas, load=load, exp_kind=exp_kind,
        )
    )
    written.append(
        write_g_vs_time(
            sim, plot_dir, igas=igas, load=load, exp_kind=exp_kind, 
        )
    )
    if skip_vtk == False:
        written.append(
            write_sa_vs_time(
                sim, exp, exp_kind=exp_kind, plot_dir=plot_dir, igas=igas, load=load
            )
        )
        written.append(
            write_sz_vs_time(
                sim, exp, exp_kind=exp_kind, plot_dir=plot_dir, igas=igas, load=load
            )
        )

        written.append(
            write_sa_vs_time_norm(
                sim, exp, exp_kind=exp_kind, plot_dir=plot_dir, igas=igas, load=load
            )
        )
        written.append(
            write_sz_vs_time_norm(
                sim, exp, exp_kind=exp_kind, plot_dir=plot_dir, igas=igas, load=load
            )
        )
        written.append(
            write_max_slip_vs_time(
                sim, plot_dir, igas=igas, load=load, exp_kind=exp_kind
            )
        )
        written.append(
            write_mean_slip_vs_time(
                sim, plot_dir, igas=igas, load=load, exp_kind=exp_kind
            )
        )

    # written.append(write_max_slip_vs_strain(sim, plot_dir, igas=igas, load=load))
    # written.append(write_mean_slip_vs_strain(sim, plot_dir, igas=igas, load=load))

    return written

from __future__ import annotations

from pathlib import Path

import numpy as np
import sys
import pandas as pd
from matplotlib.axes import Axes
from data_utils import Diagnostics

sys.path.append("../")
from plotting.common import (
    VerticalLineStyle,
    draw_vertical_lines,
    new_fig_ax,
    save_fig,
    style_axes,
    material_params_suptitle,
)

_TDOT_LINE_STYLE = VerticalLineStyle(
    color="k", linestyle="--", alpha=0.5, linewidth=1.5, zorder=0
)

GAS_COLOR = "darkred"
SOL_COLOR = "tab:blue"
ALL_COLOR = "k"
TARGET_COLOR = "tab:orange"


def _clean(y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float).copy()
    y[~np.isfinite(y)] = np.nan
    return y


def _axis_context(
    diag: Diagnostics,
    axis_kind: str,
    plot_dir: Path | str,
    exp_kind: str,
    stem: str,
) -> tuple[np.ndarray, str, Path, str]:
    axis_key = axis_kind.lower()
    plot_dir = Path(plot_dir)

    if axis_key == "time":
        if exp_kind == "creep":
            return (
                _clean(diag.sim_time / 3600),
                "Time (hrs)",
                plot_dir / "time" / f"{stem}_vs_time.png",
                "Time",
            )
        else:
            return (
                _clean(diag.sim_time),
                "Time (s)",
                plot_dir / "time" / f"{stem}_vs_time.png",
                "Time",
            )
    if axis_key == "macrostep":
        return (
            _clean(diag.macrostep),
            "Macrostep",
            plot_dir / "macrostep" / f"{stem}_vs_macrostep.png",
            "Macrostep",
        )

    raise ValueError(f"Unsupported diagnostics axis: {axis_kind}")


def _tdot_target_change_markers(
    diag: Diagnostics,
    axis_kind: str,
    exp_kind: str,
    rel_tol: float = 1e-10,
    min_target_s: float = 10.0,
) -> list[tuple[float, float]]:
    """Return markers for TDOT target changes.

    Returns a list of (x_position, new_target_value).

    We deliberately compute these markers for both time and macrostep axes.
    """

    tdot_target = diag.tdotref
    x, _, _, _ = _axis_context(diag, axis_kind, ".", exp_kind, "_unused")
    xs: list[tuple[float, float]] = []
    for i in range(1, len(tdot_target)):
        a = float(tdot_target[i - 1])
        b = float(tdot_target[i])
        if not (np.isfinite(a) and np.isfinite(b)):
            continue
        # Only mark *new* targets above a threshold to avoid clutter.
        # TDOT is interpreted as a time increment with units of seconds.
        if b <= float(min_target_s):
            continue
        if a == 0.0 and b == 0.0:
            continue
        if a == 0.0 or b == 0.0:
            xs.append((float(x[i]), b))
            continue
        if abs(b - a) / abs(a) > rel_tol:
            xs.append((float(x[i]), b))

    return xs


def filter_markers(
    x: np.ndarray,
    *,
    min_abs: float | None = None,
    min_frac: float | None = None,
    xlim: tuple[float, float] | None = None,
) -> np.ndarray:
    if x.size == 0:
        return x

    x = np.unique(x, axis=0)
    x = x[np.argsort(x[:, 0])]
    kept = [x[0]]

    if min_frac is not None:
        if xlim is None:
            xmin, xmax = float(x.min(axis=0)), float(x.max(axis=0))
        else:
            xmin, xmax = map(float, xlim)
        min_frac_abs = min_frac * max(xmax - xmin, 1e-12)
    else:
        min_frac_abs = None
    for xi in x[1:]:
        thresholds = []
        if min_abs is not None:
            thresholds.append(float(min_abs))
        if min_frac_abs is not None:
            thresholds.append(float(min_frac_abs))

        min_sep = max(thresholds) if thresholds else 0.0

        if xi[0] - kept[-1][0] >= min_sep:
            kept.append(xi)

    return np.asarray(kept, dtype=float)


def _proc_change_markers(
    diag: Diagnostics,
    axis_kind: str,
    exp_kind: str,
    rel_tol: float = 1e-10,
) -> np.ndarray:

    proc = diag.ipr
    x, _, _, _ = _axis_context(diag, axis_kind, ".", exp_kind, "_unused")
    xs = np.zeros((proc.size - 1, 2))

    for i in range(1, len(proc)):
        a = float(proc[i - 1])
        b = float(proc[i])

        if not (np.isfinite(a) and np.isfinite(b)):
            continue
        if a == 0.0 and b == 0.0:
            continue
        if a == 0.0 or b == 0.0:
            xs[i - 1] = [float(x[i]), b]

            continue
        if abs(b - a) / abs(a) > rel_tol:
            xs[i - 1] = [float(x[i]), b]

    xs = filter_markers(xs, min_frac=0.05, xlim=((np.nanmin(x), np.nanmax(x))))
    xs = xs[np.any(xs != 0.0, axis=1)]
    return xs


def _annotate_tdot_target_changes(ax, markers: np.ndarray) -> None:
    xs = [m[0] for m in markers]
    draw_vertical_lines(
        ax,
        xs,
        label="dt",
        style=_TDOT_LINE_STYLE,
    )
    if not markers:
        return

    # Text annotations: label with the new target value.
    ymin, ymax = ax.get_ylim()
    y = ymin + 0.92 * (ymax - ymin)
    for x, new_t in markers:
        if not np.isfinite(new_t):
            continue
        ax.text(
            x,
            y,
            f"{new_t:g}",
            rotation=90,
            va="top",
            ha="right",
            fontsize=8,
            color=_TDOT_LINE_STYLE.color,
            alpha=0.9,
        )


def _annotate_tdot_target_changes(ax, markers: np.ndarray) -> None:
    xs = [m[0] for m in markers]
    draw_vertical_lines(
        ax,
        xs,
        label="proc",
        style=_TDOT_LINE_STYLE,
    )
    if markers.size == 0:
        return

    # Text annotations: label with the new target value.
    ymin, ymax = ax.get_ylim()
    y = ymin + 0.92 * (ymax - ymin)
    for x, new_t in markers:
        if not np.isfinite(new_t):
            continue
        ax.text(
            x,
            y,
            f"{new_t:g}",
            rotation=90,
            va="top",
            ha="right",
            fontsize=8,
            color=_TDOT_LINE_STYLE.color,
            alpha=0.9,
        )


def write_tdot_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    x, xlabel, save_path, axis_title = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "tdot"
    )
    title = rf"$\Delta t$ vs {axis_title}"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(x, _clean(diag.tdot), linewidth=2.5, color=SOL_COLOR, label=r"$\Delta t$")
    ax.plot(
        x,
        _clean(diag.tdotref),
        linewidth=2.5,
        color=TARGET_COLOR,
        linestyle="--",
        label=r"$\Delta t_{\mathrm{target}}$",
    )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\Delta t$ (s)")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_field_errs_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = "Field Errors"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "field_errs"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(
        x, _clean(diag.erre), linewidth=2.5, color=SOL_COLOR, label=r"$\varepsilon$"
    )
    ax.plot(
        x,
        _clean(diag.errs),
        linewidth=2.5,
        color=GAS_COLOR,
        label=r"$\sigma$",
    )
    ax.plot(
        x,
        _clean(diag.errsbc),
        linewidth=2.5,
        color=ALL_COLOR,
        linestyle="--",
        label=r"$\sigma_{BC}$",
    )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Error")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_stress_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = "Stress Tracking"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "stress"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(x, _clean(diag.s33), linewidth=2.5, color=SOL_COLOR, label=r"$\sigma_{33}$")
    ax.plot(
        x,
        _clean(diag.s33_target),
        linewidth=2.5,
        color=TARGET_COLOR,
        linestyle="--",
        label=r"$\sigma_{33,\mathrm{target}}$",
    )
    ax.plot(
        x,
        _clean(diag.s33_sol),
        linewidth=2.5,
        color=ALL_COLOR,
        label=r"$\sigma_{33,\mathrm{sol}}$",
    )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Stress (MPa)")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_detfmin_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = r"Minimum $\det(\mathbf{F})$"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "detfmin"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(
        x,
        _clean(diag.detF_min_sol),
        linewidth=2.5,
        color=SOL_COLOR,
        label="Solid",
    )
    if igas != 0:
        ax.plot(
            x,
            _clean(diag.detF_min_gas),
            linewidth=2.5,
            color=GAS_COLOR,
            label="Gas",
        )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\min \det(\mathbf{F})$")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_rel_detfmin_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = r"Relative Minimum $\det(\mathbf{F})$"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "rel_detfmin"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(
        x,
        _clean(diag.detF_min_sol_rel),
        linewidth=2.5,
        color=SOL_COLOR,
        label="Solid",
    )
    if igas != 0:
        ax.plot(
            x,
            _clean(diag.detF_min_gas_rel),
            linewidth=2.5,
            color=GAS_COLOR,
            label="Gas",
        )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Relative Change")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_volinc_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = "Max Volume Increment"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "volinc"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(
        x,
        _clean(diag.volinc_max_sol),
        linewidth=2.5,
        color=SOL_COLOR,
        label="Solid",
    )
    if igas != 0:
        ax.plot(
            x,
            _clean(diag.volinc_max_gas),
            linewidth=2.5,
            color=GAS_COLOR,
            label="Gas",
        )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\max |\Delta V / V|$")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_rel_volinc_max(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = "Relative Max Volume Increment"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "rel_volinc"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(
        x,
        _clean(diag.volinc_max_sol_rel),
        linewidth=2.5,
        color=SOL_COLOR,
        label="Solid",
    )
    if igas != 0:
        ax.plot(
            x,
            _clean(diag.volinc_max_gas_rel),
            linewidth=2.5,
            color=GAS_COLOR,
            label="Gas",
        )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Relative Change")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_stretchinc_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = "Max Stretch Increment"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "stretchinc"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(
        x,
        _clean(diag.stretchinc_max_sol),
        linewidth=2.5,
        color=SOL_COLOR,
        label="Solid",
    )
    if igas != 0:
        ax.plot(
            x,
            _clean(diag.stretchinc_max_gas),
            linewidth=2.5,
            color=GAS_COLOR,
            label="Gas",
        )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Stretch Increment")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_rel_stretchinc_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = "Relative Max Stretch Increment"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "rel_stretchinc"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(
        x,
        _clean(diag.stretchinc_max_sol_rel),
        linewidth=2.5,
        color=SOL_COLOR,
        label="Solid",
    )
    if igas != 0:
        ax.plot(
            x,
            _clean(diag.stretchinc_max_gas_rel),
            linewidth=2.5,
            color=GAS_COLOR,
            label="Gas",
        )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Relative Change")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_rotang_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = "Max Rotation Increment"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "rotang"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(
        x,
        _clean(diag.rotang_max_sol),
        linewidth=2.5,
        color=SOL_COLOR,
        label="Solid",
    )
    if igas != 0:
        ax.plot(
            x,
            _clean(diag.rotang_max_gas),
            linewidth=2.5,
            color=GAS_COLOR,
            label="Gas",
        )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\max \Delta \theta$ (deg)")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_dl_step_rms_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = r"RMS $\Delta \mathbf{L}$ per Step"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "dl_step_rms"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(
        x,
        _clean(diag.dL_step_rms_sol),
        linewidth=2.5,
        color=SOL_COLOR,
        label="Solid",
    )
    if igas != 0:
        ax.plot(
            x,
            _clean(diag.dL_step_rms_gas),
            linewidth=2.5,
            color=GAS_COLOR,
            label="Gas",
        )
    ax.plot(
        x,
        _clean(diag.dL_step_rms),
        linewidth=2.5,
        color=ALL_COLOR,
        linestyle="--",
        label="All",
    )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\|\Delta \mathbf{L}\|_{\mathrm{RMS}}$")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_max_dl_step_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = r"Max $\Delta \mathbf{L}$ per Step"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "max_dl_step"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(
        x,
        _clean(diag.dL_step_max_sol),
        linewidth=2.5,
        color=SOL_COLOR,
        label="Solid",
    )
    if igas != 0:
        ax.plot(
            x,
            _clean(diag.dL_step_max_gas),
            linewidth=2.5,
            color=GAS_COLOR,
            label="Gas",
        )
    ax.plot(
        x,
        _clean(diag.dL_step_max),
        linewidth=2.5,
        color=ALL_COLOR,
        linestyle="--",
        label="All",
    )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\max \|\Delta \mathbf{L}\|$")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_dvelgradavg_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = r"Mean $\Delta \mathbf{L}$ Magnitude"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "dvelgradavg"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(
        x,
        _clean(diag.dvelgradavg_mag),
        linewidth=2.5,
        color=SOL_COLOR,
        label=r"$\|\Delta \langle \mathbf{L} \rangle\|$",
    )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\|\Delta \langle \mathbf{L} \rangle\|$")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_conditioning_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    title = "Operator Conditioning"
    x, xlabel, save_path, _ = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "conditioning"
    )
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"

    ax.plot(
        x,
        _clean(diag.cond_c066mod_max_sol),
        linewidth=2.5,
        color=SOL_COLOR,
        label=r"$c^{066}_{\mathrm{mod}}$ solid",
    )
    if igas != 0:
        ax.plot(
            x,
            _clean(diag.cond_c066mod_max_gas),
            linewidth=2.5,
            color=GAS_COLOR,
            label=r"$c^{066}_{\mathrm{mod}}$ gas",
        )
        ax.plot(
            x,
            _clean(diag.cond_cgas_max),
            linewidth=2.5,
            color=ALL_COLOR,
            linestyle="--",
            label=r"$c_{\mathrm{gas}}$",
        )
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\kappa$")
    ax.set_title(title)
    ax.legend()
    style_axes(ax)
    return save_fig(fig, save_path)


def write_iters_vs_macrostep(
    diag: Diagnostics,
    vlines: list[tuple[float, float]],
    plot_dir: Path,
    igas: int,
    load: int,
    axis_kind: str,
    exp_kind: str,
) -> Path:

    fig, ax = new_fig_ax()
    x, xlabel, save_path, axis_title = _axis_context(
        diag, axis_kind, plot_dir, exp_kind, "iters"
    )
    ylabel = "Iterations"
    title = f"Iterations per {axis_title}"
    material_params_suptitle(fig, start_dir=plot_dir, gas=igas == 2)
    if load:
        title += f" {load} MPa"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    iters = diag.iter
    ax.plot(x, iters, linewidth=2.5)
    _annotate_tdot_target_changes(ax, vlines)

    ax.set_title(title)
    style_axes(ax)
    return save_fig(fig, save_path)


def save_diagnostics_plots(
    *,
    diag: Diagnostics,
    plot_dir: Path | str,
    load: int,
    igas: int,
    exp_kind: str,
) -> list[Path]:
    """Save solver diagnostics plots."""

    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    writers = (
        write_iters_vs_macrostep,
        write_tdot_vs_macrostep,
        write_field_errs_vs_macrostep,
        write_stress_vs_macrostep,
        write_detfmin_vs_macrostep,
        write_rel_detfmin_vs_macrostep,
        write_volinc_vs_macrostep,
        write_rel_volinc_max,
        write_stretchinc_vs_macrostep,
        write_rel_stretchinc_vs_macrostep,
        write_rotang_vs_macrostep,
        write_dl_step_rms_vs_macrostep,
        write_max_dl_step_vs_macrostep,
        write_dvelgradavg_vs_macrostep,
        write_conditioning_vs_macrostep,
    )

    for axis_kind in ("macrostep", "time"):
        vlines = _proc_change_markers(diag, axis_kind, exp_kind)

        for writer in writers:
            written.append(
                writer(
                    diag,
                    vlines,
                    plot_dir,
                    igas,
                    load,
                    axis_kind,
                    exp_kind,
                )
            )

    return written

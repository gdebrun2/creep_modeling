from __future__ import annotations
from typing import TYPE_CHECKING, Literal
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
from config import (
    CREEP_CSV_PATH,
    SRJ_CSV_PATH,
    TENSILE_CSV_PATH,
    MICROSTRUCTURE_PATH,
    PHASE_GAS,
    PHASE_SOLID,
)

if TYPE_CHECKING:
    from data_utils import MicrostructureInfo


def to_true_strain(engineering_strain: np.ndarray) -> np.ndarray:
    """Convert engineering strain to true strain: ln(1 + e)."""

    return np.log1p(engineering_strain)


def to_true_stress(
    engineering_stress: np.ndarray, engineering_strain: np.ndarray
) -> np.ndarray:
    """Convert engineering stress to true stress: s_true = s_eng * (1 + e)."""

    return engineering_stress * (1.0 + engineering_strain)


def nearest_indices(query_x: np.ndarray, ref_x: np.ndarray) -> np.ndarray:
    """Return indices into ref_x closest to each query_x entry (L1 distance)."""

    if ref_x.ndim != 1 or query_x.ndim != 1:
        raise ValueError("nearest_indices expects 1D arrays")
    if len(ref_x) == 0:
        raise ValueError("ref_x must be non-empty")

    # O(N*M) but N here is typically small (macrostep count) and avoids interp
    # artifacts when time grids are irregular.
    idx = np.empty_like(query_x, dtype=int)
    for i, x in enumerate(query_x):
        idx[i] = int(np.argmin(np.abs(ref_x - x)))
    return idx


def mean_abs_error_aligned(
    *,
    sim_x: np.ndarray,
    sim_y: np.ndarray,
    exp_x: np.ndarray,
    exp_y: np.ndarray,
) -> float:
    """Mean absolute error between exp_y and sim_y sampled at exp_x locations."""

    if len(exp_x) == 0 or len(sim_x) == 0:
        return np.inf
    idx = nearest_indices(exp_x, sim_x)
    y_sim = sim_y[idx]
    r = np.abs(exp_y - y_sim) / exp_y.max()
    MAE = r.mean()
    return MAE


def mean_sq_error_aligned(
    *,
    sim_x: np.ndarray,
    sim_y: np.ndarray,
    exp_x: np.ndarray,
    exp_y: np.ndarray,
) -> float:
    """Mean squared error between exp_y and sim_y sampled at exp_x locations."""

    if len(exp_x) == 0 or len(sim_x) == 0:
        return np.inf
    idx = nearest_indices(exp_x, sim_x)
    y_sim = sim_y[idx]
    r = (exp_y - y_sim) / exp_y.max()
    SE = r**2
    MSE = SE.mean()
    return MSE


def mean_rel_err_aligned(
    *,
    sim_x: np.ndarray,
    sim_y: np.ndarray,
    exp_x: np.ndarray,
    exp_y: np.ndarray,
    eps: float = 1e-12,
) -> float:
    if len(exp_x) == 0 or len(sim_x) == 0:
        return np.inf
    idx = nearest_indices(exp_x, sim_x)
    y_sim = sim_y[idx]

    eps = np.ones_like(y_sim) * eps
    a = np.abs(y_sim - exp_y)
    b = np.maximum(np.abs(y_sim), np.abs(exp_y), eps)

    return np.mean(a / b)


def norm_max(series: np.ndarray) -> np.ndarray:
    """Normalize by max(series) with zero-protection."""

    m = float(np.max(series)) if len(series) else 0.0
    if m == 0.0:
        return np.zeros_like(series, dtype=float)
    return series / m


def node_id_to_cell_id(info: MicrostructureInfo):

    nx, ny, nz = info.nx, info.ny, info.nz
    Nnodes = (nz + 1) * (ny + 1) * (nx + 1)

    node_id = np.arange(Nnodes, dtype=np.int64).reshape((nz + 1, ny + 1, nx + 1))

    z = np.arange(nz)[:, None, None]
    y = np.arange(ny)[None, :, None]
    x = np.arange(nx)[None, None, :]

    n000 = node_id[z, y, x]
    n100 = node_id[z, y, x + 1]
    n110 = node_id[z, y + 1, x + 1]
    n010 = node_id[z, y + 1, x]
    n001 = node_id[z + 1, y, x]
    n101 = node_id[z + 1, y, x + 1]
    n111 = node_id[z + 1, y + 1, x + 1]
    n011 = node_id[z + 1, y + 1, x]

    node_to_cell = np.stack([n000, n100, n110, n010, n001, n101, n111, n011], axis=-1)

    return node_to_cell


def points_to_cells(points: np.ndarray) -> np.ndarray:

    p000 = points[0:-1, 0:-1, 0:-1]
    p100 = points[0:-1, 0:-1, 1:]
    p110 = points[0:-1, 1:, 1:]
    p010 = points[0:-1, 1:, 0:-1]
    p001 = points[1:, 0:-1, 0:-1]
    p101 = points[1:, 0:-1, 1:]
    p111 = points[1:, 1:, 1:]
    p011 = points[1:, 1:, 0:-1]

    centers = (p000 + p100 + p110 + p010 + p001 + p101 + p111 + p011) / 8.0
    # centers = centers_zyx.reshape((-1, 3))

    return centers


def get_face_idx(info: MicrostructureInfo) -> tuple[int, int, int, int]:

    nz, ny, nx = info.nz, info.ny, info.nx
    phase = info.ph_id.reshape(nz, ny, nx)
    phase_mask = phase == info.gas_phase
    ngas_xmax = phase_mask[nz // 2, ny // 2, nx // 2 :].sum()
    ngas_xmin = phase_mask[nz // 2, ny // 2, 0 : nx // 2].sum()
    ngas_ymin = phase_mask[nz // 2, 0 : ny // 2, nx // 2].sum()
    ngas_ymax = phase_mask[nz // 2, ny // 2 :, nx // 2].sum()

    ixmin = ngas_xmin
    ixmax = nx - ngas_xmax - 1
    iymin = ngas_ymin
    iymax = ny - ngas_ymax - 1

    return ixmin, ixmax, iymin, iymax


def calc_displacement(
    pos: np.ndarray,
    pos0: np.ndarray,
    info: MicrostructureInfo,
    mode: Literal["node", "cell"] = "cell",
) -> np.ndarray:

    nz, ny, nx = info.nz, info.ny, info.nx
    ixmin, ixmax, iymin, iymax = get_face_idx(info)
    if mode == "cell":
        disp = np.zeros((nz, ny, nx))
    elif mode == "node":
        disp = np.zeros((nz + 1, ny + 1, nx + 1))
        ixmax += 1
        iymax += 1

    dx_min = -(pos[:, :, ixmin, 0] - pos0[:, :, ixmin, 0])
    dx_max = pos[:, :, ixmax, 0] - pos0[:, :, ixmax, 0]
    dy_min = -(pos[:, iymin, :, 1] - pos0[:, iymin, :, 1])
    dy_max = pos[:, iymax, :, 1] - pos0[:, iymax, :, 1]

    disp[:, :, ixmin] = dx_min
    disp[:, :, ixmax] = dx_max
    disp[:, iymin, :] = dy_min
    disp[:, iymax, :] = dy_max

    # average the corners
    disp[:, iymin, ixmin] = (dx_min[:, 0] + dy_min[:, 0]) / 2
    disp[:, iymin, ixmax] = (dx_max[:, 0] + dy_min[:, -1]) / 2
    disp[:, iymax, ixmin] = (dx_min[:, -1] + dy_max[:, 0]) / 2
    disp[:, iymax, ixmax] = (dx_max[:, -1] + dy_max[:, -1]) / 2

    return disp


def calc_height(
    pos: np.ndarray, info: MicrostructureInfo, mode: Literal["node", "cell"] = "cell"
) -> np.ndarray:

    nz, ny, nx = info.nz, info.ny, info.nx
    ixmin, ixmax, iymin, iymax = get_face_idx(info)

    if mode == "cell":
        height = np.zeros((nz, ny, nx))
    elif mode == "node":
        height = np.zeros((nz + 1, ny + 1, nx + 1))
        ixmax += 1
        iymax += 1

    mean_xmin = pos[:, :, ixmin, 0].mean()
    mean_ymin = pos[:, iymin, :, 1].mean()
    mean_xmax = pos[:, :, ixmax, 0].mean()
    mean_ymax = pos[:, iymax, :, 1].mean()

    height_xmin = -(pos[:, :, ixmin, 0] - mean_xmin)
    height_ymin = -(pos[:, iymin, :, 1] - mean_ymin)
    height_xmax = pos[:, :, ixmax, 0] - mean_xmax
    height_ymax = pos[:, iymax, :, 1] - mean_ymax

    height[:, :, ixmin] = height_xmin
    height[:, :, ixmax] = height_xmax
    height[:, iymin, :] = height_ymin
    height[:, iymax, :] = height_ymax

    # average the corners
    height[:, iymin, ixmin] = (height_xmin[:, 0] + height_ymin[:, 0]) / 2
    height[:, iymin, ixmax] = (height_xmax[:, 0] + height_ymin[:, -1]) / 2
    height[:, iymax, ixmin] = (height_xmin[:, -1] + height_ymax[:, 0]) / 2
    height[:, iymax, ixmax] = (height_xmax[:, -1] + height_ymax[:, -1]) / 2

    return height


def calculate_sa(
    height: np.ndarray, info: MicrostructureInfo, mode: Literal["node", "cell"] = "cell"
):

    ixmin, ixmax, iymin, iymax = get_face_idx(info)
    if mode == "node":
        ixmax += 1
        iymax += 1

    sa_xmin = np.abs(height[:, :, ixmin]).mean()
    sa_xmax = np.abs(height[:, :, ixmax]).mean()
    sa_ymin = np.abs(height[:, iymin, :]).mean()
    sa_ymax = np.abs(height[:, iymax, :]).mean()

    return sa_xmin, sa_xmax, sa_ymin, sa_ymax


def calculate_sz(
    height: np.ndarray, info: MicrostructureInfo, mode: Literal["node", "cell"] = "cell"
):

    ixmin, ixmax, iymin, iymax = get_face_idx(info)
    if mode == "node":
        ixmax += 1
        iymax += 1

    sz_xmin = height[:, :, ixmin].max() - height[:, :, ixmin].min()
    sz_xmax = height[:, :, ixmax].max() - height[:, :, ixmax].min()
    sz_ymin = height[:, iymin, :].max() - height[:, iymin, :].min()
    sz_ymax = height[:, iymax, :].max() - height[:, iymax, :].min()

    return sz_xmin, sz_xmax, sz_ymin, sz_ymax


def calc_max_slip(
    slip: np.ndarray, info: MicrostructureInfo, mode: Literal["node", "cell"] = "cell"
):

    ixmin, ixmax, iymin, iymax = get_face_idx(info)
    if mode == "node":
        ixmax += 1
        iymax += 1

    max_slip_xmin = slip[:, :, ixmin].max()
    max_slip_xmax = slip[:, :, ixmax].max()
    max_slip_ymin = slip[:, iymin, :].max()
    max_slip_ymax = slip[:, iymax, :].max()

    return max_slip_xmin, max_slip_xmax, max_slip_ymin, max_slip_ymax


def calc_mean_slip(
    slip: np.ndarray, info: MicrostructureInfo, mode: Literal["node", "cell"] = "cell"
):

    ixmin, ixmax, iymin, iymax = get_face_idx(info)
    if mode == "node":
        ixmax += 1
        iymax += 1

    avg_slip_xmin = slip[:, :, ixmin].mean()
    avg_slip_xmax = slip[:, :, ixmax].mean()
    avg_slip_ymin = slip[:, iymin, :].mean()
    avg_slip_ymax = slip[:, iymax, :].mean()

    return avg_slip_xmin, avg_slip_xmax, avg_slip_ymin, avg_slip_ymax


def total_strain_to_plastic(series: np.ndarray, ref_index: int) -> np.ndarray:
    """Align a series by subtracting its value at `ref_index`, then drop [0:ref_index]."""

    if ref_index < 0 or ref_index >= len(series):
        raise IndexError(
            f"ref_index {ref_index} out of bounds for series length {len(series)}"
        )
    return series[ref_index:] - series[ref_index]


def convert_load(
    load_mpa: int,
    info: MicrostructureInfo,
) -> float:
    """Convert a user-requested axial load to the macroscopic scauchy33 value.

    Problem
    -------
    In this solver, the macroscopic stress boundary condition (process%scauchy)
    is enforced against the *domain-average* Cauchy stress <sigma> (all voxels).
    When a compliant dummy gas phase is present (igas=1 or igas=2), the gas
    voxels are intended to carry little traction, so the solid-only average
    stress can exceed the imposed domain-average stress.

    Convention
    ----------
    We interpret `--load` / cfg.load as the
    desired *solid-average* sigma_33 in MPa.

    This function returns the stress value that should be written into fft.in
    (scauchy33) so that, approximately:

        <sigma_33>_solid ~= load_mpa

    Conversion model (approximate)
    ------------------------------
    Let f_s be the reference solid volume fraction and assume gas voxels carry
    negligible sigma_33 compared to the solid (reasonable for a highly compliant
    gas layer). Then the domain average is:

        <sigma_33>_all ~= f_s * <sigma_33>_solid

    Therefore, to target a solid-average load of L, we impose:

        scauchy33 := f_s * L

    Parameters
    ----------
    load_mpa:
        Desired solid-average sigma_33 (MPa).
    microstructure:
        Path to the microstructure file used for the run.
    igas:
        Phase-#2 gas mode used by the solver:
        - 0: no gas behavior (no conversion)
        - 1/2: gas-like behavior enabled (conversion applied if gas voxels exist)

    Returns
    -------
    Applied macroscopic scauchy33 value (MPa) to write into fft.in.
    """

    if load_mpa is None:
        raise TypeError("load_mpa must be a float")

    solid_frac, gas_frac = info.solid_frac, info.gas_frac
    if gas_frac <= 0.0:
        return float(load_mpa)
    if solid_frac <= 0.0:
        raise ValueError(
            f"Invalid solid fraction {solid_frac} inferred from microstructure={info.path}"
        )

    return float(load_mpa) * float(solid_frac)


def fit_plane(height: np.ndarray) -> np.ndarray:
    """Remove best-fit plane z = a*x + b*y + c, ignoring NaNs in the fit."""
    rows, cols = height.shape
    X, Y = np.meshgrid(np.arange(cols), np.arange(rows))

    mask = np.isfinite(height)
    if not np.any(mask):
        return np.full_like(height, np.nan, dtype=float)

    A = np.c_[X[mask], Y[mask], np.ones(mask.sum())]
    b = height[mask]

    C, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    Z_fit = C[0] * X + C[1] * Y + C[2]

    return Z_fit


def subtract_initial_plane(height_stack: np.ndarray) -> np.ndarray:
    """
    Subtract the t=0 best-fit plane for each sample from all timesteps.

    Parameters
    ----------
    height_stack : ndarray, shape (time, sample, y, x)

    Returns
    -------
    corrected : ndarray, shape (time, sample, y, x)
    """
    corrected = height_stack.copy()

    for s in range(height_stack.shape[1]):
        z0 = height_stack[0, s]
        if np.isnan(z0).all():
            continue
        p0 = fit_plane(z0)  # (y, x)
        corrected[:, s] = height_stack[:, s] - p0[None, :, :]

    return corrected


def calc_exp_roughness(
    heights: dict[str, np.ndarray], resolution: Literal[10, 50]
) -> dict[str, np.ndarray]:

    roughness = {}
    for key, height_stack in heights.items():

        load = key.split("_")[-1]
        centered = height_stack - np.nanmean(height_stack, axis=(2, 3), keepdims=True)
        mad = np.abs(centered)
        sa_per_sample = np.nanmean(mad, axis=(2, 3))  # (time, sample)
        sa_mean = np.nanmean(sa_per_sample, axis=1)
        sa_std = np.nanstd(sa_per_sample, axis=1, ddof=1)
        sz_per_sample = np.nanmax(height_stack, axis=(2, 3)) - np.nanmin(
            height_stack, axis=(2, 3)
        )
        sz_mean = np.nanmean(sz_per_sample, axis=1)
        sz_std = np.nanstd(sz_per_sample, axis=1, ddof=1)

        roughness[f"mean_sa_{load}_{resolution}x"] = sa_mean
        roughness[f"std_sa_{load}_{resolution}x"] = sa_std
        roughness[f"mean_sz_{load}_{resolution}x"] = sz_mean
        roughness[f"std_sz_{load}_{resolution}x"] = sz_std

    return roughness

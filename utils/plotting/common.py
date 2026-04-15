from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import re

import matplotlib.pyplot as plt


DEFAULT_FIGSIZE: tuple[float, float] = (6, 4)


def style_axes(ax: plt.Axes, *, grid_alpha: float = 0.3) -> None:
    """Apply common styling to axes."""

    ax.grid(True, alpha=grid_alpha)


def new_fig_ax(
    *, figsize: tuple[float, float] = DEFAULT_FIGSIZE
) -> tuple[plt.Figure, plt.Axes]:
    """Create a new figure/axes pair with the repo default sizing."""

    fig, ax = plt.subplots(figsize=figsize, dpi=200)
    return fig, ax


def save_fig(fig: plt.Figure, path: Path, *, dpi: int = 200) -> Path:
    """Save a figure and close it.

    Returns the written file path.
    """

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Reserve extra headroom so the suptitle doesn't collide with axes titles.
    # Two-line suptitles (plastic + elastic) need more padding.
    top = 1.0
    suptitle = getattr(fig, "_suptitle", None)
    if suptitle is not None:
        try:
            if "\n" in suptitle.get_text():
                top = 0.96
        except Exception:
            pass
    fig.tight_layout(rect=(0.0, 0.0, 1.0, top))
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def format_exp_sci(x: float) -> str:
    s = f"{x:.0e}"
    # Collapse exponent: 1e-03 -> 1e-3, 1e+06 -> 1e+6
    s = re.sub(r"e([+-])0*(\d+)$", r"e\1\2", s)
    return s


def parse_numeric_prefix(line: str) -> list[float]:
    toks = line.strip().split()
    nums: list[float] = []
    for tok in toks:
        try:
            nums.append(float(tok.replace("D", "E").replace("d", "e")))
        except ValueError:
            break
    return nums


def format_thousands_mpa(x_mpa: float) -> str:
    """Format MPa as thousands of MPa with no decimals."""

    return str(int(round(float(x_mpa) / 1000.0)))


def find_cupl2_sx(
    start_dir: Path, *, gas: bool = False, max_depth: int = 6
) -> Path | None:
    cur = Path(start_dir)
    for _ in range(int(max_depth) + 1):
        if gas:
            cand = cur / "cupl2_gas.sx"
        else:
            cand = cur / "cupl2.sx"
        if cand.exists():
            return cand
        if cur.parent == cur:
            break
        cur = cur.parent
    return None


def find_cuel_sx(
    start_dir: Path, *, gas: bool = False, max_depth: int = 6
) -> Path | None:
    cur = Path(start_dir)
    for _ in range(int(max_depth) + 1):
        if gas:
            cand = cur / "cuel_gas.sx"
        else:
            cand = cur / "cuel.sx"
        if cand.exists():
            return cand
        if cur.parent == cur:
            break
        cur = cur.parent
    return None


def parse_numeric_prefix(line: str) -> list[float]:
    toks = line.strip().split()
    nums: list[float] = []
    for tok in toks:
        try:
            nums.append(float(tok.replace("D", "E").replace("d", "e")))
        except ValueError:
            break
    return nums


@lru_cache(maxsize=256)
def plastic_params_from_cupl2(cupl2_path: str) -> dict[str, float] | None:
    """Extract a small, stable subset of plastic parameters from cupl2.sx.

    We look for lines whose trailing comments list the expected parameter names:
    - modex,nsmx,nrsx,gamd0x,...
    - tau0xf,tau0xb,tau1x,thet0,thet1,...
    """

    p = Path(cupl2_path)
    if not p.exists():
        return None

    nrsx: float | None = None
    gamd0x: float | None = None
    tau0xf: float | None = None
    tau1x: float | None = None
    thet0: float | None = None
    thet1: float | None = None

    for raw in p.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        if "modex,nsmx,nrsx" in line and nrsx is None:
            nums = parse_numeric_prefix(line)
            if len(nums) >= 4:
                nrsx = float(nums[2])
                gamd0x = float(nums[3])
                continue

        if "tau0xf,tau0xb,tau1x,thet0,thet1" in line and tau0xf is None:
            nums = parse_numeric_prefix(line)
            if len(nums) >= 5:
                tau0xf = float(nums[0])
                tau1x = float(nums[2])
                thet0 = float(nums[3])
                thet1 = float(nums[4])
                continue

    vals = (nrsx, gamd0x, tau0xf, tau1x, thet0, thet1)
    if any(v is None for v in vals):
        return None

    return {
        "nrsx": float(nrsx),  # type: ignore[arg-type]
        "gamd0x": float(gamd0x),  # type: ignore[arg-type]
        "tau0xf": float(tau0xf),  # type: ignore[arg-type]
        "tau1x": float(tau1x),  # type: ignore[arg-type]
        "thet0": float(thet0),  # type: ignore[arg-type]
        "thet1": float(thet1),  # type: ignore[arg-type]
    }


@lru_cache(maxsize=256)
def elastic_params_from_cuel(cuel_path: str) -> dict[str, float] | None:
    """Extract cubic elastic constants from cuel.sx (ISO=0 format only)."""

    p = Path(cuel_path)
    if not p.exists():
        return None

    lines = [ln.strip() for ln in p.read_text().splitlines() if ln.strip()]
    if not lines:
        return None

    first = parse_numeric_prefix(lines[0])
    if not first:
        return None
    iso = int(round(first[0]))
    if iso != 0:
        return None

    rows: list[list[float]] = []
    for raw in lines[1:]:
        nums = parse_numeric_prefix(raw)
        if len(nums) >= 6:
            rows.append(nums[:6])
        if len(rows) >= 6:
            break

    if len(rows) < 6:
        return None

    return {
        "c11": float(rows[0][0]),
        "c12": float(rows[0][1]),
        "c44": float(rows[3][3]),
    }


@lru_cache(maxsize=256)
def plastic_suptitle_text(start_dir: str) -> str | None:
    cupl2 = find_cupl2_sx(Path(start_dir))
    if cupl2 is None:
        return None
    params = plastic_params_from_cupl2(str(cupl2))
    if params is None:
        return None

    n = int(round(params["nrsx"]))
    gamd0 = format_exp_sci(params["gamd0x"])
    tau0 = int(round(params["tau0xf"]))
    tau1 = int(round(params["tau1x"]))
    theta0 = int(round(params["thet0"]))
    theta1 = int(round(params["thet1"]))

    # Use mathtext for compact, readable parameter symbols.
    return (
        r"$n=%d,\ \dot{\gamma}_0=%s,\ \tau_0=%d,\ \tau_1=%d,\ \theta_0=%d,\ \theta_1=%d$"
        % (n, gamd0, tau0, tau1, theta0, theta1)
    )


@lru_cache(maxsize=256)
def elastic_suptitle_text(start_dir: str, gas: bool = False) -> str | None:
    cuel = find_cuel_sx(Path(start_dir), gas=gas)
    if cuel is None:
        return ""
    params = elastic_params_from_cuel(str(cuel))
    if params is None:
        return ""

    c11 = format_thousands_mpa(params["c11"])
    c12 = format_thousands_mpa(params["c12"])
    c44 = format_thousands_mpa(params["c44"])

    return r"$C_{11}=%s,\ C_{12}=%s,\ C_{44}=%s$" % (c11, c12, c44)


@lru_cache(maxsize=256)
def gas_plastic_suptitle_text(start_dir: str) -> str | None:
    cupl2 = find_cupl2_sx(Path(start_dir), gas=True)
    if cupl2 is None:
        return ""
    params = plastic_params_from_cupl2(str(cupl2))
    if params is None:
        return ""

    n = int(round(params["nrsx"]))
    gamd0 = format_exp_sci(params["gamd0x"])
    tau0 = int(round(params["tau0xf"]))
    tau1 = int(round(params["tau1x"]))
    theta0 = int(round(params["thet0"]))
    theta1 = int(round(params["thet1"]))

    # Use mathtext for compact, readable parameter symbols.
    return (
        r"$n=%d,\ \dot{\gamma}_0=%s,\ \tau_0=%d,\ \tau_1=%d,\ \theta_0=%d,\ \theta_1=%d$"
        % (n, gamd0, tau0, tau1, theta0, theta1)
    )


@lru_cache(maxsize=256)
def material_suptitle_text(
    start_dir: str,
    gas: bool = False,
) -> str | None:
    lines = []
    plastic = plastic_suptitle_text(start_dir)
    elastic = elastic_suptitle_text(start_dir)
    lines.append(plastic)
    lines.append(elastic)
    if gas:
        plastic_gas = gas_plastic_suptitle_text(start_dir)
        elastic_gas = elastic_suptitle_text(start_dir, gas=True)
        lines.append(plastic_gas)
        lines.append(elastic_gas)

    text = "\n".join(line for line in lines if line)

    return text


def material_params_suptitle(
    fig: plt.Figure,
    *,
    start_dir: str | Path,
    text: None | str = None,
    gas: bool = False,
) -> bool:
    """Add a material-parameter suptitle when parameter files are discoverable.

    If the parameter files can't be found or parsed, this is a no-op.
    """

    # Don't overwrite an existing suptitle.
    if getattr(fig, "_suptitle", None) is not None:
        t = fig._suptitle.get_text()  # type: ignore[attr-defined]
        if t:
            return True

    if text is None:
        text = material_suptitle_text(
            str(Path(start_dir)),
            gas=gas,
        )

    if not text:
        return False

    fig.suptitle(text, y=0.95, fontsize=11)

    return True


@dataclass(frozen=True)
class VerticalLineStyle:
    """Styling options for vertical marker lines."""

    color: str = "k"
    linestyle: str = "--"
    alpha: float = 0.5
    linewidth: float = 1.5
    zorder: int = 0


def draw_vertical_lines(
    ax: plt.Axes,
    xs: list[float],
    *,
    label: str | None = None,
    style: VerticalLineStyle = VerticalLineStyle(),
) -> None:
    """Draw vertical lines at x-positions.

    If label is provided, it is used only for the first line to avoid legend
    spam.
    """

    if not xs:
        return
    first = True
    for x in xs:
        ax.axvline(
            x,
            color=style.color,
            linestyle=style.linestyle,
            alpha=style.alpha,
            linewidth=style.linewidth,
            zorder=style.zorder,
            label=label if (label is not None and first) else None,
        )
        first = False

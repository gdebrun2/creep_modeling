from __future__ import annotations

from pathlib import Path
import os
from typing import TYPE_CHECKING, Any
from config import (
    SolidElastic,
    SolidPlastic,
    GasElastic,
    GasPlastic,
    RESULTS_DIR,
    RUN_INDEX_FILE,
    RESULTS_NEXT_ID_FILE,
)
import datetime
import numpy as np
import tempfile
import re
import json

# from utils.run import RunConfig

if TYPE_CHECKING:
    from write_fft import SimCase


def write_text(path: Path, text: str) -> None:
    """Write text to disk with a trailing newline."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if not text.endswith("\n"):
        text += "\n"
    path.write_text(text)


def write_cuel(
    cuel_path: Path,
    *,
    elastic: GasElastic | SolidElastic,
    iso: int = 0,
) -> None:
    """Write a `cuel.sx` file.

    Parameters
    ----------
    cuel_path:
        Output path.
    sigma_xx, sigma_xy, tau:
        Cubic elastic constants for ISO=0 (MPa).
    iso:
        Which format to write.
        - 0: cubic elastic constants (6x6) format used for steel.
        - 1: isotropic Young's modulus / Poisson ratio line.
    young, nu:
        Isotropic properties used only when iso=1.
    """

    if iso == 0:
        text = f"""0 					ISO
{elastic.c11} {elastic.c12} {elastic.c12} 0.0 0.0 0.0
{elastic.c12} {elastic.c11} {elastic.c12} 0.0 0.0 0.0
{elastic.c12} {elastic.c12} {elastic.c11} 0.0 0.0 0.0
0.0 0.0 0.0 {elastic.c44} 0.0 0.0
0.0 0.0 0.0 0.0 {elastic.c44} 0.0
0.0 0.0 0.0 0.0 0.0 {elastic.c44}
"""
    elif iso == 1:
        text = f"""1 					ISO
{elastic.young}   {elastic.nu}   YOUNG(MPa),NU (V+R/2)
"""
    else:
        raise ValueError(f"Unsupported iso={iso}; expected 0 or 1")

    write_text(cuel_path, text)


def write_cupl2(cupl2_path: Path, *, plastic: GasPlastic | SolidPlastic) -> None:
    """Write `cupl2.sx` (FCC {111}<110> slip) with Voce-style hardening params.

    Parameter meanings (typical EVPFFT / CP-Voce conventions)
    --------------------------------------------------------
    Rate dependence
    - gamd0x: reference shear rate \\dot{\\gamma}_0 (1/s)
    - nrsx: rate sensitivity exponent n (often m = 1/n in the literature)

    Slip resistance / hardening (Voce-type)
    - tau0xf: initial CRSS / slip resistance in forward direction (often g0 or tau0)
    - tau0xb: initial CRSS in reverse direction (captures forward/backward asymmetry)
    - tau1x: saturation increment; often g_sat = g0 + tau1
    - thet0: initial hardening rate / modulus (theta0)
    - thet1: asymptotic hardening rate / modulus (theta1), typically <= thet0

    Latent hardening
    - hselfx: self hardening interaction coefficient (diagonal)
    - hlatex: latent hardening interaction coefficient (off-diagonal)
    """

    text = f"""SLIP SYSTEMS FOR CUBIC CRYSTAL
CUBIC          icryst
1.   1.   1.   crystal axes (cdim(i))
1              nmodesx (total # of modes listed in the file)
1              nmodes  (# of modes to be used in the calculation)
1              mode(i) (label of the modes to be used)
{{111}}<110>   SLIP
1   12   {plastic.nrsx:.8f}   {plastic.gamd0x:.8e}   0.0  0                 modex,nsmx,nrsx,gamd0x,twshx,isectwx
{plastic.tau0xf}  {plastic.tau0xb}   {plastic.tau1x}   {plastic.thet0}  {plastic.thet1}  0  tau0xf,tau0xb,tau1x,thet0,thet1,hpfac
{plastic.hselfx}   {plastic.hlatex}                                               hselfx,hlatex
1    1   -1      0    1    1                        SLIP (n-b)
1    1   -1      1    0    1
1    1   -1      1   -1    0
1   -1   -1      0    1   -1
1   -1   -1      1    0    1
1   -1   -1      1    1    0
1   -1    1      0    1    1
1   -1    1      1    0   -1
1   -1    1      1    1    0
1    1    1      0    1   -1
1    1    1      1    0   -1
1    1    1      1   -1    0
"""

    write_text(cupl2_path, text)


def os_system(cmd: str) -> int:
    """Thin wrapper for running shell commands.

    This exists so we can mock or intercept calls in future tests.
    """

    return os.system(cmd)


def allocate_next_run_id() -> int:
    """Allocate the next sequential numeric run id.

    Convention
    ----------
    - `_next_id.txt` stores the *next* ID to allocate.

    This is intentionally simple for a single-user workstation workflow.
    """

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Initialize from directory scan if the file doesn't exist.
    # We scan results/*/<run_id>/ directories, excluding results/calibrations.
    if not RESULTS_NEXT_ID_FILE.exists():
        max_existing = 0
        for sim_dir in RESULTS_DIR.iterdir():
            if not sim_dir.is_dir():
                continue
            if sim_dir.name == "calibrations":
                continue
            for p in sim_dir.iterdir():
                if p.is_dir() and p.name.isdigit():
                    max_existing = max(max_existing, int(p.name))
        RESULTS_NEXT_ID_FILE.write_text(str(max_existing + 1))

    nxt = int(RESULTS_NEXT_ID_FILE.read_text().strip())
    RESULTS_NEXT_ID_FILE.write_text(str(nxt + 1))
    return nxt


def append_run_index_row(row: dict[str, Any]) -> None:
    """Append a row to the global results index (results/run_index.csv).

    This is intentionally lightweight (CSV append) for a single-user workflow.

    Notes
    -----
    - We do not include git hashes by request.
    - Missing columns are left blank.
    - New keys in `row` expand the header (existing rows remain valid CSV).
    """

    import csv

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # If the file doesn't exist yet, create it with the provided columns.
    if not RUN_INDEX_FILE.exists():
        with RUN_INDEX_FILE.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(row.keys()))
            w.writeheader()
            w.writerow(row)
        return

    # If it exists, we may need to expand header.
    with RUN_INDEX_FILE.open("r", newline="") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            header = []

    keys = list(dict.fromkeys(header + list(row.keys())))
    if keys != header:
        # Rewrite with expanded header.
        with RUN_INDEX_FILE.open("r", newline="") as f:
            old = list(csv.DictReader(f))
        with RUN_INDEX_FILE.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(old)
            w.writerow(row)
        return

    with RUN_INDEX_FILE.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writerow(row)
    return


def run_dir_from_id(
    sim_case: SimCase,
    run_id: int,
    *,
    out_dir: Path | None = None,
) -> Path:
    """Return run directory path for a known run_id.

    This function is intentionally *path-only* (plus mkdir), so callers don't
    need to thread all run parameters through it.

    Convention
    ----------
    results/<sim_case>/<run_id>/
    """

    # New convention is: results/<sim_case>/<run_id>/ for *any* sim_case.
    # Keep results/calibrations reserved for calibration outputs.
    if str(sim_case) == "calibrations":
        raise ValueError("sim_case='calibrations' is reserved")

    run_dir = RESULTS_DIR / str(sim_case) / str(int(run_id))

    if out_dir is not None:
        run_dir = out_dir

    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def register_run(
    *,
    run_id: int,
    sim_case: SimCase,
    calib_id: int | None,
    igas: int,
    nthreads: int,
    notes: str | None = None,
    microstructure: Path | None,
    fft: Path | None,
    run_dir: Path,
    load: float | None = None,
    complgas: float | None = None,
    modal: int | None = None,
    mix: int | None = None,
):
    append_run_index_row(
        {
            "run_id": run_id,
            "created_at": datetime.datetime.now().isoformat(timespec="seconds"),
            "sim_case": sim_case,
            "calib_id": calib_id,
            "igas": igas,
            "load": load,
            "complgas": complgas,
            "mix": mix,
            "modal": modal,
            "fft": fft,
            "nthreads": nthreads,
            "run_dir": run_dir,
            "microstructure": microstructure,
            "notes": notes,
        }
    )
    return None


def allocate_run_dir(
    sim_case: SimCase,
    *,
    out_dir: Path | None = None,
) -> tuple[int, Path]:
    """Allocate a new run_id and create its directory.

    Returns
    -------
    (run_id, run_dir)
    """

    run_id = allocate_next_run_id()
    return run_id, run_dir_from_id(sim_case, run_id, out_dir=out_dir)


def append_fields(
    vtk_path: str | Path,
    *,
    cell_fields: dict[str, np.ndarray] | None = None,
    point_fields: dict[str, np.ndarray] | None = None,
    dtype_str: str = "double",
    n_per_line: int = 100,
) -> None:
    """Write or replace FIELD arrays in CELL_DATA / POINT_DATA sections.

    Existing FIELD arrays with the same name are replaced within their section.
    New field names are appended to that section's FIELD block. If the section
    or FIELD block does not exist yet, it is created.
    """

    _POINTS_RE = re.compile(r"^\s*POINTS\s+(\d+)\s+")
    _CELLS_RE = re.compile(r"^\s*CELLS\s+(\d+)\s+")
    _CELL_DATA_RE = re.compile(r"^\s*CELL_DATA\s+(\d+)\s*$")
    _POINT_DATA_RE = re.compile(r"^\s*POINT_DATA\s+(\d+)\s*$")
    _FIELD_RE = re.compile(r"^\s*FIELD\s+FieldData\s+(\d+)\s*$")

    def _as_field_array(values: np.ndarray, ntuples: int) -> tuple[np.ndarray, int]:
        """
        Normalize values into (flat, ncomp) where:
        - values can be (ntuples,) or (ntuples, ncomp)
        - flat is 1D float64 of length ntuples*ncomp in cell/point-major order
        """
        v = np.asarray(values)
        if v.ndim == 1:
            if v.size != ntuples:
                raise ValueError(f"Expected ({ntuples},) got {v.shape}")
            ncomp = 1
            flat = v.astype(np.float64, copy=False)
        elif v.ndim == 2:
            if v.shape[0] != ntuples:
                raise ValueError(f"Expected ({ntuples}, ncomp) got {v.shape}")
            ncomp = int(v.shape[1])
            flat = v.astype(np.float64, copy=False).ravel(
                order="C"
            )  # tuple-major: t0 comps, t1 comps, ...
        else:
            raise ValueError(f"values must be 1D or 2D; got {v.ndim}D {v.shape}")
        return flat, ncomp

    def _write_field_array(
        fout,
        name: str,
        flat: np.ndarray,
        ncomp: int,
        ntuples: int,
        *,
        dtype_str: str = "double",
        n_per_line: int = 100,
    ) -> None:
        # FIELD header: <name> <num_components> <num_tuples> <type>
        fout.write(f"{name} {ncomp} {ntuples} {dtype_str}\n")

        for i in range(0, flat.size, n_per_line):
            chunk = flat[i : i + n_per_line]
            s = np.array2string(
                chunk,
                separator=" ",
                max_line_width=10**9,
                formatter={"float_kind": lambda x: f"{x:.8E}"},
            )
            # array2string wraps 1D arrays in brackets: "[ ... ]" -> strip them
            fout.write(s[1:-1] + "\n")

    vtk_path = Path(vtk_path)
    cell_fields = dict(cell_fields or {})
    point_fields = dict(point_fields or {})

    if not cell_fields and not point_fields:
        return

    # ---------- Pass 1: discover npoints and ncells ----------
    npoints = None
    ncells = None
    cells_from_cells_header = None

    has_cell_data = False
    has_point_data = False

    with vtk_path.open("r", encoding="utf-8", newline="") as f:
        for line in f:
            if npoints is None:
                m = _POINTS_RE.match(line)
                if m:
                    npoints = int(m.group(1))
            if cells_from_cells_header is None:
                m = _CELLS_RE.match(line)
                if m:
                    cells_from_cells_header = int(m.group(1))
            m = _CELL_DATA_RE.match(line)
            if m:
                has_cell_data = True
                ncells = int(m.group(1))
            m = _POINT_DATA_RE.match(line)
            if m:
                has_point_data = True
                # we *could* validate this equals npoints later

    if npoints is None:
        raise RuntimeError("Could not find a 'POINTS <n> <type>' line.")
    if ncells is None:
        # if there's no CELL_DATA header, fall back to CELLS header
        if cells_from_cells_header is None:
            raise RuntimeError(
                "Could not determine number of cells (no CELL_DATA and no CELLS header)."
            )
        ncells = cells_from_cells_header

    # Validate shapes up front and pre-flatten (also fixes dtype)
    cell_prepped: list[tuple[str, np.ndarray, int]] = []
    for name, arr in cell_fields.items():
        flat, ncomp = _as_field_array(arr, ncells)
        cell_prepped.append((name, flat, ncomp))

    point_prepped: list[tuple[str, np.ndarray, int]] = []
    for name, arr in point_fields.items():
        flat, ncomp = _as_field_array(arr, npoints)
        point_prepped.append((name, flat, ncomp))

    cell_field_names = {name for name, _, _ in cell_prepped}
    point_field_names = {name for name, _, _ in point_prepped}

    k_cell = len(cell_prepped)
    k_point = len(point_prepped)

    def _write_prepped_fields(
        fout,
        prepped: list[tuple[str, np.ndarray, int]],
        ntuples: int,
    ) -> None:
        for nm, flat, ncomp in prepped:
            _write_field_array(
                fout,
                nm,
                flat,
                ncomp,
                ntuples,
                dtype_str=dtype_str,
                n_per_line=n_per_line,
            )

    def _read_field_block(fin, n_arrays: int) -> list[tuple[str, list[str]]]:
        blocks: list[tuple[str, list[str]]] = []

        for _ in range(n_arrays):
            header = fin.readline()
            if header == "":
                raise RuntimeError("Unexpected EOF while reading FIELD array header.")

            parts = header.strip().split()
            if len(parts) < 4:
                raise RuntimeError(f"Malformed FIELD array header: {header.rstrip()}")

            name = parts[0]
            try:
                ncomp = int(parts[1])
                ntuples = int(parts[2])
            except ValueError as exc:
                raise RuntimeError(
                    f"Malformed FIELD array header: {header.rstrip()}"
                ) from exc

            needed = ncomp * ntuples
            n_values = 0
            block_lines = [header]
            while n_values < needed:
                data_line = fin.readline()
                if data_line == "":
                    raise RuntimeError(
                        f"Unexpected EOF while reading FIELD array '{name}'."
                    )
                block_lines.append(data_line)
                n_values += len(data_line.split())

            blocks.append((name, block_lines))

        return blocks

    # ---------- Pass 2: stream-copy + inject / replace ----------
    tmp_fd, tmp_name = tempfile.mkstemp(
        prefix=vtk_path.name + ".", suffix=".tmp", dir=str(vtk_path.parent)
    )
    os.close(tmp_fd)
    tmp_path = Path(tmp_name)

    section: str | None = None  # "cell" | "point" | None
    cell_inserted = False
    point_inserted = False

    try:
        with (
            vtk_path.open("r", encoding="utf-8", newline="") as fin,
            tmp_path.open("w", encoding="utf-8", newline="") as fout,
        ):
            while True:
                line = fin.readline()
                if line == "":
                    break

                if (
                    section == "point"
                    and (not point_inserted)
                    and _CELL_DATA_RE.match(line)
                ):
                    if k_point:
                        fout.write(f"FIELD FieldData {k_point}\n")
                        _write_prepped_fields(fout, point_prepped, npoints)
                    point_inserted = True
                    section = None

                if (
                    section == "cell"
                    and (not cell_inserted)
                    and _POINT_DATA_RE.match(line)
                ):
                    if k_cell:
                        fout.write(f"FIELD FieldData {k_cell}\n")
                        _write_prepped_fields(fout, cell_prepped, ncells)
                    cell_inserted = True
                    section = None

                m_field = _FIELD_RE.match(line)
                if m_field and section == "cell" and (not cell_inserted) and k_cell:
                    existing_blocks = _read_field_block(fin, int(m_field.group(1)))
                    retained_blocks = [
                        block
                        for name, block in existing_blocks
                        if name not in cell_field_names
                    ]
                    fout.write(f"FIELD FieldData {len(retained_blocks) + k_cell}\n")
                    for block in retained_blocks:
                        for block_line in block:
                            fout.write(block_line)
                    _write_prepped_fields(fout, cell_prepped, ncells)
                    cell_inserted = True
                    continue

                if m_field and section == "point" and (not point_inserted) and k_point:
                    existing_blocks = _read_field_block(fin, int(m_field.group(1)))
                    retained_blocks = [
                        block
                        for name, block in existing_blocks
                        if name not in point_field_names
                    ]
                    fout.write(f"FIELD FieldData {len(retained_blocks) + k_point}\n")
                    for block in retained_blocks:
                        for block_line in block:
                            fout.write(block_line)
                    _write_prepped_fields(fout, point_prepped, npoints)
                    point_inserted = True
                    continue

                fout.write(line)

                m = _CELL_DATA_RE.match(line)
                if m:
                    section = "cell"
                    continue

                m = _POINT_DATA_RE.match(line)
                if m:
                    n_decl = int(m.group(1))
                    if n_decl != npoints:
                        raise RuntimeError(
                            f"POINT_DATA {n_decl} does not match POINTS {npoints}."
                        )
                    section = "point"
                    continue

            if section == "cell" and (not cell_inserted) and k_cell:
                fout.write(f"\nFIELD FieldData {k_cell}\n")
                _write_prepped_fields(fout, cell_prepped, ncells)
                cell_inserted = True

            if section == "point" and (not point_inserted) and k_point:
                fout.write(f"\nFIELD FieldData {k_point}\n")
                _write_prepped_fields(fout, point_prepped, npoints)
                point_inserted = True

            if (not has_cell_data) and k_cell:
                fout.write(f"\nCELL_DATA {ncells}\n")
                fout.write(f"FIELD FieldData {k_cell}\n")
                _write_prepped_fields(fout, cell_prepped, ncells)

            if (not has_point_data) and k_point:
                fout.write(f"\nPOINT_DATA {npoints}\n")
                fout.write(f"FIELD FieldData {k_point}\n")
                _write_prepped_fields(fout, point_prepped, npoints)

        tmp_path.replace(vtk_path)

    finally:
        tmp_path.unlink(missing_ok=True)


def append_cell(
    vtk_path: str | Path,
    name: str,
    values: np.ndarray,
    buffer: int = 100,
) -> None:
    vtk_path = Path(vtk_path)

    tmp_fd, tmp_name = tempfile.mkstemp(
        prefix=vtk_path.name + ".", suffix=".tmp", dir=str(vtk_path.parent)
    )
    os.close(tmp_fd)
    tmp_path = Path(tmp_name)

    if values.ndim == 1:
        ncomp = 1
    elif values.ndim == 2:
        ncomp = values.shape[-1]
    else:
        raise ValueError("Error appending cell. Values must be 1D or 2D")

    patched = False
    try:
        with (
            vtk_path.open("r", encoding="utf-8", newline="") as fin,
            tmp_path.open("w", encoding="utf-8", newline="") as fout,
        ):
            for line in fin:
                if line.strip().startswith("FIELD"):
                    # if not patched:
                    n_old = int(line.split()[-1])
                    fout.write(f"FIELD FieldData {n_old + 1}\n")
                    patched = True
                else:
                    fout.write(line)

            # Append the new field
            ncells = values.size
            fout.write(f"{name} {ncomp} {ncells} double\n")
            for i in range(0, ncells, buffer):
                chunk = values[i : i + buffer]
                s = np.array2string(
                    chunk,
                    separator=" ",
                    max_line_width=10**9,
                    formatter={"float_kind": lambda x: f"{x:.8E}"},
                )
                fout.write(s[1:-1] + "\n")

        tmp_path.replace(vtk_path)

    finally:

        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)


def write_vtk_series(
    vtk_dir: str | Path, time: np.ndarray, out_name: str = "microstr.vtk.series"
) -> Path:
    vtk_dir = Path(vtk_dir)
    files = sorted(vtk_dir.glob("*.vtk"))

    if not files:
        raise FileNotFoundError(f"No .vtk files found in {vtk_dir}")

    entries = []
    for i, f in enumerate(files):
        entries.append({"name": f.name, "time": time[i]})

    meta = {"file-series-version": "1.0", "files": entries}

    out_path = vtk_dir / out_name
    out_path.write_text(json.dumps(meta, indent=2))
    print(f"Wrote: {out_path}")
    return None

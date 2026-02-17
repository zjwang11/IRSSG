#!/usr/bin/env python3
"""Generate tbbox.in from Wannier90 inputs/outputs in the current directory.

Usage:
  python tools/gen_tbbox.py --kmesh 10
"""

from __future__ import annotations

import argparse
import math
import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Tuple


@lru_cache(maxsize=None)
def read_text(path: Path) -> str:
    print(f"Reading file: {path}")
    return path.read_text(encoding="utf-8", errors="ignore")


def find_block(text: str, name: str) -> List[str]:
    pattern = re.compile(
        rf"begin\s+{re.escape(name)}\s*(.*?)\s*end\s+{re.escape(name)}",
        re.IGNORECASE | re.DOTALL,
    )
    m = pattern.search(text)
    if not m:
        return []
    lines = [ln.strip() for ln in m.group(1).strip().splitlines() if ln.strip()]
    return lines


def transpose(m: List[List[float]]) -> List[List[float]]:
    return [list(row) for row in zip(*m)]


def invert_3x3(m: List[List[float]]) -> List[List[float]]:
    a, b, c = m[0]
    d, e, f = m[1]
    g, h, i = m[2]
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if abs(det) < 1e-14:
        raise ValueError("Singular 3x3 matrix; cannot invert.")
    inv_det = 1.0 / det
    return [
        [(e * i - f * h) * inv_det, (c * h - b * i) * inv_det, (b * f - c * e) * inv_det],
        [(f * g - d * i) * inv_det, (a * i - c * g) * inv_det, (c * d - a * f) * inv_det],
        [(d * h - e * g) * inv_det, (b * g - a * h) * inv_det, (a * e - b * d) * inv_det],
    ]


def matvec(m: List[List[float]], v: List[float]) -> List[float]:
    return [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]


def wrap_frac(v: List[float]) -> List[float]:
    return [x - math.floor(x) for x in v]


def parse_spinors(text: str) -> bool:
    m = re.search(r"^\\s*spinors\\s*=\\s*([.\\w]+)", text, re.IGNORECASE | re.MULTILINE)
    if not m:
        return False
    val = m.group(1).strip().lower()
    return val in (".true.", "true", "t", ".t.")


def first_existing(candidates: List[Path]) -> Path | None:
    for p in candidates:
        if p.exists():
            return p
    return None


def detect_ispin_from_win(base: Path) -> int:
    # User-defined rule:
    # if wannier90.up.win or wannier90.1.win exists -> ISPIN=2, else ISPIN=1.
    if (base / "wannier90.up.win").exists() or (base / "wannier90.1.win").exists():
        return 2
    return 1


def find_wannier_win(base: Path, ispin: int) -> Path:
    if ispin == 2:
        win = first_existing(
            [
                base / "wannier90.up.win",
                base / "wannier90.1.win",
                base / "wannier90.dn.win",
                base / "wannier90.2.win",
            ]
        )
    else:
        win = base / "wannier90.win" if (base / "wannier90.win").exists() else None

    if win is None:
        raise SystemExit("Missing Wannier win file for inferred ISPIN mode.")
    return win


def find_matching_wout(base: Path, win: Path) -> Path:
    prefix = win.name[: -len(".win")]
    return first_existing(
        [
            base / f"{prefix}.wout",
            base / "wannier90.wout",
            base / "wannier90.up.wout",
            base / "wannier90.dn.wout",
            base / "wannier90.1.wout",
            base / "wannier90.2.wout",
        ]
    ) or (base / f"{prefix}.wout")


def find_matching_labelinfo(base: Path, win: Path) -> Path:
    prefix = win.name[: -len(".win")]
    return first_existing(
        [
            base / f"{prefix}_band.labelinfo.dat",
            base / "wannier90_band.labelinfo.dat",
            base / "wannier90.up_band.labelinfo.dat",
            base / "wannier90.dn_band.labelinfo.dat",
            base / "wannier90.1_band.labelinfo.dat",
            base / "wannier90.2_band.labelinfo.dat",
        ]
    ) or (base / f"{prefix}_band.labelinfo.dat")


def find_spinpol_hr_files(base: Path) -> Tuple[Path | None, Path | None]:
    preferred_pairs = [
        ("wannier90.up_hr.dat", "wannier90.dn_hr.dat"),
        ("wannier90.1_hr.dat", "wannier90.2_hr.dat"),
        ("sp.up_hr.dat", "sp.dn_hr.dat"),
        ("sp.1_hr.dat", "sp.2_hr.dat"),
    ]
    for up_name, dn_name in preferred_pairs:
        up = base / up_name
        dn = base / dn_name
        if up.exists() and dn.exists():
            return up, dn

    for up in sorted(base.glob("*.up_hr.dat")):
        stem = up.name[: -len(".up_hr.dat")]
        dn = base / f"{stem}.dn_hr.dat"
        if dn.exists():
            return up, dn
    for h1 in sorted(base.glob("*.1_hr.dat")):
        stem = h1.name[: -len(".1_hr.dat")]
        h2 = base / f"{stem}.2_hr.dat"
        if h2.exists():
            return h1, h2
    return None, None


def find_nonspin_hr_file(base: Path) -> Path:
    return first_existing([base / "wannier90_hr.dat", base / "sp_hr.dat"]) or (base / "wannier90_hr.dat")


def parse_projections(text: str) -> List[Tuple[str, str]]:
    lines = find_block(text, "projections")
    out: List[Tuple[str, str]] = []
    for ln in lines:
        if ":" not in ln:
            continue
        species, proj = ln.split(":", 1)
        out.append((species.strip(), proj.strip()))
    return out


def proj_to_norb(proj: str) -> int:
    fams = set()
    tokens = re.split(r"[^A-Za-z]+", proj.lower())
    for tok in tokens:
        if not tok:
            continue
        for ch in set(tok):
            if ch in ("s", "p", "d", "f"):
                fams.add(ch)
    counts = {"s": 1, "p": 3, "d": 5, "f": 7}
    if not fams:
        return 1
    return sum(counts[ch] for ch in fams)


def parse_unit_cell_cart(text: str) -> List[List[float]]:
    lines = find_block(text, "unit_cell_cart")
    if not lines:
        raise ValueError("unit_cell_cart block not found in wannier90.win")
    cell = [[float(x) for x in ln.split()] for ln in lines]
    if len(cell) != 3 or any(len(row) != 3 for row in cell):
        raise ValueError("unit_cell_cart should have 3 lines with 3 numbers each.")
    return cell


def parse_atoms_cart(text: str) -> List[Tuple[str, List[float]]]:
    lines = find_block(text, "atoms_cart")
    if not lines:
        raise ValueError("atoms_cart block not found in wannier90.win")
    atoms = []
    for ln in lines:
        parts = ln.split()
        if len(parts) < 4:
            continue
        atoms.append((parts[0], [float(parts[1]), float(parts[2]), float(parts[3])]))
    return atoms


def parse_kpoint_path(text: str) -> List[List[float]]:
    lines = find_block(text, "kpoint_path")
    nodes: List[List[float]] = []
    for ln in lines:
        parts = ln.split()
        if len(parts) < 8:
            continue
        # label1 x1 y1 z1 label2 x2 y2 z2
        x1, y1, z1 = map(float, parts[1:4])
        x2, y2, z2 = map(float, parts[5:8])
        p1 = [x1, y1, z1]
        p2 = [x2, y2, z2]
        if not nodes or any(abs(a - b) > 1e-12 for a, b in zip(nodes[-1], p1)):
            nodes.append(p1)
        nodes.append(p2)
    return nodes


def parse_labelinfo(path: Path) -> List[List[float]]:
    if not path.exists():
        return []
    nodes: List[List[float]] = []
    for ln in read_text(path).splitlines():
        parts = ln.split()
        if len(parts) < 5:
            continue
        try:
            kx, ky, kz = map(float, parts[-3:])
        except ValueError:
            continue
        nodes.append([kx, ky, kz])
    return nodes


def infer_kmesh_from_labelinfo(path: Path) -> int | None:
    if not path.exists():
        return None
    idxs = []
    for ln in read_text(path).splitlines():
        parts = ln.split()
        if len(parts) < 2:
            continue
        try:
            idxs.append(int(parts[1]))
        except ValueError:
            continue
    if len(idxs) < 2:
        return None
    steps = [idxs[i + 1] - idxs[i] for i in range(len(idxs) - 1)]
    if all(s == steps[0] for s in steps):
        return steps[0]
    return None


def parse_kmesh_from_wout(path: Path) -> int | None:
    if not path.exists():
        return None
    text = read_text(path)
    m = re.search(r"Divisions along first K-path section\\s*:\\s*(\\d+)", text)
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def parse_last_final_state_rows(path: Path) -> List[Tuple[int, List[float], float]]:
    if not path.exists():
        return []
    lines = read_text(path).splitlines()
    starts = [i for i, ln in enumerate(lines) if "Final State" in ln]
    if not starts:
        return []
    start = starts[-1]

    pat = re.compile(
        r"WF centre and spread\s+(\d+)\s+\(\s*([-+0-9.Ee]+)\s*,\s*([-+0-9.Ee]+)\s*,\s*([-+0-9.Ee]+)\s*\)\s*([-+0-9.Ee]+)"
    )
    rows: List[Tuple[int, List[float], float]] = []
    for ln in lines[start + 1 :]:
        if "Sum of centres and spreads" in ln:
            break
        m = pat.search(ln)
        if not m:
            continue
        try:
            rows.append(
                (
                    int(m.group(1)),
                    [float(m.group(2)), float(m.group(3)), float(m.group(4))],
                    float(m.group(5)),
                )
            )
        except ValueError:
            continue

    rows.sort(key=lambda x: x[0])
    return rows


def infer_spincov_from_final_state_rows(rows: List[Tuple[int, List[float], float]]) -> int | None:
    if len(rows) < 2 or len(rows) % 2 != 0:
        return None

    vals = [(c[0], c[1], c[2], spr) for _, c, spr in rows]
    n = len(vals)
    half = n // 2

    def pair_dist(a: Tuple[float, float, float, float], b: Tuple[float, float, float, float]) -> float:
        dx = a[0] - b[0]
        dy = a[1] - b[1]
        dz = a[2] - b[2]
        ds = abs(a[3] - b[3])
        return math.sqrt(dx * dx + dy * dy + dz * dz) + ds

    block_score = sum(pair_dist(vals[i], vals[i + half]) for i in range(half)) / half
    inter_score = sum(pair_dist(vals[2 * i], vals[2 * i + 1]) for i in range(half)) / half
    # spincov=1: 1up..Nup,1dn..Ndn ; spincov=2: 1up,1dn,2up,2dn,...
    return 1 if block_score <= inter_score else 2


def infer_spincov_from_wout_final_state(path: Path) -> int | None:
    return infer_spincov_from_final_state_rows(parse_last_final_state_rows(path))


def select_atoms_nearest_to_final_state_centers(
    atoms: List[Tuple[str, List[float]]],
    cell: List[List[float]],
    final_state_rows: List[Tuple[int, List[float], float]],
) -> List[Tuple[str, List[float]]]:
    if not atoms or not final_state_rows:
        return atoms

    cell_t = transpose(cell)
    inv_cell_t = invert_3x3(cell_t)
    atom_frac = [wrap_frac(matvec(inv_cell_t, cart)) for _, cart in atoms]

    selected: List[int] = []
    seen = set()
    best_shift_by_atom: Dict[int, Tuple[float, List[int]]] = {}
    for _, center_cart, _ in final_state_rows:
        center_frac = matvec(inv_cell_t, center_cart)
        best_i = -1
        best_d2 = float("inf")
        best_n = [0, 0, 0]
        for i, af in enumerate(atom_frac):
            raw = [center_frac[j] - af[j] for j in range(3)]
            nimg = [int(round(x)) for x in raw]
            dfrac = [raw[j] - nimg[j] for j in range(3)]
            dcart = matvec(cell_t, dfrac)
            d2 = dcart[0] * dcart[0] + dcart[1] * dcart[1] + dcart[2] * dcart[2]
            if d2 < best_d2:
                best_d2 = d2
                best_i = i
                best_n = nimg
        if best_i >= 0 and best_i not in seen:
            seen.add(best_i)
            selected.append(best_i)
        if best_i >= 0:
            prev = best_shift_by_atom.get(best_i)
            if prev is None or best_d2 < prev[0]:
                best_shift_by_atom[best_i] = (best_d2, best_n)

    if not selected:
        return atoms

    shifted_atoms: List[Tuple[str, List[float]]] = []
    for i in selected:
        sp, cart = atoms[i]
        _, nimg = best_shift_by_atom.get(i, (0.0, [0, 0, 0]))
        shift_cart = matvec(cell_t, [float(nimg[0]), float(nimg[1]), float(nimg[2])])
        shifted_atoms.append((sp, [cart[0] + shift_cart[0], cart[1] + shift_cart[1], cart[2] + shift_cart[2]]))
    return shifted_atoms


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate tbbox.in from Wannier90 files.")
    ap.add_argument("--orbt", type=int, default=2, help="orbt value (default: 2).")
    ap.add_argument(
        "--spincov",
        type=int,
        choices=(1, 2),
        default=None,
        help="spincov override (default: auto from the last Final State block in wannier90.wout).",
    )
    args = ap.parse_args()

    base = Path(".")
    ispin = detect_ispin_from_win(base)
    win = find_wannier_win(base, ispin)

    text = read_text(win)
    soc = parse_spinors(text)
    cell = parse_unit_cell_cart(text)
    atoms = parse_atoms_cart(text)
    wout = find_matching_wout(base, win)
    final_state_rows = parse_last_final_state_rows(wout)
    has_final_state = bool(final_state_rows)
    atoms = select_atoms_nearest_to_final_state_centers(atoms, cell, final_state_rows)
    proj_list = parse_projections(text)

    # species order and norb mapping
    species_order: List[str] = []
    norb_by_species: Dict[str, int] = {}
    for sp, proj in proj_list:
        if sp not in species_order:
            species_order.append(sp)
        norb_by_species[sp] = proj_to_norb(proj)

    # fallback if projections missing
    if not species_order:
        species_order = sorted({sp for sp, _ in atoms})
        for sp in species_order:
            norb_by_species[sp] = 1

    type_index = {sp: i + 1 for i, sp in enumerate(species_order)}

    # cart -> fractional
    cell_t = transpose(cell)
    inv_cell_t = invert_3x3(cell_t)
    frac_atoms: List[Tuple[str, List[float]]] = []
    for sp, cart in atoms:
        frac = matvec(inv_cell_t, cart)
        if not has_final_state:
            frac = wrap_frac(frac)
        frac_atoms.append((sp, frac))

    # kpoint path nodes
    labelinfo = find_matching_labelinfo(base, win)
    nodes = parse_labelinfo(labelinfo)
    if not nodes:
        nodes = parse_kpoint_path(text)
    if not nodes:
        # fallback: Gamma to X
        nodes = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]

    if args.spincov is not None:
        spincov = args.spincov
    else:
        spincov = 1 if ispin == 2 else (infer_spincov_from_final_state_rows(final_state_rows) or 1)
    kmesh = parse_kmesh_from_wout(wout) or infer_kmesh_from_labelinfo(labelinfo) or 10

    # spinpol follows inferred ISPIN from win naming.
    hr_up, hr_dn = find_spinpol_hr_files(base)
    spinpol = ispin == 2
    if spinpol and (hr_up is None or hr_dn is None):
        if (base / "wannier90.1.win").exists():
            hr_up = base / "wannier90.1_hr.dat"
            hr_dn = base / "wannier90.2_hr.dat"
        else:
            hr_up = base / "wannier90.up_hr.dat"
            hr_dn = base / "wannier90.dn_hr.dat"
    hr_name = find_nonspin_hr_file(base)

    out = Path("tbbox.in")
    with out.open("w", encoding="utf-8") as f:
        f.write(f"  spinpol = {'True' if spinpol else 'False'}\n")
        f.write(f"  soc = {'T' if soc else 'F'}\n")
        if spinpol:
            f.write(f"  hr_up_name = {hr_up.name}\n")
            f.write(f"  hr_dn_name = {hr_dn.name}\n")
        f.write(f"  hr_name = {hr_name.name}\n")
        f.write("\n")
        f.write(" proj:\n")
        f.write(f" orbt = {args.orbt}\n")
        f.write(f" spincov = {spincov}\n")
        f.write(f" ntau = {len(frac_atoms)}\n")
        for sp, frac in frac_atoms:
            t = type_index.get(sp, 1)
            norb = norb_by_species.get(sp, 1)
            f.write(
                f"  {frac[0]:.15f}  {frac[1]:.15f}  {frac[2]:.15f}  0  0  0  {t} {norb}\n"
            )
        f.write(" end projections\n\n")
        f.write(" kpoint:\n")
        f.write(f" kmesh = {kmesh}\n")
        f.write(f" Nk = {len(nodes)}\n")
        for kx, ky, kz in nodes:
            f.write(f"  {kx:.10f}  {ky:.10f}  {kz:.10f}\n")
        f.write(" end kpoint_path\n\n")
        f.write(" unit_cell:\n")
        for row in cell:
            f.write(f"    {row[0]:.15f}  {row[1]:.15f}  {row[2]:.15f}\n")
        f.write(" end unit_cell_cart\n")

    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


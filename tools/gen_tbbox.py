#!/usr/bin/env python3
"""Generate tbbox.in from Wannier90 inputs/outputs in the current directory.

Usage:
  python tools/gen_tbbox.py --kmesh 10
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, List, Tuple


def read_text(path: Path) -> str:
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
    for ln in path.read_text(encoding="utf-8", errors="ignore").splitlines():
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
    for ln in path.read_text(encoding="utf-8", errors="ignore").splitlines():
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
    text = path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"Divisions along first K-path section\\s*:\\s*(\\d+)", text)
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate tbbox.in from Wannier90 files.")
    ap.add_argument("--orbt", type=int, default=2, help="orbt value (default: 2).")
    ap.add_argument("--spincov", type=int, default=1, help="spincov value (default: 1).")
    args = ap.parse_args()

    base = Path(".")
    win = base / "wannier90.win"
    if not win.exists():
        raise SystemExit(f"Missing {win}")

    text = read_text(win)
    soc = False
    cell = parse_unit_cell_cart(text)
    atoms = parse_atoms_cart(text)
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
        frac = wrap_frac(frac)
        frac_atoms.append((sp, frac))

    # kpoint path nodes
    labelinfo = base / "wannier90_band.labelinfo.dat"
    nodes = parse_labelinfo(labelinfo)
    if not nodes:
        nodes = parse_kpoint_path(text)
    if not nodes:
        # fallback: Gamma to X
        nodes = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]

    kmesh = parse_kmesh_from_wout(base / "wannier90.wout") or infer_kmesh_from_labelinfo(labelinfo) or 10

    # spinpol: check for up/dn HR files
    hr_up = next((p for p in [base / "wannier90.up_hr.dat"] if p.exists()), None)
    hr_dn = next((p for p in [base / "wannier90.dn_hr.dat"] if p.exists()), None)
    spinpol = bool(hr_up and hr_dn)
    hr_name = base / "wannier90_hr.dat"

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
        f.write(f" spincov = {args.spincov}\n")
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

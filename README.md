# IRSSG

<p align="right">
  English | <a href="README.zh-CN.md">中文</a>
</p>

`IRSSG` is a command-line program for spin-space-group workflows in magnetic materials without spin-orbit coupling (SOC). It identifies spin space group (SSG) operations from magnetic structures, generates the `ssg.data` file used by later steps, and assigns band coirreps from either VASP plane-wave results or Wannier tight-binding Hamiltonians.

For a complete tutorial and detailed examples, please read the user manual first:

- English manual: [irssg_manual.pdf](irssg_manual.pdf)
- Chinese manual: [irssg_manual_zh.pdf](irssg_manual_zh.pdf)

## Citation

If you use `IRSSG` in your research, please cite:

Sheng Zhang, Ziyin Song, Zhong Fang, Hongming Weng, and Zhijun Wang, "IRSSG: An Open-Source Software Package for Spin Space Groups," *Computer Physics Communications*, 110190 (2026). [https://doi.org/10.1016/j.cpc.2026.110190](https://doi.org/10.1016/j.cpc.2026.110190)

## Installation

`IRSSG` requires Python 3.8-3.13. Install it from PyPI:

```bash
pip install irssg
```

Required Python packages:

- `numpy >= 1.20.1`
- `scipy >= 1.5.4`
- `spglib >= 1.15.0`
- `pymatgen >= 2021.2.8`
- `phonopy >= 2.7.0`

If dependency installation fails, install common numerical packages first and then reinstall `irssg`:

```bash
conda install numpy h5py contourpy pandas
pip install irssg
```

## Quick Commands

```bash
# Identify the spin space group and generate ssg.data
irssg -ssg -c POSCAR > ssg.out

# Assign band coirreps from VASP plane-wave output
irssg -pw > irssg.out

# Assign band coirreps from a Wannier tight-binding model
irssg -wann > irssg.out
```

If `ssg.data` is missing when running `-pw` or `-wann`, `IRSSG` automatically calls the `-ssg` workflow first.

## Input Files

### `irssg -ssg`

Use one of the following structure files:

- `POSCAR`
- `*.mcif`

For `POSCAR`, append the Cartesian magnetic moment `(mx, my, mz)` after the coordinate columns of each magnetic atom. A line with a nonzero moment contains six numbers. If the moment is zero, the moment columns may be omitted.

Example:

```text
Mn3Sn
1.0
 5.665000   0.000000   0.000000
-2.832500   4.906034   0.000000
 0.000000   0.000000   4.531000
Mn   Sn
6    2
Direct
 0.838800   0.677600   0.250000   1.5  2.5981 0
 0.161200   0.838800   0.750000  -3    0      0
 0.838800   0.161200   0.250000  -3    0      0
 0.161200   0.322400   0.750000   1.5  2.5981 0
 0.322400   0.161200   0.250000   1.5 -2.5981 0
 0.677600   0.838800   0.750000   1.5 -2.5981 0
 0.333333   0.666667   0.250000
 0.666667   0.333333   0.750000
```

### `irssg -pw`

Required files:

- `OUTCAR`
- `WAVECAR`
- `ssg.data`

For VASP noncollinear calculations without spin-orbit coupling:

```text
LNONCOLLINEAR = .TRUE.
LSORBIT = .FALSE.
```

For collinear calculations:

```text
LNONCOLLINEAR = .FALSE.
ISPIN = 2
```

### `irssg -wann`

Required files:

- `case_hr.dat`, or another Wannier90-format `*_hr.dat` Hamiltonian file
- `tbbox.in`
- `ssg.data`

Minimal `tbbox.in` layout:

```text
spinpol = False
hr_name = symm_hr.dat

proj:
orbt = 2
spincov = 1
ntau = 8
  x1 x2 x3 itau iorbit
  ...
end projections

kpoint:
kmesh = 10
Nk = 2
  0.0 0.0 0.0
  0.5 0.0 0.0
end kpoint

unit_cell:
  a1x a1y a1z
  a2x a2y a2z
  a3x a3y a3z
end unit
```

For separate up/down Wannier Hamiltonians:

```text
spinpol = True
hr_name_up = sp.up_hr.dat
hr_name_dn = sp.dn_hr.dat
```

## Options

### SSG options

```bash
irssg -ssg [SSG-related options]
```

| Option | Meaning |
| --- | --- |
| `-c FILE` | Input structure file. Default: `POSCAR`. `*.mcif` is also supported. |
| `--tolerance TOL` | Tolerance for spatial symmetry operations. Default: `1.0e-3`. |
| `--magtolerance MTOL` | Tolerance for spin/magnetic moments. Default: `1.0e-4`. |
| `--standardize` | Standardize the structure after SSG identification. |

### Band options

The `-pw` and `-wann` modes support the same band-selection options:

```bash
irssg -pw [PW-related options]
irssg -wann [Wannier-related options]
```

| Option | Meaning |
| --- | --- |
| `-nk k_start k_end` | Process only a subset of k points. |
| `-nb band_min band_max` | Restrict the band window for coirrep assignment. |
| `-tolE dE` | Upper bound on the energy difference used to group near-degenerate bands and construct coirreps. Default: `1.0e-3`. |

## Output Files

| File | Produced by | Meaning |
| --- | --- | --- |
| `ssg.out` | `irssg -ssg > ssg.out` | Screen output containing SSG type, SSG number, spin-only group `S0`, subgroup `G0`, SSG international notation, operations, and related information. |
| `POSCAR.symm` | `irssg -ssg` | Symmetrized structure generated from the input structure and identified SSG operations. |
| `ssg.data` | `irssg -ssg` | Binary file storing SSG operations and related data, used by `-pw` and `-wann`. |
| `irssg.out` | `irssg -pw > irssg.out`, `irssg -wann > irssg.out` | Screen output containing little-group information, character tables, band characters, and band coirrep labels. |
| `chart.dat` | `irssg -pw`, `irssg -wann` | High-symmetry k-point/line/plane information, little-group character tables, and compatibility relations. |
| `fort.180` | `irssg -pw`, `irssg -wann` | Characters of band representations at each k point under SSG operations. |

## POS2SSG

`POS2SSG` is a small auxiliary tool, similar in spirit to `POS2MSG`, for generating highest-symmetry magnetic configurations by assigning magnetic moments to magnetic atoms. Each generated configuration corresponds to a specific SSG number.

## Default Behavior

Running `irssg` without a mode option first runs the `-ssg` workflow. If both `OUTCAR` and `WAVECAR` are present in the current directory, it then continues with the `-pw` workflow.

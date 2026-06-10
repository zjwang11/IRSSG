# IRSSG SSG Output Reference

Use this reference when parsing or explaining `irssg -ssg` / `MOM2SSG` output, or when deciding whether the downstream `-pw` / `-wann` band-corep modes are valid for the files in the working directory.

## Primary Command

```bash
irssg -ssg -c INPUT_FILE --tolerance 1e-3 --magtolerance 1e-4 > ssg.out
```

The bundled wrapper runs the same stage, writes `ssg_summary.json`, and then fetches HSP API records into `hsp_group_info.json` unless `--skip-hsp-api` is passed.

Do not use bare `irssg` for identification-only work. Per the IRSSG manual, bare `irssg` first behaves like `irssg -ssg`, then automatically enters `-pw` if `OUTCAR` and `WAVECAR` exist.

## Installation Notes

IRSSG supports Python 3.8-3.13. Check whether an installed package is already usable:

```bash
command -v irssg
python -c "import irssg.ssg.MOM2SSG; print('irssg ok')"
```

Install from PyPI when no local checkout is required:

```bash
python -m pip install irssg
```

If dependency or compilation errors occur, preinstall common numerical packages, then retry:

```bash
conda install numpy h5py contourpy pandas
python -m pip install irssg
```

For the local IRSSG checkout used in this environment, install the packaged publish tree:

```bash
python -m pip install /data/home/szhang/soft_sz/irssg_0702/publish
```

If the active Python environment is not writable:

```bash
python -m pip install --user irssg
export PATH="$HOME/.local/bin:$PATH"
```

To force the wrapper to use a local package source:

```bash
python skills/irssg-hsp-symmetry/scripts/run_irssg_ssg.py INPUT_FILE \
  --irssg-install-source /path/to/irssg/publish \
  --output-dir ssg-analysis --overwrite
```

## Key Fields

- `The SSG number`: canonical SSG identifier, e.g. a number with suffix indicating magnetic class.
- `The SSG international symbol`: formatted SSG label for reporting.
- `I/II/III ... So`: magnetic configuration class and spin-only group.
- `P (spin part of Go)`: spin point group part of `Go`.
- `H (lattice part of Go)`: lattice space-group part of `Go`.
- `# Number` under `Spin space group operations`: number of SSG operations.
- `Atomic space group`: nonmagnetic atomic SG detected from positions only.
- `N_ASG/N_SSG`: ratio between atomic space-group and spin-space-group operation counts in the reported convention.
- `The MSG number`: OG/BNS/SSG setting identifiers when MSG mapping is available.
- `The MSG international symbol`: MSG label.
- `# Number` under `Magnetic space group operations`: number of MSG operations.
- `hsp_group_identifiers` in `ssg_summary.json`: parsed identifiers used for API lookup:
  - `space_group_number`
  - `magnetic_group_og_number`
  - `magnetic_group_bns_number`
  - `spin_group_number`
- `hsp_api_status` in `ssg_summary.json`: `ok`, `partial_or_failed`, or `skipped_no_group_identifiers`.

## Generated Artifacts

- `ssg.out`: full stdout log; preserve for provenance and operation tables.
- `ssg.data`: IRSSG handoff file for later `irssg -pw` or `irssg -wann` character/corep analysis.
- `msg.data`: MSG-mode handoff file when generated.
- `ssgop.npy`: NumPy pickle of symmetry operation data.
- `POSCAR.symm`: symmetrized input structure with magnetic data.
- `POSCAR.ssg_primitive`: SSG primitive cell, written when the input POSCAR cell is not SSG primitive.
- `ssg_summary.json`: wrapper-generated parsed summary.
- `hsp_group_info.json`: wrapper-generated API payload containing complete HSP records:
  - `space_group`: full `/space-groups/{number}` JSON, including operations, Wyckoff positions, and high-symmetry k-vectors.
  - `magnetic_group`: full `/magnetic-groups/{og_or_bns}` JSON, including MSG operations, Wyckoff positions, and high-symmetry k-vectors.
  - `spin_group`: full `/spin-groups/{ssg_number}` JSON, including spin operations, Wyckoff positions, and high-symmetry k-vectors.
  - `requests`: URL/status/error metadata for each API request.
  - `web_links`: frontend URLs for the matching SG/MSG/SSG pages.

## Band/Corep Modes

Use these only after SSG identification has generated `ssg.data`, and only when the user is asking for band characters, irreducible projective representations, or coirreps from electronic-structure data.

`irssg -pw` is the VASP plane-wave workflow:

```bash
irssg -pw -nk 5 20 -nb 10 30 -tolE 1.0e-4 > irssg.out
```

Required files are `OUTCAR`, `WAVECAR`, and `ssg.data` in the current directory. The VASP calculation should be without SOC. For noncollinear workflows use `LNONCOLLINEAR=.TRUE.` and `LSORBIT=.FALSE.`. For collinear spin-separated workflows, `LNONCOLLINEAR=.FALSE.` and `ISPIN=2` can be used; if a k-little-group operation flips spin, IRSSG reports `spin = up & down`, otherwise it may output separate `spin = up` and `spin = down` blocks.

`irssg -wann` is the Wannier tight-binding workflow:

```bash
irssg -wann -nb 10 30 > irssg.out
```

Required files are `tbbox.in`, the referenced Wannier90-format `*_hr.dat` file(s), and `ssg.data`. In `tbbox.in`, `spinpol = False` uses one `hr_name`, while `spinpol = True` requires `hr_name_up` and `hr_name_dn`. The file also defines `unit_cell`, `kpoint`, and `proj` blocks used by IRSSG.

Both `-pw` and `-wann` support:

- `-nk k_start k_end`: 1-based k-point range.
- `-nb band_min band_max`: band range.
- `-tolE tol`: energy tolerance for degeneracy/coirrep assignment.

Both modes write:

- `irssg.out`: SSG operations, k-point names/coordinates, unitary and anti-unitary little-group elements, character tables, phases, coirreps, band characters, and band representation labels.
- `chart.dat`: high-symmetry point/line/plane names, little-group elements, character tables, and compatibility relations.
- `fort.180`: band wavefunction characters under SSG operations.

## Reporting Rules

Report values from `ssg_summary.json` first. If a field is missing but present in `ssg.out`, quote or paraphrase the exact relevant line and note that the wrapper did not parse it.

For complete symmetry data, report from `hsp_group_info.json`; do not reconstruct operations, Wyckoff positions, or high-symmetry k-vectors by hand from `ssg.out` unless the API record is missing.

For high-symmetry k-points, read the existing API payload fields:

- `space_group.k_points`
- `magnetic_group.k_points`
- `spin_group.k_points`

For each k-point, report `multiplicity`, `little_cogroup`, `coordinates_primitive`, and `coordinates_conventional` exactly as provided. `multiplicity` follows the HSP/API cell convention and is not required to equal the number of primitive coordinate strings in `coordinates_primitive`.

Do not add a new wrapper output file for each requested subgroup of fields; the wrapper's contract is to fetch complete HSP records once, then the agent extracts the requested fields.

Do not report an SSG number when IRSSG prints `The SG number` only; that means the moments are all below `--magtolerance` or missing.

Treat `ssg.data` as the required checkpoint before band-character or corep workflows. Do not run `-pw` or `-wann` until this file exists in the working directory alongside the required electronic-structure inputs.

When the user asks only for HSP group information, high-symmetry k-points, Wyckoff positions, operations, or database representation payloads, do not run `-pw` or `-wann`; use `hsp_group_info.json` and HSP API records instead.

## Common Failure Checks

- Input path: ensure the wrapper copied the intended mPOSCAR or `.mcif` into the output directory.
- Magnetic moments: ensure POSCAR rows include moment triplets and moments exceed `--magtolerance`.
- Tolerances: try looser `--tolerance` for slightly distorted coordinates; try looser `--magtolerance` for noisy moments.
- Source checkout imports: pass `--repo-root /path/to/irssg` only when explicitly testing an IRSSG source checkout. The default workflow must not scan parent directories for `src_mom2ssg/MOM2SSG.py`.
- Installed-package imports: run from an isolated output directory so a local source-tree `irssg/` folder does not shadow the installed `irssg.ssg` package.
- Bare command behavior: if an agent ran `irssg` without `-ssg`, re-run with explicit `-ssg` unless the user intentionally requested `-pw` auto-entry.
- Plane-wave inputs: `-pw` requires `OUTCAR`, `WAVECAR`, and `ssg.data`; missing any of these means the workflow is not ready.
- Wannier inputs: `-wann` requires `tbbox.in`, valid `hr_name`/`hr_name_up`/`hr_name_dn` files, and `ssg.data`.
- HSP API auth: set `HSP_API_KEY` or pass `--hsp-api-key-file`. The default header is `Authorization: Bearer <key>`.
- HSP API base URL: default is `https://cmpdc.iphy.ac.cn/hsp/api/v1`; override with `--hsp-api-base-url` for local development or a different deployment.
- API fetch failures are recorded in `hsp_group_info.json.requests`. They are nonfatal unless `--require-hsp-api` is used.

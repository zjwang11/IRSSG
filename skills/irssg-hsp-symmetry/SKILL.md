---
name: irssg-hsp-symmetry
description: Determine material space, magnetic, and spin space group information from magnetic POSCAR-like mPOSCAR files or magnetic CIF/mcif files using IRSSG, then fetch full related HSP web/API records for the identified SG, MSG, and SSG. Use when Codex needs to run or guide IRSSG-based symmetry identification, extract SG/MSG/SSG numbers and symbols, output complete symmetry operations, Wyckoff positions, high-symmetry k-vectors, generated ssg.data/msg.data artifacts, or prepare downstream spin-group workflows from mPOSCAR or mcif input.
---

# HSP Symmetry Skill

## Purpose

Use IRSSG to identify a material's SG/MSG/SSG from an mPOSCAR/POSCAR-like file or `.mcif`, then fetch complete Huairou Symmetry Platform (HSP) records for those groups. Treat IRSSG as the source of truth for group identification and HSP API records as the source of truth for operations, Wyckoff positions, high-symmetry k-vectors, and database-side representation payloads.

## Inputs

Accept either:

- mPOSCAR/POSCAR-like files whose magnetic atom rows include Cartesian magnetic moments after fractional/Cartesian positions, typically `x y z mx my mz`. Rows for atoms with zero moment may contain only `x y z`.
- `.mcif` files containing magnetic structure data readable by IRSSG.

Confirm that magnetic moments are present and nonzero. If IRSSG reports only a nonmagnetic SG, report that no SSG was identified because the input is nonmagnetic under the selected `--magtolerance`.

## Core Rules

- Use the bundled wrapper for normal material analysis. It runs explicit `irssg -ssg`, parses SG/MSG/SSG identifiers, and fetches HSP API records.
- Do not run bare `irssg`. The IRSSG manual says bare `irssg` first runs `-ssg`, then automatically enters `-pw` if `OUTCAR` and `WAVECAR` exist.
- Do not search or import from an IRSSG source tree unless the user explicitly asks to debug an IRSSG checkout. Use `--repo-root` only for that case.
- Do not invent or recompute HSP k-points, Wyckoff positions, operations, or representation payloads when `hsp_group_info.json` contains the API record.
- Run `irssg -pw` or `irssg -wann` only when the user asks for band characters/coreps from electronic-structure data and the required files exist.

## Main Workflow

1. Choose tolerances. Start with `--tolerance 1e-3` and `--magtolerance 1e-4` unless the user provides values or the structure is numerically noisy.
2. Provide HSP auth through `HSP_API_KEY`, `HSP_WEB_API_KEY`, `CMPDC_API_KEY`, `--hsp-api-key`, or `--hsp-api-key-file` when needed.
3. Run the wrapper from the skill directory:

```bash
python skills/irssg-hsp-symmetry/scripts/run_irssg_ssg.py INPUT_FILE --output-dir ssg-analysis --overwrite
```

4. Inspect `ssg-analysis/ssg_summary.json` first, then `ssg-analysis/hsp_group_info.json`, then `ssg-analysis/ssg.out` for debugging or full IRSSG operation tables.
5. Report identifiers, symbols, API status, generated files, and full-record paths. For high-symmetry k-points, read `space_group.k_points`, `magnetic_group.k_points`, and `spin_group.k_points` from `hsp_group_info.json`.

The wrapper auto-installs IRSSG if neither the `irssg` command nor `irssg.ssg.MOM2SSG` is available. It prefers the local package root when this skill lives under the IRSSG publish tree, otherwise it installs `irssg` from pip. For a specific local package, pass `--irssg-install-source /path/to/irssg/publish` or set `IRSSG_INSTALL_SOURCE`.

If the environment must not auto-install, pass `--no-auto-install-irssg`. If auto-install fails, follow the installation notes in `references/irssg-output.md`.

Use `--require-hsp-api` when missing HSP records should make the wrapper fail. Use `--skip-hsp-api` only when the user wants IRSSG identification without HSP records.

## Output Contract

- `ssg_summary.json`: parsed IRSSG identifiers, symbols, operation counts, generated files, attempted commands, and HSP API status.
- `hsp_group_info.json`: complete fetched HSP records plus request metadata and frontend links.
- `ssg.out` and `ssg.err`: raw IRSSG stdout/stderr for provenance and triage.
- `ssg.data`: binary handoff for later `irssg -pw` or `irssg -wann`.
- `POSCAR.symm`, `POSCAR.ssg_primitive`, `msg.data`, and `ssgop.npy`: IRSSG-generated artifacts when produced.

Prefer citing identifiers and run status from `ssg_summary.json`, and detailed group data from `hsp_group_info.json`. If parsing misses a field but `ssg.out` contains it, cite the relevant line and mark the JSON field as unavailable instead of guessing.

For high-symmetry k-points, report `multiplicity`, `little_cogroup`, and coordinate lists exactly as HSP returns them. Do not compare `multiplicity` against the number of primitive coordinate strings as a validity check.

## HSP API

The wrapper fetches these records after successful IRSSG identification unless `--skip-hsp-api` is passed:

- Space group: `GET /space-groups/{number}` using the leading number from `Atomic space group`.
- Magnetic group: `GET /magnetic-groups/{og_or_bns}` using OG first from `The MSG number`, then BNS if OG is absent.
- Spin group: `GET /spin-groups/{ssg_number}` using `The SSG number`.

Defaults:

- API base: `https://cmpdc.iphy.ac.cn/hsp/api/v1`
- Web links: `https://cmpdc.iphy.ac.cn/hsp`
- API key header: `Authorization: Bearer <key>`

Use `--hsp-api-key-header`, `--hsp-api-key-prefix`, or `--hsp-api-key-query-param` only when the deployment expects a different API-key convention.

## Direct IRSSG Commands

Run direct commands only when the wrapper is unnecessary or the user explicitly asks for IRSSG-native output:

```bash
irssg -ssg -c INPUT_FILE --tolerance 1e-3 --magtolerance 1e-4 > ssg.out
```

For band/corep workflows, first generate or preserve `ssg.data`, then run the appropriate mode only in a directory with the required electronic-structure files:

```bash
irssg -pw -nk k_start k_end -nb band_min band_max -tolE 1.0e-4 > irssg.out
irssg -wann -nb band_min band_max > irssg.out
```

## Reference

Read `references/irssg-output.md` when you need exact IRSSG field meanings, installation triage, generated artifact meanings, `-pw` / `-wann` details, or failure checks.

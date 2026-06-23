---
name: hsp-symmetry
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
python skills/hsp-symmetry/scripts/run_irssg_ssg.py INPUT_FILE --output-dir ssg-analysis --overwrite
```

4. Inspect `ssg-analysis/ssg_summary.json` first, then `ssg-analysis/hsp_group_info.json`, then `ssg-analysis/ssg.out` for debugging or full IRSSG operation tables.
5. Report identifiers, symbols, API status, generated files, and full-record paths. For high-symmetry k-points, read `space_group.k_points`, `magnetic_group.k_points`, and `spin_group.k_points` from `hsp_group_info.json`.

The wrapper auto-installs IRSSG if neither the `irssg` command nor `irssg.ssg.MOM2SSG` is available. It prefers the local package root when this skill lives under the IRSSG publish tree, otherwise it installs `irssg` from pip. For a specific local package, pass `--irssg-install-source /path/to/irssg/publish` or set `IRSSG_INSTALL_SOURCE`.

If the environment must not auto-install, pass `--no-auto-install-irssg`. If auto-install fails, follow the installation notes in `references/irssg-output.md`.

Use `--require-hsp-api` when missing HSP records should make the wrapper fail. Use `--skip-hsp-api` only when the user wants IRSSG identification without HSP records.

## Output Contract

Strictly follow this output contract in the final answer. Do not replace requested
tables with a file path, summary, count, or "see attached/generated file" response.
Generated files may be cited only as provenance or as an additional artifact after
the requested data has been printed directly in the final answer.

- `ssg_summary.json`: parsed IRSSG identifiers, symbols, operation counts, generated files, attempted commands, and HSP API status.
- `hsp_group_info.json`: complete fetched HSP records plus request metadata and frontend links.
- `ssg.out` and `ssg.err`: raw IRSSG stdout/stderr for provenance and triage.
- `ssg.data`: binary handoff for later `irssg -pw` or `irssg -wann`.
- `POSCAR.symm`, `POSCAR.ssg_primitive`, `msg.data`, and `ssgop.npy`: IRSSG-generated artifacts when produced.

Prefer citing identifiers and run status from `ssg_summary.json`, and detailed group data from `hsp_group_info.json`. If parsing misses a field but `ssg.out` contains it, cite the relevant line and mark the JSON field as unavailable instead of guessing.

When the user asks to output, list, report, summarize, or return identifiers, group data, operations, Wyckoff positions, or k-points, present the requested structured data as Markdown tables directly in the final answer. High-symmetry k-point output must be tabular, including every coordinate string from the relevant HSP `k_points` record. When a k-point has multiple coordinate strings, use one table row per coordinate and add a coordinate index column, instead of combining coordinates in one cell. Never answer a k-point coordinate request with only labels, counts, a summary, or a path to a generated file.

Wyckoff-position output must include coordinates by default. For `space_group.wyckoff_positions`, `magnetic_group.wyckoff_positions`, and `spin_group.wyckoff_positions`, report `label`, `multiplicity`, `site_symmetry` when present, coordinate index, and each coordinate string from HSP. When a Wyckoff position has multiple coordinate strings, use one table row per coordinate and add a coordinate index column, instead of combining coordinates in one cell. Do not answer a Wyckoff-position request with only a summary of labels, multiplicities, or site symmetries unless the user explicitly asks for a summary-only response.

For space-group and magnetic-group high-symmetry k-point tables, do not output HSP's internal `letter` field as a user-visible column. Use the HSP `label` field as the k-point name, then report `multiplicity`, `little_cogroup` when present, coordinate index, and coordinates. Apply the same omission to spin-group output unless the user explicitly asks for the internal letter field.

Do not put HTML tags or renderer-specific markup such as `<br>` in final output. Render coordinate values in the user-visible `(k | spin)` or `(position | spin)` form, including outer parentheses; if HSP returns a coordinate string without outer parentheses, add them for display without changing the internal values. In pipe-delimited Markdown tables, encode only the coordinate separator as `&#124;` in the Markdown source so it renders visibly as `|` without splitting the table. Do not show backslash-escaped pipes such as `\|`. Do not use inline code to hide this problem. The rendered coordinate must look exactly like `(1/2, 1/2, w | 0, 0, 0)`.

For all user-visible coordinate strings, including Wyckoff-position coordinates and high-symmetry k-point coordinates, render standard numeric coefficients in exact symbolic form even when the coefficient is attached to a variable or magnetic moment. Do not display `0.866025m1`, `-0.866025mx`, `0.5x`, or `-0.5my`; display them as `(sqrt(3)/2)m1`, `-(sqrt(3)/2)mx`, `(1/2)x`, and `-(1/2)my`. Apply the same rule to standalone entries and algebraic expressions while otherwise preserving the HSP coordinate expression.

Before producing the final answer, omit any displayed column whose values would all be `null`, `None`, empty, or unavailable for the displayed rows. Never print a `null` column or the literal `null` in final user-facing output. If a field is absent only for some rows, omit that field for those rows instead of printing `null`.

For symmetry operation tables, preserve the HSP operation data but render user-visible numeric matrix and vector entries in exact symbolic form whenever they correspond to standard crystallographic values. Do not display rounded decimal approximations such as `0.5`, `-0.5`, `0.866025`, `-0.866025`, `0.8660254037844388`, or small floating-point noise like `6.123233995736766e-17`; display them as `1/2`, `-1/2`, `sqrt(3)/2`, `-sqrt(3)/2`, and `0` respectively. Apply this to space, spin, parent, translation, axis, SU(2), coordinate strings, and any other operation matrices/vectors shown to the user. If a value cannot be confidently mapped to an exact standard value, keep the HSP value and state that it is shown as returned by the API.

When reporting symmetry operations for a material run, output the complete operation set for each requested group. For SSG operations, use the IRSSG `Spin space group operations` table in `ssg.out` and include every operation up to `ssg_operation_count` / the table `# Number` value. For MSG operations, include every operation up to `msg_operation_count` / the MSG table `# Number` value. Treat `Indices of MSG operations within the list of SSG operations` only as a cross-reference between MSG and SSG rows; do not use it to filter or truncate the SSG table unless the user explicitly asks for only that subset. If `hsp_group_info.json` contains more spin-group operation rows than the IRSSG material SSG operation count, explain the distinction and prefer the IRSSG material operation table for the user's POSCAR-specific operation output.

For high-symmetry k-points, report `multiplicity`, `little_cogroup` when present, and coordinate lists from HSP without changing the coordinate values except for display-only exact-symbol rendering and outer-parenthesis normalization as described above. Do not compare `multiplicity` against the number of primitive coordinate strings as a validity check.

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

# IRSSG Unified Build

This repository contains a reorganized IRSSG project that merges all sources
into a single, unified directory structure and supports a one‑shot build.

## Directory Layout

```
unified_irssg/
├── lib/           # Library sources (originally lib_ssg)
├── src/           # Main program sources (originally src_irssg)
├── Makefile       # Unified build file
├── libIRSSG.a     # Generated static library
├── irssg          # Generated executable
└── README.md      # This document
```

## Build Instructions

1) Load the Intel OneAPI environment

   ```bash
   module load oneapi22.3
   ```

2) Enter the project directory

   ```bash
   cd unified_irssg
   ```

3) Build

   ```bash
   make
   ```

4) Clean build artifacts

   ```bash
   make clean
   ```

## Build Pipeline

The new Makefile automatically:
1. Compiles all library sources under `lib/`
2. Creates the static library `libIRSSG.a`
3. Compiles all main program sources under `src/`
4. Links and generates the final executable `irssg`

## Benefits

- Simplified workflow: previously two stages (lib_ssg → src_irssg), now a single step.
- Unified directory structure for easier maintenance.
- Original compiler flags and dependencies are preserved.

## Single Executable with Runtime Mode Selection

Two main programs are merged into a single executable `irssg`. Select the mode at runtime via flags:

- Default: plain PW path
- With `--wann` or `--mode wann`: use the Wannier path

Examples

```
./irssg -nk 1 10                      # PW path
./irssg --wann -nk 1 10               # Wannier path
./irssg --mode wann -tolE 1e-3        # Equivalent Wannier path
```

Notes

- All original flags remain the same. The new flags `--wann` / `--mode` are ignored by the legacy parser and do not affect argument parsing.
- In the source, `src_wann` modules `comms/init/wave_data` were renamed (with `_wann` suffix) to avoid name collisions with the PW path. Both flows are dispatched by a unified `main` through `driver` subroutines.

## MOM2SSG Integration (`-ssg`)

The Python package MOM2SSG is bundled in the installation. Usage:

- Run the Fortran main flow (default):
  - `irssg [original flags]`
- Call MOM2SSG (Python):
  - `irssg -ssg [MOM2SSG flags]`

Notes

- `-ssg` forwards subsequent arguments to `MOM2SSG.MOM2SSG:main`, equivalent to running `python -m MOM2SSG.MOM2SSG ...`.
- MOM2SSG depends on `numpy/scipy/spglib/pymatgen/phonopy`. Please ensure these dependencies are available before using `-ssg`.

## Local pip Install (Build from Source)

Local pip install is supported; the Fortran binary is compiled and the `irssg` command is installed automatically.

Requirements

- Intel ifort/ifx compiler and MKL (Makefile uses ifort with `-qmkl`)
- Python ≥ 3.8

Recommended for offline/HPC clusters: disable isolation to avoid network downloads during build.

```
python -m pip install . --no-build-isolation
```

Notes

- During build, `make USE_ABS_DATA_PATH=0` is invoked so the binary and `kLittleGroups` data are bundled into the Python package.
- At runtime, the `irssg` entry point sets `IRVSPDATA` to the package data directory automatically; no manual configuration is needed.
- To repackage a platform wheel locally: `python -m pip wheel . --no-build-isolation`. The resulting wheel includes the compiled binary.


# irssg_construction

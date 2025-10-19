#!/usr/bin/env bash
set -euo pipefail

# Clean all build and packaging artifacts to restore repo to pre-build state.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
cd "${repo_root}"

echo "Cleaning build and packaging artifacts in: ${repo_root}" >&2

# Use Makefile clean if present
if [[ -f Makefile ]]; then
  make clean || true
fi

# Remove Python packaging artifacts and wheels
rm -rf build dist dist_* *.egg-info irssg.egg-info irssg_pkg.egg-info || true

# Remove leftover binaries and archives at root
rm -f irssg libIRSSG.a *.o *.mod || true

# Remove library objects, modules, and archive
rm -f lib/*.o lib/*.mod lib/libIRSSG.a || true

# Remove PW and Wann objects/modules
rm -f src_pw/*.o src_pw/*.mod || true
rm -f src_wann/*.o src_wann/*.mod || true

# Remove Fortran module output directory
rm -rf build/mod || true

echo "Cleanup complete." >&2

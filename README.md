# POS2SSG

Convert POSCAR to Spin Space Group (SSG) with magnetic moments.

POS2SSG is a Python package that analyzes crystal structures and generates magnetic configurations corresponding to specific spin space groups (SSG). It is designed to work with VASP POSCAR files and can automatically identify magnetic atoms and generate appropriate magnetic moments based on spin space group symmetry.

## Features

- **Automatic SSG Detection**: Find possible spin space groups for a given crystal structure
- **Magnetic Moment Generation**: Generate magnetic moments based on SSG symmetry
- **Multiple Magnetic Atom Support**: Handle structures with multiple types of magnetic atoms
- **Custom Magnetic Configuration**: Specify custom magnetic atoms and moments
- **VASP Integration**: Works seamlessly with VASP POSCAR files
- **Flexible Output**: Generate SPOSCAR files with magnetic moments

### Dependencies

- Python >= 3.8
- numpy >= 1.20.0
- scipy >= 1.7.0
- pymatgen >= 2022.0.0
- spglib >= 1.16.0
- irssg >= 1.0.1

## Installation

Install from PyPI (recommended):

```bash
pip install pos2ssg
```

Install from source (editable for development):

```bash
cd POS2SSG
pip install -e .
```

## Usage

### Command Line Interface

```bash
# Find possible SSG numbers for a structure
pos2ssg -c POSCAR

# Find SSG numbers with specific Ik value
pos2ssg -c POSCAR -Ik 2

# Generate SPOSCAR for specific SSG
pos2ssg -c POSCAR -ssg 43.2.4.1.P

# Use custom magnetic atoms and moments
pos2ssg -c POSCAR -ssg 43.2.4.1.P --magatom Fe Cu --magmom 2.5 1.8

# Specify output file
pos2ssg -c POSCAR -ssg 43.2.4.1.P --output my_structure.vasp

# Enable verbose output
pos2ssg -c POSCAR -ssg 43.2.4.1.P --verbose

# Adjust tolerance for symmetry analysis
pos2ssg -c POSCAR -ssg 43.2.4.1.P --tolerance 1e-3
```

## Command Line Options

- `-c, --poscar`: Input POSCAR file (default: POSCAR)
- `-Ik`: Filter by specific Ik value
- `-ssg`: Target spin space group number (e.g., 43.2.4.1.P)
- `--tolerance`: Tolerance for symmetry analysis (default: 1e-2)
- `--magatom`: List of magnetic atom types (e.g., Fe Cu)
- `--magmom`: List of magnetic moments (e.g., 3.0 2.0)
- `--output, -o`: Output file name (default: SPOSCAR)
- `--verbose, -v`: Enable verbose output

### Basic Usage

```bash
# Find all possible SSG numbers
pos2ssg -c POSCAR

# Find SSG numbers for Ik=2
pos2ssg -c POSCAR -Ik 2
```

### Generate Magnetic Structure

```bash
# Generate SPOSCAR for SSG 43.2.4.1.P
pos2ssg -c POSCAR -ssg 43.2.4.1.P

# Use custom magnetic moments
pos2ssg -c POSCAR -ssg 43.2.4.1.P --magatom Fe Cu --magmom 2.5 1.8

# Save to custom file
pos2ssg -c POSCAR -ssg 43.2.4.1.P --output magnetic_structure.vasp
```

## Output

The program generates a SPOSCAR file containing:
- The original crystal structure
- Magnetic moments for all magnetic atoms
- Proper symmetry operations for the target SSG

## Theory

POS2SSG is based on the theory of spin space groups, which describe the symmetry of magnetic crystals. The program:

1. Analyzes the crystal structure to determine the space group
2. Identifies magnetic atoms and their Wyckoff positions
3. Applies spin space group symmetry operations to generate magnetic moments
4. Ensures the magnetic configuration is consistent with the target SSG

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.


## Support

For questions and support, please open an issue on GitHub or contact the development team.

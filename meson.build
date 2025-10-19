project('irssg-pkg', 'fortran',
  version : '0.1.1',
  meson_version : '>=0.64.0',
  default_options : [
    'buildtype=release'
  ]
)

# Fortran compiler and Python
fortran_compiler = meson.get_compiler('fortran')
py_mod = import('python')
py3 = py_mod.find_installation('python3', pure: false)
py3_dep = py3.dependency()

# Filesystem helper and package layout switch (support both irssg and legacy irssg_pkg)
fs = import('fs')
new_pkg = fs.exists('irssg/__init__.py')
pkg_name = 'irssg'
if not new_pkg
  pkg_name = 'irssg_pkg'
endif

# Build flags and macros for gfortran compatibility
add_project_arguments([
  '-cpp',
  '-ffree-line-length-none',
  '-std=legacy',
  '-DIRVSPDATA',
], language: 'fortran')

# Fortran sources
lib_sources = files(
  'lib/lib_params.f90',
  'lib/lib_comms.f90',
  'lib/mathlib.f90',
  'lib/lib_bilbao.f90',
  'lib/lib_chrct.f90',
  'lib/invreal33.f90',
  'lib/invmati.f90',
  'lib/kgroup.f90',
  'lib/irrep_ssg.f90',
  'lib/get_ssg.f90',
  'lib/comprel.f90',
  'lib/linear_rep.f90',
)

pw_sources = files(
  'src_pw/comms.f90',
  'src_pw/init.f90',
  'src_pw/wave_data_pw.f90',
  'src_pw/main.f90',
)

wann_sources = files(
  'src_wann/file_util.f90',
  'src_wann/comms.f90',
  'src_wann/init.f90',
  'src_wann/wave_data.f90',
  'src_wann/main.f90',
)

main_prog = files('src_pw/unified_main.f90')

# Link LAPACK/BLAS on Linux
link_args = []
if host_machine.system() == 'linux'
  link_args += ['-L/usr/lib64', '-llapack', '-lblas']
endif

# Build and install executable into the Python package bin dir
irssg_exe = executable('irssg',
  lib_sources + pw_sources + wann_sources + main_prog,
  link_args: link_args,
  install: true,
  install_dir: py3.get_install_dir() / pkg_name / 'bin',
)

# Install Python package sources under site-packages (CLI and __init__)
if new_pkg
  py3.install_sources(files('irssg/__init__.py', 'irssg/cli.py'), subdir: 'irssg')
else
  # Install legacy package
  py3.install_sources(files('irssg_pkg/__init__.py', 'irssg_pkg/cli.py'), subdir: 'irssg_pkg')
  # Also install an irssg shim that forwards to irssg_pkg
  py3.install_sources(files('irssg_shim/__init__.py', 'irssg_shim/cli.py'), subdir: 'irssg')
endif

# Install SSG Python modules into <pkg_name>/ssg
py3.install_sources(files(
  'src_mom2ssg/__init__.py',
  'src_mom2ssg/MOM2SSG.py',
  'src_mom2ssg/small_func.py',
  'src_mom2ssg/SG_utils.py',
  'src_mom2ssg/SSGLabel.py',
  'src_mom2ssg/SG_isomorphism.py',
  'src_mom2ssg/lirssg.py',
  'src_mom2ssg/poscar_io.py',
  'src_mom2ssg/load_ssgdata.py',
  'src_mom2ssg/wyckoff.py',
  'src_mom2ssg/find_ssg_operation.py',
  'src_mom2ssg/get_ssg_number.py',
  'src_mom2ssg/std2ssg.py',
  'src_mom2ssg/spintrans.py',
  'src_mom2ssg/eqvpg2label.py',
), subdir: pkg_name + '/ssg')

# Additional SSG Python modules required by std2ssg and others
py3.install_sources(files(
  'src_mom2ssg/pos2abr.py',
  'src_mom2ssg/find_magprim_unit.py',
  'src_mom2ssg/getcell.py',
), subdir: pkg_name + '/ssg')

# Install SSG data directories
install_subdir('src_mom2ssg/ssg_data', install_dir: py3.get_install_dir() / pkg_name / 'ssg' / 'ssg_data', strip_directory: true)
install_subdir('src_mom2ssg/smith_form', install_dir: py3.get_install_dir() / pkg_name / 'ssg' / 'smith_form', strip_directory: true)

# Install data directory
install_subdir('lib/kLittleGroups', install_dir: py3.get_install_dir() / pkg_name / 'data', strip_directory: false)

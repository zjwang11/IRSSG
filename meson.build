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

# Add tolerant Fortran flags and preprocessor macro
add_project_arguments([
  '-cpp',
  '-ffree-line-length-none',
  '-std=legacy',
  '-DIRVSPDATA',
], language: 'fortran')

# Detect if building under cibuildwheel
is_cibw = run_command('python', '-c', 'import os; print(os.environ.get("CIBUILDWHEEL", "0"))', check: false).stdout().strip() == '1'

# Collect Fortran sources
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

# Link arguments: under manylinux, link to system LAPACK/BLAS
link_args = []
if host_machine.system() == 'linux'
  # Prefer manual linking in CIBW environment
  link_args += ['-L/usr/lib64', '-llapack', '-lblas', '-Wl,-rpath,/usr/lib64']
endif

# Build executable
irssg_exe = executable('irssg',
  lib_sources + pw_sources + wann_sources + main_prog,
  link_args: link_args,
  install: true,
  install_dir: py3.get_install_dir() / 'irssg_pkg' / 'bin',
)

# Install Python package sources
install_subdir('irssg_pkg', install_dir: py3.get_install_dir())

# Install SSG Python module under irssg_pkg/ssg
install_subdir('src_mom2ssg', install_dir: py3.get_install_dir() / 'irssg_pkg' / 'ssg', strip_directory: true)

# Install data directory
install_subdir('lib/kLittleGroups', install_dir: py3.get_install_dir() / 'irssg_pkg' / 'data', strip_directory: false)

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
import subprocess
import shutil
import pathlib
import os


class build_py(_build_py):
    def run(self):
        # Build the Fortran binary and data via Makefile
        root = pathlib.Path(__file__).parent.resolve()
        env = os.environ.copy()

        # Prefer runtime-resolved data path in lib_bilbao (IRVSPDATA)
        subprocess.check_call(["make", "clean"], cwd=root, env=env)
        subprocess.check_call(["make", "USE_ABS_DATA_PATH=0"], cwd=root, env=env)

        # Prepare package layout in build_lib
        build_lib = pathlib.Path(self.build_lib)
        pkg_dir = build_lib / "irssg_pkg"
        (pkg_dir / "bin").mkdir(parents=True, exist_ok=True)
        (pkg_dir / "data").mkdir(parents=True, exist_ok=True)

        # Copy binary
        shutil.copy2(root / "irssg", pkg_dir / "bin" / "irssg")

        # Copy data directory (kLittleGroups)
        data_src = root / "lib" / "kLittleGroups"
        if data_src.exists():
            shutil.copytree(data_src, pkg_dir / "data" / "kLittleGroups", dirs_exist_ok=True)

        # Vendor MOM2SSG code under irssg_pkg/ssg so imports use the irssg namespace
        mom_pkg_src = root / "src_mom2ssg"
        mom_pkg_dst = pkg_dir / "ssg"
        if mom_pkg_src.exists():
            shutil.copytree(mom_pkg_src, mom_pkg_dst, dirs_exist_ok=True)

        super().run()


cmdclass = {"build_py": build_py}
try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

    class bdist_wheel(_bdist_wheel):
        def finalize_options(self):
            super().finalize_options()
            # Mark wheel as non-pure so it gets a platform tag
            self.root_is_pure = False

    cmdclass["bdist_wheel"] = bdist_wheel
except Exception:
    pass

setup(cmdclass=cmdclass)

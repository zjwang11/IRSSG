import os
import sys
import subprocess

try:
    from importlib.resources import files
except ImportError:  # Python <3.9 fallback
    from importlib_resources import files  # type: ignore


def _base_paths():
    base = files("irssg")
    data_dir = base.joinpath("data")
    bin_path = base.joinpath("bin", "irssg")
    return base, data_dir, bin_path


def _has_ssg_data():
    return os.path.exists("ssg.data")


def _require_file(path: str, err: str) -> int:
    if not os.path.exists(path):
        sys.stderr.write(err + "\n")
        return 1
    return 0


def _has_pw_inputs():
    return os.path.exists("OUTCAR") and os.path.exists("WAVECAR")


def _run_ssg(ssg_args):
    # Run integrated SSG Python tool under irssg namespace
    cmd = [sys.executable, "-m", "irssg.ssg.MOM2SSG"] + ssg_args
    return subprocess.call(cmd)


def _run_pw(bin_path, pw_args):
    cmd = [str(bin_path)] + pw_args
    env = os.environ.copy()
    if not env.get("IRSSG_KEEP_LD_LIBRARY_PATH"):
        env.pop("LD_LIBRARY_PATH", None)
        _debug("LD_LIBRARY_PATH removed for PW child process")
    return subprocess.call(cmd, env=env)


def _run_wann(bin_path, wann_args):
    cmd = [str(bin_path), "--wann"] + wann_args
    env = os.environ.copy()
    if not env.get("IRSSG_KEEP_LD_LIBRARY_PATH"):
        env.pop("LD_LIBRARY_PATH", None)
        _debug("LD_LIBRARY_PATH removed for WANN child process")
    return subprocess.call(cmd, env=env)

def _run_deplicate(extra_args=None):
    # Run post-processing script after successful Fortran run
    args = [] if extra_args is None else list(extra_args)
    cmd = [sys.executable, "-m", "irssg.deplicate_rep"] + args
    return subprocess.call(cmd)

def _debug(msg: str) -> None:
    if os.environ.get("IRSSG_DEBUG"):
        sys.stderr.write(f"[irssg-debug] {msg}\n")

def main() -> int:
    # Resolve packaged data path and binary
    _, data_dir, bin_path = _base_paths()
    _debug(f"cwd={os.getcwd()}")
    _debug(f"bin_path={bin_path}")
    _debug(f"data_dir={data_dir}")

    # Ensure IRVSPDATA points to directory containing kLittleGroups/
    os.environ.setdefault("IRVSPDATA", str(data_dir))

    # Parse and classify arguments by option names with optional mode tags
    argv = sys.argv[1:]
    modes = {"ssg": [], "pw": [], "wann": []}
    present = {"ssg": False, "pw": False, "wann": False}
    current = None

    # Known option specifications
    ssg_spec = {
        "-c": 1,
        "--standardize": 0,
        "--tolerance": 1,
        "--magtolerance": 1,
    }
    pw_spec = {
        "-nk": 2,
        "-nb": 2,
        "-tolE": 1,
    }

    def fallback_pw_mode():
        # If user explicitly chose -wann and not -pw, prefer wann bucket for PW-like opts
        return "wann" if (present["wann"] and not present["pw"]) else "pw"

    i = 0
    while i < len(argv):
        a = argv[i]
        # Mode tags
        if a in ("-ssg", "--ssg"):
            present["ssg"] = True
            current = "ssg"
            i += 1
            continue
        if a in ("-pw", "--pw"):
            present["pw"] = True
            current = "pw"
            i += 1
            continue
        if a in ("-wann", "--wann"):
            present["wann"] = True
            current = "wann"
            i += 1
            continue

        # Option classification
        if a.startswith("-"):
            # Known SSG option
            if a in ssg_spec:
                argc = ssg_spec[a]
                modes["ssg"].append(a)
                for _ in range(argc):
                    if i + 1 < len(argv):
                        modes["ssg"].append(argv[i + 1])
                        i += 1
                i += 1
                continue
            # Known PW/Wann option
            if a in pw_spec:
                argc = pw_spec[a]
                bucket = current if current in ("pw", "wann") else fallback_pw_mode()
                modes[bucket].append(a)
                for _ in range(argc):
                    if i + 1 < len(argv):
                        modes[bucket].append(argv[i + 1])
                        i += 1
                i += 1
                continue

            # Unknown option: route to current mode if set, else fallback PW; attach one value if present
            bucket = current if current in ("ssg", "pw", "wann") else fallback_pw_mode()
            modes[bucket].append(a)
            if i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                modes[bucket].append(argv[i + 1])
                i += 1
            i += 1
            continue

        # Positional argument: route to current mode if set, else fallback PW
        bucket = current if current in ("ssg", "pw", "wann") else fallback_pw_mode()
        modes[bucket].append(a)
        i += 1

    # No explicit mode flags: default to PW flow
    if not any(present.values()):
        # Ensure ssg.data exists; if not, run SSG with any provided SSG opts (none by default)
        _debug(f"has_ssg_data={_has_ssg_data()}")
        if not _has_ssg_data():
            _debug(f"run ssg: {sys.executable} -m irssg.ssg.MOM2SSG {' '.join(modes['ssg'])}")
            rc = _run_ssg(modes["ssg"])  # may be empty
            _debug(f"ssg rc={rc}")
            if rc != 0:
                return rc
        # Require WAVECAR and OUTCAR
        _debug(f"has_pw_inputs={_has_pw_inputs()}")
        if not _has_pw_inputs():
            sys.stderr.write("PW mode requires WAVECAR and OUTCAR in current directory.\n")
            return 1
        _debug(f"run pw: {bin_path} {' '.join(modes['pw'])}")
        rc = _run_pw(bin_path, modes["pw"])
        _debug(f"pw rc={rc}")
        if rc != 0:
            return rc
        return _run_deplicate([])

    # Only SSG
    if present["ssg"] and not present["pw"] and not present["wann"]:
        _debug(f"run ssg only: {sys.executable} -m irssg.ssg.MOM2SSG {' '.join(modes['ssg'])}")
        return _run_ssg(modes["ssg"])

    # If both -pw and -wann are specified, this is invalid
    if present["pw"] and present["wann"]:
        sys.stderr.write("Cannot specify both -pw and -wann in one run.\n")
        return 2

    # -ssg with -pw: run SSG first, then PW
    if present["ssg"] and present["pw"] and not present["wann"]:
        _debug(f"run ssg: {sys.executable} -m irssg.ssg.MOM2SSG {' '.join(modes['ssg'])}")
        rc = _run_ssg(modes["ssg"])  # run SSG unconditionally first
        if rc != 0:
            return rc
        # Require WAVECAR and OUTCAR
        _debug(f"has_pw_inputs={_has_pw_inputs()}")
        if not _has_pw_inputs():
            sys.stderr.write("PW mode requires WAVECAR and OUTCAR in current directory.\n")
            return 1
        _debug(f"run pw: {bin_path} {' '.join(modes['pw'])}")
        rc = _run_pw(bin_path, modes["pw"])
        _debug(f"pw rc={rc}")
        if rc != 0:
            return rc
        return _run_deplicate([])

    # -ssg with -wann: run SSG first, then Wann
    if present["ssg"] and present["wann"] and not present["pw"]:
        _debug(f"run ssg: {sys.executable} -m irssg.ssg.MOM2SSG {' '.join(modes['ssg'])}")
        rc = _run_ssg(modes["ssg"])  # run SSG unconditionally first
        if rc != 0:
            return rc
        # Require tbbox.in
        if _require_file("tbbox.in", "WANN mode requires tbbox.in in current directory.") != 0:
            return 1
        _debug(f"run wann: {bin_path} --wann {' '.join(modes['wann'])}")
        rc = _run_wann(bin_path, modes["wann"])
        _debug(f"wann rc={rc}")
        if rc != 0:
            return rc
        return _run_deplicate([])

    # PW flow
    if present["pw"] and not present["wann"]:
        # Ensure ssg.data exists (run SSG first if missing)
        _debug(f"has_ssg_data={_has_ssg_data()}")
        if not _has_ssg_data():
            _debug(f"run ssg: {sys.executable} -m irssg.ssg.MOM2SSG {' '.join(modes['ssg'])}")
            rc = _run_ssg(modes["ssg"])  # may be empty
            _debug(f"ssg rc={rc}")
            if rc != 0:
                return rc
        # Require WAVECAR and OUTCAR
        _debug(f"has_pw_inputs={_has_pw_inputs()}")
        if not _has_pw_inputs():
            sys.stderr.write("PW mode requires WAVECAR and OUTCAR in current directory.\n")
            return 1
        _debug(f"run pw: {bin_path} {' '.join(modes['pw'])}")
        rc = _run_pw(bin_path, modes["pw"])
        _debug(f"pw rc={rc}")
        if rc != 0:
            return rc
        return _run_deplicate([])

    # WANN flow
    if present["wann"] and not present["pw"]:
        # Ensure ssg.data exists (run SSG first if missing)
        _debug(f"has_ssg_data={_has_ssg_data()}")
        if not _has_ssg_data():
            _debug(f"run ssg: {sys.executable} -m irssg.ssg.MOM2SSG {' '.join(modes['ssg'])}")
            rc = _run_ssg(modes["ssg"])  # may be empty
            _debug(f"ssg rc={rc}")
            if rc != 0:
                return rc
        # Require tbbox.in
        if _require_file("tbbox.in", "WANN mode requires tbbox.in in current directory.") != 0:
            return 1
        _debug(f"run wann: {bin_path} --wann {' '.join(modes['wann'])}")
        rc = _run_wann(bin_path, modes["wann"])
        _debug(f"wann rc={rc}")
        if rc != 0:
            return rc
        return _run_deplicate([])

    # If both -pw and -wann are specified, this is ambiguous.
    sys.stderr.write("Cannot specify both -pw and -wann in one run.\n")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())

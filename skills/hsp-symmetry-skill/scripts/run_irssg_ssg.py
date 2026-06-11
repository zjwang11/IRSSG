#!/usr/bin/env python3
"""Run IRSSG SSG identification in an isolated directory and summarize output."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode, urlsplit, urlunsplit
from urllib.request import Request, urlopen


KNOWN_OUTPUTS = {
    "ssg.out",
    "ssg.err",
    "ssg_summary.json",
    "hsp_group_info.json",
    "ssg.data",
    "msg.data",
    "ssgop.npy",
    "POSCAR.symm",
    "POSCAR.ssg_primitive",
}


def bundled_irssg_package_root() -> Optional[Path]:
    """Return the IRSSG package root when this skill is inside the publish tree."""
    root = Path(__file__).resolve().parents[3]
    if (root / "pyproject.toml").is_file() and (root / "irssg").is_dir():
        return root
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run IRSSG -ssg/MOM2SSG on an mPOSCAR or mcif input and write a JSON summary."
    )
    parser.add_argument("input_file", help="mPOSCAR/POSCAR-like file or .mcif file")
    parser.add_argument("--output-dir", default="ssg-analysis", help="Directory for IRSSG outputs")
    parser.add_argument("--tolerance", type=float, default=1e-3, help="Spatial symmetry tolerance")
    parser.add_argument("--magtolerance", type=float, default=1e-4, help="Magnetic moment tolerance")
    parser.add_argument("--standardize", action="store_true", help="Pass --standardize to MOM2SSG/IRSSG")
    parser.add_argument("--repo-root", help="IRSSG source checkout containing src_mom2ssg/MOM2SSG.py")
    parser.add_argument("--irssg-cmd", default="irssg", help="IRSSG command name or path")
    parser.add_argument(
        "--irssg-install-source",
        default=os.environ.get("IRSSG_INSTALL_SOURCE"),
        help=(
            "Path to the installable IRSSG package root or wheel used when irssg is not installed. "
            "Defaults to the publish tree containing this skill, when available."
        ),
    )
    parser.add_argument(
        "--no-auto-install-irssg",
        action="store_true",
        help="Do not run python -m pip install automatically when IRSSG is missing.",
    )
    parser.add_argument(
        "--hsp-api-base-url",
        default=os.environ.get("HSP_API_BASE_URL", "https://cmpdc.iphy.ac.cn/hsp/api/v1"),
        help="HSP API base URL, e.g. https://cmpdc.iphy.ac.cn/hsp/api/v1",
    )
    parser.add_argument(
        "--hsp-web-base-url",
        default=os.environ.get("HSP_WEB_BASE_URL", "https://cmpdc.iphy.ac.cn/hsp"),
        help="HSP frontend base URL used only for human-readable links",
    )
    parser.add_argument("--hsp-api-key", help="API key/token for HSP requests")
    parser.add_argument("--hsp-api-key-file", help="File containing the HSP API key/token")
    parser.add_argument(
        "--hsp-api-key-env",
        default="HSP_API_KEY,HSP_WEB_API_KEY,CMPDC_API_KEY",
        help="Comma-separated environment variable names to read an HSP API key from",
    )
    parser.add_argument(
        "--hsp-api-key-header",
        default="Authorization",
        help="Header name used for the API key. Default: Authorization",
    )
    parser.add_argument(
        "--hsp-api-key-prefix",
        default="Bearer",
        help="Prefix prepended to the API key header value. Use '' for a raw key.",
    )
    parser.add_argument(
        "--hsp-api-key-query-param",
        default="",
        help="Optional query parameter name for APIs that expect the key in the URL instead of a header.",
    )
    parser.add_argument("--hsp-timeout", type=float, default=30.0, help="Timeout in seconds for HSP API requests")
    parser.add_argument("--skip-hsp-api", action="store_true", help="Do not fetch group details from the HSP API")
    parser.add_argument(
        "--require-hsp-api",
        action="store_true",
        help="Return a nonzero exit code if any expected HSP API group detail request fails",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite known files in output directory")
    return parser.parse_args()


def resolve_explicit_repo_root(explicit: Optional[str]) -> Optional[Path]:
    if not explicit:
        return None

    root = Path(explicit).expanduser().resolve()
    mom2ssg_path = root / "src_mom2ssg" / "MOM2SSG.py"
    if not mom2ssg_path.is_file():
        raise SystemExit(f"--repo-root does not contain src_mom2ssg/MOM2SSG.py: {root}")
    return root


def irssg_command_available(irssg_cmd: str) -> bool:
    if os.path.sep in irssg_cmd:
        return Path(irssg_cmd).expanduser().is_file()
    return shutil.which(irssg_cmd) is not None


def irssg_module_available() -> bool:
    try:
        return importlib.util.find_spec("irssg.ssg.MOM2SSG") is not None
    except (AttributeError, ImportError, ModuleNotFoundError, ValueError):
        return False


def shell_command(command: List[str]) -> str:
    return " ".join(shlex_quote(part) for part in command)


def shlex_quote(value: str) -> str:
    if re.match(r"^[A-Za-z0-9_@%+=:,./-]+$", value):
        return value
    return "'" + value.replace("'", "'\"'\"'") + "'"


def resolve_irssg_install_source(explicit: Optional[str]) -> Optional[Path]:
    if explicit:
        return Path(explicit).expanduser().resolve()
    return bundled_irssg_package_root()


def ensure_irssg_installed(args: argparse.Namespace, repo_root: Optional[Path]) -> None:
    if repo_root is not None:
        return
    if irssg_command_available(args.irssg_cmd) or irssg_module_available():
        return

    install_source = resolve_irssg_install_source(args.irssg_install_source)
    install_hint = (
        "Install IRSSG with:\n"
        "  python -m pip install irssg\n"
        "If dependency compilation fails, install common numerical dependencies first, for example:\n"
        "  conda install numpy h5py contourpy pandas\n"
        "  python -m pip install irssg\n"
        "For this local checkout you can also install from the publish tree:\n"
        "  python -m pip install /path/to/irssg/publish\n"
        "  python -m pip install /data/home/szhang/soft_sz/irssg_0702/publish\n"
        "Then verify with:\n"
        "  irssg --help\n"
        "  python -c \"import irssg.ssg.MOM2SSG; print('irssg ok')\""
    )

    if args.no_auto_install_irssg:
        raise SystemExit(f"IRSSG is not installed and --no-auto-install-irssg was passed.\n{install_hint}")
    package_source = str(install_source) if install_source is not None and install_source.exists() else "irssg"
    command = [sys.executable, "-m", "pip", "install", package_source]
    print(f"IRSSG is not installed; installing with: {shell_command(command)}", file=sys.stderr)
    proc = subprocess.run(command, text=True, capture_output=True)
    if proc.returncode != 0:
        raise SystemExit(
            "IRSSG automatic installation failed.\n"
            f"Command: {shell_command(command)}\n"
            f"stdout:\n{proc.stdout[-4000:]}\n"
            f"stderr:\n{proc.stderr[-4000:]}\n"
            f"{install_hint}"
        )

    importlib.invalidate_caches()
    if not (irssg_command_available(args.irssg_cmd) or irssg_module_available()):
        raise SystemExit(
            "IRSSG installation command completed, but the wrapper still cannot find the irssg command "
            "or irssg.ssg.MOM2SSG module. If pip installed console scripts outside PATH, add that bin "
            f"directory to PATH or pass --irssg-cmd explicitly.\n{install_hint}"
        )


def clean_output_dir(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in KNOWN_OUTPUTS:
        path = output_dir / name
        if path.is_file() or path.is_symlink():
            path.unlink()
        elif path.is_dir():
            shutil.rmtree(path)


def prepare_output_dir(output_dir: Path, overwrite: bool) -> None:
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise SystemExit(f"Output directory is not empty: {output_dir}. Use --overwrite or choose another directory.")
    if overwrite:
        clean_output_dir(output_dir)
    else:
        output_dir.mkdir(parents=True, exist_ok=True)


def base_mom2ssg_args(local_input: str, args: argparse.Namespace) -> List[str]:
    cmd_args = [
        "-c",
        local_input,
        "--tolerance",
        str(args.tolerance),
        "--magtolerance",
        str(args.magtolerance),
    ]
    if args.standardize:
        cmd_args.append("--standardize")
    return cmd_args


def build_attempts(args: argparse.Namespace, repo_root: Optional[Path], local_input: str) -> List[Tuple[str, List[str], Dict[str, str]]]:
    attempts: List[Tuple[str, List[str], Dict[str, str]]] = []

    module_args = base_mom2ssg_args(local_input, args)
    if repo_root is not None:
        env = os.environ.copy()
        env["PYTHONPATH"] = str(repo_root) + os.pathsep + env.get("PYTHONPATH", "")
        attempts.append(("source src_mom2ssg.MOM2SSG", [sys.executable, "-m", "src_mom2ssg.MOM2SSG", *module_args], env))

    irssg_path = shutil.which(args.irssg_cmd) if os.path.sep not in args.irssg_cmd else str(Path(args.irssg_cmd).expanduser())
    if irssg_path and Path(irssg_path).exists():
        attempts.append(("irssg command", [irssg_path, "-ssg", *module_args], os.environ.copy()))

    attempts.append(("installed irssg.ssg.MOM2SSG", [sys.executable, "-m", "irssg.ssg.MOM2SSG", *module_args], os.environ.copy()))
    return attempts


def run_attempts(attempts: Iterable[Tuple[str, List[str], Dict[str, str]]], cwd: Path) -> Tuple[str, List[str], str, str, List[Dict[str, object]]]:
    records: List[Dict[str, object]] = []
    last_stdout = ""
    last_stderr = ""

    for label, command, env in attempts:
        try:
            proc = subprocess.run(command, cwd=cwd, env=env, text=True, capture_output=True)
        except FileNotFoundError as exc:
            records.append(
                {
                    "label": label,
                    "command": command,
                    "returncode": 127,
                    "stderr_tail": str(exc),
                }
            )
            continue
        last_stdout = proc.stdout
        last_stderr = proc.stderr
        records.append(
            {
                "label": label,
                "command": command,
                "returncode": proc.returncode,
                "stderr_tail": "\n".join(proc.stderr.splitlines()[-20:]),
            }
        )
        if proc.returncode == 0:
            return label, command, proc.stdout, proc.stderr, records

    return "", [], last_stdout, last_stderr, records


def first_match(pattern: str, text: str) -> Optional[str]:
    match = re.search(pattern, text, re.MULTILINE)
    return match.group(1).strip() if match else None


def parse_atomic_space_group_number(value: object) -> Optional[str]:
    if not isinstance(value, str):
        return None
    match = re.match(r"\s*(\d{1,3})\b", value)
    return match.group(1) if match else None


def parse_msg_numbers(value: object) -> Dict[str, Optional[str]]:
    result: Dict[str, Optional[str]] = {"og": None, "bns": None, "ssg_setting": None}
    if not isinstance(value, str):
        return result

    og_match = re.search(r"([0-9]+\.[0-9]+\.[0-9]+)\s*\(OG setting\)", value)
    bns_match = re.search(r"([0-9]+\.[0-9]+)\s*\(BNS setting\)", value)
    ssg_match = re.search(r"([0-9]+(?:\.[0-9]+){3}(?:\.[LP])?)\s*\(SSG setting\)", value)
    if og_match:
        result["og"] = og_match.group(1)
    if bns_match:
        result["bns"] = bns_match.group(1)
    if ssg_match:
        result["ssg_setting"] = ssg_match.group(1)
    return result


def normalize_base_url(value: str) -> str:
    return value.strip().rstrip("/")


def append_query_param(url: str, key: str, value: str) -> str:
    parts = urlsplit(url)
    query = parts.query
    suffix = urlencode({key: value})
    query = f"{query}&{suffix}" if query else suffix
    return urlunsplit((parts.scheme, parts.netloc, parts.path, query, parts.fragment))


def read_hsp_api_key(args: argparse.Namespace) -> Optional[str]:
    if args.hsp_api_key:
        return args.hsp_api_key.strip()

    if args.hsp_api_key_file:
        key_path = Path(args.hsp_api_key_file).expanduser()
        if key_path.is_file():
            return key_path.read_text(encoding="utf-8").strip()

    for env_name in (item.strip() for item in args.hsp_api_key_env.split(",")):
        if not env_name:
            continue
        value = os.environ.get(env_name)
        if value:
            return value.strip()
    return None


def build_hsp_headers(args: argparse.Namespace, api_key: Optional[str]) -> Dict[str, str]:
    headers = {"Accept": "application/json"}
    if api_key and args.hsp_api_key_header and not args.hsp_api_key_query_param:
        prefix = args.hsp_api_key_prefix.strip()
        headers[args.hsp_api_key_header] = f"{prefix} {api_key}".strip() if prefix else api_key
    return headers


def fetch_json(url: str, headers: Dict[str, str], timeout: float) -> Dict[str, object]:
    request = Request(url, headers=headers, method="GET")
    try:
        with urlopen(request, timeout=timeout) as response:
            raw_body = response.read()
            text = raw_body.decode("utf-8")
            payload = json.loads(text) if text else None
            return {
                "ok": True,
                "status_code": response.status,
                "url": url,
                "data": payload,
                "error": None,
            }
    except HTTPError as exc:
        try:
            body = exc.read().decode("utf-8", errors="replace")
        except Exception:
            body = ""
        return {
            "ok": False,
            "status_code": exc.code,
            "url": url,
            "data": None,
            "error": body[:2000] or str(exc),
        }
    except (URLError, TimeoutError, json.JSONDecodeError, OSError) as exc:
        return {
            "ok": False,
            "status_code": None,
            "url": url,
            "data": None,
            "error": str(exc),
        }


def hsp_detail_url(base_url: str, *parts: str) -> str:
    encoded = "/".join(quote(part.strip("/"), safe="") for part in parts)
    return f"{normalize_base_url(base_url)}/{encoded}"


def build_hsp_identifiers(summary: Dict[str, object]) -> Dict[str, Optional[str]]:
    msg_numbers = parse_msg_numbers(summary.get("msg_number"))
    ssg_number = summary.get("ssg_number")
    return {
        "space_group_number": parse_atomic_space_group_number(summary.get("atomic_space_group"))
        or parse_atomic_space_group_number(summary.get("nonmagnetic_space_group")),
        "magnetic_group_og_number": msg_numbers["og"],
        "magnetic_group_bns_number": msg_numbers["bns"],
        "magnetic_group_ssg_setting": msg_numbers["ssg_setting"],
        "spin_group_number": ssg_number.strip() if isinstance(ssg_number, str) and ssg_number.strip() else None,
    }


def fetch_hsp_group_info(summary: Dict[str, object], output_dir: Path, args: argparse.Namespace) -> Dict[str, object]:
    identifiers = build_hsp_identifiers(summary)
    api_base_url = normalize_base_url(args.hsp_api_base_url)
    web_base_url = normalize_base_url(args.hsp_web_base_url)
    api_key = read_hsp_api_key(args)
    headers = build_hsp_headers(args, api_key)

    def api_url(*parts: str) -> str:
        url = hsp_detail_url(api_base_url, *parts)
        if api_key and args.hsp_api_key_query_param:
            return append_query_param(url, args.hsp_api_key_query_param, api_key)
        return url

    requests: Dict[str, object] = {}

    space_group_number = identifiers["space_group_number"]
    if space_group_number:
        requests["space_group"] = fetch_json(api_url("space-groups", space_group_number), headers, args.hsp_timeout)

    magnetic_group_number = identifiers["magnetic_group_og_number"] or identifiers["magnetic_group_bns_number"]
    if magnetic_group_number:
        requests["magnetic_group"] = fetch_json(
            api_url("magnetic-groups", magnetic_group_number), headers, args.hsp_timeout
        )

    spin_group_number = identifiers["spin_group_number"]
    if spin_group_number:
        requests["spin_group"] = fetch_json(api_url("spin-groups", spin_group_number), headers, args.hsp_timeout)

    web_links = {
        "space_group": hsp_detail_url(web_base_url, "space-groups", space_group_number) if space_group_number else None,
        "magnetic_group": hsp_detail_url(
            web_base_url, "magnetic-groups", magnetic_group_number.split(".", 1)[0], magnetic_group_number
        )
        if magnetic_group_number
        else None,
        "spin_group": hsp_detail_url(web_base_url, "spin-groups", spin_group_number.split(".", 1)[0], spin_group_number)
        if spin_group_number
        else None,
    }

    group_info: Dict[str, object] = {
        "status": "ok" if requests and all(bool(item.get("ok")) for item in requests.values()) else "partial_or_failed",
        "api_base_url": api_base_url,
        "web_base_url": web_base_url,
        "used_api_key": bool(api_key),
        "api_key_header": args.hsp_api_key_header if api_key and not args.hsp_api_key_query_param else None,
        "api_key_query_param": args.hsp_api_key_query_param if api_key and args.hsp_api_key_query_param else None,
        "identifiers": identifiers,
        "web_links": web_links,
        "requests": {
            key: {
                "ok": value.get("ok"),
                "status_code": value.get("status_code"),
                "url": value.get("url"),
                "error": value.get("error"),
            }
            for key, value in requests.items()
            if isinstance(value, dict)
        },
        "space_group": requests.get("space_group", {}).get("data") if isinstance(requests.get("space_group"), dict) else None,
        "magnetic_group": requests.get("magnetic_group", {}).get("data")
        if isinstance(requests.get("magnetic_group"), dict)
        else None,
        "spin_group": requests.get("spin_group", {}).get("data") if isinstance(requests.get("spin_group"), dict) else None,
    }

    if not requests:
        group_info["status"] = "skipped_no_group_identifiers"

    output_path = output_dir / "hsp_group_info.json"
    output_path.write_text(json.dumps(group_info, indent=2, ensure_ascii=False), encoding="utf-8")
    return group_info


def parse_operation_count(text: str, section_header: str) -> Optional[int]:
    idx = text.find(section_header)
    if idx < 0:
        return None
    match = re.search(r"^# Number:\s*(\d+)\s*$", text[idx:], re.MULTILINE)
    return int(match.group(1)) if match else None


def parse_spin_only_group(text: str) -> Optional[str]:
    for pattern in (
        r"^I:\s+Collinear SSG.*?So:\s*(.+)$",
        r"^II:\s+Coplanar SSG.*?So:\s*(.+)$",
        r"^III:\s+Noncoplanar SSG.*?So:\s*(.+)$",
    ):
        value = first_match(pattern, text)
        if value:
            return value
    return None


def parse_magnetic_class(text: str) -> Optional[str]:
    if re.search(r"^I:\s+Collinear SSG", text, re.MULTILINE):
        return "collinear"
    if re.search(r"^II:\s+Coplanar SSG", text, re.MULTILINE):
        return "coplanar"
    if re.search(r"^III:\s+Noncoplanar SSG", text, re.MULTILINE):
        return "noncoplanar"
    return None


def parse_summary(text: str, output_dir: Path, input_name: str, method: str, command: List[str], attempts: List[Dict[str, object]]) -> Dict[str, object]:
    generated = {path.name for path in output_dir.iterdir() if path.name in KNOWN_OUTPUTS or path.name.startswith("POSCAR.")}
    generated.add("ssg_summary.json")
    summary: Dict[str, object] = {
        "input_file": input_name,
        "run_method": method,
        "command": command,
        "ssg_number": first_match(r"^The SSG number:\s*(.+)$", text),
        "ssg_international_symbol": first_match(r"^The SSG international symbol:\s*(.+)$", text),
        "magnetic_configuration": parse_magnetic_class(text),
        "spin_only_group": parse_spin_only_group(text),
        "spin_part_of_go": first_match(r"^P \(spin part of Go\):\s*(.+)$", text),
        "lattice_part_of_go": first_match(r"^H \(lattice part of Go\):\s*(.+)$", text),
        "ssg_operation_count": parse_operation_count(text, "Spin space group operations"),
        "atomic_space_group": first_match(r"^Atomic space group:\s*(.+)$", text),
        "n_asg_over_n_ssg": first_match(r"^N_ASG/N_SSG\s*=\s*(.+)$", text),
        "msg_number": first_match(r"^The MSG number:\s*(.+)$", text),
        "msg_international_symbol": first_match(r"^The MSG international symbol:\s*(.+)$", text),
        "msg_operation_count": parse_operation_count(text, "Magnetic space group operations"),
        "ssg_url": first_match(r"^(https://cmpdc\.iphy\.ac\.cn/ssg/ssgs/\S+)\s*$", text),
        "generated_files": sorted(generated),
        "attempts": attempts,
    }

    if summary["ssg_number"] is None:
        sg = first_match(r"^The SG number:\s*(.+)$", text)
        if sg:
            summary["nonmagnetic_space_group"] = sg
            summary["status"] = "nonmagnetic"
        else:
            summary["status"] = "parsed_without_ssg_number"
    else:
        summary["status"] = "ok"
    return summary


def main() -> int:
    args = parse_args()
    input_path = Path(args.input_file).expanduser().resolve()
    if not input_path.is_file():
        raise SystemExit(f"Input file not found: {input_path}")

    output_dir = Path(args.output_dir).expanduser().resolve()
    prepare_output_dir(output_dir, args.overwrite)

    local_input = input_path.name
    local_input_path = output_dir / local_input
    if local_input_path.resolve() != input_path:
        shutil.copy2(input_path, local_input_path)

    repo_root = resolve_explicit_repo_root(args.repo_root)
    ensure_irssg_installed(args, repo_root)
    attempts = build_attempts(args, repo_root, local_input)
    method, command, stdout, stderr, records = run_attempts(attempts, output_dir)

    (output_dir / "ssg.out").write_text(stdout, encoding="utf-8")
    (output_dir / "ssg.err").write_text(stderr, encoding="utf-8")

    if not method:
        failure = {
            "status": "failed",
            "input_file": str(input_path),
            "output_dir": str(output_dir),
            "attempts": records,
        }
        (output_dir / "ssg_summary.json").write_text(json.dumps(failure, indent=2), encoding="utf-8")
        print(f"IRSSG SSG identification failed. See {output_dir / 'ssg.err'}", file=sys.stderr)
        return 1

    summary = parse_summary(stdout, output_dir, str(input_path), method, command, records)
    summary["output_dir"] = str(output_dir)

    if not args.skip_hsp_api:
        group_info = fetch_hsp_group_info(summary, output_dir, args)
        summary["hsp_group_info_file"] = str(output_dir / "hsp_group_info.json")
        summary["hsp_group_identifiers"] = group_info.get("identifiers")
        summary["hsp_api_status"] = group_info.get("status")
        generated_files = set(summary.get("generated_files", []))
        generated_files.add("hsp_group_info.json")
        summary["generated_files"] = sorted(generated_files)

    summary_path = output_dir / "ssg_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    if args.require_hsp_api and summary.get("hsp_api_status") != "ok":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

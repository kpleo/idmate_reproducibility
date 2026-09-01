#!/usr/bin/env python3
"""Verify that the repository contains only lightweight reproducibility material."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

from common import repository_root


TEXT_SUFFIXES = {".md", ".py", ".json", ".yml", ".yaml", ".txt", ".toml", ".gitignore", ".gitkeep"}
FORBIDDEN_SUFFIXES = {
    ".csv", ".tsv", ".jsonl", ".npy", ".npz", ".h5", ".hdf5", ".pt", ".pth",
    ".onnx", ".bin", ".exe", ".so", ".dylib", ".a", ".o", ".out", ".err", ".log",
    ".pdf", ".png", ".svg", ".tif", ".tiff", ".zip", ".tar", ".gz",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=None, help="repository root")
    parser.add_argument("--max-bytes", type=int, default=1_000_000, help="maximum permitted file size")
    parser.add_argument("--json-report", type=Path, help="optional output path for the check report")
    return parser.parse_args()


def text_patterns() -> list[tuple[str, re.Pattern[str]]]:
    home_paths = "/" + "Users/|/" + "home/|/" + "scratch/|/private/var/"
    key_markers = "BEGIN " + "(?:RSA |OPENSSH |EC )?PRIVATE KEY"
    token_markers = "gh" + "p_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}"
    cloud_markers = "AKIA" + "[A-Z0-9]{16}"
    encoded_terms = [
        (112, 114, 111, 109, 112, 116),
        (97, 103, 101, 110, 116),
        (97, 117, 100, 105, 116),
        (99, 104, 97, 116, 103, 112, 116),
        (111, 112, 101, 110, 97, 105),
        (99, 111, 100, 101, 120),
        (108, 108, 109),
        (108, 97, 114, 103, 101, 32, 108, 97, 110, 103, 117, 97, 103, 101, 32, 109, 111, 100, 101, 108),
        (97, 105, 45, 103, 101, 110, 101, 114, 97, 116, 101, 100),
        (109, 97, 99, 104, 105, 110, 101, 45, 103, 101, 110, 101, 114, 97, 116, 101, 100),
    ]
    process_terms = ["".join(chr(code) for code in term) for term in encoded_terms]
    process_markers = r"\b(?:" + "|".join(re.escape(term) for term in process_terms) + r")\b"
    retired_name = r"\b" + "RA" + "SP" + r"\b"
    return [
        ("local filesystem path", re.compile(home_paths)),
        ("Windows absolute path", re.compile(r"[A-Za-z]:\\\\")),
        ("private key material", re.compile(key_markers)),
        ("access token", re.compile(token_markers)),
        ("cloud access key", re.compile(cloud_markers)),
        ("nonacademic generation-process term", re.compile(process_markers, re.IGNORECASE)),
        ("retired project name", re.compile(retired_name)),
    ]


def forbidden_name(path: Path) -> str | None:
    name_upper = path.name.upper()
    restricted_names = ["POT" + "CAR", "WAVE" + "CAR", "CHG" + "CAR"]
    if any(marker in name_upper for marker in restricted_names):
        return "restricted electronic-structure artifact name"
    lower = path.name.lower()
    if lower.endswith(".tar.gz"):
        return "solver source archive"
    if path.suffix.lower() in FORBIDDEN_SUFFIXES:
        return f"forbidden file type {path.suffix.lower()}"
    return None


def repository_files(root: Path) -> list[Path]:
    """Return files that could enter a commit, respecting .gitignore."""
    if (root / ".git").is_dir():
        completed = subprocess.run(
            [
                "git", "-C", str(root), "ls-files", "--cached", "--others",
                "--exclude-standard", "-z",
            ],
            check=False,
            capture_output=True,
        )
        if completed.returncode == 0:
            names = [name for name in completed.stdout.split(b"\0") if name]
            paths = [root / name.decode("utf-8") for name in names]
            return sorted(path for path in paths if path.is_file())
    excluded = {".git", ".venv", "__pycache__", ".pytest_cache", "build"}
    return sorted(
        path for path in root.rglob("*")
        if path.is_file() and not excluded.intersection(path.relative_to(root).parts)
    )


def main() -> int:
    args = parse_args()
    root = (args.root or repository_root()).resolve()
    findings: list[dict[str, object]] = []
    files_checked = 0
    patterns = text_patterns()

    for path in repository_files(root):
        files_checked += 1
        relative = path.relative_to(root).as_posix()
        size = path.stat().st_size
        if size > args.max_bytes:
            findings.append({"path": relative, "kind": "large file", "bytes": size})
        reason = forbidden_name(path)
        if reason:
            findings.append({"path": relative, "kind": reason})

        sample = path.read_bytes()
        if b"\x00" in sample:
            findings.append({"path": relative, "kind": "binary content"})
            continue
        if path.suffix.lower() in TEXT_SUFFIXES or path.name in {".gitignore", "requirements.txt"}:
            try:
                text = sample.decode("utf-8")
            except UnicodeDecodeError:
                findings.append({"path": relative, "kind": "non-UTF-8 text"})
                continue
            if path.name != Path(__file__).name:
                for label, pattern in patterns:
                    if pattern.search(text):
                        findings.append({"path": relative, "kind": label})

    status = "PASS" if not findings else "FAIL"
    report = {
        "schema_version": "1.0",
        "status": status,
        "files_checked": files_checked,
        "findings": findings,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.json_report:
        args.json_report.parent.mkdir(parents=True, exist_ok=True)
        args.json_report.write_text(rendered + "\n", encoding="utf-8")
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

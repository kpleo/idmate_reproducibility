#!/usr/bin/env python3
"""Validate supplied records and regenerate all mapped figures and tables."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from common import repository_root


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True, help="directory containing derived inputs")
    parser.add_argument("--output-root", type=Path, default=Path("build"), help="destination root")
    parser.add_argument("--manifest", type=Path, default=None, help="alternate manifest JSON")
    return parser.parse_args()


def run(command: list[str], *, root: Path) -> int:
    completed = subprocess.run(command, cwd=root, check=False)
    return completed.returncode


def main() -> int:
    args = parse_args()
    root = repository_root()
    manifest_args = ["--manifest", str(args.manifest)] if args.manifest else []
    commands = [
        [
            sys.executable, "scripts/validate_inputs.py", "--data-root", str(args.data_root),
            *manifest_args,
        ],
        [
            sys.executable, "scripts/generate_figures.py", "--data-root", str(args.data_root),
            "--output-dir", str(args.output_root / "figures"), *manifest_args,
        ],
        [
            sys.executable, "scripts/generate_tables.py", "--data-root", str(args.data_root),
            "--output-dir", str(args.output_root / "tables"), *manifest_args,
        ],
    ]
    for command in commands:
        return_code = run(command, root=root)
        if return_code:
            return return_code
    print(f"PASS: outputs written under {args.output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

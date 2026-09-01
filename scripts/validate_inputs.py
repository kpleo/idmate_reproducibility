#!/usr/bin/env python3
"""Validate author-supplied derived records against the manuscript manifest."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import (
    ReproductionError,
    load_manifest,
    repository_root,
    safe_path,
    selected_artifacts,
    validate_csv,
    validate_table_csv,
    validate_worklaw_json,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=None, help="manifest JSON; defaults to repository manifest.json")
    parser.add_argument("--data-root", type=Path, default=None, help="directory containing derived inputs")
    parser.add_argument("--only", nargs="*", help="artifact identifiers to validate; default is all")
    parser.add_argument("--allow-missing", action="store_true", help="report absent non-distributed inputs without returning a nonzero status")
    parser.add_argument("--json-report", type=Path, help="optional path for the validation report")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = repository_root()
    manifest_path = args.manifest or root / "manifest.json"
    try:
        manifest = load_manifest(manifest_path)
        data_root = args.data_root or root / manifest["default_data_root"]
        figures, tables = selected_artifacts(manifest, args.only)
    except ReproductionError as exc:
        print(f"ERROR: {exc}")
        return 1

    declarations: dict[str, tuple[str, object]] = {}
    for figure in figures:
        for source in figure["inputs"]:
            declarations[source["path"]] = (source["format"], source)
    for table in tables:
        for block in table["blocks"]:
            declarations[block["input"]] = ("table_csv", block)

    files: list[dict[str, object]] = []
    missing_count = 0
    error_count = 0
    for relative in sorted(declarations):
        kind, declaration = declarations[relative]
        try:
            path = safe_path(data_root, relative)
            if not path.is_file():
                files.append({"path": relative, "status": "missing"})
                missing_count += 1
                continue
            if kind == "csv":
                rows = validate_csv(path, declaration["columns"])
            elif kind == "table_csv":
                rows = validate_table_csv(path, declaration["columns"])
            elif kind == "json" and relative == "fig3_worklaw_summary.json":
                validate_worklaw_json(path)
                rows = None
            else:
                raise ReproductionError(f"no validator is defined for {relative}")
            record: dict[str, object] = {"path": relative, "status": "valid"}
            if rows is not None:
                record["rows"] = rows
            files.append(record)
        except ReproductionError as exc:
            files.append({"path": relative, "status": "invalid", "message": str(exc)})
            error_count += 1

    if error_count:
        status = "ERROR"
    elif missing_count:
        status = "INCOMPLETE"
    else:
        status = "PASS"
    report = {
        "schema_version": "1.0",
        "status": status,
        "files_checked": len(files),
        "missing": missing_count,
        "invalid": error_count,
        "files": files,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.json_report:
        args.json_report.parent.mkdir(parents=True, exist_ok=True)
        args.json_report.write_text(rendered + "\n", encoding="utf-8")

    if error_count:
        return 1
    if missing_count and not args.allow_missing:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


#!/usr/bin/env python3
"""Generate manuscript LaTeX table fragments from declared derived CSV files."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import (
    ReproductionError,
    latex_escape,
    load_manifest,
    read_csv_rows,
    repository_root,
    safe_path,
    validate_table_csv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=None, help="manifest JSON; defaults to repository manifest.json")
    parser.add_argument("--data-root", type=Path, default=None, help="directory containing derived table CSV files")
    parser.add_argument("--output-dir", type=Path, default=None, help="destination for LaTeX fragments")
    parser.add_argument("--only", action="append", help="one table identifier; repeat to select more than one")
    parser.add_argument("--list", action="store_true", help="list table identifiers and exit")
    return parser.parse_args()


def make_tabular(rows: list[dict[str, str]], columns: list[dict[str, object]]) -> list[str]:
    alignment = "l" + "r" * (len(columns) - 1)
    headings = " & ".join(latex_escape(str(column["heading"])) for column in columns)
    lines = [f"\\begin{{tabular}}{{@{{}}{alignment}@{{}}}}", "\\toprule", headings + " \\\\", "\\midrule"]
    for row in rows:
        values = [latex_escape(row[str(column["key"])].strip()) for column in columns]
        lines.append(" & ".join(values) + " \\\\")
    lines.extend(["\\bottomrule", "\\end{tabular}"])
    return lines


def render_table(table: dict[str, object], data_root: Path) -> str:
    lines = [
        "% Generated from author-supplied derived records; do not edit values here.",
        "\\begin{table*}[t]",
        "\\centering",
        f"\\caption{{{latex_escape(str(table['display']))}.}}",
        f"\\label{{{table['id']}}}",
    ]
    for block_index, block in enumerate(table["blocks"]):
        path = safe_path(data_root, str(block["input"]))
        if not path.is_file():
            raise ReproductionError(
                f"required derived input is not distributed: {block['input']}; "
                "supply it under --data-root"
            )
        validate_table_csv(path, block["columns"])
        rows = read_csv_rows(path)
        if block_index:
            lines.append("\\vspace{4pt}")
        lines.extend(make_tabular(rows, block["columns"]))
    lines.extend(["\\end{table*}", ""])
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    root = repository_root()
    try:
        manifest = load_manifest(args.manifest or root / "manifest.json")
        tables = manifest["tables"]
        if args.list:
            for table in tables:
                print(f"{table['id']}\t{table['display']}")
            return 0
        requested = set(args.only or [])
        known = {table["id"] for table in tables}
        unknown = sorted(requested - known)
        if unknown:
            raise ReproductionError(f"unknown table identifier(s): {', '.join(unknown)}")
        selected = [table for table in tables if not requested or table["id"] in requested]
        data_root = args.data_root or root / manifest["default_data_root"]
        output_dir = args.output_dir or root / "outputs" / "tables"
        output_dir.mkdir(parents=True, exist_ok=True)
        for table in selected:
            rendered = render_table(table, data_root)
            destination = output_dir / table["output"]
            destination.write_text(rendered, encoding="utf-8")
            print(f"wrote {destination}")
    except ReproductionError as exc:
        print(f"INCOMPLETE: {exc}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

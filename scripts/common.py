"""Shared, standard-library helpers for the reproduction commands."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable


class ReproductionError(RuntimeError):
    """Raised when a declared input or repository contract is invalid."""


def repository_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_manifest(path: Path) -> dict[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ReproductionError(f"manifest not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ReproductionError(f"manifest is not valid JSON: {path}: {exc}") from exc
    if document.get("schema_version") != "1.0":
        raise ReproductionError(
            f"unsupported manifest schema: {document.get('schema_version')!r}"
        )
    return document


def safe_path(root: Path, relative: str) -> Path:
    requested = Path(relative)
    if requested.is_absolute() or ".." in requested.parts:
        raise ReproductionError(f"unsafe relative input path: {relative!r}")
    root_resolved = root.resolve()
    result = (root_resolved / requested).resolve()
    if result != root_resolved and root_resolved not in result.parents:
        raise ReproductionError(f"input path escapes data root: {relative!r}")
    return result


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ReproductionError(f"CSV has no header: {path}")
            rows = list(reader)
    except UnicodeDecodeError as exc:
        raise ReproductionError(f"CSV is not UTF-8 text: {path}") from exc
    if not rows:
        raise ReproductionError(f"CSV has no data rows: {path}")
    return rows


def read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except UnicodeDecodeError as exc:
        raise ReproductionError(f"JSON is not UTF-8 text: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ReproductionError(f"invalid JSON: {path}: {exc}") from exc


def parse_value(value: str, declared_type: str, *, nullable: bool, context: str) -> Any:
    if value == "":
        if nullable:
            return None
        raise ReproductionError(f"{context}: empty value is not permitted")
    try:
        if declared_type == "string":
            return value
        if declared_type == "number":
            number = float(value)
            if not math.isfinite(number):
                raise ValueError("non-finite")
            return number
        if declared_type == "integer":
            number = int(value)
            return number
        if declared_type == "boolean":
            normalized = value.strip().lower()
            if normalized in {"true", "1", "yes"}:
                return True
            if normalized in {"false", "0", "no"}:
                return False
            raise ValueError("expected true/false")
    except ValueError as exc:
        raise ReproductionError(
            f"{context}: {value!r} is not a valid {declared_type}"
        ) from exc
    raise ReproductionError(f"{context}: unsupported declared type {declared_type!r}")


def validate_csv(path: Path, columns: list[dict[str, Any]]) -> int:
    rows = read_csv_rows(path)
    header = set(rows[0])
    expected = {column["name"] for column in columns}
    missing = sorted(expected - header)
    if missing:
        raise ReproductionError(f"{path}: missing CSV columns: {', '.join(missing)}")
    for row_number, row in enumerate(rows, start=2):
        for column in columns:
            name = column["name"]
            parse_value(
                row.get(name, ""),
                column["type"],
                nullable=bool(column.get("nullable", False)),
                context=f"{path}:{row_number}:{name}",
            )
    return len(rows)


def validate_table_csv(path: Path, columns: list[dict[str, Any]]) -> int:
    normalized = [
        {
            "name": column["key"],
            "type": column["type"],
            "nullable": column.get("nullable", False),
        }
        for column in columns
    ]
    return validate_csv(path, normalized)


def validate_worklaw_json(path: Path) -> None:
    document = read_json(path)
    if not isinstance(document, dict):
        raise ReproductionError(f"{path}: top-level JSON value must be an object")
    required = {"models", "categories", "counts"}
    missing = sorted(required - set(document))
    if missing:
        raise ReproductionError(f"{path}: missing JSON keys: {', '.join(missing)}")

    models = document["models"]
    if not isinstance(models, list) or not models:
        raise ReproductionError(f"{path}: models must be a non-empty array")
    model_keys = {"name", "improvement_fraction", "ci_low", "ci_high", "criterion"}
    for index, model in enumerate(models):
        if not isinstance(model, dict) or not model_keys <= set(model):
            raise ReproductionError(f"{path}: models[{index}] has an invalid structure")
        for key in model_keys - {"name"}:
            if not isinstance(model[key], (int, float)) or not math.isfinite(float(model[key])):
                raise ReproductionError(f"{path}: models[{index}].{key} must be finite")

    categories = document["categories"]
    if not isinstance(categories, list) or not categories:
        raise ReproductionError(f"{path}: categories must be a non-empty array")
    for index, category in enumerate(categories):
        if not isinstance(category, dict) or not {"name", "effect_fraction", "status"} <= set(category):
            raise ReproductionError(f"{path}: categories[{index}] has an invalid structure")
        effect = category["effect_fraction"]
        if effect is not None and (
            not isinstance(effect, (int, float)) or not math.isfinite(float(effect))
        ):
            raise ReproductionError(f"{path}: categories[{index}].effect_fraction must be finite or null")
        if category["status"] not in {"evaluated", "insufficient_data"}:
            raise ReproductionError(f"{path}: categories[{index}].status is invalid")

    counts = document["counts"]
    count_keys = {"trajectories", "converged_trajectories", "paired_steps", "parent_lineages"}
    if not isinstance(counts, dict) or not count_keys <= set(counts):
        raise ReproductionError(f"{path}: counts has an invalid structure")
    for key in count_keys:
        if not isinstance(counts[key], int) or counts[key] < 0:
            raise ReproductionError(f"{path}: counts.{key} must be a non-negative integer")


def selected_artifacts(
    manifest: dict[str, Any], selected: Iterable[str] | None
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    requested = set(selected or [])
    figures = manifest["figures"]
    tables = manifest["tables"]
    known = {item["id"] for item in figures + tables}
    unknown = sorted(requested - known)
    if unknown:
        raise ReproductionError(f"unknown artifact identifier(s): {', '.join(unknown)}")
    if not requested:
        return figures, tables
    return (
        [item for item in figures if item["id"] in requested],
        [item for item in tables if item["id"] in requested],
    )


def latex_escape(value: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(character, character) for character in value)

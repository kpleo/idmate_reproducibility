"""Reconstruct the historical-subspace comparison table from supplied scalar data."""

import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
LABELS = {"si": "diamond Si", "graphene": "graphene", "al": "fcc Al"}


def number(value, scientific=False):
    if value is None:
        return "--"
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("Nonfinite value")
    if value == 0:
        return "$0.000$"
    if not scientific and 1e-2 <= abs(value) < 1e3:
        return "${:.3f}$".format(value)
    mantissa, exponent = "{:.3e}".format(value).split("e")
    return "$" + mantissa + r"\times10^{" + str(int(exponent)) + "}$"


def validate(data):
    rows = data["rows"]
    if [row["system"] for row in rows] != ["si", "graphene", "al"]:
        raise ValueError("The table must retain every prescribed material in order")
    if data["reference_maps"] != 6 or any(row["maps"] != 2 for row in rows):
        raise ValueError("Reference count mismatch")
    for row in rows:
        if row["stop_required"] or row["oracle"]:
            raise ValueError("Stopped attempt is not a completed comparison")
        if row["decision"] == "accept":
            if row["window_bound"] is None or row["window_distance"] is None:
                raise ValueError("An accepted point needs window quantities")
        elif row["decision"] == "abstain":
            if not row["reason"] or row["window_bound"] is not None:
                raise ValueError("An unavailable window must remain unavailable")
        else:
            raise ValueError("Unexpected recorded decision")
        for key in ("full_distance", "density_distance", "next_mixed_density_distance"):
            if not math.isfinite(row[key]) or row[key] < 0:
                raise ValueError("Invalid distance")
        if row["next_mixed_density_distance"] is not None:
            residual = row["next_mixed_density_distance"] - row["alpha"] * row["density_distance"]
            if abs(residual) > 1e-12:
                raise ValueError("Mixing identity mismatch")
    return rows


def mesh(values):
    if len(set(values)) == 1:
        return "$" + str(values[0]) + "^3$"
    return "$" + r"\times".join(str(value) for value in values) + "$"


def table_rows(data):
    rows = validate(data)
    main_rows, settings_rows, diagnostic_rows = [], [], []
    for row in rows:
        fields = [LABELS[row["system"]], row["decision"], number(row["window_bound"]),
                  number(row["window_distance"]), number(row["full_distance"]),
                  number(row["density_distance"])]
        main_rows.append(" & ".join(fields) + r" \\")
        setting = data["settings"][row["system"]]
        kmesh = [len({point["frac"][axis] for point in setting["kpoints"]})
                 for axis in range(3)]
        if kmesh[0] * kmesh[1] * kmesh[2] != len(setting["kpoints"]):
            raise ValueError("Expected a complete product k mesh")
        atoms = sum(len(species["positions_frac"]) for species in setting["species"])
        fields = [LABELS[row["system"]], str(atoms),
                  str(int(setting["target_electrons"])), number(setting["cutoff_ha"]),
                  mesh(kmesh), mesh(setting["grid_dims"]),
                  "$(%.3f,%.3f)$" % (setting["tau_ha"], setting["linear_mixing"])]
        settings_rows.append(" & ".join(fields) + r" \\")
        fields = [LABELS[row["system"]], number(row["delta_F_full_ha"], True),
                  number(row["next_mixed_density_distance"]),
                  number(row["mixing_identity_error"], True)]
        diagnostic_rows.append(" & ".join(fields) + r" \\")
    return {"material_rows.tex": main_rows, "settings_rows.tex": settings_rows,
            "diagnostic_rows.tex": diagnostic_rows}


def main(output=None):
    data = json.loads((HERE / "materials.json").read_text())
    tables = table_rows(data)
    output = Path(output) if output is not None else HERE.parent / "build" / "material_tables"
    output.mkdir(parents=True, exist_ok=True)
    for name, lines in tables.items():
        (output / name).write_text("\n".join(lines) + "\n")
    fields = ["system", "decision", "reason", "window_bound", "window_distance",
              "full_distance", "density_distance", "next_mixed_density_distance",
              "delta_F_full_ha", "delta_F_window_ha"]
    with (output / "materials.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data["rows"])
    print("\n".join(tables["material_rows.tex"]))


if __name__ == "__main__":
    main()

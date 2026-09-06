"""Portable numerical readers shared by the five figure builders."""

import csv
from bisect import bisect_right
import json
import math
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RELEASE = HERE.parent
DATA = RELEASE / "data"
OUTPUT = Path(os.environ.get("IDMATE_FIGURE_OUTPUT", RELEASE / "build" / "figures"))
ONLINE = "online_low_cutoff"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, label, atol=2e-12):
    require(math.isclose(actual, expected, rel_tol=0, abs_tol=atol), label)


def ecdf(values):
    values = np.sort(values)
    return values, np.arange(1, len(values) + 1) / len(values)


def read_json(name):
    return json.loads((DATA / name).read_text())


def read_csv(name):
    with (DATA / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def load_response():
    arrays = {}
    counts = {"e0": (10000, 114), "e3_graphene": (164, 0),
              "e3_al": (1614, 0), "gk_graphene": (310, 2), "gk_al": (400, 0)}
    for name, (count, missing) in counts.items():
        rows = read_csv(name + ".csv")
        require(len(rows) == count, "Distance-series row count changed")
        values = []
        undefined = 0
        for row in rows:
            distance, bound, allowance = (float(row[key]) for key in ("distance", "bound", "allowance"))
            require(all(math.isfinite(x) and x >= 0 for x in (distance, bound, allowance)),
                    "Nonfinite or negative distance input")
            require(distance - bound - allowance <= 0, "Recorded numerical allowance exceeded")
            if row["ratio"] == "":
                undefined += 1
            else:
                ratio = float(row["ratio"])
                require(bound > 0 and math.isclose(ratio, distance / bound, rel_tol=1e-13),
                        "Stored ratio differs from distance divided by bound")
                values.append(ratio)
        require(undefined == missing, "Undefined ratios must remain unplotted")
        arrays[name] = np.array(values)
    shells = read_csv("shells.csv")
    retained = [[int(row["shell"]), float(row["share"]), float(row["cumulative"])]
                for row in shells if row["retained"] == "True"]
    tail = [[int(row["shell"]), float(row["share"]), float(row["cumulative"])]
            for row in shells if int(row["shell"]) <= 40]
    return arrays, retained, {}, tail, {}


def load_rows():
    rows = read_csv("screen_scope.csv")
    require(len(rows) == 10, "Expected ten historical candidates")
    for row in rows:
        for key in ("window_error", "full_error", "window_bound", "window_over_bound"):
            row[key] = float(row[key])
            require(math.isfinite(row[key]) and row[key] > 0, "Invalid log coordinate")
        require(math.isclose(row["window_over_bound"], row["window_error"] / row["window_bound"],
                             rel_tol=1e-14), "Historical ratio mismatch")
    return rows


def bm3_energy(volume, e0, v0, b0, bp):
    eta = (v0 / np.asarray(volume)) ** (2.0 / 3.0)
    return e0 + 9.0 * v0 * b0 / 16.0 * ((eta - 1.0) ** 3 * bp + (eta - 1.0) ** 2 * (6.0 - 4.0 * eta))


def load_physics():
    data = read_json("silicon_physics.json")
    data["bm3"] = bm3_energy
    for record in data["eos"].values():
        for key in ("volumes", "energies"):
            record[key] = np.array(record[key])
        fit = record["fit"]
        curve = bm3_energy(record["volumes"], fit["E0_eV"], fit["V0_A3"], fit["B0_eV_A3"], fit["B_prime"])
        record["shifted"] = (record["energies"] - fit["E0_eV"]) * 500
        record["residual_meV_atom"] = (record["energies"] - curve) * 500
        np.testing.assert_allclose(np.sqrt(np.mean((record["energies"] - curve) ** 2)),
                                   fit["rms_residual_eV"], rtol=0, atol=2e-12)
    candidate, reference = data["eos"]["idmate"], data["eos"]["vasp"]
    cf, rf = candidate["fit"], reference["fit"]
    lower = max(candidate["volumes"].min() / cf["V0_A3"], reference["volumes"].min() / rf["V0_A3"])
    upper = min(candidate["volumes"].max() / cf["V0_A3"], reference["volumes"].max() / rf["V0_A3"])
    require(lower < upper, "Nonoverlapping scaled energy-volume ranges")
    np.testing.assert_allclose([lower, upper], data["eos_shape"]["scaled_volume_range"], rtol=0, atol=2e-12)
    shape_differences = []
    for index in range(41):
        scaled = lower + (upper - lower) * index / 40
        cenergy = bm3_energy(scaled * cf["V0_A3"], cf["E0_eV"], cf["V0_A3"], cf["B0_eV_A3"], cf["B_prime"]) - cf["E0_eV"]
        renergy = bm3_energy(scaled * rf["V0_A3"], rf["E0_eV"], rf["V0_A3"], rf["B0_eV_A3"], rf["B_prime"]) - rf["E0_eV"]
        shape_differences.append((cenergy - renergy) * 500)
    np.testing.assert_allclose(np.sqrt(np.mean(np.array(shape_differences) ** 2)),
                               data["eos_shape"]["shape_rms_meV_per_atom"], rtol=0, atol=2e-12)
    np.testing.assert_allclose(max(abs(value) for value in shape_differences),
                               data["eos_shape"]["shape_max_meV_per_atom"], rtol=0, atol=2e-12)
    bands = data["bands"]
    for key in ("cx", "rx", "candidate", "reference", "candidate_grid_ang_inv",
                "reference_grid_ang_inv", "reference_all"):
        bands[key] = np.array(bands[key])
    # Interpolate on the physical reciprocal-distance grid, not the display grid.
    def interpolate(x, y, point):
        tolerance = 1e-10
        require(x[0] - tolerance <= point <= x[-1] + tolerance, "Refusing band extrapolation")
        if math.isclose(point, x[0], rel_tol=0, abs_tol=tolerance):
            return y[0]
        if math.isclose(point, x[-1], rel_tol=0, abs_tol=tolerance):
            return y[-1]
        upper = bisect_right(x, point)
        lower = upper - 1
        weight = (point - x[lower]) / (x[upper] - x[lower])
        return y[lower] + weight * (y[upper] - y[lower])
    paired = np.array([[interpolate(bands["reference_grid_ang_inv"], band, x)
                        for x in bands["candidate_grid_ang_inv"]]
                       for band in bands["reference_all"]])
    bands["difference"] = bands["candidate"] - paired
    np.testing.assert_allclose(np.abs(bands["difference"]).mean(), bands["metrics"]["mae_eV"], rtol=0, atol=1e-14)
    np.testing.assert_allclose(np.abs(bands["difference"]).max(), bands["metrics"]["max_error_eV"], rtol=0, atol=1e-14)
    return data


def load_supplement():
    data = read_json("supplement.json")
    data["radii"] = {float(key): value for key, value in data["radii"].items()}
    data["ak"], data["dp"] = np.array(data["ak"]), np.array(data["dp"])
    slope, intercept = np.polyfit(np.log(data["ak"]), np.log(data["dp"]), 1)
    np.testing.assert_allclose([slope, intercept], [data["slope"], data["intercept"]], rtol=0, atol=1e-12)
    data["slope"], data["intercept"] = slope, intercept
    return data

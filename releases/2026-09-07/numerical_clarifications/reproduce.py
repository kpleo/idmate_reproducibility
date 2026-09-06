#!/usr/bin/env python3
"""Reconstruct the supplied IDMate numerical comparisons using Python only."""

import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
MATERIALS = ("Si", "graphene", "Al")
CONFIGURATIONS = ("exploratory", "selection", "frozen")
FAMILIES = ("low_cutoff", "target_subspace", "coarse_k")
INPUTS = (
    "historical_candidates.csv", "material_candidates.csv",
    "excluded_candidates.csv", "shell_profiles.csv",
    "shell_probe_scattering.csv", "shell_configurations.json",
    "h2_reuse.csv", "source_coverage.json",
)
THRESHOLD = 0.05


def require(condition, message):
    if not condition:
        raise ValueError(message)


def number(value, nonnegative=True):
    require(not isinstance(value, bool), "Boolean is not a numerical observation")
    parsed = float(value)
    require(math.isfinite(parsed), "Nonfinite numerical observation")
    require(not nonnegative or parsed >= 0.0, "Negative nonnegative observation")
    return parsed


def integer(value):
    parsed = number(value)
    require(parsed.is_integer(), "Noninteger count or index")
    return int(parsed)


def boolean(value):
    require(value in ("true", "false"), "Boolean must be true or false")
    return value == "true"


def table(directory, filename, numeric=(), integers=(), booleans=(), signed=(), optional=()):
    with (directory / filename).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames and len(set(reader.fieldnames)) == len(reader.fieldnames),
                "Missing or duplicate column names")
        rows = list(reader)
    require(rows, "Empty numerical table")
    for row in rows:
        require(None not in row and all(value is not None for value in row.values()),
                "Uneven numerical table")
        for key in numeric:
            row[key] = number(row[key])
        for key in signed:
            row[key] = number(row[key], nonnegative=False)
        for key in integers:
            row[key] = integer(row[key])
        for key in booleans:
            row[key] = boolean(row[key])
        for key in optional:
            row[key] = number(row[key]) if row[key] else None
    return rows


def unique(rows, keys):
    identities = [tuple(row[key] for key in keys) for row in rows]
    require(len(set(identities)) == len(identities), "Duplicate observation identity")


def verify_inputs(directory):
    manifest = json.loads((directory / "input_manifest.json").read_text())
    require(set(manifest) == set(INPUTS), "Unexpected input manifest entries")
    for name in INPUTS:
        payload = (directory / name).read_bytes()
        entry = manifest[name]
        require(len(payload) == entry["bytes"], "Input size mismatch: " + name)
        require(hashlib.sha256(payload).hexdigest() == entry["sha256"],
                "Input checksum mismatch: " + name)


def load_inputs(directory=HERE):
    directory = Path(directory)
    verify_inputs(directory)
    distances = ("raw_window_bound", "window_distance", "full_distance")
    history = table(directory, "historical_candidates.csv",
                    numeric=distances + ("trace_deviation",), integers=("iteration",))
    candidates = table(directory, "material_candidates.csv", numeric=distances,
                       integers=("iteration", "observation_order"))
    excluded = table(directory, "excluded_candidates.csv", numeric=("full_distance",),
                     integers=("iteration", "observation_order"))
    shells = table(directory, "shell_profiles.csv", integers=("shell_key",),
                   numeric=("radius_bohr_inv", "energy_share", "cumulative_share"))
    scattering = table(directory, "shell_probe_scattering.csv",
                       integers=("input_shell_key", "columns"), numeric=("mean_self_shell_fraction",))
    flags = ("retained_window_open", "occupied_window_open", "residual_converged",
             "density_within_tolerance", "orthogonality_within_tolerance",
             "provisional_reuse_passed")
    h2 = table(directory, "h2_reuse.csv", integers=("iteration", "anchor_iteration"),
               numeric=("potential_change_linf_ha", "adiabatic_number_squared", "projector_distance"),
               signed=("anchor_retained_edge_ha", "anchor_guard_ha", "retained_margin_ha",
                       "occupied_margin_ha", "retained_candidate_separation_ha",
                       "occupied_candidate_separation_ha"), booleans=flags,
               optional=("retained_angle_bound",))
    for rows in (history, candidates, excluded):
        unique(rows, ("candidate_id",))
        require(all(r["material"] in MATERIALS and r["proposal_family"] in FAMILIES for r in rows),
                "Unknown material or proposal family")
    combined = candidates + excluded
    unique(combined, ("candidate_id",))
    unique(combined, ("run_id", "iteration"))
    require(sorted(r["observation_order"] for r in combined) == list(range(1, 67)),
            "Expected 66 ordered observations")
    require(all(r["decision"] in ("accept", "reject") for r in candidates),
            "Unexpected binary decision")
    require(all(r["decision"] == "abstain" for r in excluded), "Unexpected exclusion")
    candidates.sort(key=lambda r: r["observation_order"])
    unique(shells, ("configuration", "shell_key"))
    unique(scattering, ("configuration", "input_shell_key"))
    unique(h2, ("iteration",))
    configs = json.loads((directory / "shell_configurations.json").read_text())
    require({c["configuration"] for c in configs} == set(CONFIGURATIONS) and len(configs) == 3,
            "Expected three shell configurations")
    return {"history": history, "candidates": candidates, "excluded": excluded,
            "shells": shells, "scattering": scattering, "configs": configs, "h2": h2}


def auc_pairs(accepted, rejected):
    accepted = [number(value) for value in accepted]
    rejected = [number(value) for value in rejected]
    require(accepted and rejected, "AUC requires both decision groups")
    concordant = sum(r > a for a in accepted for r in rejected)
    tied = sum(r == a for a in accepted for r in rejected)
    pairs = len(accepted) * len(rejected)
    return {"accepted": len(accepted), "rejected": len(rejected), "pairs": pairs,
            "concordant": concordant, "tied": tied,
            "discordant": pairs - concordant - tied,
            "auc": (concordant + 0.5 * tied) / pairs}


def pooled_pairs(cells):
    totals = {key: sum(cell[key] for cell in cells)
              for key in ("pairs", "concordant", "tied", "discordant")}
    require(totals["pairs"] > 0, "No pairs to pool")
    totals["auc"] = (totals["concordant"] + 0.5 * totals["tied"]) / totals["pairs"]
    return totals


def quantile(ordered, q):
    require(ordered and 0.0 <= q <= 1.0, "Invalid quantile request")
    p = (len(ordered) - 1) * q
    lo, hi = math.floor(p), math.ceil(p)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (1.0 - (p - lo)) + ordered[hi] * (p - lo)


def candidate_bootstrap(accepted, rejected, resamples=10000, seed=20260829):
    resamples = integer(resamples)
    require(resamples > 0, "Bootstrap requires positive resample count")
    auc_pairs(accepted, rejected)
    rng = random.Random(seed)
    estimates = []
    na, nr = len(accepted), len(rejected)
    for _ in range(resamples):
        sample_a = [accepted[rng.randrange(na)] for _ in range(na)]
        sample_r = [rejected[rng.randrange(nr)] for _ in range(nr)]
        estimates.append(auc_pairs(sample_a, sample_r)["auc"])
    estimates.sort()
    return {"resamples": resamples, "seed": seed,
            "ci95_percentile": [quantile(estimates, 0.025), quantile(estimates, 0.975)],
            "resampling_unit": "candidate, independently within each decision group",
            "material_or_trajectory_cluster_resampling": False}


def history_analysis(rows):
    comparisons = []
    for row in rows:
        bound = row["raw_window_bound"]
        require(bound > 0, "Historical bound must be positive for its ratio")
        distance = row["window_distance"]
        excess = max(0.0, distance - bound)
        comparisons.append({**row, "window_over_raw_bound": distance / bound,
                            "positive_window_excess": excess,
                            "trace_covers_positive_excess": excess <= row["trace_deviation"],
                            "window_above_criterion": distance > THRESHOLD})
    return {"count": len(rows), "window_criterion": THRESHOLD,
            "raw_bound_exceedances": sum(r["positive_window_excess"] > 0 for r in comparisons),
            "window_criterion_exceedances": sum(r["window_above_criterion"] for r in comparisons),
            "full_distance_above_window_criterion": sum(r["full_distance"] > THRESHOLD for r in rows),
            "trace_covers_positive_excess_count": sum(r["trace_covers_positive_excess"] for r in comparisons),
            "largest_ratio_candidate": max(comparisons, key=lambda r: r["window_over_raw_bound"])["candidate_id"],
            "largest_positive_excess_candidate": max(comparisons, key=lambda r: r["positive_window_excess"])["candidate_id"],
            "comparisons": comparisons}


def material_analysis(rows, excluded):
    accepted = [r for r in rows if r["decision"] == "accept"]
    rejected = [r for r in rows if r["decision"] == "reject"]
    values_a = [r["full_distance"] for r in accepted]
    values_r = [r["full_distance"] for r in rejected]
    cells = []
    for material_a in MATERIALS:
        for material_r in MATERIALS:
            values_am = [r["full_distance"] for r in accepted if r["material"] == material_a]
            values_rm = [r["full_distance"] for r in rejected if r["material"] == material_r]
            cells.append({"accepted_material": material_a, "rejected_material": material_r,
                          **auc_pairs(values_am, values_rm)})
    same = [c for c in cells if c["accepted_material"] == c["rejected_material"]]
    different = [c for c in cells if c["accepted_material"] != c["rejected_material"]]
    total = auc_pairs(values_a, values_r)
    within, between = pooled_pairs(same), pooled_pairs(different)
    require(within["pairs"] + between["pairs"] == total["pairs"], "Pair decomposition failed")
    require(within["concordant"] + between["concordant"] == total["concordant"],
            "Concordance decomposition failed")
    deletion = {}
    for material in MATERIALS:
        deletion[material] = auc_pairs(
            [r["full_distance"] for r in accepted if r["material"] != material],
            [r["full_distance"] for r in rejected if r["material"] != material])
    decisions = {}
    family_decisions = {}
    for material in MATERIALS:
        decisions[material] = dict(Counter(r["decision"] for r in rows + excluded
                                            if r["material"] == material))
    for family in FAMILIES:
        family_decisions[family] = dict(Counter(r["decision"] for r in rows + excluded
                                                if r["proposal_family"] == family))
    return {"included": len(rows), "excluded_abstentions": len(excluded),
            "decisions_by_material": decisions, "decisions_by_proposal_family": family_decisions,
            "pooled": total, "same_material": within, "cross_material": between,
            "same_material_pair_fraction": within["pairs"] / total["pairs"],
            "equal_material_mean_auc": math.fsum(c["auc"] for c in same) / len(same),
            "per_material": {c["accepted_material"]: c for c in same},
            "material_pair_cells": cells, "deletion_of_material": deletion,
            "candidate_bootstrap": candidate_bootstrap(values_a, values_r)}


def regression_slope(xs, ys):
    require(len(xs) == len(ys) and len(xs) >= 2, "Regression needs paired observations")
    mx, my = math.fsum(xs) / len(xs), math.fsum(ys) / len(ys)
    denominator = math.fsum((x - mx) ** 2 for x in xs)
    require(denominator > 0, "Constant regression predictor")
    return math.fsum((x - mx) * (y - my) for x, y in zip(xs, ys)) / denominator


def shell_analysis(rows, configs, scattering):
    require({r["configuration"] for r in rows} == set(CONFIGURATIONS), "Unknown shell configuration")
    summaries = {}
    for config in configs:
        name = config["configuration"]
        group = sorted((r for r in rows if r["configuration"] == name), key=lambda r: r["shell_key"])
        require(len(group) == config["printed_shells"], "Printed shell count mismatch")
        require(len(group) < config["total_shells_reported"], "Expected partially printed profile")
        cumulative = 0.0
        running_share = 0.0
        for row in group:
            require(row["shell_key"] > 0 and 0 <= row["energy_share"] <= 1, "Invalid shell share")
            require(cumulative <= row["cumulative_share"] <= 1, "Nonmonotone cumulative share")
            exact_radius = 2 * math.pi / config["lattice_constant_bohr"] * math.sqrt(row["shell_key"])
            require(abs(exact_radius - row["radius_bohr_inv"]) <= 5.01e-5,
                    "Inconsistent printed shell radius")
            running_share += row["energy_share"]
            require(abs(running_share - row["cumulative_share"]) <= 7e-7,
                    "Printed shares disagree beyond their available precision")
            cumulative = row["cumulative_share"]
        crossings = {}
        for target in (0.9, 0.99):
            matching = [r for r in group if r["cumulative_share"] >= target]
            require(matching, "Unreached shell coverage target")
            first = matching[0]
            crossings[f"{target:.2f}"] = {k: first[k] for k in
                                         ("shell_key", "radius_bohr_inv", "cumulative_share")}
        selected = {str(r["shell_key"]): r["cumulative_share"] for r in group
                    if r["shell_key"] in (11, 24, 27)}
        require(len(selected) == 3, "Missing selected shell")
        summaries[name] = {"derivative_columns": config["derivative_columns"],
                           "printed_shells": len(group), "total_shells_reported": config["total_shells_reported"],
                           "crossings": crossings, "cumulative_at_selected_shells": selected}
        if name == "exploratory":
            fit = [r for r in group if r["energy_share"] > 1e-8]
            xs = [2 * math.pi / config["lattice_constant_bohr"] * math.sqrt(r["shell_key"]) for r in fit]
            ys = [math.log10(r["energy_share"]) for r in fit]
            summaries[name]["decay_fit"] = {"shells": len(fit), "share_threshold": 1e-8,
                "slope_decades_per_bohr_inverse": regression_slope(xs, ys)}
    require(all(r["configuration"] == "exploratory" and 0 <= r["mean_self_shell_fraction"] <= 1
                for r in scattering), "Invalid self-shell observation")
    shares = [r["mean_self_shell_fraction"] for r in scattering]
    return {"configurations": summaries, "exploratory_self_shell": {
        "input_shell_groups": len(scattering), "columns": sum(r["columns"] for r in scattering),
        "minimum_mean_fraction": min(shares), "maximum_mean_fraction": max(shares),
        "monotonically_decreasing_with_input_shell": all(b <= a for a, b in zip(shares, shares[1:]))}}


def reuse_decision(row):
    return (row["residual_converged"] and row["retained_margin_ha"] > 0
            and row["retained_candidate_separation_ha"] > 0 and row["retained_angle_bound"] is not None
            and row["occupied_margin_ha"] > 0 and row["occupied_candidate_separation_ha"] > 0
            and row["density_within_tolerance"] and row["orthogonality_within_tolerance"])


def h2_analysis(rows):
    rows = sorted(rows, key=lambda r: r["iteration"])
    require([r["iteration"] for r in rows] == list(range(2, 62)), "Expected consecutive H2 attempts")
    open_iterations, passed_iterations, open_failed = [], [], []
    margin_errors, ratios = [], []
    for row in rows:
        computed = row["anchor_guard_ha"] - row["anchor_retained_edge_ha"] - 2 * row["potential_change_linf_ha"]
        margin_errors.append(abs(computed - row["retained_margin_ha"]))
        require(margin_errors[-1] <= 3e-12, "Margin inconsistent with rounded anchor energies")
        require(row["retained_window_open"] == (row["retained_margin_ha"] > 0), "Retained margin decision differs")
        require(row["occupied_window_open"] == (row["occupied_margin_ha"] > 0), "Occupied margin decision differs")
        passed = reuse_decision(row)
        require(passed == row["provisional_reuse_passed"], "Composite reuse decision differs")
        if row["retained_window_open"]:
            open_iterations.append(row["iteration"])
            if not passed:
                open_failed.append(row)
        if passed:
            passed_iterations.append(row["iteration"])
        require(row["adiabatic_number_squared"] > 0 and row["projector_distance"] > 0,
                "Positive parameter and distance needed for the perturbative comparison")
        ratios.append(row["projector_distance"] / math.sqrt(row["adiabatic_number_squared"]))
    errors = [abs(ratio - 1) for ratio in ratios]
    errors_over_distance = [abs(1 - 1 / ratio) for ratio in ratios]
    xs = [math.log10(math.sqrt(r["adiabatic_number_squared"])) for r in rows]
    ys = [math.log10(r["projector_distance"]) for r in rows]
    return {"attempts": len(rows), "retained_window_open": len(open_iterations),
            "occupied_window_open": sum(r["occupied_window_open"] for r in rows),
            "provisional_reuse_passed": len(passed_iterations),
            "retained_open_iterations": open_iterations, "reuse_passed_iterations": passed_iterations,
            "retained_open_but_reuse_failed": len(open_failed),
            "open_failures_with_unconverged_residual": sum(not r["residual_converged"] for r in open_failed),
            "maximum_margin_reconstruction_difference_ha": max(margin_errors),
            "separate_perturbative_comparison": {"within_relative_error_0p10": sum(e < 0.1 for e in errors),
                "relative_error_0p10_denominator": "A",
                "mean_distance_over_A": math.fsum(ratios) / len(ratios),
                "mean_relative_error_over_A": math.fsum(errors) / len(errors),
                "mean_relative_error_over_distance": math.fsum(errors_over_distance) / len(errors_over_distance),
                "maximum_relative_error_over_A": max(errors),
                "log_distance_vs_log_A_slope": regression_slope(xs, ys)}}


def calculate(data):
    return {"format_version": 1,
            "historical_candidates": history_analysis(data["history"]),
            "material_separation": material_analysis(data["candidates"], data["excluded"]),
            "silicon_shells": shell_analysis(data["shells"], data["configs"], data["scattering"]),
            "h2_reuse": h2_analysis(data["h2"])}


def fixed(value):
    return f"{value:.3f}"


def scientific(value):
    return f"{value:.3e}"


def ratio_display(value):
    return scientific(value) if 0 < abs(value) < 0.001 else fixed(value)


def equivalent(actual, expected):
    """Allow tiny library-level differences in derived floating-point values."""
    if type(actual) is not type(expected):
        return False
    if isinstance(actual, dict):
        return actual.keys() == expected.keys() and all(equivalent(actual[k], expected[k]) for k in actual)
    if isinstance(actual, list):
        return len(actual) == len(expected) and all(equivalent(a, b) for a, b in zip(actual, expected))
    if isinstance(actual, float):
        return (math.isfinite(actual) and math.isfinite(expected)
                and math.isclose(actual, expected, rel_tol=2e-13, abs_tol=1e-300))
    return actual == expected


def markdown(results):
    history = results["historical_candidates"]
    material = results["material_separation"]
    shells = results["silicon_shells"]
    h2 = results["h2_reuse"]
    lines = ["# IDMate numerical comparisons", "", "## Historical window bounds", "",
        f"{history['raw_bound_exceedances']}/{history['count']} distances exceed the raw bound; "
        f"{history['window_criterion_exceedances']}/{history['count']} exceed the {fixed(THRESHOLD)} window criterion.",
        "Every positive excess is smaller than that row's recorded trace deviation.",
        "This empirical comparison is not a reconstruction of a trace-corrected bound.", "",
        "| Candidate | Raw bound | Window distance | Full distance | Trace deviation | Window / raw bound |",
        "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for row in history["comparisons"]:
        values = [scientific(row[k]) for k in ("raw_window_bound", "window_distance", "full_distance", "trace_deviation")]
        lines.append("| " + " | ".join([row["candidate_id"], *values, ratio_display(row["window_over_raw_bound"])]) + " |")
    lines.extend(["", "## Candidate separation", "",
        f"{material['pooled']['accepted']} accepted, {material['pooled']['rejected']} rejected; "
        f"{material['excluded_abstentions']} abstentions excluded.", "",
        "| Pair set | Concordant | Tied | Total pairs | AUC |",
        "| --- | ---: | ---: | ---: | ---: |"])
    for label, group in [("All", material["pooled"]), ("Same material", material["same_material"]),
                         ("Cross material", material["cross_material"])]:
        lines.append(f"| {label} | {group['concordant']} | {group['tied']} | {group['pairs']} | {fixed(group['auc'])} |")
    lines.extend(["", "| Accepted material | Rejected material | Concordant / total | AUC |",
                  "| --- | --- | ---: | ---: |"])
    for cell in material["material_pair_cells"]:
        lines.append(f"| {cell['accepted_material']} | {cell['rejected_material']} | "
                     f"{cell['concordant']} / {cell['pairs']} | {fixed(cell['auc'])} |")
    lines.extend(["", "| Material | Within-material AUC | AUC after deleting this material |",
                  "| --- | ---: | ---: |"])
    for name in MATERIALS:
        lines.append(f"| {name} | {fixed(material['per_material'][name]['auc'])} | "
                     f"{fixed(material['deletion_of_material'][name]['auc'])} |")
    lo, hi = material["candidate_bootstrap"]["ci95_percentile"]
    lines.extend(["", f"Equal-material mean AUC: {fixed(material['equal_material_mean_auc'])}.",
        f"Candidate-level bootstrap 95% interval: [{fixed(lo)}, {fixed(hi)}].",
        "The interval does not use material or trajectory clusters. Deletion is a sensitivity",
        "calculation, not held-out prediction. Decision and proposal family are confounded.", "",
        "| Proposal family | Accept | Reject | Abstain |", "| --- | ---: | ---: | ---: |"])
    for family in FAMILIES:
        counts = material["decisions_by_proposal_family"][family]
        lines.append(f"| {family} | {counts.get('accept', 0)} | {counts.get('reject', 0)} | {counts.get('abstain', 0)} |")
    lines.extend(["", "## Silicon shell profiles", "",
        "The three calculations differ in settings and number of derivative columns.", "",
        "| Configuration | Columns | Share at S=11 | Share at S=24 | Share at S=27 | 90% radius | 99% radius |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"])
    for name in CONFIGURATIONS:
        group = shells["configurations"][name]
        shares = [fixed(group["cumulative_at_selected_shells"][s]) for s in ("11", "24", "27")]
        radii = [fixed(group["crossings"][target]["radius_bohr_inv"]) for target in ("0.90", "0.99")]
        lines.append("| " + " | ".join([name, str(group["derivative_columns"]), *shares, *radii]) + " |")
    fit = shells["configurations"]["exploratory"]["decay_fit"]
    self_shell = shells["exploratory_self_shell"]
    lines.extend(["", "Radii are in inverse bohr. Rounded printed profiles do not reproduce full-matrix closure.",
        f"Exploratory log-share slope: {fixed(fit['slope_decades_per_bohr_inverse'])} decades per inverse bohr "
        f"({fit['shells']} shells with share above {scientific(fit['share_threshold'])}).",
        f"Exploratory mean self-shell fractions span {fixed(self_shell['minimum_mean_fraction'])} to "
        f"{fixed(self_shell['maximum_mean_fraction'])}; they do not decrease monotonically.", "",
        "## H2 reuse conditions", "",
        "| Condition | Number of attempts |", "| --- | ---: |",
        f"| Retained-state spectral margin positive | {h2['retained_window_open']} / {h2['attempts']} |",
        f"| Occupied-state spectral margin positive | {h2['occupied_window_open']} / {h2['attempts']} |",
        f"| Full provisional reuse test passed | {h2['provisional_reuse_passed']} / {h2['attempts']} |",
        f"| Separate first-order comparison within 10% | {h2['separate_perturbative_comparison']['within_relative_error_0p10']} / {h2['attempts']} |", "",
        f"Retained windows are open at iterations {h2['retained_open_iterations'][0]}-{h2['retained_open_iterations'][-1]}; "
        f"reuse passes occur at {h2['reuse_passed_iterations'][0]}-{h2['reuse_passed_iterations'][-1]}.",
        f"The {h2['retained_open_but_reuse_failed']} open-window failures all have unconverged candidate residuals.",
        "The perturbative comparison is distinct from the provisional reuse test and is not a decision condition.", ""])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Compare fresh results with saved files without writing")
    args = parser.parse_args()
    results = calculate(load_inputs())
    text = markdown(results)
    if args.check:
        saved = json.loads((HERE / "results.json").read_text())
        require(equivalent(saved, results), "Saved numerical results differ")
        require((HERE / "results.md").read_text() == text, "Saved displayed results differ")
        print("Input checksums and independently recomputed numerical results agree.")
    else:
        (HERE / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
        (HERE / "results.md").write_text(text)
        print("Wrote results.json and results.md.")


if __name__ == "__main__":
    main()

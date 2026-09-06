"""Numerical and portability tests for the supplied comparisons."""

import copy
import json
import shutil
import subprocess
import sys
import tempfile
import unittest

import reproduce as r


class NumericalResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = r.load_inputs()
        cls.results = r.calculate(cls.data)

    def test_historical_counts_and_distinct_criteria(self):
        history = self.results["historical_candidates"]
        self.assertEqual(history["count"], 10)
        self.assertEqual(history["raw_bound_exceedances"], 6)
        self.assertEqual(history["window_criterion_exceedances"], 0)
        self.assertEqual(history["full_distance_above_window_criterion"], 7)
        self.assertEqual(history["trace_covers_positive_excess_count"], 10)

    def test_same_row_trace_comparison(self):
        rows = {row["candidate_id"]: row for row in self.results["historical_candidates"]["comparisons"]}
        silicon = rows["Si_2"]
        self.assertAlmostEqual(silicon["positive_window_excess"], 7.011837845739798e-10, places=24)
        self.assertEqual(silicon["trace_deviation"], 2.3283943662022466e-9)
        self.assertLess(silicon["positive_window_excess"], silicon["trace_deviation"])

    def test_largest_ratio_and_largest_excess_are_different_rows(self):
        history = self.results["historical_candidates"]
        self.assertEqual(history["largest_ratio_candidate"], "Al_2")
        self.assertEqual(history["largest_positive_excess_candidate"], "Al_5")
        rows = {row["candidate_id"]: row for row in history["comparisons"]}
        self.assertEqual(rows["Al_2"]["window_over_raw_bound"], 5.41425895208607)
        self.assertEqual(rows["Al_5"]["positive_window_excess"], 1.4293377219514513e-9)

    def test_ratio_is_computed_before_display_rounding(self):
        row = next(row for row in self.data["history"] if row["candidate_id"] == "Al_2")
        numerator = r.scientific(row["window_distance"])
        denominator = r.scientific(row["raw_window_bound"])
        self.assertEqual((numerator, denominator), ("1.751e-09", "3.234e-10"))
        self.assertEqual(r.fixed(float(numerator) / float(denominator)), "5.414")
        self.assertEqual(r.fixed(row["window_distance"] / row["raw_window_bound"]), "5.414")

    def test_binary_cohort_and_exclusions(self):
        material = self.results["material_separation"]
        self.assertEqual((material["included"], material["excluded_abstentions"]), (58, 8))
        self.assertEqual((material["pooled"]["accepted"], material["pooled"]["rejected"]), (40, 18))
        self.assertTrue(all(row["material"] == "graphene" for row in self.data["excluded"]))
        self.assertEqual(material["decisions_by_material"], {
            "Si": {"accept": 16, "reject": 6},
            "graphene": {"accept": 8, "reject": 6, "abstain": 8},
            "Al": {"accept": 16, "reject": 6}})

    def test_pooled_pairs(self):
        pooled = self.results["material_separation"]["pooled"]
        self.assertEqual((pooled["concordant"], pooled["tied"], pooled["discordant"], pooled["pairs"]),
                         (500, 0, 220, 720))
        self.assertEqual(pooled["auc"], 500 / 720)

    def test_pair_decomposition(self):
        material = self.results["material_separation"]
        same, cross = material["same_material"], material["cross_material"]
        self.assertEqual((same["concordant"], same["pairs"], same["auc"]), (156, 240, 0.65))
        self.assertEqual((cross["concordant"], cross["pairs"], cross["auc"]), (344, 480, 344 / 480))
        self.assertAlmostEqual((same["auc"] + 2 * cross["auc"]) / 3, material["pooled"]["auc"], places=15)
        self.assertEqual(material["same_material_pair_fraction"], 1 / 3)

    def test_all_nine_material_pair_cells(self):
        cells = self.results["material_separation"]["material_pair_cells"]
        expected = ((48, 96), (80, 96), (96, 96), (0, 48), (12, 48), (24, 48),
                    (48, 96), (96, 96), (96, 96))
        self.assertEqual([(c["concordant"], c["pairs"]) for c in cells], list(expected))
        self.assertTrue(all(cell["tied"] == 0 for cell in cells))

    def test_per_material_and_equal_material_summary(self):
        material = self.results["material_separation"]
        self.assertEqual([material["per_material"][name]["auc"] for name in r.MATERIALS], [0.5, 0.25, 1.0])
        self.assertEqual(material["equal_material_mean_auc"], 1.75 / 3)
        self.assertNotEqual(material["equal_material_mean_auc"], material["same_material"]["auc"])

    def test_material_deletions(self):
        deletion = self.results["material_separation"]["deletion_of_material"]
        self.assertEqual([deletion[name]["auc"] for name in r.MATERIALS], [0.7916666666666666, 0.75, 0.4861111111111111])
        self.assertEqual([deletion[name]["pairs"] for name in r.MATERIALS], [288, 384, 288])

    def test_candidate_bootstrap_reproduces_saved_interval(self):
        boot = self.results["material_separation"]["candidate_bootstrap"]
        self.assertEqual(boot["ci95_percentile"], [0.5444444444444444, 0.8305555555555556])
        self.assertEqual((boot["seed"], boot["resamples"]), (20260829, 10000))
        self.assertFalse(boot["material_or_trajectory_cluster_resampling"])

    def test_proposal_family_confounding(self):
        family = self.results["material_separation"]["decisions_by_proposal_family"]
        self.assertEqual(family, {"low_cutoff": {"accept": 22, "abstain": 2},
                                  "target_subspace": {"accept": 18, "abstain": 6},
                                  "coarse_k": {"reject": 18}})

    def test_historical_observations_are_overlapping_not_extra_samples(self):
        multi = {(row["material"], row["iteration"]): row for row in self.data["candidates"]
                 if row["run_id"].endswith("_multiple") and row["decision"] == "accept"}
        self.assertEqual(len(multi), 10)
        for row in self.data["history"]:
            other = multi[row["material"], row["iteration"]]
            for key in ("raw_window_bound", "window_distance", "full_distance", "proposal_family"):
                self.assertEqual(row[key], other[key])

    def test_shell_values_and_configuration_specific_crossings(self):
        configurations = self.results["silicon_shells"]["configurations"]
        expected = {"exploratory": (128, 0.927397, 0.991931, 0.996323, 24),
                    "selection": (136, 0.908607, 0.987312, 0.992773, 27),
                    "frozen": (168, 0.904053, 0.982448, 0.992685, 27)}
        for name, (columns, s11, s24, s27, crossing) in expected.items():
            group = configurations[name]
            self.assertEqual(group["derivative_columns"], columns)
            self.assertEqual(group["cumulative_at_selected_shells"], {"11": s11, "24": s24, "27": s27})
            self.assertEqual(group["crossings"]["0.90"]["shell_key"], 11)
            self.assertEqual(group["crossings"]["0.99"]["shell_key"], crossing)

    def test_printed_shells_do_not_cover_all_reported_shells(self):
        configs = self.results["silicon_shells"]["configurations"]
        self.assertEqual([(configs[name]["printed_shells"], configs[name]["total_shells_reported"])
                          for name in r.CONFIGURATIONS], [(37, 99), (142, 715), (142, 715)])

    def test_exploratory_decay_fit_uses_34_shells(self):
        fit = self.results["silicon_shells"]["configurations"]["exploratory"]["decay_fit"]
        self.assertEqual(fit["shells"], 34)
        self.assertEqual(r.fixed(fit["slope_decades_per_bohr_inverse"]), "-1.529")

    def test_corrected_self_shell_fractions(self):
        scattering = self.results["silicon_shells"]["exploratory_self_shell"]
        self.assertEqual((scattering["input_shell_groups"], scattering["columns"]), (9, 128))
        self.assertEqual(scattering["minimum_mean_fraction"], 0.140539)
        self.assertEqual(scattering["maximum_mean_fraction"], 0.589853)
        self.assertFalse(scattering["monotonically_decreasing_with_input_shell"])

    def test_h2_distinct_conditions(self):
        h2 = self.results["h2_reuse"]
        self.assertEqual((h2["attempts"], h2["retained_window_open"], h2["provisional_reuse_passed"]), (60, 54, 29))
        self.assertEqual(h2["occupied_window_open"], 60)
        self.assertEqual(h2["retained_open_iterations"], list(range(8, 62)))
        self.assertEqual(h2["reuse_passed_iterations"], list(range(33, 62)))
        self.assertEqual(h2["retained_open_but_reuse_failed"], 25)
        self.assertEqual(h2["open_failures_with_unconverged_residual"], 25)

    def test_h2_perturbative_comparison_is_not_reuse_pass_count(self):
        h2 = self.results["h2_reuse"]
        self.assertEqual(h2["separate_perturbative_comparison"]["within_relative_error_0p10"], 60)
        self.assertNotEqual(h2["provisional_reuse_passed"], 60)
        self.assertLess(h2["maximum_margin_reconstruction_difference_ha"], 3e-12)

    def test_h2_relative_error_denominators_are_distinct(self):
        comparison = self.results["h2_reuse"]["separate_perturbative_comparison"]
        self.assertEqual(comparison["relative_error_0p10_denominator"], "A")
        self.assertAlmostEqual(comparison["mean_relative_error_over_A"], 0.005450149650073216, places=15)
        self.assertAlmostEqual(comparison["mean_relative_error_over_distance"], 0.005399816582960773, places=15)

    def test_each_reuse_requirement_is_necessary(self):
        row = next(row for row in self.data["h2"] if row["provisional_reuse_passed"])
        changes = {"residual_converged": False, "retained_margin_ha": 0.0,
                   "retained_candidate_separation_ha": 0.0, "retained_angle_bound": None,
                   "occupied_margin_ha": 0.0, "occupied_candidate_separation_ha": 0.0,
                   "density_within_tolerance": False, "orthogonality_within_tolerance": False}
        for key, value in changes.items():
            with self.subTest(condition=key):
                modified = dict(row, **{key: value})
                self.assertFalse(r.reuse_decision(modified))

    def test_h2_inconsistent_record_is_rejected(self):
        rows = copy.deepcopy(self.data["h2"])
        rows[0]["provisional_reuse_passed"] = True
        with self.assertRaisesRegex(ValueError, "Composite reuse decision differs"):
            r.h2_analysis(rows)

    def test_invalid_shell_radius_is_rejected(self):
        rows = copy.deepcopy(self.data["shells"])
        rows[0]["radius_bohr_inv"] += 0.01
        with self.assertRaisesRegex(ValueError, "shell radius"):
            r.shell_analysis(rows, self.data["configs"], self.data["scattering"])

    def test_previous_summary_checks_are_independently_reproduced(self):
        coverage = json.loads((r.HERE / "source_coverage.json").read_text())
        original = coverage["material_candidates"]["previous_summary_comparison"]
        result = self.results["material_separation"]
        self.assertEqual(original["auc_full"], result["pooled"]["auc"])
        self.assertEqual(original["discordant_pairs"], result["pooled"]["discordant"])
        self.assertEqual(original["bootstrap"]["ci95_percentile"], result["candidate_bootstrap"]["ci95_percentile"])
        self.assertEqual(sum(row["matching_numeric_fields"] for row in coverage["historical_candidates"]["comparisons"]), 40)


class PrimitiveChecks(unittest.TestCase):
    def test_ties_receive_half_weight(self):
        result = r.auc_pairs([1, 2], [2, 3])
        self.assertEqual((result["concordant"], result["tied"], result["pairs"], result["auc"]), (3, 1, 4, 0.875))

    def test_all_ties_and_orientation(self):
        self.assertEqual(r.auc_pairs([2, 2], [2])["auc"], 0.5)
        self.assertEqual(r.auc_pairs([1], [2])["auc"], 1.0)
        self.assertEqual(r.auc_pairs([2], [1])["auc"], 0.0)

    def test_permutation_does_not_change_auc(self):
        self.assertEqual(r.auc_pairs([3, 1, 2], [5, 2]), r.auc_pairs([2, 3, 1], [2, 5]))

    def test_empty_decision_groups_are_rejected(self):
        for accepted, rejected in (([], [1]), ([1], []), ([], [])):
            with self.assertRaises(ValueError):
                r.auc_pairs(accepted, rejected)

    def test_nonfinite_and_negative_values_are_rejected(self):
        for value in (float("nan"), float("inf"), -1, True):
            with self.subTest(value=value), self.assertRaises(ValueError):
                r.auc_pairs([value], [1])

    def test_integer_and_boolean_validation(self):
        for value in (1.5, True, -1):
            with self.assertRaises(ValueError):
                r.integer(value)
        self.assertTrue(r.boolean("true"))
        self.assertFalse(r.boolean("false"))
        with self.assertRaises(ValueError):
            r.boolean("yes")

    def test_duplicate_identities_are_rejected(self):
        with self.assertRaises(ValueError):
            r.unique([{"candidate_id": "a"}, {"candidate_id": "a"}], ("candidate_id",))

    def test_quantile_interpolation(self):
        self.assertEqual(r.quantile([0, 10], 0.25), 2.5)
        self.assertEqual(r.quantile([7], 0.99), 7)
        self.assertEqual(r.quantile([0, 10], 1.0), 10)

    def test_small_bootstrap_is_repeatable(self):
        self.assertEqual(r.candidate_bootstrap([1, 2], [2, 3], 50, 8),
                         r.candidate_bootstrap([1, 2], [2, 3], 50, 8))
        with self.assertRaises(ValueError):
            r.candidate_bootstrap([1], [2], 0)

    def test_three_decimal_mantissa_does_not_zero_small_values(self):
        self.assertEqual(r.scientific(1.8906878275335727e-13), "1.891e-13")
        self.assertEqual(r.ratio_display(4.534767628888115e-8), "4.535e-08")
        self.assertEqual(r.fixed(0.6944444444444444), "0.694")

    def test_result_comparison_only_allows_roundoff(self):
        self.assertTrue(r.equivalent({"count": 10, "value": [0.5]}, {"count": 10, "value": [0.5 + 1e-15]}))
        self.assertFalse(r.equivalent(10, 11))
        self.assertFalse(r.equivalent(10, 10.0))
        self.assertFalse(r.equivalent(1e-9, 1.001e-9))
        self.assertFalse(r.equivalent([0.5], []))
        self.assertFalse(r.equivalent(float("nan"), float("nan")))


class PortableExecution(unittest.TestCase):
    def test_input_change_is_detected(self):
        with tempfile.TemporaryDirectory(prefix="input_check_", dir=r.HERE) as directory:
            root = r.Path(directory)
            for name in (*r.INPUTS, "input_manifest.json"):
                shutil.copyfile(r.HERE / name, root / name)
            path = root / "historical_candidates.csv"
            path.write_text(path.read_text().replace("Si_2,", "Si_9,", 1))
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                r.verify_inputs(root)

    def test_detached_execution_and_read_only_check(self):
        with tempfile.TemporaryDirectory(prefix="portable_", dir=r.HERE) as directory:
            root = r.Path(directory)
            package = root / "numerical_inputs"
            package.mkdir()
            for name in (*r.INPUTS, "input_manifest.json", "reproduce.py"):
                shutil.copyfile(r.HERE / name, package / name)
            generated = subprocess.run([sys.executable, "-B", str(package / "reproduce.py")],
                                       cwd=root, capture_output=True, text=True, timeout=60)
            self.assertEqual(generated.returncode, 0, generated.stderr)
            for name in ("results.json", "results.md"):
                self.assertEqual((package / name).read_bytes(), (r.HERE / name).read_bytes())
            before = {path.name: (path.read_bytes(), path.stat().st_mtime_ns) for path in package.iterdir()}
            checked = subprocess.run([sys.executable, "-B", str(package / "reproduce.py"), "--check"],
                                     cwd=root, capture_output=True, text=True, timeout=60)
            self.assertEqual(checked.returncode, 0, checked.stderr)
            after = {path.name: (path.read_bytes(), path.stat().st_mtime_ns) for path in package.iterdir()}
            self.assertEqual(before, after)


if __name__ == "__main__":
    unittest.main()

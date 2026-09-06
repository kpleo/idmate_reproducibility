"""Checks for the supplied material tables and their reconstruction."""

import copy
import json
import unittest

import reproduce


class MaterialTables(unittest.TestCase):
    def setUp(self):
        self.data = json.loads((reproduce.HERE / "materials.json").read_text())

    def test_all_materials_in_each_table(self):
        tables = reproduce.table_rows(self.data)
        self.assertEqual(len(tables), 3)
        self.assertTrue(all(len(rows) == 3 for rows in tables.values()))

    def test_missing_is_not_zero(self):
        self.assertEqual(reproduce.number(None), "--")
        self.assertEqual(reproduce.number(0), "$0.000$")
        self.assertEqual(reproduce.table_rows(self.data)["material_rows.tex"][1],
                         r"graphene & abstain & -- & -- & $0.095$ & $0.060$ \\")

    def test_small_values_remain_resolved(self):
        self.assertEqual(reproduce.number(3.4891214588154704e-3),
                         r"$3.489\times10^{-3}$")

    def test_nonfinite_is_rejected(self):
        for value in (float("nan"), float("inf"), -float("inf")):
            with self.assertRaises(ValueError):
                reproduce.number(value)

    def test_dropped_or_reordered_material_is_rejected(self):
        for rows in (self.data["rows"][:2], self.data["rows"][::-1]):
            altered = copy.deepcopy(self.data)
            altered["rows"] = rows
            with self.assertRaises(ValueError):
                reproduce.validate(altered)

    def test_wrong_map_count_is_rejected(self):
        self.data["reference_maps"] = 5
        with self.assertRaises(ValueError):
            reproduce.validate(self.data)

    def test_stopped_result_is_rejected(self):
        self.data["rows"][0]["stop_required"] = True
        with self.assertRaises(ValueError):
            reproduce.validate(self.data)

    def test_target_informed_result_is_rejected(self):
        self.data["rows"][0]["oracle"] = True
        with self.assertRaises(ValueError):
            reproduce.validate(self.data)

    def test_abstention_cannot_be_filled_with_zero(self):
        self.data["rows"][1]["window_bound"] = 0.0
        with self.assertRaises(ValueError):
            reproduce.validate(self.data)

    def test_mixing_identity_is_checked(self):
        self.data["rows"][0]["next_mixed_density_distance"] += 1e-4
        with self.assertRaises(ValueError):
            reproduce.validate(self.data)

    def test_reference_crosscheck_contains_every_material(self):
        report = json.loads((reproduce.HERE / "numerical_crosscheck.json").read_text())
        self.assertEqual({point["material"] for point in report["materials"]},
                         {"si", "graphene", "al"})
        self.assertEqual(sum(point["comparison_count"] for point in report["materials"]), 1135)
        self.assertTrue(all(point["failed_comparisons"] == 0 for point in report["materials"]))


if __name__ == "__main__":
    unittest.main()

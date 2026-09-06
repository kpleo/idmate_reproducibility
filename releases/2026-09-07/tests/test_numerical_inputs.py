"""Focused reconstruction and invalid-input checks for the supplied figures."""

import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "figures"))
import common
from crystal_motifs import geometry


class NumericalInputs(unittest.TestCase):
    def test_five_ratio_series_and_missing_values(self):
        arrays, retained, _, tail, _ = common.load_response()
        self.assertEqual({key:len(value) for key, value in arrays.items()},
                         {"e0": 9886, "e3_graphene": 164, "e3_al": 1614,
                          "gk_graphene": 308, "gk_al": 400})
        self.assertEqual(len(retained), 10)
        self.assertEqual(len(tail), 14)
        self.assertAlmostEqual(sum(row[1] for row in retained), 0.992685, places=6)

    def test_all_ten_candidates_and_six_uncorrected_exceedances(self):
        rows = common.load_rows()
        self.assertEqual(sum(row["window_over_bound"] > 1 for row in rows), 6)
        self.assertAlmostEqual(max(row["window_over_bound"] for row in rows), 5.41425895208607)

    def test_corrupted_ratio_is_rejected(self):
        rows = common.read_csv("screen_scope.csv")
        rows[0]["window_over_bound"] = "1.1"
        with patch.object(common, "read_csv", return_value=rows):
            with self.assertRaises(ValueError):
                common.load_rows()

    def test_nan_is_rejected(self):
        rows = common.read_csv("screen_scope.csv")
        rows[0]["full_error"] = "nan"
        with patch.object(common, "read_csv", return_value=rows):
            with self.assertRaises(ValueError):
                common.load_rows()

    def test_missing_historical_row_is_rejected(self):
        with patch.object(common, "read_csv", return_value=common.read_csv("screen_scope.csv")[:-1]):
            with self.assertRaises(ValueError):
                common.load_rows()

    def test_bands_recompute_all_155_differences(self):
        data = common.load_physics()
        self.assertEqual(data["bands"]["difference"].shape, (5, 31))
        self.assertEqual(data["bands"]["reference_all"].shape, (5, 78))
        self.assertEqual(data["bands"]["reference"].shape, (5, 39))
        self.assertAlmostEqual(np.abs(data["bands"]["difference"]).mean(),
                               data["bands"]["metrics"]["mae_eV"], places=13)

    def test_energy_fit_reconstruction(self):
        data = common.load_physics()
        for record in data["eos"].values():
            self.assertEqual(len(record["volumes"]), 5)
            self.assertAlmostEqual(np.max(np.abs(record["residual_meV_atom"])) / 500,
                                   record["fit"]["max_residual_eV"], places=12)

    def test_crystal_coordination(self):
        from crystal_motifs import periodic_neighbours
        for key, expected in (("si_scope", 4), ("graphene_scope", 3), ("al_scope", 12), ("si_bands", 4)):
            cell, points, bonds, frame, edges, nearest = geometry(key)
            self.assertEqual(periodic_neighbours(cell)[0], [expected] * len(cell["fractional"]))
            self.assertGreater(nearest, 0)

    def test_projector_slope_recomputed(self):
        data = common.load_supplement()
        self.assertEqual(len(data["ak"]), 60)
        self.assertAlmostEqual(data["slope"], 0.9990746935, places=9)

    def test_stored_interval_not_reinterpreted(self):
        data = common.load_supplement()["verdict"]
        self.assertLess(data["bootstrap"]["ci_lower"], 0)
        self.assertGreater(data["bootstrap"]["ci_upper"], 0)
        margin = data["g2_gate_evaluation"]["criteria"]["margin"]
        self.assertLess(margin["improvement"], margin["threshold"])


if __name__ == "__main__":
    unittest.main()

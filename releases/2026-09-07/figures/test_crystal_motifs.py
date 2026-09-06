import unittest
import numpy as np
from crystal_motifs import camera, geometry, periodic_neighbours, read_cell


class CrystalTests(unittest.TestCase):
    def test_camera_is_orthonormal(self):
        for azimuth, elevation in [(-32, 22), (-10, 55), (0, 0)]:
            rotation = camera(azimuth, elevation)
            np.testing.assert_allclose(rotation @ rotation.T, np.eye(3), atol=1e-14)
            self.assertAlmostEqual(np.linalg.det(rotation), 1)

    def test_coordination(self):
        for key, expected in [("si_scope", 4), ("si_bands", 4), ("graphene_scope", 3), ("al_scope", 12)]:
            cell = read_cell(key)
            self.assertEqual(periodic_neighbours(cell)[0], [expected]*len(cell["fractional"]))

    def test_nearest_lengths(self):
        for key, factor in [("si_scope", np.sqrt(3)/4), ("si_bands", np.sqrt(3)/4), ("al_scope", 1/np.sqrt(2))]:
            cell, _, _, _, _, nearest = geometry(key)
            self.assertAlmostEqual(nearest, cell["conventional_a"]*factor, places=10)
        cell, _, _, _, _, nearest = geometry("graphene_scope")
        self.assertAlmostEqual(nearest, cell["lattice"][0][0]/np.sqrt(3), places=8)

    def test_periodic_display_crops(self):
        for key, expected in [("si_scope", 18), ("si_bands", 18), ("al_scope", 14), ("graphene_scope", 24)]:
            _, points, bonds, _, _, nearest = geometry(key)
            self.assertEqual(len(points), expected)
            self.assertTrue(np.isfinite(points).all())
            self.assertEqual(len(points), len(np.unique(points, axis=0)))
            for i, j in bonds:
                self.assertAlmostEqual(np.linalg.norm(points[i]-points[j]), nearest, places=8)


if __name__ == "__main__":
    unittest.main(verbosity=2)

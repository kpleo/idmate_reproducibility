"""Numbers quoted in the manuscript, recomputed from the derived inputs of this release."""
import csv, json, math, os, unittest
import numpy as np
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')


def load(name):
    with open(os.path.join(DATA, name)) as f:
        return json.load(f)


class SavedHamiltonianDecomposition(unittest.TestCase):
    """Sec. IV C: 296 distinct candidates derived from the three saved Hamiltonians."""

    @classmethod
    def setUpClass(cls):
        o = load('decomposition_results.json')
        cls.rows = [r for m in o for kind in ('extra_sweep', 'tau_sweep') for r in o[m][kind]
                    if not r['duplicate_of_buffer_sweep']]

    def test_count(self):
        self.assertEqual(len(self.rows), 296)

    def test_leading_order_accuracy(self):
        rel = np.array([abs(r['pyth'] - r['d_full']) / r['d_full'] for r in self.rows])
        rel2 = np.array([abs(r['pyth2'] - r['d_full']) / r['d_full'] for r in self.rows])
        self.assertAlmostEqual(float(np.median(rel)), 1.78e-4, delta=0.01e-4)
        self.assertAlmostEqual(float(rel.max()), 0.1334, delta=1e-4)
        self.assertAlmostEqual(float(rel2.max()), 0.0170, delta=1e-4)
        self.assertEqual(int((rel > 0.01).sum()), 4)

    def test_indicator_and_bound(self):
        cold = [r for r in self.rows if r['tau'] <= 0.06]
        ii = np.array([r['indicator'] / r['d_full'] for r in cold])
        self.assertEqual(len(cold), 284)
        self.assertAlmostEqual(float(ii.min()), 1.0145, delta=1e-4)
        self.assertAlmostEqual(float(ii.max()), 3.467, delta=1e-3)
        self.assertTrue(all(r['d_full'] <= r['rigorous'] * (1 + 1e-12) for r in self.rows))


class InLoopRerun(unittest.TestCase):
    """Sec. IV D: every candidate presented in the SCF tests, rerun with saved Hamiltonians."""

    @classmethod
    def setUpClass(cls):
        cls.rows = load('inloop_results.json')
        cls.ritz = [r for r in cls.rows if r['ritz']]

    def test_counts(self):
        self.assertEqual(len(self.rows), 66)
        self.assertEqual(len(self.ritz), 48)
        self.assertEqual(sum(r['decision'] == 'accept' for r in self.rows), 40)
        self.assertEqual(sum(r['decision'] == 'reject' for r in self.rows), 18)
        self.assertTrue(all(not r['ritz'] for r in self.rows if r['decision'] == 'reject'))

    def test_reproduction(self):
        rel = max(abs(r['d_full_reproduced'] - r['d_full_recorded']) / r['d_full_recorded'] for r in self.rows)
        self.assertLess(rel, 4e-13)
        self.assertLess(max(r['ref_vs_exact_frame'] for r in self.rows), 5e-13)
        self.assertLess(max(r['ritz_offdiag_rel'] for r in self.ritz), 2e-15)
        self.assertLess(max(r['identity_residual'] for r in self.ritz), 3e-15)

    def test_window_misses_bad_accepts(self):
        acc = [r for r in self.ritz if r['decision'] == 'accept']
        self.assertEqual(sum(r['d_full'] > 0.05 for r in acc), 26)
        self.assertLessEqual(max(r['d_window'] for r in acc), 1.8e-9)

    def test_indicator_never_underestimates(self):
        ratios = np.array([r['indicator'] / r['d_full'] for r in self.ritz])
        self.assertTrue((ratios >= 1).all())
        self.assertAlmostEqual(float(ratios.min()), 1.183, delta=1e-3)
        self.assertAlmostEqual(float(ratios.max()), 11.66, delta=1e-2)
        self.assertAlmostEqual(float(np.median(ratios)), 1.537, delta=1e-3)

    def test_indicator_decisions(self):
        acc = [r for r in self.ritz if r['decision'] == 'accept']
        bad = [r for r in acc if r['d_full'] > 0.05]
        good = [r for r in acc if r['d_full'] <= 0.05]
        self.assertEqual(sum(r['indicator'] > 0.05 for r in bad), 26)
        self.assertEqual(sum(r['indicator'] <= 0.05 for r in good), 14)
        self.assertEqual(sum((r['indicator'] <= 0.05) == (r['d_full'] <= 0.05) for r in self.ritz), 48)
        self.assertAlmostEqual(min(r['indicator'] for r in bad), 0.1181, delta=1e-4)
        self.assertAlmostEqual(max(r['indicator'] for r in good), 0.04461, delta=1e-5)

    def test_rigorous_bound(self):
        ratio = [r['rigorous_computable'] / r['d_full'] for r in self.ritz]
        self.assertAlmostEqual(min(ratio), 30.57, delta=0.01)
        self.assertAlmostEqual(max(ratio), 330.0, delta=0.1)
        self.assertEqual(sum(r['rigorous_computable'] <= 0.05 for r in self.ritz), 7)

    def test_low_cutoff_coupling(self):
        low = [r for r in self.ritz if r['kind'] == 'low_cutoff_continuation']
        lead = [(r['pyth'] - r['d_full']) / r['d_full'] for r in low]
        self.assertEqual(len(low), 24)
        self.assertLess(max(lead), 0)
        self.assertAlmostEqual(min(lead), -0.2114, delta=1e-4)
        share = [r['blk_off'] ** 2 / r['d_full'] ** 2 for r in low]
        self.assertAlmostEqual(min(share), 0.623, delta=1e-3)


class RerunReproduction(unittest.TestCase):
    """Sec. IV D and Supplemental Sec. S6: the rerun compared with the original records."""

    @classmethod
    def setUpClass(cls):
        cls.p = load('rerun_provenance.json')
        cls.ref = cls.p['comparison_with_original']['reference_configurations']

    def test_runs(self):
        runs = self.p['runs']
        self.assertEqual(len(runs), 51)
        self.assertTrue(all(r['exit_code'] == 0 for r in runs))
        self.assertEqual(sum(r['insertion'] == 'single' for r in runs), 48)

    def test_single_insertion_bitwise(self):
        s = self.p['comparison_with_original']['single_insertion']
        self.assertEqual((s['runs'], s['arms_compared'], s['arms_bitwise_identical']), (48, 96, 96))

    def test_reference_configurations(self):
        rs = [r['reference_scoring'] for r in self.ref]
        ig = [r['integrated'] for r in self.ref]
        self.assertEqual([a['iterations']['rerun'] for a in rs], [45, 69, 44])
        for a in rs + ig:
            self.assertTrue(a['decisions_equal'])
            self.assertEqual(a['iterations']['rerun'], a['iterations']['original'])
            self.assertFalse(a['density_trajectory_bitwise'])
        self.assertLess(max(a['max_relative_difference']['full_state_distance'] for a in rs), 2e-11)
        self.assertLess(max(a['max_relative_difference']['full_state_distance'] for a in ig), 1.3e-8)
        self.assertAlmostEqual(max(a['max_relative_difference']['window_bound'] for a in rs), 0.138, delta=1e-3)
        self.assertLessEqual(max(a['free_energy_ha']['abs_difference'] for a in rs + ig), 4.5e-15)
        self.assertTrue(all(r['arms_identical_within_rerun'] for r in self.ref))


class CertificateVerification(unittest.TestCase):
    """Fig. 1(c): 12 372 defined distance/bound ratios, none above one."""

    def test_ratios(self):
        total, worst = 0, 0.0
        for name in ('e0', 'e3_graphene', 'e3_al', 'gk_graphene', 'gk_al'):
            with open(os.path.join(DATA, name + '.csv')) as f:
                r = [float(x['ratio']) for x in csv.DictReader(f) if x['ratio'] not in ('', None)]
            total += len(r); worst = max(worst, max(r))
        self.assertEqual(total, 12372)
        self.assertLessEqual(worst, 1.0)


if __name__ == '__main__':
    unittest.main()

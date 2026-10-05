#!/usr/bin/env python3
"""Number ledger for the in-loop rerun (Sec. IV D and SM): every statistic quoted from inloop_results.json."""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 0.05
rows = json.load(open(os.path.join(HERE, 'inloop_results.json')))
FAM = {'low_cutoff_continuation': 'low_cutoff', 'converged_subspace_reuse': 'oracle', 'coarse_k_mesh_upsample': 'coarse_k'}
L = {}
L['n_presented'] = len(rows)
L['by_decision'] = {d: sum(1 for r in rows if r['decision'] == d) for d in ('accept', 'reject', 'abstain')}
L['by_family'] = {f: sum(1 for r in rows if FAM[r['kind']] == f) for f in ('low_cutoff', 'oracle', 'coarse_k')}
rec = np.array([abs(r['d_full_reproduced'] - r['d_full_recorded']) / r['d_full_recorded'] for r in rows])
L['reproduction'] = dict(max_rel_dfull=float(rec.max()), median_rel_dfull=float(np.median(rec)),
                         max_ref_vs_exact_frame=float(max(r['ref_vs_exact_frame'] for r in rows)),
                         max_herm=float(max(r['herm_max'] for r in rows)))
R = [r for r in rows if r['ritz']]
L['n_ritz'] = len(R)
L['ritz_structure'] = dict(max_orthonormality=float(max(r['orthonormality'] for r in R)),
                           max_offdiag_rel=float(max(r['ritz_offdiag_rel'] for r in R)),
                           max_eig_dev=float(max(r['ritz_vs_recorded_eigs'] for r in R)),
                           max_candidate_reconstruction=float(max(r['candidate_reconstruction'] for r in R)),
                           max_identity_residual=float(max(r['identity_residual'] for r in R)))


def stats(v):
    v = np.asarray(v, float)
    return dict(n=int(v.size), min=float(v.min()), median=float(np.median(v)), max=float(v.max()))


for name, sel in (('all', R), ('low_cutoff', [r for r in R if FAM[r['kind']] == 'low_cutoff']),
                  ('oracle', [r for r in R if FAM[r['kind']] == 'oracle']),
                  ('accepted', [r for r in R if r['decision'] == 'accept'])):
    d = np.array([r['d_full'] for r in sel])
    L[name] = dict(
        d_full=stats(d),
        lead_signed=stats([(r['pyth'] - r['d_full']) / r['d_full'] for r in sel]),
        lead_abs=stats([abs(r['pyth'] - r['d_full']) / r['d_full'] for r in sel]),
        I1_ratio=stats([r['indicator'] / r['d_full'] for r in sel]),
        I1b_ratio=stats([r['indicator2'] / r['d_full'] for r in sel]),
        bound_ratio=stats([r['rigorous'] / r['d_full'] for r in sel]),
        bound_comp_ratio=stats([r['rigorous_computable'] / r['d_full'] for r in sel]),
        coupling_share=stats([r['blk_off'] ** 2 / r['d_full'] ** 2 for r in sel]),
        diag_over_d=stats([math.hypot(r['blk_P'], r['blk_Q']) / r['d_full'] for r in sel]),
        eP_eQ_max=float(max(max(r['eP'], r['eQ']) for r in sel)),
        BE_over_tol=stats([r['bound_cpl'] / TOL for r in sel]),
        underestimates=int(sum(1 for r in sel if r['indicator'] < r['d_full'])),
        lead2_abs=stats([abs(r['pyth2'] - r['d_full']) / r['d_full'] for r in sel]))
acc = [r for r in R if r['decision'] == 'accept']
bad = [r for r in acc if r['d_full'] > TOL]
good = [r for r in acc if r['d_full'] <= TOL]
L['screen_I1'] = dict(accepted=len(acc), bad=len(bad), good=len(good),
                      bad_flagged=sum(1 for r in bad if r['indicator'] > TOL),
                      good_kept=sum(1 for r in good if r['indicator'] <= TOL),
                      good_lost=[(r['run_id'], r['iteration'], r['d_full'], r['indicator']) for r in good if r['indicator'] > TOL],
                      abstained=[(r['run_id'], r['iteration'], r['d_full'], r['indicator']) for r in R if r['decision'] == 'abstain'],
                      window_bad_accepts=len(bad))
L['screen_rigorous'] = dict(accepted_by_bound=sum(1 for r in R if r['rigorous_computable'] <= TOL))
L['worst_lead'] = sorted([(abs(r['pyth'] - r['d_full']) / r['d_full'], r['run_id'], r['iteration'], r['kind'],
                           r['d_full'], r['pyth'], r['indicator']) for r in R], reverse=True)[:6]
L['min_I1_ratio'] = sorted([(r['indicator'] / r['d_full'], r['run_id'], r['iteration'], r['kind'], r['d_full'])
                            for r in R])[:5]
json.dump(L, open(os.path.join(HERE, 'inloop_ledger.json'), 'w'), indent=1, default=float)
print(json.dumps(L, indent=1, default=float))

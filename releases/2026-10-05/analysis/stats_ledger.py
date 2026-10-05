"""Number ledger for the revised manuscript: every statistic quoted from the decomposition sweeps."""
import json, math, os
import numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
o = json.load(open('decomposition_results.json'))
L = {}
rows = []
for m in o:
    for kind in ('extra_sweep', 'tau_sweep'):
        for r in o[m][kind]:
            r = dict(r); r['_m'] = m; r['_k'] = kind
            if not r['duplicate_of_buffer_sweep']:
                rows.append(r)
L['n_distinct'] = len(rows)
L['n_checks'] = sum(len(o[m]['extra_sweep']) + len(o[m]['tau_sweep']) for m in o)
for m in o:
    es = o[m]['extra_sweep']
    L[m] = dict(tau=o[m]['tau'], b_min=es[0]['extra'], b_max=es[-1]['extra'], n_buffer=len(es),
                n=sorted(set(es[0]['n'])), nk=o[m]['nk'],
                BE_min=min(r['bound_cpl'] for r in es), BE_max=max(r['bound_cpl'] for r in es),
                d_max=max(r['d_full'] for r in es), d_min=min(r['d_full'] for r in es))
    L[m]['BE_ratio'] = L[m]['BE_max'] / L[m]['BE_min']; L[m]['d_ratio'] = L[m]['d_max'] / L[m]['d_min']
    w = es[2]
    L[m]['b2'] = {k: w[k] for k in ('d_full', 'd_dk', 'blk_P', 'blk_Q', 'blk_off', 'eP', 'eQ', 'indicator', 'indicator2',
                                    'bound_cpl', 'rigorous', 'eps_rho', 'pyth', 'pyth2', 'm', 'mu_c', 'mu0', 'mu_star', 'gapQ')}
    L[m]['b2']['diag_sq'] = w['blk_P'] ** 2 + w['blk_Q'] ** 2
    L[m]['b2']['off_minus_dk_sq'] = w['blk_off'] ** 2 - w['d_dk'] ** 2
    L[m]['b2']['rel_dev_pyth'] = abs(w['pyth'] - w['d_full']) / w['d_full']
    L[m]['b2']['rel_dev_off_vs_dk'] = abs(w['blk_off'] - w['d_dk']) / w['blk_off']
    ts = o[m]['tau_sweep']
    L[m]['tau_eQ_over_dE'] = [(r['tau'], r['eQ'] / r['d_dk']) for r in ts]
    L[m]['tau_I_over_d'] = [(r['tau'], r['indicator'] / r['d_full']) for r in ts]
    L[m]['eQest_over_eQ'] = [(r['tau'], r['eQ_est'] / r['eQ']) for r in ts if r['eQ'] > 1e-3 * r['d_full']]
rel = np.array([abs(r['pyth'] - r['d_full']) / r['d_full'] for r in rows])
rel2 = np.array([abs(r['pyth2'] - r['d_full']) / r['d_full'] for r in rows])
i1 = np.array([r['indicator'] / r['d_full'] for r in rows]); i2 = np.array([r['indicator2'] / r['d_full'] for r in rows])
bd = np.array([r['d_full'] / r['rigorous'] for r in rows])
L['quad'] = dict(median=float(np.median(rel)), p95=float(np.percentile(rel, 95)), max=float(rel.max()),
                 n_gt_1pct=int((rel > 0.01).sum()), n_gt_0p4pct=int((rel > 0.004).sum()),
                 outliers=[(r['_m'], r['_k'], r['extra'], r['tau'], float(a), float(b)) for r, a, b in zip(rows, rel, rel2) if a > 0.004])
L['quad2'] = dict(median=float(np.median(rel2)), max=float(rel2.max()), n_gt_1pct=int((rel2 > 0.01).sum()),
                  argmax=[(r['_m'], r['_k'], r['extra'], r['tau']) for r, b in zip(rows, rel2) if b == rel2.max()])
work = np.array([(r['_k'] == 'extra_sweep') or (r['tau'] <= 0.054) for r in rows])
L['ind'] = dict(min=float(i1.min()), median=float(np.median(i1)), max=float(i1.max()),
                work_n=int(work.sum()), work_min=float(i1[work].min()), work_median=float(np.median(i1[work])), work_max=float(i1[work].max()),
                hot_max_by_mat={m: float(max(i for i, r in zip(i1, rows) if r['_m'] == m and not ((r['_k'] == 'extra_sweep') or (r['tau'] <= 0.054)))) for m in o},
                argmin=[(r['_m'], r['_k'], r['extra'], r['tau']) for r, a in zip(rows, i1) if a == i1.min()],
                ind2_min=float(i2.min()), ind2_maxchange=float((i2 / i1 - 1).max()), ind2_n_changed_gt1pct=int(((i2 / i1 - 1) > 0.01).sum()))
L['bound'] = dict(violations=int((bd > 1 + 1e-12).sum()), max_ratio=float(bd.max()),
                  argmax=[(r['_m'], r['_k'], r['extra'], r['tau']) for r, a in zip(rows, bd) if a == bd.max()],
                  n_full_bound_below_tol=int(sum(r['rigorous'] < 0.05 for r in rows)),
                  full_bound_below_tol=[(r['_m'], r['tau']) for r in rows if r['rigorous'] < 0.05],
                  BE_below_tol=[(r['_m'], r['_k'], r['extra'], r['tau'], r['bound_cpl']) for r in rows if r['bound_cpl'] < 0.05])
json.dump(L, open('stats_ledger.json', 'w'), indent=1, default=float)
if __name__ == '__main__':
    print(json.dumps({k: L[k] for k in ('n_distinct', 'n_checks', 'quad', 'quad2', 'ind', 'bound')}, indent=1, default=float))
    for m in o:
        x = L[m]; print(m, 'b %d..%d (%d), n %s; BE %.3g-%.3g (%.2fx); d %.3g-%.3g (%.0fx)' % (x['b_min'], x['b_max'], x['n_buffer'], x['n'], x['BE_min'], x['BE_max'], x['BE_ratio'], x['d_min'], x['d_max'], x['d_ratio']))
        b = x['b2']; print('   b2: d %.4e dk %.4e P %.3e Q %.3e off %.4e | diag^2 %.3e off^2-dk^2 %.3e | dev %.2e offdev %.2e | I %.4e I2 %.4e BE %.4g eps %.4f m %s' % (
            b['d_full'], b['d_dk'], b['blk_P'], b['blk_Q'], b['blk_off'], b['diag_sq'], b['off_minus_dk_sq'], b['rel_dev_pyth'], b['rel_dev_off_vs_dk'], b['indicator'], b['indicator2'], b['bound_cpl'], b['eps_rho'], sorted(set(b['m']))))
        print('   eQ/dE vs tau:', ' '.join('%.3g:%.2g' % t for t in x['tau_eQ_over_dE'][-6:]))
        print('   eQest/eQ:', ' '.join('%.3g:%.1f' % t for t in x['eQest_over_eQ']))

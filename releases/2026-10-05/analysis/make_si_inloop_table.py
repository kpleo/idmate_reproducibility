#!/usr/bin/env python3
"""Rows of SM Table S-inloop: every candidate presented in the SCF tests, rerun with saved Hamiltonians."""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(HERE, 'inloop_results.json')))
FAM = {'low_cutoff_continuation': 'low cutoff', 'converged_subspace_reuse': 'oracle',
       'coarse_k_mesh_upsample': r'coarse $k$'}
DEC = {'accept': r'$\checkmark$', 'reject': r'$\times$', 'abstain': r'$\varnothing$'}
MORD = {'Si': 0, 'graphene': 1, 'Al': 2}


def sci(x, nd=2):
    """three significant digits (trailing zeros kept) between 0.01 and 1000, otherwise a x 10^e."""
    if x is None:
        return '--'
    if x == 0:
        return '0'
    if 1e-2 <= abs(x) < 1e3:
        dec = max(0, nd - int(np.floor(np.log10(abs(x)))))
        return f'{x:.{dec}f}'
    e = int(np.floor(np.log10(abs(x))))
    return f'${x / 10 ** e:.{nd - 1}f}\\times10^{{{e}}}$'


def run_label(r):
    if 'multiple' in r['run_id']:
        return 'multi'
    return 'single'


out = []
key = lambda r: (MORD[r['material']], 0 if 'multiple' in r['run_id'] else 1, r['iteration'], r['kind'])
for r in sorted(rows, key=key):
    dw = r['d_window']
    cells = [r['material'], run_label(r), str(r['iteration']), FAM[r['kind']], DEC[r['decision']],
             sci(dw) if dw is not None else '--', sci(r['d_full_recorded'])]
    if r['ritz']:
        cells += [sci(r['pyth']), sci(r['indicator']), sci(r['rigorous_computable'])]
    else:
        cells += ['--', '--', '--']
    out.append(' & '.join(cells) + r' \\')
open(os.path.join(HERE, '..', 'manuscript', 'tabS_inloop_rows.tex'), 'w').write('\n'.join(out) + '\n')
print(len(out), 'rows')


# ---------------------------------------------------------------- summary by family (SM Table S-rerun-summary)
def rng(v, fmt):
    v = np.asarray(v, float)
    a, b = fmt(v.min()), fmt(v.max())
    return a if a == b else f'{a}--{b}'


def pct(x):
    return f'{x:+.2f}' if abs(x) < 1 else f'{x:+.1f}'


summ = []
for kind, lab in (('low_cutoff_continuation', 'low cutoff'), ('converged_subspace_reuse', 'oracle')):
    for dec, dlab in (('accept', 'accepted'), ('abstain', 'abstained')):
        S = [r for r in rows if r['kind'] == kind and r['decision'] == dec]
        if not S:
            continue
        lead = [100 * (r['pyth'] - r['d_full']) / r['d_full'] for r in S]
        la, lb = min(lead), max(lead)
        lead_s = f'${pct(la)}$' if abs(la - lb) < 1e-9 else f'${pct(la)}$ to ${pct(lb)}$'
        cells = [f'{lab} ({dlab})', str(len(S)),
                 rng([r['d_full'] for r in S], lambda x: sci(x)),
                 lead_s,
                 rng([100 * r['blk_off'] ** 2 / r['d_full'] ** 2 for r in S], lambda x: f'{x:.1f}' if x < 99.95 else '100'),
                 rng([r['indicator'] / r['d_full'] for r in S], lambda x: f'{x:.3g}'),
                 rng([r['rigorous_computable'] / r['d_full'] for r in S], lambda x: f'{x:.0f}'),
                 str(sum((r['indicator'] <= 0.05) == (r['d_full'] <= 0.05) for r in S))]
        summ.append(' & '.join(cells) + r' \\')
open(os.path.join(HERE, '..', 'manuscript', 'tabS_inloop_summary_rows.tex'), 'w').write('\n'.join(summ) + '\n')
print('\n'.join(summ))

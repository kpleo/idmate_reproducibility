"""Rebuild the LaTeX rows of the Supplemental tables from the derived inputs in ../data and compare them with
the rows typeset in the manuscript (tables/expected/).

Tables: S1 (buffer sweep) and S2 (temperature sweep) of the saved-Hamiltonian decomposition, and the summary and
per-candidate tables of the in-loop rerun with saved Hamiltonians."""
import json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, '..', 'build', 'tables')
EXP = os.path.join(HERE, 'expected')
LABEL = {'si': 'Si', 'graphene': 'graphene', 'al': 'Al'}


# ---------------------------------------------------------------- S1, S2 (saved-Hamiltonian sweeps)
def num(x, small=1e-30):
    if x == 0 or abs(x) < small:
        return '$<10^{-30}$'
    if 0.01 <= abs(x) < 1000:
        digits = max(0, 3 - int(math.floor(math.log10(abs(x)))) - 1)
        return f'{x:.{digits}f}'
    e = int(math.floor(math.log10(abs(x)))); m = x / 10 ** e
    return f'${m:.2f}\\times10^{{{e}}}$'


def pct(x):
    return f'${100 * x:.2f}$' if 100 * x >= 0.005 else '$<0.01$'


def ratio(x):
    return f'${x:.3g}$' if x >= 100 else f'${x:.2f}$'


def share(x):
    if x >= 0.01:
        return f'${x:.3f}$'
    e = int(math.floor(math.log10(abs(x)))); m = x / 10 ** e
    return f'${m:.1f}\\times10^{{{e}}}$'


def sweep_tables(o):
    BUF = {'si': [0, 2, 4, 8, 16, 32, 50], 'graphene': [0, 2, 4, 8, 14, 15, 16, 17, 32, 64, 128],
           'al': [0, 2, 4, 8, 16, 32, 64]}
    s1 = []
    for key in ['si', 'graphene', 'al']:
        rows = {r['extra']: r for r in o[key]['extra_sweep']}
        for i, b in enumerate(BUF[key]):
            r = rows[b]; d = r['d_full']
            dev1 = abs(r['pyth'] - d) / d; dev2 = abs(r['pyth2'] - d) / d
            cells = [LABEL[key] if i == 0 else '', str(b), str(r['m'][0]), num(d), num(r['d_dk']), num(r['blk_P']),
                     pct(dev1), pct(dev2), f"${r['indicator'] / d:.2f}$", f"${r['rigorous'] / d:.0f}$", num(r['eps_rho'])]
            s1.append(' & '.join(cells) + ' \\\\')
        if key != 'al':
            s1.append('\\midrule')
    TAUS = [0.002, 0.00746, 0.02, 0.0386, 0.0746, 0.1036, 0.1439, 0.2]
    s2 = []
    for key in ['si', 'graphene', 'al']:
        ts = o[key]['tau_sweep']
        for i, t in enumerate(TAUS):
            r = min(ts, key=lambda r: abs(r['tau'] - t)); d = r['d_full']
            dev1 = abs(r['pyth'] - d) / d
            cells = [LABEL[key] if i == 0 else '', f"${r['tau']:.4f}$", num(d), num(r['d_dk']), num(r['eQ']),
                     num(r['eP']) if r['eP'] > 1e-12 else '$<10^{-12}$', share(r['blk_off'] ** 2 / d ** 2), pct(dev1),
                     ratio(r['indicator'] / d), f"${r['rigorous'] / d:.1f}$"]
            s2.append(' & '.join(cells) + ' \\\\')
        if key != 'al':
            s2.append('\\midrule')
    return s1, s2


# ---------------------------------------------------------------- in-loop rerun tables
FAM = {'low_cutoff_continuation': 'low cutoff', 'converged_subspace_reuse': 'oracle',
       'coarse_k_mesh_upsample': r'coarse $k$'}
DEC = {'accept': r'$\checkmark$', 'reject': r'$\times$', 'abstain': r'$\varnothing$'}
MORD = {'Si': 0, 'graphene': 1, 'Al': 2}


def sci(x, nd=2):
    if x is None:
        return '--'
    if x == 0:
        return '0'
    if 1e-2 <= abs(x) < 1e3:
        dec = max(0, nd - int(np.floor(np.log10(abs(x)))))
        return f'{x:.{dec}f}'
    e = int(np.floor(np.log10(abs(x))))
    return f'${x / 10 ** e:.{nd - 1}f}\\times10^{{{e}}}$'


def inloop_tables(rows):
    key = lambda r: (MORD[r['material']], 0 if 'multiple' in r['run_id'] else 1, r['iteration'], r['kind'])
    full = []
    for r in sorted(rows, key=key):
        dw = r['d_window']
        cells = [r['material'], 'multi' if 'multiple' in r['run_id'] else 'single', str(r['iteration']), FAM[r['kind']],
                 DEC[r['decision']], sci(dw) if dw is not None else '--', sci(r['d_full_recorded'])]
        cells += ([sci(r['pyth']), sci(r['indicator']), sci(r['rigorous_computable'])] if r['ritz'] else ['--'] * 3)
        full.append(' & '.join(cells) + r' \\')

    def rng(v, fmt):
        v = np.asarray(v, float)
        a, b = fmt(v.min()), fmt(v.max())
        return a if a == b else f'{a}--{b}'

    def pc(x):
        return f'{x:+.2f}' if abs(x) < 1 else f'{x:+.1f}'

    summ = []
    for kind, lab in (('low_cutoff_continuation', 'low cutoff'), ('converged_subspace_reuse', 'oracle')):
        for dec, dlab in (('accept', 'accepted'), ('abstain', 'abstained')):
            S = [r for r in rows if r['kind'] == kind and r['decision'] == dec]
            if not S:
                continue
            lead = [100 * (r['pyth'] - r['d_full']) / r['d_full'] for r in S]
            la, lb = min(lead), max(lead)
            lead_s = f'${pc(la)}$' if abs(la - lb) < 1e-9 else f'${pc(la)}$ to ${pc(lb)}$'
            cells = [f'{lab} ({dlab})', str(len(S)), rng([r['d_full'] for r in S], sci), lead_s,
                     rng([100 * r['blk_off'] ** 2 / r['d_full'] ** 2 for r in S],
                         lambda x: f'{x:.1f}' if x < 99.95 else '100'),
                     rng([r['indicator'] / r['d_full'] for r in S], lambda x: f'{x:.3g}'),
                     rng([r['rigorous_computable'] / r['d_full'] for r in S], lambda x: f'{x:.0f}'),
                     str(sum((r['indicator'] <= 0.05) == (r['d_full'] <= 0.05) for r in S))]
            summ.append(' & '.join(cells) + r' \\')
    return full, summ


def main():
    os.makedirs(OUT, exist_ok=True)
    o = json.load(open(os.path.join(DATA, 'decomposition_results.json')))
    rows = json.load(open(os.path.join(DATA, 'inloop_results.json')))
    s1, s2 = sweep_tables(o)
    full, summ = inloop_tables(rows)
    out = {'tabS_sweep_rows.tex': s1, 'tabS_tau_rows.tex': s2, 'tabS_inloop_rows.tex': full,
           'tabS_inloop_summary_rows.tex': summ}
    status = {}
    for name, lines in out.items():
        text = '\n'.join(lines) + '\n'
        open(os.path.join(OUT, name), 'w').write(text)
        status[name] = text == open(os.path.join(EXP, name)).read()
    print(json.dumps(status, indent=1))
    return 0 if all(status.values()) else 1


if __name__ == '__main__':
    sys.exit(main())

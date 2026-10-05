"""Write the LaTeX rows of SI Tables S1 (buffer sweep) and S2 (temperature sweep)
from decomposition_results.json."""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
o = json.load(open(os.path.join(HERE, 'decomposition_results.json')))
OUT = os.path.join(HERE, '..', 'manuscript')
LABEL = {'si': 'Si', 'graphene': 'graphene', 'al': 'Al'}


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


BUF = {'si': [0, 2, 4, 8, 16, 32, 50], 'graphene': [0, 2, 4, 8, 14, 15, 16, 17, 32, 64, 128], 'al': [0, 2, 4, 8, 16, 32, 64]}
lines = []
for key in ['si', 'graphene', 'al']:
    rows = {r['extra']: r for r in o[key]['extra_sweep']}
    for i, b in enumerate(BUF[key]):
        r = rows[b]; d = r['d_full']
        dev1 = abs(r['pyth'] - d) / d; dev2 = abs(r['pyth2'] - d) / d
        cells = [LABEL[key] if i == 0 else '', str(b), str(r['m'][0]), num(d), num(r['d_dk']), num(r['blk_P']),
                 pct(dev1), pct(dev2), f"${r['indicator'] / d:.2f}$", f"${r['rigorous'] / d:.0f}$", num(r['eps_rho'])]
        lines.append(' & '.join(cells) + ' \\\\')
    if key != 'al':
        lines.append('\\midrule')
open(os.path.join(OUT, 'tabS_sweep_rows.tex'), 'w').write('\n'.join(lines) + '\n')

TAUS = [0.002, 0.00746, 0.02, 0.0386, 0.0746, 0.1036, 0.1439, 0.2]
lines = []
for key in ['si', 'graphene', 'al']:
    ts = o[key]['tau_sweep']
    for i, t in enumerate(TAUS):
        r = min(ts, key=lambda r: abs(r['tau'] - t)); d = r['d_full']
        dev1 = abs(r['pyth'] - d) / d
        cells = [LABEL[key] if i == 0 else '', f"${r['tau']:.4f}$", num(d), num(r['d_dk']), num(r['eQ']),
                 num(r['eP']) if r['eP'] > 1e-12 else '$<10^{-12}$', share(r['blk_off'] ** 2 / d ** 2), pct(dev1),
                 ratio(r['indicator'] / d), f"${r['rigorous'] / d:.1f}$"]
        lines.append(' & '.join(cells) + ' \\\\')
    if key != 'al':
        lines.append('\\midrule')
open(os.path.join(OUT, 'tabS_tau_rows.tex'), 'w').write('\n'.join(lines) + '\n')
print(open(os.path.join(OUT, 'tabS_sweep_rows.tex')).read())
print(open(os.path.join(OUT, 'tabS_tau_rows.tex')).read())

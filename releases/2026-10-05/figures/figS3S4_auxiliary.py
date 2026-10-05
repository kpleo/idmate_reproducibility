"""Supplemental Figs. S3 and S4 restyled from the stored records (supplement.json).
S3 | Fixed-point-preserving interventions move the intrinsic margin; mixers change only amplification.
S4 | The recorded-counter work regression missed its prespecified target.
Drawn at the printed SM width (6.5 in) in the letter style of prb_style.py."""
import json, os, sys
import numpy as np
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prb_style as S

D = json.load(open(os.path.join(S.DATA, 'supplement.json')))
OUT = S.OUT
C1, C2, C3 = S.BLUE, S.GREEN, S.ORANGE

# ------------------------------------------------------------------ Fig. S3
W, H = 6.5, 2.55
fig, L = S.figure(W, H)
y0, h = 0.48, 1.62
a = L.ax(0.60, y0, 1.55, h); b = L.ax(2.72, y0, 1.55, h); c = L.ax(5.05, y0, 1.38, h)
L.head(0.0, y0 + h + 0.2, 'a', 'Jacobian update error')
L.head(2.30, y0 + h + 0.2, 'b', 'Intrinsic margin')
L.head(4.40, y0 + h + 0.2, 'c', 'Mixer amplification')
rungs = D['h2']['rank_one_rungs']
fr = [r['fraction'] for r in rungs]
by = {r['lambda']: r for r in D['coupling']}
lams = sorted(by, key=abs)
series = [(r'H$_2$/Grid4', fr, [r['rank_one_update_error'] for r in rungs], 'o', C1),
          ('Plane-wave model', D['pilot']['fixed_fractions'], D['pilot']['jacobian_update_errors'], 's', C2),
          ('HGH Si', [0.2, 0.4], [by[x]['rank_one_update_error'] for x in lams[1:]], 'D', C3)]
for name, x, y, mk, col in series:
    a.plot(x, y, '-', color=col, lw=S.LW_DATA)
    a.plot(x, y, mk, color=col, ms=4.2, label=name, **S.MKW)
a.set_yscale('log'); a.set_xlim(0.06, 0.44); a.set_ylim(2e-7, 4e-5); a.set_xticks([0.1, 0.2, 0.3, 0.4])
S.xlabel(a, r'$|\lambda|/|\lambda_c|$'); S.ylabel(a, 'Rank-one update error', labelpad=1.0)
a.legend(loc='upper left', fontsize=S.FS, borderaxespad=0.1, handlelength=1.0)
h2b = rungs[0]['m_fp'] - rungs[0]['delta_m_fp']
hx = [0.0] + [f / fr[-1] for f in fr]; hy = [h2b] + [r['m_fp'] for r in rungs]
b.plot(hx, hy, '-', color=C1, lw=S.LW_DATA)
b.plot(hx, hy, 'o', color=C1, ms=4.0, label=r'H$_2$/Grid4', **S.MKW)
styles = [('s', '-', '$1^3$'), ('D', (0, (3, 1.5)), '$2^3$'), ('^', (0, (1, 1.2)), '$4^3$')]
for mesh, (mk, ls, lab) in zip(D['k2']['meshes'], styles):
    cr = mesh['coupling_rungs']
    x = [abs(r['lambda']) / abs(cr[-1]['lambda']) for r in cr]; y = [r['fixed_point_margin'] for r in cr]
    b.plot(x, y, ls=ls, color=C2, lw=S.LW_DATA)
    kw = dict(mfc='white', mec=C2, mew=0.9) if mk == 'D' else dict(mfc=C2, mec='white', mew=0.35)
    b.plot(x, y, mk, ms=3.8, label=f'Weighted-$k$, {lab}', **kw)
b.set_xlim(-0.05, 1.05); b.set_ylim(0.5, 1.04); b.set_xticks([0, 0.5, 1])
S.xlabel(b, r'$|\lambda|/|\lambda_\mathrm{max}|$'); S.ylabel(b, r'$m_\mathrm{fp}=\sigma_\mathrm{min}(I-J)$', labelpad=1.0)
b.legend(loc='lower left', fontsize=S.FS, borderaxespad=0.1, handlelength=1.8)
names = ['unit', 'scalar_0.25', 'scalar_0.70', 'baseline_inverse_control']
labels = ['Unit', 'Scalar 0.25', 'Scalar 0.70', 'Exact inverse\n(control)']
mks = ['o', 's', 'D']; faces = ['white', C2, C1]; edges = [S.INK2, 'white', 'white']
for row, name in zip([3, 2, 1, 0], names):
    rkeys = sorted(D['radii'], key=lambda s: abs(float(s)))
    x = [D['radii'][k][name] for k in rkeys]; y = [row + dy for dy in (0.17, 0, -0.17)]
    c.plot(x, y, color=S.GRID, lw=S.LW_GUIDE, zorder=1)
    for v, p, mk, fc, ec in zip(x, y, mks, faces, edges):
        c.plot(v, p, mk, ms=4.0, mfc=fc, mec=ec, mew=0.9 if fc == 'white' else 0.35, zorder=3)
c.axvline(1, color=S.INK, lw=S.LW_GUIDE, ls=(0, (3, 2)))
S.ann(c, 0.985, -0.55, r'$\rho=1$', ha='right', va='bottom')
c.set_yticks([3, 2, 1, 0]); c.set_yticklabels(labels); c.tick_params(axis='y', length=0)
for t in c.get_yticklabels():
    t.set_color(S.INK)
c.set_xlim(-0.05, 1.12); c.set_ylim(-0.6, 3.6); c.set_xticks([0, 0.5, 1])
S.xlabel(c, r'$\rho(M_B)$, HGH Si')
c.legend(handles=[Line2D([], [], ls='', marker=m, mfc=f, mec=e if f != 'white' else S.INK2, mew=0.9 if f == 'white' else 0.35,
                         ms=4.0) for m, f, e in zip(mks, faces, edges)],
         labels=[r'$\lambda=0$', r'$0.2\lambda_c$', r'$0.4\lambda_c$'], loc='upper left', fontsize=S.FS,
         bbox_to_anchor=(0.0, 1.04), borderaxespad=0.1, handletextpad=0.2)
S.save(fig, os.path.join(OUT, 'figS3_interventions'))

# ------------------------------------------------------------------ Fig. S4
W, H = 6.5, 2.6
fig, L = S.figure(W, H)
y0, h = 0.48, 1.66
a = L.ax(0.18, y0, 2.15, h); b = L.ax(2.86, y0, 1.45, h); c = L.ax(5.0, y0, 1.42, h)
L.head(0.0, y0 + h + 0.2, 'a', 'Held-out error reduction')
L.head(2.52, y0 + h + 0.2, 'b', 'Per category')
L.head(4.44, y0 + h + 0.2, 'c', 'Projector geometry')
v = D['verdict']; gate = v['g2_gate_evaluation']['criteria']
est = gate['margin']['improvement'] * 100; target = gate['margin']['threshold'] * 100
lo, hi = (v['bootstrap'][k] * 100 for k in ('ci_lower', 'ci_upper'))
a.axvline(0, color=S.INK3, lw=S.LW_GUIDE)
a.axvline(target, color=S.FAIL, lw=S.LW_GUIDE, ls=(0, (3, 2)))
a.plot([lo, hi], [0.6, 0.6], color=S.INK2, lw=S.LW_DATA)
a.plot([lo, hi], [0.6, 0.6], ls='', marker='|', ms=8, mew=1.2, color=S.INK2)
a.plot([est], [0.6], 'o', ms=5.2, color=S.BLUE, zorder=3, **S.MKW)
S.callout(a, est - 4, 0.70, f'Corrected fit {est:+.2f}%', ha='right', color=S.BLUE, fontsize=S.FA)
S.ann(a, target + 2, 0.94, f'Prespecified\ntarget {target:.0f}%', color=S.FAIL, va='top')
S.ann(a, -185, 0.47, f'95% hierarchical interval\n[{lo:.1f}%, +{hi:.1f}%]'.replace('-', '−'), va='top', ha='left')
n = v['analysis']
S.ann(a, -185, 0.29, f"{n['rows_admitted_paired']:,} paired SCF steps,\n{n['trajectories']} trajectories, "
      f"{n['trajectories_effective']} distinct series,\n{len(n['parents'])} parent lineages", va='top', ha='left')
a.set_xlim(-190, 60); a.set_ylim(0, 1); a.set_yticks([]); a.spines['left'].set_visible(False)
S.xlabel(a, 'Error reduction vs best baseline (%)')
per = gate['category_reversal']['per_category_improvement']
b.axvline(0, color=S.INK3, lw=S.LW_GUIDE, zorder=0); b.axvline(-5, color=S.FAIL, lw=S.LW_GUIDE, ls=(0, (3, 2)), zorder=0)
for pos, name in zip(range(7, -1, -1), [f'C{i}' for i in range(1, 9)]):
    if name not in per:
        S.ann(b, 2, pos, 'No eligible steps', va='center')
        continue
    val = per[name] * 100; rev = val < -5
    col = S.FAIL if rev else S.BLUE
    b.plot([0, val], [pos, pos], color=col, lw=S.LW_DATA)
    b.plot(val, pos, 's' if rev else 'o', ms=4.2, color=col, **S.MKW)
    b.text(68, pos, f'{val:+.1f}'.replace('-', '−'), ha='right', va='center', fontsize=S.FS, color=S.INK2)
b.set_yticks(range(7, -1, -1)); b.set_yticklabels([f'C{i}' for i in range(1, 9)]); b.tick_params(axis='y', length=0)
for t in b.get_yticklabels():
    t.set_color(S.INK)
b.set_xlim(-25, 70); b.set_ylim(-0.7, 7.6); S.xlabel(b, 'Reduction (%)')
ak, dp = np.array(D['ak']), np.array(D['dp'])
c.plot(ak, dp, 'o', ms=3.4, mfc='none', mec=S.BLUE, mew=0.8, zorder=3)
xs = np.geomspace(ak.min(), ak.max(), 100)
c.plot(xs, np.exp(D['intercept']) * xs ** D['slope'], color=S.INK2, lw=S.LW_GUIDE)
c.set_xscale('log'); c.set_yscale('log')
S.xlabel(c, r'Adiabaticity $A_k$'); S.ylabel(c, r'Projector distance $d_P$', labelpad=1.0)
S.callout(c, 0.04, 0.96, f"Slope {D['slope']:.3f}", transform=c.transAxes, ha='left', va='top', fontsize=S.FA)
S.ann(c, 0.04, 0.80, '60 pairs, one trajectory', transform=c.transAxes, ha='left', va='top')
S.save(fig, os.path.join(OUT, 'figS4_worklaw'))
print('ok')

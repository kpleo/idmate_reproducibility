"""Fig. 3 | The coupling block carries the full-state error of Rayleigh-Ritz candidates.
Saved IDMate Hamiltonians H2 and preceding-map frames U1 (three Hamiltonians, 296 distinct candidates).
(a-c) Error versus number of buffer bands kept beyond the occupied set (every buffer size).
(d) Share of the squared error carried by the coupling block versus temperature (exact identity).
(e) Indicator relative to the measured error.
(f) Accuracy of the leading-order decomposition, with and without the second-order near-Fermi shift.
Drawn at printed size in the letter style of prb_style.py."""
import json, os, sys
import numpy as np
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prb_style as S

res = json.load(open(os.path.join(S.DATA, 'decomposition_results.json')))
t2 = {r['material']: r for r in json.load(open(os.path.join(S.DATA, 'table2_decomposition.json')))}
MATS = [('si', 'Si'), ('graphene', 'graphene'), ('al', 'Al')]
HOT = 0.06          # tau above which the complement is thermally populated (tau-sweep points only)


def sci(x):
    e = int(np.floor(np.log10(abs(x))))
    return f'${x / 10 ** e:.1f}\\times10^{{{e}}}$'


def distinct(key):
    return [r for r in res[key]['extra_sweep']] + [r for r in res[key]['tau_sweep'] if not r['duplicate_of_buffer_sweep']]


W, H = S.W2, 4.70
fig, L = S.figure(W, H)
pw, ph = 1.93, 1.46
xs = [0.52, 2.79, 5.02]
yt, yb = 2.70, 0.44
top = [L.ax(x, yt, pw, ph) for x in xs]
axd, axe, axf = (L.ax(x, yb, pw, ph) for x in xs)
for ax, x, (key, mat) in zip(top, xs, MATS):
    pass
L.head(0.0, yt + ph + 0.2, 'a', r'Diamond Si, $\tau=0.02$ Ha')
L.head(xs[1] - 0.42, yt + ph + 0.2, 'b', r'Graphene, $\tau=0.01$ Ha')
L.head(xs[2] - 0.42, yt + ph + 0.2, 'c', r'fcc Al, $\tau=0.02$ Ha')
L.head(0.0, yb + ph + 0.2, 'd', 'Coupling share of the error')
L.head(xs[1] - 0.42, yb + ph + 0.2, 'e', 'Indicator over measured error')
L.head(xs[2] - 0.42, yb + ph + 0.2, 'f', 'Accuracy of the decomposition')

# ---------------------------------------------------------------- (a)-(c) buffer sweeps
for ax, (key, mat) in zip(top, MATS):
    rows = res[key]['extra_sweep']
    x = np.array([r['extra'] for r in rows])
    d = np.array([r['d_full'] for r in rows]); py = np.array([r['pyth'] for r in rows])
    ind = np.array([r['indicator'] for r in rows]); gb = np.array([r['rigorous'] for r in rows])
    diag = np.sqrt(np.array([r['blk_P'] ** 2 + r['blk_Q'] ** 2 for r in rows]))
    ax.plot(x, gb, color=S.BOUND, lw=S.LW_DATA, ls=(0, (4, 2)))
    ax.plot(x, ind, color=S.IND, lw=S.LW_DATA, ls=(0, (1.2, 1.2)), zorder=3)
    ax.plot(x, py, color=S.LEAD, lw=S.LW_DATA, zorder=2)
    ax.plot(x, diag, color=S.DIAG, lw=S.LW_DATA)
    ax.plot(x, d, 'o', ms=2.6, color=S.MEAS, zorder=4, **S.MKW)
    ax.axvline(2, color=S.GRID, lw=S.LW_GUIDE, zorder=0)
    if t2[key]['B_W'] is not None:
        ya = 1.5e-5 if key == 'si' else 2.0e-4
        ax.annotate('', xy=(2, 1.05e-6), xytext=(2, ya * 3.0),
                    arrowprops=dict(arrowstyle='-|>', lw=0.9, color=S.WINDOW, mutation_scale=7))
        S.ann(ax, 2 + 0.05 * x.max(), ya, f"Window bound off scale:\n$B_W=${sci(t2[key]['B_W'])}",
              color=S.WINDOW, va='center')
    else:
        S.ann(ax, 2 + 0.05 * x.max(), 6.0e-6, 'Window test abstains\n' + r'(no level within $\pm18\tau$)',
              color=S.WINDOW, va='center')
    ax.set_yscale('log'); ax.set_ylim(1.0e-6, 2e4)
    ax.set_xlim(-0.03 * x.max() - 1, x.max() + 0.03 * x.max() + 1)
    ax.set_yticks([1e-6, 1e-4, 1e-2, 1e0, 1e2, 1e4])
    S.xlabel(ax, r'Buffer bands $b$')
    print(key, 'n range', min(rows[0]['n']), max(rows[0]['n']))
S.ylabel(top[0], 'Embedded distance', labelpad=1.0)
for ax in top[1:]:
    ax.tick_params(labelleft=False)
handles = [Line2D([], [], ls='', marker='o', ms=3.4, color=S.MEAS, mec='white', mew=0.35, label=r'Measured $d_\mathrm{full}$'),
           Line2D([], [], color=S.LEAD, lw=S.LW_DATA, label='Leading order'),
           Line2D([], [], color=S.IND, lw=S.LW_DATA, ls=(0, (1.2, 1.2)), label=r'Indicator $I_1$'),
           Line2D([], [], color=S.BOUND, lw=S.LW_DATA, ls=(0, (4, 2)), label='Rigorous bound'),
           Line2D([], [], color=S.DIAG, lw=S.LW_DATA, label='Diagonal blocks')]
top[0].legend(handles=handles, loc='upper right', ncol=2, fontsize=S.FS, handlelength=1.6, columnspacing=0.8,
              borderaxespad=0.1, bbox_to_anchor=(1.02, 1.03))
# graphene near-Fermi second-order case
g = {r['extra']: r for r in res['graphene']['extra_sweep']}[16]
top[1].annotate('', xy=(16, g['d_full'] * 1.25), xytext=(30, 0.9),
                arrowprops=dict(arrowstyle='-', lw=0.6, color=S.INK2, shrinkA=0.5, shrinkB=1.5))
S.ann(top[1], 31, 0.9, f"$b=16$: first order {100 * abs(g['pyth'] - g['d_full']) / g['d_full']:.0f}% low,\n"
      f"second order {100 * abs(g['pyth2'] - g['d_full']) / g['d_full']:.2f}%", ha='left', va='center', color=S.INK)

# ---------------------------------------------------------------- (d) coupling share versus temperature
for key, mat in MATS:
    rows = res[key]['tau_sweep']
    tau = np.array([r['tau'] for r in rows])
    share = np.array([r['blk_off'] ** 2 / r['d_full'] ** 2 for r in rows])
    axd.plot(tau, share, '-', color=S.MAT[mat], lw=S.LW_DATA, zorder=2)
    axd.plot(tau, share, S.MAT_MK[mat], color=S.MAT[mat], ms=3.4, zorder=3, label=S.MAT_LABEL[mat], **S.MKW)
for tw in sorted({res[k]['tau'] for k, _ in MATS}):
    axd.axvline(tw, color=S.GRID, lw=S.LW_GUIDE, zorder=0)
S.ann(axd, 0.0141, 0.42, r'Working $\tau$', rotation=90, va='bottom', ha='center')
axd.axhline(0.5, color=S.INK3, lw=S.LW_THIN, ls=(0, (3, 2)), zorder=0)
axd.set_xscale('log'); axd.set_xlim(1.8e-3, 0.23); axd.set_ylim(-0.03, 1.06)
S.xlabel(axd, r'Temperature $\tau$ (Ha)'); S.ylabel(axd, r'$\delta_E^2/d_\mathrm{full}^2$', labelpad=1.0)
axd.legend(loc='lower left', fontsize=S.FS, handlelength=1.0, borderaxespad=0.1)
S.ann(axd, 0.04, 0.88, r'Exact block identity, $b=2$', transform=axd.transAxes, ha='left', va='top', fontsize=S.FS)

# ---------------------------------------------------------------- (e) indicator / measured
axe.axhspan(0.6, 1.0, color=S.PALE, lw=0, zorder=0)
axe.axhline(1.0, color=S.INK, lw=S.LW_GUIDE)
for key, mat in MATS:
    rows = distinct(key)
    d = np.array([r['d_full'] for r in rows]); hot = np.array([r['tau'] > HOT for r in rows])
    ii = np.array([r['indicator'] / r['d_full'] for r in rows])
    axe.plot(d[~hot], ii[~hot], S.MAT_MK[mat], ms=3.6, mfc='white', mec=S.MAT[mat], mew=0.75, zorder=3)
    axe.plot(d[hot], ii[hot], S.MAT_MK[mat], ms=3.6, color=S.MAT[mat], zorder=3, **S.MKW)
axe.set_xscale('log'); axe.set_yscale('log'); axe.set_xlim(3e-4, 0.3); axe.set_ylim(0.6, 60)
S.xlabel(axe, r'Measured $d_\mathrm{full}$'); S.ylabel(axe, r'$I_1/d_\mathrm{full}$', labelpad=1.0)
S.ann(axe, 0.03, 0.04, 'Underestimate', transform=axe.transAxes, ha='left', va='bottom')
allr = [r for k, _ in MATS for r in distinct(k)]
ii = np.array([r['indicator'] / r['d_full'] for r in allr]); cold = np.array([r['tau'] <= HOT for r in allr])
bb = np.array([r['d_full'] / r['rigorous'] for r in allr])
S.callout(axe, 0.03, 0.97, f"{cold.sum()} cases, $\\tau\\leq0.054$ Ha:\nratio {ii[cold].min():.3f}–{ii[cold].max():.1f}, "
          f"median {np.median(ii[cold]):.2f}", transform=axe.transAxes, ha='left', va='top', fontsize=S.FA)
he = [Line2D([], [], ls='', marker='o', mfc='none', mec=S.INK2, ms=3.2, mew=0.8, label=r'$\tau\leq0.054$ Ha'),
      Line2D([], [], ls='', marker='o', mfc=S.INK2, mec='white', mew=0.35, ms=3.6, label=r'$\tau\geq0.075$ Ha')]
axe.legend(handles=he, loc='upper left', fontsize=S.FS, handletextpad=0.3, bbox_to_anchor=(0.0, 0.74),
           borderaxespad=0.1)

# ---------------------------------------------------------------- (f) accuracy of the decomposition
rel = np.sort([abs(r['pyth'] - r['d_full']) / r['d_full'] for r in allr])
rel2 = np.sort([abs(r['pyth2'] - r['d_full']) / r['d_full'] for r in allr])
yy = np.arange(1, len(rel) + 1) / len(rel)
axf.step(rel, yy, where='post', color=S.LEAD, lw=S.LW_DATA, label='Leading order')
axf.step(rel2, yy, where='post', color=S.INK, lw=1.1, ls=(0, (3, 1.5)), label='+ second order')
axf.axvline(0.01, color=S.GRID, lw=S.LW_GUIDE, zorder=0)
axf.set_xscale('log'); axf.set_xlim(1e-8, 0.3); axf.set_ylim(0, 1.03)
axf.set_xticks([1e-8, 1e-6, 1e-4, 1e-2])
S.xlabel(axf, r'Relative deviation'); S.ylabel(axf, 'Cumulative fraction', labelpad=1.0)
axf.annotate('', xy=(rel.max(), 0.995), xytext=(rel.max(), 0.52),
             arrowprops=dict(arrowstyle='-', lw=0.6, color=S.INK2, shrinkB=1.0))
S.ann(axf, 0.2, 0.42, f'Max {rel.max() * 100:.0f}%\n(graphene, $b=16$)', ha='right', va='center')
S.callout(axf, 1.4e-8, 0.97, f'Median {np.median(rel) * 100:.2f}%\n{(rel > 0.01).sum()} of {len(rel)} above 1%',
          va='top', ha='left', fontsize=S.FA)
axf.legend(loc='lower right', fontsize=S.FS, handlelength=1.6, borderaxespad=0.1)

out = os.path.join(S.OUT, 'fig3_decomposition')
S.save(fig, out)
print('saved', out, 'n', len(allr), 'median rel', np.median(rel), 'max rel', rel.max(), 'max rel2', rel2.max(),
      'I1 cold range', ii[cold].min(), ii[cold].max(), 'I1 all max', ii.max(), 'viol', (bb > 1).sum())

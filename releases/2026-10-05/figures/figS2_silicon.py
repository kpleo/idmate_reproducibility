"""Fig. S2 | Silicon equation of state and bands: IDMate (HGH-LDA) versus VASP (PAW-LDA).
Shape-level comparison (different operators, smearing and meshes; see text).
Drawn at the printed SM width (6.5 in) in the letter style of prb_style.py."""
import json, os, sys
import numpy as np
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prb_style as S

d = json.load(open(os.path.join(S.DATA, 'silicon_physics.json')))


def bm3(V, E0, V0, B0, Bp):
    eta = (V0 / V) ** (2.0 / 3.0)
    return E0 + 9.0 * V0 * B0 / 16.0 * ((eta - 1) ** 3 * Bp + (eta - 1) ** 2 * (6 - 4 * eta))


W, H = 6.5, 3.55
fig, L = S.figure(W, H)
xa, xb, wa, wb = 0.62, 3.78, 2.45, 2.55
y_lo, h_lo, y_hi, h_hi = 0.42, 0.66, 1.20, 1.88
a = L.ax(xa, y_hi, wa, h_hi); c = L.ax(xa, y_lo, wa, h_lo, sharex=a)
b = L.ax(xb, y_hi, wb, h_hi); e = L.ax(xb, y_lo, wb, h_lo, sharex=b)
L.head(0.0, y_hi + h_hi + 0.2, 'a', 'Equation of state')
L.head(3.16, y_hi + h_hi + 0.2, 'b', r'Bands along L$-\Gamma-$X')
S.ann(c, 0.02, 1.0, 'c', transform=c.transAxes, ha='left', va='top', fontsize=S.FL, fontweight='bold', color=S.INK)
S.ann(e, 0.02, 1.0, 'd', transform=e.transAxes, ha='left', va='top', fontsize=S.FL, fontweight='bold', color=S.INK)

STY = {'idmate': dict(color=S.BLUE, marker='o', ls='-', filled=True, label='IDMate (HGH-LDA)'),
       'vasp': dict(color=S.INK, marker='s', ls=(0, (4, 2)), filled=False, label='VASP (PAW-LDA)')}


def mk(st):
    return dict(mfc=st['color'], mec='white', mew=0.35) if st['filled'] else dict(mfc='white', mec=st['color'], mew=0.9)


for key in ('idmate', 'vasp'):
    rec = d['eos'][key]; fit = rec['fit']; st = STY[key]
    V = np.array(rec['volumes']); E = np.array(rec['energies'])
    Vg = np.linspace(V.min() - 0.3, V.max() + 0.3, 400)
    Ef = bm3(Vg, fit['E0_eV'], fit['V0_A3'], fit['B0_eV_A3'], fit['B_prime'])
    a.plot(Vg, (Ef - fit['E0_eV']) * 1000 / 2, color=st['color'], ls=st['ls'], lw=S.LW_DATA)
    a.plot(V, (E - fit['E0_eV']) * 1000 / 2, ls='', marker=st['marker'], ms=4.4, zorder=4, **mk(st))
    res = (E - bm3(V, fit['E0_eV'], fit['V0_A3'], fit['B0_eV_A3'], fit['B_prime'])) * 1000 / 2
    c.plot(V, res, ls='', marker=st['marker'], ms=4.2, zorder=4, **mk(st))
fi, fv = d['eos']['idmate']['fit'], d['eos']['vasp']['fit']
tx = dict(transform=a.transAxes, va='top', fontsize=S.FA)
COLS = {'idmate': (0.50, 0.565), 'vasp': (0.765, 0.81)}      # (header text x, value column centre)
for k, (xt, xv) in COLS.items():
    st = STY[k]
    a.text(xt, 0.97, 'IDMate' if k == 'idmate' else 'VASP', ha='left', color=st['color'], fontweight='bold', **tx)
    a.plot([xt - 0.075, xt - 0.015], [0.945, 0.945], transform=a.transAxes, color=st['color'], ls=st['ls'],
           lw=S.LW_DATA, clip_on=False)
    a.plot([xt - 0.045], [0.945], transform=a.transAxes, ls='', marker=st['marker'], ms=4.2, clip_on=False, **mk(st))
for (lab, key, fmt, y) in ((r'$a_0$ (Å)', 'a0_A', '.3f', 0.87), (r'$B_0$ (GPa)', 'B0_GPa', '.1f', 0.77)):
    a.text(0.40, y, lab, ha='right', color=S.INK2, **tx)
    a.text(COLS['idmate'][1], y, format(fi[key], fmt), ha='center', color=S.BLUE, **tx)
    a.text(COLS['vasp'][1], y, format(fv[key], fmt), ha='center', color=S.INK, **tx)
S.ann(a, 0.5, 0.42, f"Shape RMS {d['eos_shape']['shape_rms_meV_per_atom']:.2f} meV/atom\n(41 scaled volumes)",
      transform=a.transAxes, ha='center', va='center')
S.ylabel(a, r'$E-E_0$ (meV/atom)', labelpad=1.0); a.tick_params(labelbottom=False)
a.set_ylim(-8, 235)
c.axhline(0, color=S.INK3, lw=S.LW_GUIDE, ls=(0, (1, 1.5)))
c.set_ylim(-0.35, 0.35); c.set_yticks([-0.2, 0, 0.2])
S.xlabel(c, r'Cell volume (Å$^3$)'); S.ylabel(c, 'Residual\n(meV/atom)', labelpad=1.0)
B = d['bands']
cx = np.array(B['cx']); rx = np.array(B['rx'])
cg = np.array(B['candidate_grid_ang_inv']); rg = np.array(B['reference_grid_ang_inv'])
diffs = []
for n, (cand, ref, refall) in enumerate(zip(B['candidate'], B['reference'], B['reference_all'])):
    b.plot(cx, cand, color=S.BLUE, lw=S.LW_DATA)
    b.plot(rx, ref, color=S.INK, lw=1.0, ls=(0, (4, 2)))
    diffs.append(np.array(cand) - np.interp(cg, rg, refall))
diffs = np.array(diffs)
b.axvline(1.0, color=S.GRID, lw=S.LW_GUIDE, zorder=0); e.axvline(1.0, color=S.GRID, lw=S.LW_GUIDE, zorder=0)
b.axhline(0, color=S.INK3, lw=S.LW_GUIDE, ls=(0, (1, 1.5)))
S.ylabel(b, r'$\varepsilon-\varepsilon_\mathrm{VBM}$ (eV)', labelpad=1.0); b.tick_params(labelbottom=False)
S.ann(b, 0.03, 0.96, r'$a=5.431$ Å, five bands nearest the gap', transform=b.transAxes, va='top', ha='left')
b.legend(handles=[Line2D([], [], color=S.BLUE, lw=S.LW_DATA, label='IDMate'),
                  Line2D([], [], color=S.INK, lw=1.0, ls=(0, (4, 2)), label='VASP')], loc='upper left',
         bbox_to_anchor=(0.0, 0.47), fontsize=S.FS, handlelength=2.4, borderaxespad=0.1)
b.set_ylim(-8.4, 5.6)
lo, hi = diffs.min(axis=0), diffs.max(axis=0)
e.fill_between(cx, lo, hi, color=S.BLUE, alpha=0.18, lw=0)
e.plot(cx, np.mean(np.abs(diffs), axis=0), color=S.BLUE, lw=S.LW_DATA)
e.axhline(0, color=S.INK3, lw=S.LW_GUIDE, ls=(0, (1, 1.5)))
e.set_xticks([0, 1, 2]); e.set_xticklabels(['L', r'$\Gamma$', 'X']); e.set_xlim(0, 2)
for t in e.get_xticklabels():
    t.set_fontweight('bold')
S.ylabel(e, r'$\Delta\varepsilon$ (eV)', labelpad=1.0); e.set_ylim(-0.02, 0.14); e.set_yticks([0, 0.1])
S.ann(e, 0.98, 0.94, f"MAE {np.mean(np.abs(diffs)):.3f} eV, max {np.max(np.abs(diffs)):.3f} eV",
      transform=e.transAxes, va='top', ha='right')
out = os.path.join(S.OUT, 'figS2_silicon')
S.save(fig, out)
print('MAE', np.mean(np.abs(diffs)), 'max', np.max(np.abs(diffs)), 'n', diffs.size, 'stored', B['metrics']['mae_eV'],
      B['metrics']['max_error_eV'])

"""Fig. 2 | In the SCF tests the window screen cannot see the full-state error; the computable indicator can.
(a) Window distance versus full-state distance for the 58 screened in-loop candidates (+8 abstentions).
(b) The same candidates' full-state distance against the indicator I1, evaluated on the Hamiltonians saved in the
    2026-10-05 rerun (Rayleigh-Ritz families; the coarse-k candidates are not compressed states and are omitted).
(c) Full-state distance versus insertion iteration: the error follows the candidate subspace.
Drawn at printed size in the letter style of prb_style.py."""
import csv, json, os, sys
import numpy as np
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prb_style as S

rows = list(csv.DictReader(open(os.path.join(S.DATA, 'material_candidates.csv'))))
abst = list(csv.DictReader(open(os.path.join(S.DATA, 'excluded_candidates.csv'))))
INL = json.load(open(os.path.join(S.DATA, 'inloop_results.json')))
FAMK = {'low_cutoff_continuation': 'low_cutoff', 'converged_subspace_reuse': 'target_subspace',
        'coarse_k_mesh_upsample': 'coarse_k'}
MK = {'low_cutoff': 'o', 'target_subspace': 'D', 'coarse_k': 's'}
MS = {'low_cutoff': 4.2, 'target_subspace': 3.6, 'coarse_k': 3.9}
FAM = {'low_cutoff': 'Low cutoff', 'target_subspace': 'Oracle subspace', 'coarse_k': r'Coarse $k$'}
LS = {'low_cutoff': '-', 'target_subspace': (0, (3.2, 1.6)), 'coarse_k': (0, (1.1, 1.3))}
TOL = 0.05

W, H = S.W2, 4.42
fig, L = S.figure(W, H)
ya, ha = 2.50, 1.50
axa = L.ax(0.56, ya, 2.36, ha)
axs = L.ax(0.56 + 2.36 + 0.06, ya, 0.24, ha)
axb = L.ax(4.18, ya, 2.36, ha)
yc, hc = 0.42, 1.12
xc = [0.56, 2.72, 4.88]
axc = [L.ax(x, yc, 1.95, hc) for x in xc]
L.head(0.0, ya + ha + 0.2, 'a', 'What the window screen certifies')
L.head(3.62, ya + ha + 0.2, 'b', r'What the indicator $I_1$ sees')
L.head(0.0, yc + hc + 0.38, 'c', 'Full-state error follows the candidate subspace')


def face(acc, col):
    return dict(mfc=col if acc else 'white', mec='white' if acc else col, mew=0.35 if acc else 0.9)


# ---------------------------------------------------------------- (a) window vs full
axa.set_xscale('log'); axa.set_yscale('log')
axa.axvspan(1e-13, TOL, color=S.PALE, lw=0, zorder=0)
axa.axvline(TOL, color=S.INK, lw=S.LW_GUIDE, ls=(0, (3, 2)))
axa.axhline(TOL, color=S.INK, lw=S.LW_GUIDE, ls=(0, (3, 2)))
xx = np.array([1.5e-2, 3]); axa.plot(xx, xx, color=S.INK3, lw=S.LW_GUIDE, zorder=1)
for r in rows:
    x, y = float(r['window_distance']), float(r['full_distance'])
    fam = r['proposal_family']; acc = r['decision'] == 'accept'
    axa.plot(x, y, MK[fam], ms=MS[fam], zorder=3, **face(acc, S.MAT[r['material']]))
rng = np.random.default_rng(3)
for r in abst:
    fam = r['proposal_family']
    axs.plot(float(np.clip(0.5 + 0.2 * rng.standard_normal(), 0.15, 0.85)), float(r['full_distance']), MK[fam],
             ms=MS[fam], **face(False, S.MAT[r['material']]))
for a in (axa, axs):
    a.set_yscale('log'); a.set_ylim(1.5e-5, 4)
axa.set_xlim(2e-12, 4)
axa.set_xticks([1e-10, 1e-7, 1e-4, 1e-1])
S.xlabel(axa, r'Window distance $d_W$'); S.ylabel(axa, r'Full-state distance $d_\mathrm{full}$', labelpad=1.0)
axs.set_xlim(0, 1); axs.set_xticks([0.5]); axs.set_xticklabels(['Abst.'])
axs.tick_params(axis='y', which='both', left=False, labelleft=False); axs.spines['left'].set_visible(False)
nbad = sum(1 for r in rows if r['decision'] == 'accept' and float(r['full_distance']) > TOL)
nacc = sum(1 for r in rows if r['decision'] == 'accept')
S.callout(axa, 3e-8, 0.30, f'{nbad} of {nacc} accepted\nhave $d_\\mathrm{{full}}>0.05$', ha='left', va='center',
          color=S.FAIL)
S.ann(axa, 4e-9, 2.5, f'Accepted: {nacc}', ha='left', va='center', color=S.INK)
S.ann(axa, 0.55, 2.2e-2, f'Rejected:\n{len(rows) - nacc}', ha='center', va='top', color=S.INK)
S.ann(axa, 0.041, 2.2e-5, 'Tolerance', rotation=90, ha='right', va='bottom')
S.ann(axa, 1.0e-10 / 1.8, 2.2e-3, 'Si', ha='right', va='center', color=S.MAT['Si'], fontweight='bold')
S.ann(axa, 1.0e-9 * 2.0, 2.2e-3, 'Al', ha='left', va='center', color=S.MAT['Al'], fontweight='bold')
S.ann(axa, 2.5e-12, 2.5, 'Graphene', ha='left', va='center', color=S.MAT['graphene'], fontweight='bold')

# ---------------------------------------------------------------- (b) indicator vs full (in-loop rerun)
axb.set_xscale('log'); axb.set_yscale('log')
lo, hi = 1.5e-5, 12
axb.fill_between([lo, hi], [lo, hi], [hi, hi], color='#fbe3e3', lw=0, zorder=0)      # y > x: I1 would underestimate
axb.plot([lo, hi], [lo, hi], color=S.INK3, lw=S.LW_GUIDE, zorder=1)
axb.axvline(TOL, color=S.INK, lw=S.LW_GUIDE, ls=(0, (3, 2)))
axb.axhline(TOL, color=S.INK, lw=S.LW_GUIDE, ls=(0, (3, 2)))
ritz = [r for r in INL if r['ritz']]
for r in ritz:
    fam = FAMK[r['kind']]; acc = r['decision'] == 'accept'
    axb.plot(r['indicator'], r['d_full'], MK[fam], ms=MS[fam], zorder=3, **face(acc, S.MAT[r['material']]))
axb.set_xlim(lo, hi); axb.set_ylim(lo, 4)
axb.set_xticks([1e-4, 1e-2, 1]); axb.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
S.xlabel(axb, r'Indicator $I_1$'); S.ylabel(axb, r'Full-state distance $d_\mathrm{full}$', labelpad=1.0)
acc_r = [r for r in ritz if r['decision'] == 'accept']
bad = [r for r in acc_r if r['d_full'] > TOL]
caught = sum(1 for r in bad if r['indicator'] > TOL)
kept = sum(1 for r in acc_r if r['d_full'] <= TOL and r['indicator'] <= TOL)
good = sum(1 for r in acc_r if r['d_full'] <= TOL)
under = sum(1 for r in ritz if r['indicator'] < r['d_full'])
S.ann(axb, 2.2e-5, 2.2, r'$I_1<d_\mathrm{full}$' + f' (none of {len(ritz)})', ha='left', va='top', color=S.FAIL)
S.ann(axb, 2.2e-3, 3.4e-3, r'$I_1=d_\mathrm{full}$', ha='center', va='bottom', rotation=45, rotation_mode='anchor')
S.callout(axb, 9.5, 2.2e-5, f'{caught}/{len(bad)} bad accepts flagged\n{kept}/{good} good accepts kept',
          ha='right', va='bottom', fontsize=S.FA, bbox=dict(boxstyle='square,pad=0.12', fc='white', ec='none'))
assert under == 0, under

# legends (families, materials, decision) shared by (a) and (b)
h1 = [Line2D([], [], ls='', marker=MK[k], ms=MS[k], mfc=S.INK2, mec='white', mew=0.35, label=FAM[k]) for k in MK]
h2 = [Line2D([], [], ls='', marker='o', ms=4.2, mfc=S.MAT[k], mec='white', mew=0.35, label=S.MAT_LABEL[k])
      for k in S.MAT]
h3 = [Line2D([], [], ls='', marker='o', ms=4.2, mfc=S.INK2, mec='white', mew=0.35, label='Accepted'),
      Line2D([], [], ls='', marker='o', ms=4.2, mfc='white', mec=S.INK2, mew=0.9, label='Rejected / abstained')]
axa.legend(handles=h1 + h3, loc='lower left', fontsize=S.FS, ncol=1, handletextpad=0.25,
           bbox_to_anchor=(0.335, -0.01), borderaxespad=0.1, labelspacing=0.22)

# ---------------------------------------------------------------- (c) error versus insertion iteration
allr = rows + [dict(r, window_distance='nan') for r in abst]
for a, mat in zip(axc, ['Si', 'graphene', 'Al']):
    a.axhline(TOL, color=S.INK, lw=S.LW_GUIDE, ls=(0, (3, 2)), zorder=0)
    for fam in ['low_cutoff', 'coarse_k', 'target_subspace']:
        sel = [r for r in allr if r['material'] == mat and r['proposal_family'] == fam]
        single = sorted([r for r in sel if 'single' in r['run_id']], key=lambda r: int(r['iteration']))
        a.plot([int(r['iteration']) for r in single], [float(r['full_distance']) for r in single], color=S.MAT[mat],
               lw=1.0, ls=LS[fam], zorder=1)
        for r in single:
            a.plot(int(r['iteration']), float(r['full_distance']), MK[fam], ms=MS[fam] - 0.4, zorder=3,
                   **face(r['decision'] == 'accept', S.MAT[mat]))
        for r in sel:
            if 'single' in r['run_id']:
                continue
            a.plot(int(r['iteration']) + 0.45, float(r['full_distance']), MK[fam], ms=MS[fam] - 1.4, zorder=2,
                   alpha=0.85, **face(r['decision'] == 'accept', S.MAT[mat]))
    a.set_yscale('log'); a.set_xlim(0.5, 19.5); a.set_xticks([2, 6, 10, 14, 18]); a.set_ylim(1.5e-5, 4)
    a.set_yticks([1e-4, 1e-2, 1])
    a.text(0.5, 1.03, S.MAT_LABEL[mat], transform=a.transAxes, ha='center', va='bottom', fontsize=7.4,
           fontweight='bold', color=S.MAT[mat])
    if a is not axc[0]:
        a.tick_params(labelleft=False)
S.ylabel(axc[0], r'$d_\mathrm{full}$', labelpad=1.0)
S.xlabel(axc[1], 'SCF iteration at which the candidate is inserted')
hb = [Line2D([], [], color=S.INK2, lw=1.0, ls=LS[k], marker=MK[k], ms=MS[k] - 0.6, mfc=S.INK2, mec='white',
             mew=0.35, label=FAM[k]) for k in ['low_cutoff', 'coarse_k', 'target_subspace']]
axc[0].legend(handles=hb, loc='lower left', fontsize=S.FS, handlelength=2.2, borderaxespad=0.1,
              bbox_to_anchor=(-0.02, -0.02))
out = os.path.join(S.OUT, 'fig2_inloop')
S.save(fig, out)
print('saved', out, 'accepted with dfull>0.05:', nbad, '| in-loop Ritz candidates', len(ritz), 'accepted', len(acc_r),
      'bad', len(bad), 'caught', caught, 'good', good, 'kept', kept, 'underestimates', under)

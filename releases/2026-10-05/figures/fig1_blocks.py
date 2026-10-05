"""Fig. 1 | A window certificate is sharp, but it sees only one block of the error.
(a) |H| for silicon (k-point 1) in the basis [Ritz vectors of P | eigenvectors of QHQ].
(b) |D_c - D*| in the same basis: the error sits in the coupling blocks.
(c) ECDFs of distance/bound for the matrix and material verification sets.
Drawn at printed size in the letter style of prb_style.py."""
import csv, os, sys
import numpy as np
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from matplotlib.patches import Rectangle
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prb_style as S

blk = S.load_blocks('si')
m, n = int(blk['m']), int(blk['n'])


def ratios(name):
    rows = list(csv.DictReader(open(os.path.join(S.DATA, name + '.csv'))))
    r = np.array([float(x['ratio']) for x in rows if x['ratio'] not in ('', None)])
    return np.sort(r), len(rows)


def sci(x, digits=1):
    e = int(np.floor(np.log10(abs(x))))
    return f'${x / 10 ** e:.{digits}f}\\times10^{{{e}}}$'


W, H = S.W2, 2.62
fig, L = S.figure(W, H)
y0, s = 0.43, 1.66
axa = L.ax(0.47, y0, s, s)
caa = L.ax(0.47 + s + 0.06, y0, 0.065, s)
axb = L.ax(2.92, y0, s, s)
cab = L.ax(2.92 + s + 0.06, y0, 0.065, s)
axc = L.ax(5.50, y0, 1.42, s)
yh = y0 + s + 0.2
L.head(0.0, yh, 'a', r'Si Hamiltonian, $[P\,|\,Q]$ basis')
L.head(2.50, yh, 'b', r'Error $\Delta D=D_c-D^\star$')
L.head(5.02, yh, 'c', 'Window certificate')

greys = LinearSegmentedColormap.from_list('g', ['#FFFFFF', '#C9C8C4', '#52514E', '#0B0B0B'])
warm = LinearSegmentedColormap.from_list('w', ['#FFFFFF', '#D9D3F0', S.LEAD, '#1E1650'])


def heat(ax, cax, M, cmap, vmin, vmax, cbar_label):
    im = ax.imshow(np.maximum(M, vmin), norm=LogNorm(vmin=vmin, vmax=vmax), cmap=cmap, interpolation='nearest',
                   extent=(0, n, n, 0))
    ax.axhline(m, color=S.INK2, lw=S.LW_THIN, ls=(0, (2, 1.5)))
    ax.axvline(m, color=S.INK2, lw=S.LW_THIN, ls=(0, (2, 1.5)))
    for (x0, yy, w, h) in [(0, m, m, n - m), (m, 0, n - m, m)]:
        ax.add_patch(Rectangle((x0, yy), w, h, fill=False, ec=S.LEAD, lw=1.1))
    ax.set_xticks([0, m, n]); ax.set_yticks([0, m, n])
    S.frame_box(ax)
    cb = fig.colorbar(im, cax=cax)
    cb.outline.set_linewidth(S.LW_AX); cb.outline.set_edgecolor(S.INK2)
    cb.ax.tick_params(labelsize=S.FTK, length=2.2, width=S.LW_AX)
    if cbar_label:
        cb.set_label(S.boldmath(cbar_label), fontsize=S.FAX, fontweight='bold', labelpad=1.5)
    return im


heat(axa, caa, blk['H_abs'], greys, 1e-4, 10, r'$|H_{ij}|$ (Ha)')
S.xlabel(axa, 'Basis index'); S.ylabel(axa, 'Basis index', labelpad=1.0)
S.ann(axa, m / 2, m / 2, r'$H_P$', ha='center', va='center', color=S.INK, fontsize=S.FK)
S.ann(axa, (m + n) / 2 + 4, (m + n) / 2 + 13, r'$QHQ$', ha='center', va='center', color=S.INK, fontsize=S.FK)
S.ann(axa, (m + n) / 2, m / 2, r'$PHQ$', ha='center', va='center', color=S.LEAD, fontsize=S.FK, fontweight='bold')
S.ann(axa, m / 2, (m + n) / 2 + 6, r'$QHP$', ha='center', va='center', color=S.LEAD, fontsize=S.FK,
      fontweight='bold', rotation=90)

heat(axb, cab, blk['dD_abs'], warm, 1e-6, 0.3, r'$|\Delta D_{ij}|$')
S.xlabel(axb, 'Basis index')
PP = np.linalg.norm(blk['dD_abs'][:m, :m]); QQ = np.linalg.norm(blk['dD_abs'][m:, m:])
PQ = np.sqrt(2) * np.linalg.norm(blk['dD_abs'][m:, :m])
S.callout(axb, (m + n) / 2 + 2, m / 2 + 1, 'Coupling blocks: ' + f'{PQ:.3f}', ha='center', va='center',
          color=S.LEAD, fontsize=S.FA + 0.2)
S.ann(axb, (m + n) / 2 + 3, (m + n) / 2 + 12, 'Diagonal blocks\n$P$: ' + sci(PP) + '\n$Q$: ' + sci(QQ),
      ha='center', va='center', color=S.INK)

sets = [('e0', 'Random', S.INK3, '-'),
        ('e3_graphene', 'Graphene, Dirac', S.MAT['graphene'], '-'),
        ('e3_al', 'Al', S.MAT['Al'], '-'),
        ('gk_graphene', r'Graphene, shared $\mu$', S.MAT['graphene'], (0, (3.2, 1.6))),
        ('gk_al', r'Al, shared $\mu$', S.MAT['Al'], (0, (3.2, 1.6)))]
ntot, rmax = 0, 0.0
for key, lab, col, ls in sets:
    r, nrows = ratios(key); ntot += len(r); rmax = max(rmax, r.max())
    y = np.arange(1, len(r) + 1) / len(r)
    axc.step(np.r_[0, r], np.r_[0, y], where='post', color=col, lw=S.LW_DATA if ls == '-' else 1.1, ls=ls,
             label=lab)
    print(key, len(r))
axc.axvline(1.0, color=S.INK, lw=S.LW_GUIDE, ls=(0, (3, 2)))
S.ann(axc, 1.015, 0.36, 'Bound', rotation=90, ha='left', va='center')
axc.set_xlim(0, 1.1); axc.set_ylim(0, 1.03)
axc.set_xticks([0, 0.5, 1.0]); axc.set_yticks([0, 0.5, 1])
S.xlabel(axc, r'Distance / bound'); S.ylabel(axc, 'Cumulative fraction', labelpad=1.0)
axc.legend(loc='upper left', fontsize=S.FS, handlelength=1.6, borderaxespad=0.1, bbox_to_anchor=(-0.02, 1.04))
S.callout(axc, 0.97, 0.05, f'{ntot:,} ratios\nmax {rmax:.3f}', transform=axc.transAxes, ha='right', va='bottom')
out = os.path.join(S.OUT, 'fig1_certificate')
S.save(fig, out)
print('saved', out, 'PP', PP, 'QQ', QQ, 'PQ(both)', PQ, 'ntot', ntot, 'rmax', rmax)

"""Fig. S1 | In all three materials the full-state error of the preceding-map Ritz candidate occupies the
coupling blocks. |D_c - D*| at the first k point in the basis [Ritz vectors of P | eigenvectors of QHQ].
Drawn at the printed SM width (6.5 in) in the letter style of prb_style.py."""
import os, sys
import numpy as np
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from matplotlib.patches import Rectangle
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prb_style as S

warm = LinearSegmentedColormap.from_list('w', ['#FFFFFF', '#D9D3F0', S.LEAD, '#1E1650'])
W, H = 6.5, 2.55
fig, L = S.figure(W, H)
s, y0 = 1.55, 0.43
xs = [0.47, 2.36, 4.25]
axes = [L.ax(x, y0, s, s) for x in xs]
cax = L.ax(4.25 + s + 0.08, y0, 0.07, s)
names = [('si', 'Si'), ('graphene', 'graphene'), ('al', 'Al')]


def sci(x):
    e = int(np.floor(np.log10(abs(x))))
    mant = x / 10 ** e
    return f'${mant:.0f}\\times10^{{{e}}}$' if e < -2 else f'{x:.3f}'

for letter, x, ax, (key, mat) in zip('abc', xs, axes, names):
    b = S.load_blocks(key)
    m, n, M = int(b['m']), int(b['n']), b['dD_abs']
    im = ax.imshow(np.maximum(M, 1e-7), norm=LogNorm(vmin=1e-7, vmax=0.3), cmap=warm, interpolation='nearest',
                   extent=(0, n, n, 0))
    ax.axhline(m, color=S.INK2, lw=S.LW_THIN, ls=(0, (2, 1.5))); ax.axvline(m, color=S.INK2, lw=S.LW_THIN, ls=(0, (2, 1.5)))
    for (x0, yy, w, h) in [(0, m, m, n - m), (m, 0, n - m, m)]:
        ax.add_patch(Rectangle((x0, yy), w, h, fill=False, ec=S.LEAD, lw=1.1))
    S.frame_box(ax)
    ax.set_xticks([0, m, n]); ax.set_yticks([0, m, n])
    S.xlabel(ax, 'Basis index')
    PP = np.linalg.norm(M[:m, :m]); QQ = np.linalg.norm(M[m:, m:]); PQ = np.sqrt(2) * np.linalg.norm(M[m:, :m])
    L.head(x - 0.47, y0 + s + 0.2, letter, f"{S.MAT_LABEL[mat]}, $m={m}$, $n={n}$")
    S.callout(ax, 0.97 * n, 0.58 * n, f'Coupling {PQ:.2g}', ha='right', va='center', color=S.LEAD, fontsize=S.FA + 0.2)
    S.ann(ax, 0.97 * n, 0.95 * n, 'Diagonal ' + sci(PP) + ', ' + sci(QQ), ha='right', va='bottom', color=S.INK)
S.ylabel(axes[0], 'Basis index', labelpad=1.0)
cb = fig.colorbar(im, cax=cax)
cb.outline.set_linewidth(S.LW_AX); cb.outline.set_edgecolor(S.INK2)
cb.ax.tick_params(labelsize=S.FTK, length=2.2, width=S.LW_AX)
cb.set_label(S.boldmath(r'$|(D_c-D^\star)_{ij}|$'), fontsize=S.FAX, fontweight='bold', labelpad=1.5)
S.save(fig, os.path.join(S.OUT, 'figS1_blocks'))
print('ok')

"""Figure style for the IDMate PRB figures (letter style).

Printed sizes (figures are drawn at their printed width, no rescaling, no cropping):
- Liberation Sans for text and math; bold panel letter (10 pt) and bold sentence-case panel title (8.5 pt) on one
  baseline above each panel;
- bold sentence-case axis titles (8 pt) with bold math (variables bold italic); tick labels 7 pt;
- legends and annotations 6.8 pt (annotations in secondary ink with a white halo); key callouts bold 7.6 pt;
- data lines 1.3 pt, guides 0.85 pt, thin lines 0.55 pt; axes and major ticks 0.7 pt;
- filled markers with a thin white ring; open markers with a visible stroke; no grids; top and right spines open on
  line plots; frameless legends.
Widths: PRB double column 7.0 in (178 mm), single column 3.4 in (86 mm).
"""
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch

# ---------------------------------------------------------------- palette
INK, INK2, INK3, GRID = '#0b0b0b', '#52514e', '#8a8984', '#e4e3df'
PALE = '#f4f3ef'
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
BLUE, ORANGE, GREEN, AMBER, PINK, DKGREEN, VIOLET, RED = C
# materials (categorical, paper-wide)
MAT = {'Si': BLUE, 'graphene': GREEN, 'Al': ORANGE}
MAT_LABEL = {'Si': 'Diamond Si', 'graphene': 'Graphene', 'Al': 'fcc Al'}
MAT_MK = {'Si': 'o', 'graphene': 's', 'Al': '^'}
# quantities (paper-wide; never used for materials inside one panel)
MEAS = INK          # measured full-state distance (the reference data)
LEAD = VIOLET       # leading-order decomposition / coupling term
IND = AMBER         # computable indicator I1
BOUND = INK3        # rigorous gap-free bound
DIAG = PINK         # exact diagonal blocks
WINDOW = DKGREEN    # window (compressed) quantities
FAIL = RED          # underestimates, failures

# ---------------------------------------------------------------- printed sizes (pt)
FL, FT, FAX, FTK, FLEG, FA, FK, FS = 10.0, 8.5, 8.0, 7.0, 6.8, 6.8, 7.6, 6.6
LW_DATA, LW_GUIDE, LW_THIN, LW_AX = 1.3, 0.85, 0.55, 0.7
W2, W1 = 7.0, 3.4                     # double / single column width (in)
MKW = dict(mec='white', mew=0.35)     # filled markers: thin white ring

RC = {'font.family': 'Liberation Sans', 'font.size': 7.2, 'axes.labelsize': FAX, 'axes.labelweight': 'bold',
      'axes.titlesize': FT, 'legend.fontsize': FLEG, 'xtick.labelsize': FTK, 'ytick.labelsize': FTK,
      'lines.markersize': 4.0, 'lines.linewidth': LW_DATA, 'axes.linewidth': LW_AX, 'xtick.major.width': LW_AX,
      'ytick.major.width': LW_AX, 'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5, 'xtick.major.size': 3.0,
      'ytick.major.size': 3.0, 'xtick.minor.size': 1.7, 'ytick.minor.size': 1.7, 'xtick.major.pad': 1.8,
      'ytick.major.pad': 1.6, 'axes.labelpad': 2.0, 'legend.handletextpad': 0.45, 'legend.labelspacing': 0.28,
      'legend.borderaxespad': 0.25, 'legend.handlelength': 1.5, 'legend.columnspacing': 1.0, 'legend.frameon': False,
      'mathtext.fontset': 'custom', 'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic',
      'mathtext.bf': 'Liberation Sans:bold', 'mathtext.bfit': 'Liberation Sans:bold:italic',
      'mathtext.sf': 'Liberation Sans', 'mathtext.fallback': 'stixsans', 'axes.edgecolor': INK2,
      'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2, 'xtick.labelcolor': INK2,
      'ytick.labelcolor': INK2, 'text.color': INK, 'axes.spines.top': False, 'axes.spines.right': False,
      'axes.grid': False, 'axes.titlepad': 3.0, 'figure.dpi': 150, 'savefig.dpi': 600, 'pdf.fonttype': 42,
      'ps.fonttype': 42, 'svg.fonttype': 'none', 'savefig.bbox': None, 'axes.unicode_minus': True,
      'image.cmap': 'cividis'}


def use():
    mpl.rcParams.update(RC)
    return 'Liberation Sans'


# ---------------------------------------------------------------- bold math and sentence case (from figstyle.py)
GREEK = set('alpha beta gamma delta epsilon varepsilon zeta eta theta vartheta iota kappa lambda mu nu xi pi varpi rho '
            'varrho sigma varsigma tau upsilon phi varphi chi psi omega Gamma Delta Theta Lambda Xi Pi Sigma Upsilon '
            'Phi Psi Omega'.split())
FUNCS = set('ln log exp sin cos tan min max arg det Re Im sinh cosh tanh'.split())
UPRIGHT = {'rm', 'mathrm', 'mathbf', 'mathsf', 'text'}


def _group(s, i):
    if i < len(s) and s[i] == '{':
        d = 0
        for j in range(i, len(s)):
            d += {'{': 1, '}': -1}.get(s[j], 0)
            if d == 0:
                return s[i + 1:j], j + 1
        return s[i + 1:], len(s)
    if i < len(s) and s[i] == '\\':
        m = re.match(r'\\([A-Za-z]+|.)', s[i:])
        return m.group(0), i + len(m.group(0))
    return (s[i], i + 1) if i < len(s) else ('', i)


def _bm(s):
    """Bold version of a mathtext expression: variables bold italic; digits, operators and upright text bold."""
    out, i = [], 0
    while i < len(s):
        c = s[i]
        if c == '\\':
            m = re.match(r'\\([A-Za-z]+|.)', s[i:])
            name = m.group(1)
            i += len(m.group(0))
            if name in GREEK:
                out.append((r'\mathbf{\%s}' if name[0].isupper() else r'\mathbfit{\%s}') % name)
            elif name in FUNCS:
                pre = r'\,' if out and out[-1].endswith('}') and not out[-1].endswith(r'\,') else ''
                out.append(pre + r'\mathbf{%s}\,' % name)
            elif name in UPRIGHT or name in ('mathbfit', 'mathit'):
                arg, i = _group(s, i)
                if name in ('mathit', 'mathbfit'):
                    out.append(r'\mathbfit{%s}' % arg)
                else:
                    out.append(r'\mathbf{%s}' % arg)
            elif name in ('frac', 'dfrac'):
                a, i = _group(s, i)
                b, i = _group(s, i)
                out.append(r'\frac{%s}{%s}' % (_bm(a), _bm(b)))
            elif name == 'sqrt':
                a, i = _group(s, i)
                out.append(r'\sqrt{%s}' % _bm(a))
            else:
                out.append(m.group(0))
        elif c in '^_':
            a, i = _group(s, i + 1)
            out.append(c + '{' + _bm(a) + '}')
        elif c == '{':
            a, i = _group(s, i)
            out.append('{' + _bm(a) + '}')
        elif c.isalpha() and c.isascii():
            out.append(r'\mathbfit{%s}' % c)
            i += 1
        elif c.isdigit() or c == '.':
            j = i
            while j < len(s) and (s[j].isdigit() or s[j] == '.'):
                j += 1
            out.append(r'\mathbf{%s}' % s[i:j])
            i = j
        elif c in '()[]|+-=/,<>!:;*':
            out.append(r'\mathbf{%s}' % c)
            i += 1
        else:
            out.append(c)
            i += 1
    return ''.join(out)


def boldmath(t):
    """Bold math in every $...$ segment of a label (axis titles are bold, math included)."""
    t = re.sub(r'([_^])(\\[A-Za-z]+\{[^{}]*\})', r'\1{\2}', t)     # d_\mathrm{full} -> d_{\mathrm{full}}
    parts = re.split(r'(?<!\\)\$', t)
    return ''.join('$' + _bm(p) + '$' if k % 2 else p for k, p in enumerate(parts))


def xlabel(ax, t, **kw):
    ax.set_xlabel(boldmath(t), fontweight='bold', **kw)


def ylabel(ax, t, **kw):
    ax.set_ylabel(boldmath(t), fontweight='bold', **kw)


# ---------------------------------------------------------------- layout in inches, panel heads, annotations
class Lay:
    """Inch-based placement on a figure of W x H inches (origin bottom left)."""

    def __init__(self, fig):
        self.fig = fig
        self.W, self.H = fig.get_size_inches()

    def ax(self, x, y, w, h, **kw):
        return self.fig.add_axes([x / self.W, y / self.H, w / self.W, h / self.H], **kw)

    def text(self, x, y, s, **kw):
        return self.fig.text(x / self.W, y / self.H, s, **kw)

    def head(self, x, y, letter, title=None):
        """Panel letter (bold, 10 pt) and title (bold, 8.5 pt) on one baseline at height y (in)."""
        self.text(x, y, letter, fontsize=FL, fontweight='bold', va='baseline', ha='left', color=INK)
        if title:
            self.text(x + 0.16, y, boldmath(title), fontsize=FT, fontweight='bold', va='baseline', ha='left',
                      color=INK)

    def arrow(self, p0, p1, color=INK, lw=0.8):
        self.fig.add_artist(FancyArrowPatch((p0[0] / self.W, p0[1] / self.H), (p1[0] / self.W, p1[1] / self.H),
                                            transform=self.fig.transFigure,
                                            arrowstyle='-|>,head_length=3.0,head_width=1.6', color=color, lw=lw,
                                            shrinkA=0, shrinkB=0))


def figure(w, h):
    use()
    fig = plt.figure(figsize=(w, h))
    return fig, Lay(fig)


def halo(t, color='white', lw=2.0):
    t.set_path_effects([pe.withStroke(linewidth=lw, foreground=color)])
    return t


def ann(ax, x, y, s, **kw):
    """Annotation: 6.8 pt, secondary ink, white halo."""
    kw.setdefault('fontsize', FA)
    kw.setdefault('color', INK2)
    kw.setdefault('zorder', 50)
    return halo(ax.text(x, y, s, **kw), 'white', 2.0)


def callout(ax, x, y, s, color=INK, **kw):
    """Key callout: bold 7.6 pt with a white halo; math in s is made bold."""
    kw.setdefault('fontsize', FK)
    kw.setdefault('zorder', 60)
    return halo(ax.text(x, y, boldmath(s), color=color, fontweight='bold', **kw), 'white', 2.2)


def frame_box(ax):
    """Closed frame for image-like panels (heat maps)."""
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(LW_AX)
        sp.set_color(INK2)


def save(fig, path):
    fig.savefig(f'{path}.pdf', bbox_inches=None)
    fig.savefig(f'{path}.png', dpi=300, bbox_inches=None)
    plt.close(fig)


# ---------------------------------------------------------------- release paths (reproducibility package)
import json as _json
import os as _os
import numpy as _np
DATA = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'data')
OUT = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'build', 'figures')
_os.makedirs(OUT, exist_ok=True)


def load_blocks(material):
    """Block matrices of Fig. 1 and Fig. S1 (|H| and |D_c - D*| in the [P | Q] eigenbases), as numpy arrays."""
    d = _json.load(open(_os.path.join(DATA, f'blocks_{material}.json')))
    return {'m': d['m'], 'n': d['n'], 'H_abs': _np.array(d['H_abs_ha']), 'dD_abs': _np.array(d['dD_abs'])}

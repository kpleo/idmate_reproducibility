from __future__ import annotations
from itertools import combinations
import math
import numpy as np
from matplotlib import colors, font_manager
from matplotlib.lines import Line2D
from matplotlib.text import Text
from style import COLORS, MM, WIDTH_MM, plt, save_figure, use_style
from crystal_motifs import CELLS, draw_motif, read_cell
from common import close
HEIGHT_MM = 138
TRACKS_MM = {'a': (16, 43, 64, 60), 'b': (99, 43, 80, 60), 'c': (16, 12, 64, 22), 'd': (99, 12, 80, 22)}

def draw(data):
    use_style()
    plt.rcParams.update({'svg.hashsalt': 'idmate-fig3-physics-2026-09-05', 'figure.dpi': 160, 'path.simplify': False})
    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM))
    axes = {key: fig.add_axes((x / WIDTH_MM, y / HEIGHT_MM, w / WIDTH_MM, h / HEIGHT_MM)) for (key, (x, y, w, h)) in TRACKS_MM.items()}
    for (key, ax) in axes.items():
        (x, y, w, h) = TRACKS_MM[key]
        fig.text((x - 12) / WIDTH_MM, (y + h + 0.4) / HEIGHT_MM, key, fontsize=10, va='top', weight='bold')
        ax.tick_params(top=False, right=False, pad=2.5)
    handles = [Line2D([], [], color=COLORS['blue'], lw=1.3, marker='o', markersize=3.6, label='IDMate (HGH)'), Line2D([], [], color=COLORS['ink'], lw=1.0, ls=(0, (4, 2)), marker='s', markerfacecolor='white', markersize=3.6, label='VASP-LDA (PAW)')]
    fig.legend(handles=handles, loc='center', bbox_to_anchor=(0.72, 122 / HEIGHT_MM), ncol=1, handlelength=2.9, labelspacing=1.0, borderaxespad=0)
    motif = draw_motif(fig, (9, 107, 31, 29), 'si_bands', COLORS['blue'])
    fig.text(43 / WIDTH_MM, 122 / HEIGHT_MM, 'Diamond Si', fontsize=10, va='center')
    fig._crystal_record = motif
    (a, b, c, d) = (axes[key] for key in 'abcd')
    for (key, color, marker, linestyle) in (('idmate', COLORS['blue'], 'o', '-'), ('vasp', COLORS['ink'], 's', (0, (4, 2)))):
        record = data['eos'][key]
        fit = record['fit']
        volume = np.linspace(record['volumes'].min(), record['volumes'].max(), 401)
        curve = data['bm3'](volume, fit['E0_eV'], fit['V0_A3'], fit['B0_eV_A3'], fit['B_prime'])
        a.plot(volume, (curve - fit['E0_eV']) * 500, color=color, ls=linestyle, lw=1.2)
        style = {'color': color, 'marker': marker, 'ms': 4.0, 'mfc': color if key == 'idmate' else 'white', 'mew': 0.9}
        a.plot(record['volumes'], record['shifted'], ls='none', **style)
        c.plot(record['volumes'], record['residual_meV_atom'], ls='none', **style)
    a.set(xlim=(32.3, 48.7), ylim=(-7, 237), yticks=[0, 50, 100, 150, 200])
    a.set_ylabel('$E-E_0$ (meV/atom)', labelpad=5)
    for ax in (a, c):
        ax.set_xlim(32.3, 48.7)
        ax.set_xticks([34, 38, 42, 46])
    a.tick_params(labelbottom=False)
    for (text, x) in (('IDMate', 0.52), ('VASP', 0.8)):
        a.text(x, 0.91, text, transform=a.transAxes, ha='center', fontsize=7)
    for (label, field, format_, y) in (('$a_0$ (Å)', 'a0_A', '.3f', 0.81), ('$B_0$ (GPa)', 'B0_GPa', '.1f', 0.72)):
        a.text(0.15, y, label, transform=a.transAxes, ha='left', fontsize=7.2)
        for (key, x) in (('idmate', 0.52), ('vasp', 0.8)):
            a.text(x, y, format(data['eos'][key]['fit'][field], format_), transform=a.transAxes, ha='center', fontsize=7.2)
    a.text(0.5, 0.56, f"Shape RMS  {data['eos_shape']['shape_rms_meV_per_atom']:.3f} meV/atom", transform=a.transAxes, fontsize=7.1, ha='center')
    c.axhline(0, color=COLORS['ink'], lw=0.6, ls=':', zorder=0)
    c.set(ylim=(-0.3, 0.3), yticks=[-0.2, 0, 0.2])
    c.set_xlabel('Cell volume (Å$^3$)', labelpad=3)
    c.set_ylabel('$E-E_{\\mathrm{BM3}}$\n(meV/atom)', labelpad=4)
    fig.text(48 / WIDTH_MM, 37.5 / HEIGHT_MM, '5 points per method', fontsize=7, ha='center')
    bands = data['bands']
    for (candidate, reference) in zip(bands['candidate'], bands['reference']):
        b.plot(bands['cx'], candidate, color=COLORS['blue'], lw=1.2)
        b.plot(bands['rx'], reference, color=COLORS['ink'], ls=(0, (4, 2)), lw=0.85)
    for ax in (b, d):
        ax.set_xlim(-0.025, 2.025)
        ax.set_xticks([0, 1, 2], ['L', 'Γ', 'X'])
        ax.axvline(1, color=COLORS['reference'], lw=0.55, zorder=0)
    b.axhline(0, color=COLORS['ink'], lw=0.7, ls=':', zorder=0)
    b.set(ylim=(-8.4, 5.5), yticks=[-8, -6, -4, -2, 0, 2, 4])
    b.tick_params(labelbottom=False)
    b.set_ylabel('$\\varepsilon-\\varepsilon_{\\mathrm{VBM}}$ (eV)', labelpad=4)
    b.text(0.035, 0.94, 'Si, $a=5.431$ Å', transform=b.transAxes, fontsize=7.4)
    b.text(0.965, 0.94, '5 bands', transform=b.transAxes, fontsize=7.1, ha='right')
    difference_styles = [(COLORS['blue'], (0, (1, 1)), ''), (COLORS['blue'], (0, (4, 2)), ''), (COLORS['green'], '-', ''), (COLORS['orange'], (0, (4, 2)), ''), (COLORS['orange'], '-', '')]
    for (n, (delta, (color, linestyle, marker))) in enumerate(zip(bands['difference'], difference_styles), 2):
        d.plot(bands['cx'], delta, color=color, lw=1.05, ls=linestyle, marker=marker, label=str(n))
    d.axhline(0, color=COLORS['ink'], lw=0.65, ls=':', zorder=0)
    d.set(ylim=(-0.018, 0.145), yticks=[0, 0.05, 0.1])
    d.set_yticklabels(['0', '0.05', '0.10'])
    d.set_ylabel('$\\Delta\\varepsilon_n$ (eV)', labelpad=4)
    d.set_xlabel('Band path', labelpad=3)
    d.legend(loc='upper center', ncol=5, handlelength=1.5, handletextpad=0.4, columnspacing=0.85, borderaxespad=0.25, fontsize=6.7)
    metrics = bands['metrics']
    fig.text(139 / WIDTH_MM, 37.5 / HEIGHT_MM, 'Mean $|\\Delta\\varepsilon_n|$ ' + f"{metrics['mae_eV']:.3f} eV; max {metrics['max_error_eV']:.3f} eV", ha='center', fontsize=7.2)
    return (fig, axes)

def check_visuals(fig, axes, data):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = fig.bbox
    text_boxes = []
    for text in fig.findobj(Text):
        if not text.get_visible() or not text.get_text():
            continue
        extent = text.get_window_extent(renderer)
        assert bounds.contains(extent.x0, extent.y0) and bounds.contains(extent.x1, extent.y1), text.get_text()
        assert colors.to_rgba(text.get_color()) == colors.to_rgba(COLORS['ink']), text.get_text()
        assert text.get_fontweight() == 'bold', text.get_text()
        text_boxes.append((text.get_text(), extent))
    overlaps = [(left, right) for ((left, a), (right, b)) in combinations(text_boxes, 2) if min(a.x1, b.x1) - max(a.x0, b.x0) > 1 and min(a.y1, b.y1) - max(a.y0, b.y0) > 1]
    assert not overlaps, overlaps
    for (key, ax) in axes.items():
        assert all((spine.get_visible() for spine in ax.spines.values()))
        assert not ax.get_title()
        for line in ax.lines:
            if line.get_transform() != ax.transData:
                continue
            (x, y) = line.get_data()
            if len(x) == 0:
                continue
            assert np.isfinite(x).all() and np.isfinite(y).all()
            assert np.min(x) >= ax.get_xlim()[0] and np.max(x) <= ax.get_xlim()[1], key
            assert np.min(y) >= ax.get_ylim()[0] and np.max(y) <= ax.get_ylim()[1], key
    assert len(axes['b'].lines) == 12
    assert sum((len(line.get_xdata()) for line in axes['d'].lines if len(line.get_xdata()) == 31)) == 155
    for (upper, lower) in (('a', 'c'), ('b', 'd')):
        close(axes[upper].get_position().x0, axes[lower].get_position().x0, 'left alignment')
        close(axes[upper].get_position().x1, axes[lower].get_position().x1, 'right alignment')
        assert axes[upper].get_xlim() == axes[lower].get_xlim()
    return {'black_bold_text_artists': len(text_boxes), 'text_bbox_overlaps': overlaps, 'closed_axes': 4, 'tracks_mm': TRACKS_MM, 'no_data_clipping': True, 'fonts': [plt.rcParams['font.sans-serif'][0], plt.rcParams['mathtext.rm']]}

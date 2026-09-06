import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator
from style import COLORS, MM, panel_label, plt, save_figure, use_style
from crystal_motifs import CELLS, draw_motif
from common import ONLINE, load_rows, OUTPUT

def main():
    rows = load_rows()
    use_style()
    colors = {'Si': COLORS['blue'], 'graphene': COLORS['orange'], 'Al': COLORS['green']}
    height = 116
    fig = plt.figure(figsize=(183 * MM, height * MM))
    left = fig.add_axes([0.135, 0.17 * 90 / height, 0.455, 0.625 * 90 / height])
    right = fig.add_axes([0.705, 0.17 * 90 / height, 0.205, 0.625 * 90 / height], sharey=left)
    structures = []
    for (x, key, name, color, az, el) in [(7, 'si_scope', 'Diamond Si', colors['Si'], -32, 22), (67, 'graphene_scope', 'Graphene', colors['graphene'], -8, 52), (126, 'al_scope', 'FCC Al', colors['Al'], -32, 22)]:
        structures.append(draw_motif(fig, (x, 92, 25, 23), key, color, az, el))
        fig.text((x + 28) / 183, 103.5 / height, name, fontsize=8.2, va='center')
    for ax in (left, right):
        ax.set_ylim(9.6, -0.6)
        ax.set_yticks(range(10))
        ax.yaxis.set_minor_locator(NullLocator())
    left.set_yticklabels([f"{('Graphene' if r['system'] == 'graphene' else r['system'])}  {r['iteration']}" for r in rows])
    left.set_xscale('log')
    left.set_xlim(3e-12, 6)
    left.set_xticks([1e-11, 1e-08, 1e-05, 0.01, 1])
    left.xaxis.set_minor_locator(NullLocator())
    left.set_xlabel('Density-matrix distance')
    left.set_ylabel('System / iteration', labelpad=4)
    right.tick_params(axis='y', left=False, labelleft=False)
    right.set_xlim(-0.18, 5.9)
    right.set_xticks([0, 1, 2, 3, 4, 5])
    right.set_xlabel('Uncorrected ratio, $d_W/B_0$', labelpad=6)
    right.axvline(1, color=COLORS['ink'], lw=0.85, ls=(0, (3, 2)), zorder=1)
    (markers, labels) = ([], [])
    for (index, row) in enumerate(rows):
        color = colors[row['system']]
        shape = 'o' if row['proposal'] == ONLINE else 'D'
        (low, high) = (row['window_error'], row['full_error'])
        assert high > low > 0
        left.plot([low, high], [index, index], color=color, lw=0.9, zorder=2)
        for (x, face) in ((low, 'white'), (high, color)):
            markers.append(left.plot(x, index, marker=shape, linestyle='none', ms=4.5, mec=color, mew=1.0, mfc=face, zorder=3)[0])
        ratio = row['window_over_bound']
        markers.append(right.plot(ratio, index, marker=shape, linestyle='none', ms=4.5, mec=color, mew=1.0, mfc=color, zorder=3)[0])
        if ratio < 0.0001:
            (mantissa, exponent) = f'{ratio:.2e}'.split('e')
            value = f'${mantissa}\\times10^{{{int(exponent)}}}$'
        else:
            value = f'{ratio:.3f}'
        labels.append(right.text(1.37, index, value, transform=right.get_yaxis_transform(), ha='right', va='center', fontsize=7.0, clip_on=False))
    panel_label(left, 'a', x=-0.25)
    panel_label(right, 'b', x=-0.22)
    left.text(0.03, 1.025, 'Open: $d_W$', transform=left.transAxes, va='bottom')
    left.text(0.66, 1.025, 'Filled: $d_{\\mathrm{full}}$', transform=left.transAxes, va='bottom')
    right.text(0.52, 1.025, 'Above 1: 6 of 10', transform=right.transAxes, ha='center', va='bottom', fontsize=7.1)
    legend = fig.legend(handles=[Line2D([], [], marker='o', linestyle='none', mec=COLORS['ink'], mfc='white', ms=4.5, label='Low cutoff (online)'), Line2D([], [], marker='D', linestyle='none', mec=COLORS['ink'], mfc='white', ms=4.5, label='Target-converged subspace (oracle)')], loc='upper center', bbox_to_anchor=(0.55, 0.994 * 90 / height), ncol=2, handlelength=1, handletextpad=0.6, columnspacing=2.5)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label_boxes = [label.get_window_extent(renderer) for label in labels]
    assert not any((a.overlaps(b) for (i, a) in enumerate(label_boxes) for b in label_boxes[:i]))
    assert all((box.x1 < fig.bbox.x1 and box.x0 > right.bbox.x1 for box in label_boxes))
    assert not any((legend.get_window_extent(renderer).overlaps(ax.bbox) for ax in (left, right)))
    assert len(markers) == 30
    assert sum((row['window_over_bound'] > 1 for row in rows)) == 6
    for marker in markers:
        xy = marker.axes.transData.transform((marker.get_xdata()[0], marker.get_ydata()[0]))
        assert marker.axes.bbox.contains(*xy)
    target = OUTPUT / 'fig_screen_scope'
    save_figure(fig, target)

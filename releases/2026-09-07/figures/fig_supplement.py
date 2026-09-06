from __future__ import annotations
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatterMathtext, LogLocator, NullFormatter
import numpy as np
from style import COLORS as C, MM, panel_label, save_figure, use_style

def axis_mm(fig, left, bottom, width, height):
    (fw, fh) = fig.get_size_inches() / MM
    return fig.add_axes([left / fw, bottom / fh, width / fw, height / fh])

def log_axis(ax, which):
    axis = getattr(ax, which + 'axis')
    getattr(ax, 'set_' + which + 'scale')('log')
    axis.set_major_locator(LogLocator(base=10, numticks=4))
    axis.set_major_formatter(LogFormatterMathtext())
    axis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis=which, labelsize=7.2)

def record_points(ledger, ax, name, x, y, artist):
    xy = np.column_stack([x, y]).astype(float)
    if hasattr(artist, 'get_offsets'):
        actual = np.asarray(artist.get_offsets())
    else:
        actual = np.column_stack(artist.get_data())
    np.testing.assert_array_equal(actual, xy)
    assert np.isfinite(xy).all()
    ledger.append({'series': name, 'panel': ax.get_label(), 'points': xy.tolist()})

def fig_s1(data):
    fig = plt.figure(figsize=(183 * MM, 84 * MM))
    axes = [axis_mm(fig, 15, 17, 43, 60), axis_mm(fig, 77, 17, 42, 60), axis_mm(fig, 143, 17, 37, 60)]
    (a, b, c) = axes
    ledger = []
    for (ax, letter) in zip(axes, 'abc'):
        ax.set_label(letter)
        panel_label(ax, letter, x=-0.02, y=1.025)
    rungs = data['h2']['rank_one_rungs']
    fractions = [r['fraction'] for r in rungs]
    by_lam = {r['lambda']: r for r in data['coupling']}
    lams = sorted(by_lam, key=abs)
    series = [('H2/Grid4', fractions, [r['rank_one_update_error'] for r in rungs], 'o', C['blue']), ('PW model', data['pilot']['fixed_fractions'], data['pilot']['jacobian_update_errors'], 's', C['green']), ('HGH-Si', [0.2, 0.4], [by_lam[x]['rank_one_update_error'] for x in lams[1:]], 'D', C['orange'])]
    for (name, x, y, marker, color) in series:
        (line,) = a.plot(x, y, marker=marker, color=color, ms=4, zorder=3)
        record_points(ledger, a, name, x, y, line)
    a.text(0.26, 1.65e-05, 'PW model', ha='center', fontsize=7.2)
    a.text(0.26, 2.3e-06, 'H$_2$/Grid4', ha='center', fontsize=7.2)
    a.text(0.3, 7e-07, 'HGH-Si', ha='center', fontsize=7.2)
    log_axis(a, 'y')
    a.set(xlim=(0.065, 0.435), ylim=(2e-07, 3e-05), xticks=[0.1, 0.2, 0.4], xlabel='Coupling $|\\lambda|/|\\lambda_c|$', ylabel='Rank-one update error')
    h2_base = rungs[0]['m_fp'] - rungs[0]['delta_m_fp']
    hx = [0.0] + [f / fractions[-1] for f in fractions]
    hy = [h2_base] + [r['m_fp'] for r in rungs]
    (line,) = b.plot(hx, hy, marker='o', ms=3.6, color=C['blue'])
    record_points(ledger, b, 'H2/Grid4', hx, hy, line)
    mesh_styles = [('s', '-', 3.5), ('D', '--', 5.0), ('^', ':', 3.2)]
    for (mesh, (marker, ls, size)) in zip(data['k2']['meshes'], mesh_styles):
        cr = mesh['coupling_rungs']
        x = [abs(r['lambda']) / abs(cr[-1]['lambda']) for r in cr]
        y = [r['fixed_point_margin'] for r in cr]
        (line,) = b.plot(x, y, marker=marker, ls=ls, ms=size, color=C['green'], mfc='white' if marker == 'D' else C['green'], mew=0.8)
        record_points(ledger, b, mesh['mesh'], x, y, line)
    label_slots = [('k mesh $2^3$', 0.902, 1), ('k mesh $4^3$', 0.84, 2), ('k mesh $1^3$', 0.775, 0)]
    for (label, slot, index) in label_slots:
        value = data['k2']['meshes'][index]['coupling_rungs'][-1]['fixed_point_margin']
        b.annotate(label, xy=(1, value), xytext=(1.1, slot), fontsize=7.2, ha='left', va='center', arrowprops=dict(arrowstyle='-', color=C['ink'], lw=0.65, shrinkA=1, shrinkB=3))
    b.text(1.1, hy[-1], 'H$_2$', fontsize=7.2, va='center')
    b.set(xlim=(-0.06, 1.68), ylim=(0.51, 1.035), xticks=[0, 0.5, 1], yticks=[0.6, 0.8, 1.0], xlabel='Coupling $|\\lambda|/|\\lambda_{\\max}|$', ylabel='Intrinsic margin $m_{\\mathrm{fp}}$')
    b.text(0.035, 0.06, '$\\lambda_{\\max}=0.4\\lambda_c$', transform=b.transAxes, fontsize=7.4, va='bottom')
    mixer_names = ['unit', 'scalar_0.25', 'scalar_0.70', 'baseline_inverse_control']
    markers = ['o', 's', 'D']
    faces = ['white', C['green'], C['blue']]
    edges = [C['reference'], C['green'], C['blue']]
    offsets = [0.15, 0, -0.15]
    for (row, name) in zip([3, 2, 1, 0], mixer_names):
        x = [data['radii'][lam][name] for lam in lams]
        y = [row + dy for dy in offsets]
        c.plot(x, y, color='#B5B5B5', lw=0.8, zorder=1)
        for (value, pos, marker, face, edge, lam) in zip(x, y, markers, faces, edges, lams):
            point = c.scatter([value], [pos], s=17, marker=marker, facecolor=face, edgecolor=edge, lw=0.85, zorder=3)
            record_points(ledger, c, f'{name}; lambda={lam}', [value], [pos], point)
    c.axvline(1, ls='--', color=C['boundary'], lw=0.9, zorder=0)
    c.text(0.98, -0.5, '$\\rho=1$', fontsize=7.2, ha='right', va='bottom')
    handles = [Line2D([], [], marker=m, mfc=f, mec=e, ms=4, ls='none') for (m, f, e) in zip(markers, faces, edges)]
    c.legend(handles, ['$0$', '$0.2\\lambda_c$', '$0.4\\lambda_c$'], loc='upper left', bbox_to_anchor=(0.015, 0.99), fontsize=7.2, borderaxespad=0, labelspacing=0.25, handletextpad=0.4)
    c.set(xlim=(-0.055, 1.15), ylim=(-0.62, 4.7), xticks=[0, 0.5, 1], yticks=[3, 2, 1, 0], yticklabels=['Unit', 'Scalar 0.25', 'Scalar 0.70', 'Exact inverse'], xlabel='Spectral radius $\\rho(M_B)$')
    c.tick_params(axis='y', labelsize=7)
    return (fig, ledger)

def fig_s2(data):
    fig = plt.figure(figsize=(183 * MM, 98 * MM))
    a = axis_mm(fig, 11, 18, 51, 72)
    b = axis_mm(fig, 80, 18, 43, 72)
    c = axis_mm(fig, 141, 18, 36, 72)
    ledger = []
    for (ax, letter) in zip((a, b, c), 'abc'):
        ax.set_label(letter)
        panel_label(ax, letter, x=-0.02, y=1.025)
    v = data['verdict']
    gate = v['g2_gate_evaluation']['criteria']
    est = gate['margin']['improvement'] * 100
    target = gate['margin']['threshold'] * 100
    (lo, hi) = (v['bootstrap'][k] * 100 for k in ('ci_lower', 'ci_upper'))
    a.axvline(0, ymin=0.48, ymax=0.72, color=C['reference'], lw=0.8)
    a.axvline(target, ymin=0.48, ymax=0.72, color=C['boundary'], ls='--', lw=1.0)
    (interval,) = a.plot([lo, hi], [0.61, 0.61], color=C['reference'], lw=1.7)
    a.plot([lo, hi], [0.61, 0.61], ls='none', marker='|', ms=9, mew=1.0, color=C['reference'])
    (point,) = a.plot([est], [0.61], marker='o', color=C['blue'], ms=5.0, ls='none', zorder=3)
    record_points(ledger, a, 'Corrected estimate (%)', [est], [0.61], point)
    record_points(ledger, a, 'Stored 95% interval endpoints (%)', [lo, hi], [0.61, 0.61], interval)
    a.text(0.06, 0.94, f'Corrected {est:+.2f}%', transform=a.transAxes, va='top', fontsize=8.3)
    a.text(0.06, 0.84, f'{target:.0f}% target: not met', transform=a.transAxes, va='top', fontsize=7.6)
    a.text(0.5, 0.45, '95% hierarchical interval\n' + f'[{lo:.1f}%, +{hi:.1f}%]'.replace('-', '−'), transform=a.transAxes, va='top', ha='center', fontsize=7.2, linespacing=1.4)
    n = v['analysis']
    a.text(0.06, 0.265, f"{n['rows_admitted_paired']:,} paired SCF steps\n{n['trajectories']} trajectories\n{n['trajectories_effective']} unique diagnostic series\n{len(n['parents'])} parent lineages", transform=a.transAxes, va='top', fontsize=7.2, linespacing=1.45)
    a.set(xlim=(-190, 60), ylim=(0, 1), xticks=[-150, -100, -50, 0, 50], yticks=[], xlabel='Recorded-counter\nerror reduction (%)')
    per_cat = gate['category_reversal']['per_category_improvement']
    b.axvline(0, ymin=0.075, color=C['reference'], lw=0.7, zorder=0)
    b.axvline(-5, ymin=0.075, color=C['boundary'], lw=0.9, ls='--', zorder=0)
    for (pos, name) in zip(range(7, -1, -1), [f'C{i}' for i in range(1, 9)]):
        if name not in per_cat:
            assert name == 'C6'
            b.text(0.08, pos, 'No eligible samples\nfor fitted comparison', transform=b.get_yaxis_transform(), va='center', fontsize=7.0, linespacing=1.15, bbox=dict(facecolor='white', edgecolor='none', pad=1.2))
            continue
        value = per_cat[name] * 100
        reversal = value < gate['category_reversal']['threshold'] * 100
        color = C['orange'] if reversal else C['blue']
        b.plot([0, value], [pos, pos], color=color, lw=1)
        point = b.scatter([value], [pos], s=20, marker='s' if reversal else 'o', color=color, zorder=3)
        record_points(ledger, b, name, [value], [pos], point)
        b.text(0.97, pos, f'{value:+.1f}'.replace('-', '−'), transform=b.get_yaxis_transform(), ha='right', va='center', fontsize=7.0)
    b.text(0.045, 0.03, 'Reversal < −5%', transform=b.transAxes, fontsize=7.0)
    b.set(xlim=(-25, 74), ylim=(-0.8, 7.65), xticks=[-20, 0, 20, 40, 60], yticks=list(range(7, -1, -1)), yticklabels=[f'C{i}' for i in range(1, 9)], xlabel='Category\nerror reduction (%)')
    (ak, dp) = (data['ak'], data['dp'])
    raw = c.scatter(ak, dp, s=9, facecolor='white', edgecolor=C['blue'], lw=0.65, zorder=3)
    record_points(ledger, c, '60 positive pairs, one trajectory', ak, dp, raw)
    xs = np.geomspace(ak.min(), ak.max(), 128)
    c.plot(xs, np.exp(data['intercept']) * xs ** data['slope'], color=C['reference'], lw=1, zorder=2)
    for which in ('x', 'y'):
        log_axis(c, which)
    c.set_xlim(ak.min() / 2.5, ak.max() * 2.5)
    c.set_ylim(dp.min() / 2.5, dp.max() * 12)
    c.set_xlabel('Adiabaticity $A_k$')
    c.set_ylabel('Projector distance $d_P$')
    c.text(0.045, 0.965, f"Slope {data['slope']:.3f}\n60 pairs; 1 trajectory", transform=c.transAxes, fontsize=7.2, va='top', linespacing=1.45)
    return (fig, ledger)

def figure_checks(fig, ledger):
    fig.canvas.draw()
    for text in fig.findobj(mpl.text.Text):
        floor = 8.0 if '$' in text.get_text() else 7.4
        text.set_fontsize(max(floor, text.get_fontsize()))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = []
    for ax in fig.axes:
        assert not ax.get_title()
        assert all((s.get_visible() for s in ax.spines.values()))
        for row in [s for s in ledger if s['panel'] == ax.get_label()]:
            pts = np.array(row['points'])
            (lo_x, hi_x) = ax.get_xlim()
            (lo_y, hi_y) = ax.get_ylim()
            assert np.all((pts[:, 0] > lo_x) & (pts[:, 0] < hi_x))
            assert np.all((pts[:, 1] > lo_y) & (pts[:, 1] < hi_y))
        texts.extend(ax.texts)
        texts.extend([ax.xaxis.label, ax.yaxis.label])
        for (axis, limits) in [(ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())]:
            texts.extend((t.label1 for t in axis.get_major_ticks() if limits[0] <= t.get_loc() <= limits[1]))
        if ax.get_legend():
            texts.extend(ax.get_legend().get_texts())
    visible = [t for t in texts if t.get_visible() and t.get_text()]
    overlaps = []
    for (i, text) in enumerate(visible):
        bbox = mpl.text.Text.get_window_extent(text, renderer)
        assert fig.bbox.contains(bbox.x0, bbox.y0) and fig.bbox.contains(bbox.x1, bbox.y1), text.get_text()
        assert text.get_fontweight() in ('bold', 700), text.get_text()
        assert mpl.colors.to_hex(text.get_color()) == C['ink'], text.get_text()
        for other in visible[i + 1:]:
            if bbox.overlaps(mpl.text.Text.get_window_extent(other, renderer)):
                overlaps.append([text.get_text(), other.get_text()])
    assert not overlaps, overlaps
    positions = [a.get_position().bounds for a in fig.axes]
    np.testing.assert_allclose([p[1] for p in positions], positions[0][1], atol=1e-12)
    np.testing.assert_allclose([p[3] for p in positions], positions[0][3], atol=1e-12)
    return {'size_mm': (fig.get_size_inches() / MM).tolist(), 'axes': [{'panel': a.get_label(), 'xlim': list(a.get_xlim()), 'ylim': list(a.get_ylim()), 'xscale': a.get_xscale(), 'yscale': a.get_yscale(), 'closed_spines': True} for a in fig.axes], 'text_count': len(visible), 'text_overlaps': overlaps, 'minimum_text_artist_pt': min((t.get_fontsize() for t in visible)), 'all_points_within_axes': True, 'raw_series': ledger}

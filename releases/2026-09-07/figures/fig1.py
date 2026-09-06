from __future__ import annotations
import hashlib
import numpy as np
from style import COLORS as C, MM, WIDTH_MM, save_figure, use_style
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.ticker import NullLocator
from common import load_response, require, ecdf, OUTPUT
HEIGHT_MM = 156

def array_hash(values):
    return hashlib.sha256(np.asarray(values, dtype='<f8').tobytes()).hexdigest()

def schematic(fig):
    ax = fig.add_axes([0, 0, 1, 1], label='schematic')
    ax.set(xlim=(0, WIDTH_MM), ylim=(0, HEIGHT_MM))
    ax.set_axis_off()
    ax.set_zorder(-1)
    (boxes, texts, routes) = ({}, [], [])

    def box(name, x, y, w, h, color=C['reference'], dashed=False):
        patch = Rectangle((x, y), w, h, facecolor='white', edgecolor=color, linewidth=0.85, linestyle=(0, (3, 2)) if dashed else '-', zorder=2)
        ax.add_patch(patch)
        boxes[name] = patch

    def txt(x, y, label, size=7.3, owner=None, ha='center', **kwargs):
        obj = ax.text(x, y, label, fontsize=size, ha=ha, va='center', color=C['ink'], fontweight='bold', zorder=5, **kwargs)
        texts.append((owner, obj))
        return obj

    def route(name, points, color=C['ink']):
        for (p, q) in zip(points, points[1:]):
            require(p[0] == q[0] or p[1] == q[1], name + ' orthogonal connector')
        if len(points) > 2:
            ax.plot(*np.array(points[:-1]).T, color=color, lw=0.9, solid_capstyle='butt', solid_joinstyle='miter', zorder=1)
        ax.annotate('', xy=points[-1], xytext=points[-2], zorder=1, arrowprops=dict(arrowstyle='-|>', color=color, lw=0.9, shrinkA=0, shrinkB=0, mutation_scale=7))
        routes.append({'name': name, 'points_mm': points})
    txt(2.0, 152.5, 'a', 10, ha='left')
    box('state', 8, 139, 39, 12)
    txt(27.5, 147.5, 'Current state / snapshot', owner='state')
    txt(27.5, 142.5, '$S_t=(H_t,D_t,M_t)$', 8.4, owner='state')
    box('sources', 8, 106, 39, 29, C['blue'])
    txt(27.5, 131.8, 'Low cutoff: current potential', 6.9, owner='sources')
    ax.plot([10, 45], [128.4, 128.4], color='#D0D0D0', lw=0.5, zorder=3)
    txt(27.5, 124.8, 'Target-converged Ritz', 7.2, owner='sources')
    txt(27.5, 120.5, 'oracle/control', 7.2, owner='sources')
    ax.plot([10, 45], [117.2, 117.2], color='#D0D0D0', lw=0.5, zorder=3)
    txt(27.5, 113.8, 'Coarse k', 7.2, owner='sources')
    txt(27.5, 109.5, 'Auxiliary + fine frame', 6.9, owner='sources')
    box('candidate', 8, 92, 39, 10, C['blue'])
    txt(27.5, 99, 'Candidate window blocks', 7.1, owner='candidate')
    txt(27.5, 94.8, '$P_k,\\ D_{c,k},\\ N_W$', 8.4, owner='candidate')
    route('state_to_sources', [(27.5, 139), (27.5, 135)])
    route('sources_to_candidate', [(27.5, 106), (27.5, 102)], C['blue'])
    box('residual', 57, 122, 71, 29, C['blue'])
    txt(92.5, 147.6, 'Fixed-H window residual', 7.7, owner='residual')
    txt(92.5, 142.0, '$H_{P,k}=P_k^{\\dagger}H_{t,k}P_k$', 8.7, owner='residual')
    txt(92.5, 135.1, '$g_k=H_{P,k}+\\tau\\ln[D_{c,k}(I_k-D_{c,k})^{-1}]$', 8.0, owner='residual')
    txt(92.5, 127.6, '$(\\Pi_0g)_k=g_k-cI_k,\\quad r=\\|\\Pi_0g\\|_w$', 8.7, owner='residual')
    route('snapshot_H_to_compression', [(47, 144.5), (57, 144.5)])
    txt(52, 147.2, '$H_t$', 7.8)
    route('candidate_to_residual', [(47, 97), (52, 97), (52, 135.1), (57, 135.1)], C['blue'])
    box('theorem', 57, 92, 71, 26, C['reference'], dashed=True)
    txt(92.5, 114.5, 'Exact-trace theorem', 7.5, owner='theorem')
    txt(92.5, 108.4, '$\\sum_k w_k\\,\\mathrm{Tr}\\,D_{c,k}=N_W$', 9.0, owner='theorem')
    txt(92.5, 99.0, '$\\Delta_D\\leq\\frac{r}{4\\tau},\\qquad\\Delta_{\\mathcal{F}}\\leq\\frac{r^2}{8\\tau}$', 10.0, owner='theorem')
    box('gate', 138, 125, 41, 26, C['orange'])
    txt(158.5, 147.6, 'Finite numerical gate', 7.5, owner='gate')
    txt(158.5, 141.1, '$|\\sum_k w_k\\,\\mathrm{Tr}\\,D_{c,k}-N_W|\\leq\\varepsilon_N$', 7.4, owner='gate')
    txt(158.5, 134.9, '$b_D=r/(4\\tau)\\leq\\varepsilon_D$', 8.5, owner='gate')
    txt(158.5, 128.6, 'Interiority + density checks', 6.9, owner='gate')
    route('residual_to_gate', [(128, 135.1), (138, 135.1)])
    txt(133, 138.0, '$r,D_c,N_W$', 7.2)
    box('accept', 138, 111, 41, 10, C['green'])
    txt(158.5, 118.0, 'Pass: accept candidate', 7.2, owner='accept')
    txt(158.5, 113.6, 'Commit $D_c$', 7.3, owner='accept')
    box('recovery', 138, 96, 41, 11)
    txt(158.5, 103.8, 'Reject/abstain: reference', 7.1, owner='recovery')
    txt(158.5, 99.0, 'Restore $S_t$; evaluate map', 7.1, owner='recovery')
    box('mixer', 138, 83, 41, 8)
    txt(158.5, 87, 'Mixer $\\to S_{t+1}$', 8.0, owner='mixer')
    route('gate_pass', [(152, 125), (152, 121)], C['green'])
    route('gate_reject_or_abstain', [(179, 138.0), (181, 138.0), (181, 101.5), (179, 101.5)])
    route('accepted_update_to_mixer', [(138, 116), (134, 116), (134, 87), (138, 87)], C['green'])
    route('reference_update_to_mixer', [(158.5, 96), (158.5, 91)])
    txt(8, 85.8, 'Independent response interventions: Supplemental Fig. S1', 7.0, ha='left')
    return (ax, boxes, texts, routes)

def quantitative_axes(fig, x, y, w, h, letter):
    ax = fig.add_axes([x / WIDTH_MM, y / HEIGHT_MM, w / WIDTH_MM, h / HEIGHT_MM], label=letter)
    for spine in ax.spines.values():
        spine.set_visible(True)
    ax.text(-0.15, 1.015, letter, transform=ax.transAxes, fontsize=10, va='bottom', ha='left', fontweight='bold')
    return ax

def draw_ecdf(ax, series, lines, metadata):
    for (name, values, color, linestyle) in series:
        (xs, ys) = ecdf(values)
        (drawn,) = ax.step(xs, ys, where='post', color=color, lw=1.2, ls=linestyle, label=name, zorder=3)
        require(np.array_equal(drawn.get_xdata(), np.sort(values)), name + ' all ECDF values plotted')
        require(np.array_equal(drawn.get_ydata(), np.arange(1, len(values) + 1) / len(values)), name + ' exact ECDF normalization')
        metadata[name] = {'point_count': len(values), 'sorted_x_sha256': array_hash(xs), 'y_sha256': array_hash(ys), 'where': 'post'}
    ax.axvline(1, color=C['boundary'], lw=0.9, ls=(0, (3, 2)), zorder=1)
    ax.set(xlim=(0, 1.06), ylim=(0, 1.04), xticks=[0, 0.25, 0.5, 0.75, 1], yticks=[0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(['0', '0.25', '0.5', '0.75', '1'])
    ax.set_yticklabels(['0', '0.25', '0.5', '0.75', '1'])
    ax.set_xlabel('Distance / raw bound', labelpad=2.0, fontsize=7.5)
    ax.set_ylabel('ECDF', labelpad=3, fontsize=7.5)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    for (i, (label, color, style)) in enumerate(lines):
        yy = 0.93 - 0.135 * i
        if color:
            ax.plot([0.045, 0.115], [yy, yy], transform=ax.transAxes, color=color, ls=style, lw=1.3, zorder=4)
            xx = 0.14
        else:
            xx = 0.045
        ax.text(xx, yy, label, transform=ax.transAxes, fontsize=7.0, va='center', color=C['ink'], zorder=4)

def draw_shells(ax, retained, summary, tail, metadata):
    (keys, shares, _) = np.array(retained).T
    (cumulative_keys, _, cumulative) = np.array(tail).T
    stems = ax.vlines(keys, 0, shares, color=C['reference'], lw=1.4, zorder=2)
    ax.plot(keys, shares, ls='none', marker='s', ms=2.4, color=C['reference'], zorder=3)
    (curve,) = ax.step(cumulative_keys, cumulative, where='post', color=C['blue'], lw=1.3, zorder=3)
    ax.axhline(0.9, color='#AAAAAA', ls=(0, (2, 2)), lw=0.6, zorder=1)
    ax.axhline(0.99, color='#AAAAAA', ls=(0, (2, 2)), lw=0.6, zorder=1)
    ax.axvline(27, ymin=0.86, color=C['orange'], lw=0.9, ls=(0, (3, 2)), zorder=1)
    ax.plot([11, 27], [0.904053, 0.992685], ls='none', marker='o', ms=3, mfc='white', mec=C['blue'], mew=0.9, zorder=4)
    ax.set(xlim=(0, 41), ylim=(-0.04, 1.04), xticks=[3, 11, 19, 27, 40], yticks=[0, 0.5, 1])
    ax.set_yticklabels(['0', '0.5', '1'])
    ax.set_xlabel('$|G|^2$ shell (units $(2\\pi/a)^2$)', fontsize=7.5, labelpad=2.0)
    ax.set_ylabel('Derivative-energy share', fontsize=7.5, labelpad=3)
    ax.text(0.035, 0.78, 'HGH-Si', transform=ax.transAxes, fontsize=7.0, va='center')
    ax.plot([0.27, 0.34], [0.56, 0.56], transform=ax.transAxes, color=C['blue'], lw=1.3)
    ax.text(0.37, 0.56, 'Cumulative', transform=ax.transAxes, fontsize=7.0, va='center')
    ax.plot([0.27, 0.34], [0.43, 0.43], transform=ax.transAxes, color=C['reference'], lw=1.4, marker='s', markevery=[0], ms=2.4)
    ax.text(0.37, 0.43, 'Shell weight', transform=ax.transAxes, fontsize=7.0, va='center')
    ax.text(0.97, 0.26, '168 columns; S = 27', transform=ax.transAxes, fontsize=7.0, ha='right', va='center')
    omitted = 1.0 - sum((s[1] for s in retained))
    ax.text(0.97, 0.13, f'Omitted at S: {100 * omitted:.3f}%', transform=ax.transAxes, fontsize=7.0, ha='right', va='center')
    ax.annotate('90.4%', xy=(11, 0.904053), xytext=(8.0, 0.69), fontsize=7.0, ha='center', arrowprops=dict(arrowstyle='-', lw=0.6, color=C['ink']))
    ax.annotate('99.3%', xy=(27, 0.992685), xytext=(34, 0.76), fontsize=7.0, ha='center', arrowprops=dict(arrowstyle='-', lw=0.6, color=C['ink']))
    require(np.array_equal(curve.get_xdata(), cumulative_keys), 'all cumulative shell x values')
    require(np.array_equal(curve.get_ydata(), cumulative), 'all cumulative shell y values')
    require(np.array_equal([seg[-1, 1] for seg in stems.get_segments()], shares), 'all shell shares')
    metadata['shells'] = {'retained_stems': len(shares), 'cumulative_points': len(cumulative), 'omitted_share_label': omitted, 'shares_sha256': array_hash(shares), 'cumulative_sha256': array_hash(cumulative)}

def check_layout(fig, schematic_ax, boxes, texts, routes, quantitative):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    text_boxes = []
    for (owner, obj) in texts:
        bb = obj.get_window_extent(renderer)
        if owner:
            pb = boxes[owner].get_window_extent(renderer)
            require(pb.x0 + 1 <= bb.x0 and bb.x1 <= pb.x1 - 1 and (pb.y0 + 1 <= bb.y0) and (bb.y1 <= pb.y1 - 1), f'text outside {owner}: {obj.get_text()}')
        text_boxes.append((obj.get_text(), bb))
    collisions = []
    for (i, (label, box)) in enumerate(text_boxes):
        for (other_label, other) in text_boxes[i + 1:]:
            if box.overlaps(other):
                collisions.append((label, other_label))
    require(not collisions, f'schematic text overlaps: {collisions}')
    all_visible = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
    outside = []
    for obj in all_visible:
        bb = obj.get_window_extent(renderer)
        if bb.x0 < -0.5 or bb.y0 < -0.5 or bb.x1 > fig.bbox.x1 + 0.5 or (bb.y1 > fig.bbox.y1 + 0.5):
            outside.append(obj.get_text())
    require(not outside, 'text outside canvas: ' + repr(outside))
    require(all((all((s.get_visible() for s in ax.spines.values())) for ax in quantitative)), 'four closed quantitative panels')

    def at_box_edge(point):
        (x, y) = point
        return any(((np.isclose(x, patch.get_x()) or np.isclose(x, patch.get_x() + patch.get_width())) and patch.get_y() <= y <= patch.get_y() + patch.get_height() or ((np.isclose(y, patch.get_y()) or np.isclose(y, patch.get_y() + patch.get_height())) and patch.get_x() <= x <= patch.get_x() + patch.get_width()) for patch in boxes.values()))
    require(all((at_box_edge(route['points_mm'][0]) and at_box_edge(route['points_mm'][-1]) for route in routes)), 'all connector endpoints terminate at module edges')
    return {'schematic_text_collisions': collisions, 'texts_outside_canvas': outside, 'schematic_texts_inside_modules': True, 'all_connector_endpoints_on_box_edges': True, 'four_closed_quantitative_panels': True, 'orthogonal_connectors': routes, 'all_text_artist_count': len(all_visible), 'non_math_font_family': plt.rcParams['font.sans-serif'], 'non_math_font_weight': plt.rcParams['font.weight'], 'math_rm': plt.rcParams['mathtext.rm'], 'math_it': plt.rcParams['mathtext.it']}

def main():
    (arrays, retained, shell_summary, tail, checks) = load_response()
    use_style()
    plt.rcParams.update({'mathtext.fallback': 'stix', 'path.simplify': False, 'svg.hashsalt': 'idmate-main-fig1-ncs-2026-09-05'})
    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM), dpi=200)
    (ax_a, boxes, texts, routes) = schematic(fig)
    ax_b = quantitative_axes(fig, 15, 48, 72, 29, 'b')
    ax_c = quantitative_axes(fig, 107, 48, 72, 29, 'c')
    ax_d = quantitative_axes(fig, 15, 10, 72, 29, 'd')
    ax_e = quantitative_axes(fig, 107, 10, 72, 29, 'e')
    plotted = {}
    draw_ecdf(ax_b, [('e0', arrays['e0'], C['blue'], '-')], [('Seeded matrices: n = 10,000', None, '-'), ('9,886 defined; 114 undefined', None, '-')], plotted)
    draw_ecdf(ax_c, [('e3_graphene', arrays['e3_graphene'], C['blue'], '-'), ('e3_al', arrays['e3_al'], C['green'], (0, (4, 2)))], [('Graphene: n = 164', C['blue'], '-'), ('fcc Al: n = 1,614', C['green'], (0, (4, 2)))], plotted)
    maximum = float(arrays['e3_graphene'].max())
    ax_c.plot(maximum, 1, marker='D', ms=3.5, color=C['orange'], zorder=6)
    ax_c.annotate('Dirac contact\n1.000 (rounded)', xy=(maximum, 1), xytext=(0.32, 0.63), fontsize=7.0, color=C['ink'], va='center', arrowprops=dict(arrowstyle='-', color=C['ink'], lw=0.6))
    draw_ecdf(ax_d, [('gk_graphene', arrays['gk_graphene'], C['blue'], '-'), ('gk_al', arrays['gk_al'], C['green'], (0, (4, 2)))], [('Shared $\\mu$', None, '-'), ('Graphene: 308/310 defined', C['blue'], '-'), ('fcc Al: 400/400 defined', C['green'], (0, (4, 2)))], plotted)
    draw_shells(ax_e, retained, shell_summary, tail, plotted)
    for ax in (ax_b, ax_c, ax_d, ax_e):
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')
            label.set_color(C['ink'])
    checks['layout'] = check_layout(fig, ax_a, boxes, texts, routes, (ax_b, ax_c, ax_d, ax_e))
    checks['plotted_series'] = plotted
    save_figure(fig, OUTPUT / 'fig1_screen')
    plt.close(fig)

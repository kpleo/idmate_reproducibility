#!/usr/bin/env python3
"""Generate manuscript figures from author-supplied derived records."""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

from common import (
    ReproductionError,
    load_manifest,
    parse_value,
    read_csv_rows,
    read_json,
    repository_root,
    safe_path,
    validate_csv,
    validate_worklaw_json,
)


COLORS = {
    "blue": "#0072B2",
    "green": "#009E73",
    "orange": "#E69F00",
    "red": "#D55E00",
    "purple": "#9B5AA0",
    "gray": "#4D4D4D",
    "light_blue": "#DEEBF7",
    "light_orange": "#FCF4D7",
    "light_gray": "#F4F4F4",
    "ink": "#202020",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path, default=None,
        help="manifest JSON; defaults to repository manifest.json",
    )
    parser.add_argument(
        "--data-root", type=Path, default=None,
        help="directory containing the derived numerical inputs",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="destination for PDF and SVG files",
    )
    parser.add_argument(
        "--only", action="append",
        help="one figure identifier; repeat to select more than one",
    )
    parser.add_argument("--list", action="store_true", help="list figure identifiers and exit")
    parser.add_argument(
        "--schematic-only", action="store_true",
        help="render only the data-independent Fig. 1a method schematic",
    )
    return parser.parse_args()


def plotting_modules() -> tuple[Any, Any, Any]:
    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyBboxPatch
    except ImportError as exc:
        raise ReproductionError(
            "Matplotlib is required; install requirements.txt before rendering"
        ) from exc
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "mathtext.fontset": "stix",
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "font.size": 7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "axes.grid": False,
            "legend.frameon": False,
            "figure.dpi": 180,
        }
    )
    return mpl, plt, FancyBboxPatch


def save_figure(fig: Any, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg"):
        destination = output_dir / f"{stem}.{suffix}"
        fig.savefig(destination, bbox_inches="tight")
        print(f"wrote {destination}")


def close_axes(ax: Any) -> None:
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.85)
        spine.set_color(COLORS["ink"])


def bold_axis_text(ax: Any) -> None:
    ax.xaxis.label.set_fontweight("bold")
    ax.yaxis.label.set_fontweight("bold")
    for label in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        label.set_fontweight("bold")


def numeric_rows(path: Path, declaration: dict[str, Any]) -> list[dict[str, Any]]:
    validate_csv(path, declaration["columns"])
    columns = {column["name"]: column for column in declaration["columns"]}
    parsed: list[dict[str, Any]] = []
    for row_number, row in enumerate(read_csv_rows(path), start=2):
        record: dict[str, Any] = {}
        for name, column in columns.items():
            record[name] = parse_value(
                row.get(name, ""),
                column["type"],
                nullable=bool(column.get("nullable", False)),
                context=f"{path}:{row_number}:{name}",
            )
        parsed.append(record)
    return parsed


def load_figure_inputs(
    figure: dict[str, Any], data_root: Path
) -> dict[str, Any]:
    loaded: dict[str, Any] = {}
    for declaration in figure["inputs"]:
        relative = declaration["path"]
        path = safe_path(data_root, relative)
        if not path.is_file():
            raise ReproductionError(
                f"required derived input is not distributed: {relative}; "
                "supply it under --data-root"
            )
        if declaration["format"] == "csv":
            loaded[relative] = numeric_rows(path, declaration)
        elif declaration["format"] == "json" and figure["id"] == "fig3_worklaw":
            validate_worklaw_json(path)
            loaded[relative] = read_json(path)
        else:
            raise ReproductionError(f"no figure reader is defined for {relative}")
    return loaded


def ecdf(values: list[float]) -> tuple[list[float], list[float]]:
    ordered = sorted(values)
    if not ordered:
        return [], []
    n = len(ordered)
    return ordered, [(index + 1) / n for index in range(n)]


def method_schematic(ax: Any, FancyBboxPatch: Any) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.06)
    ax.axis("off")

    def box(
        x: float,
        y: float,
        w: float,
        h: float,
        label: str,
        fill: str,
        *,
        edge: str | None = None,
        fontsize: float = 6.2,
    ) -> None:
        patch = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="square,pad=0",
            facecolor=fill, edgecolor=edge or COLORS["ink"], linewidth=1.0,
            zorder=3,
        )
        ax.add_patch(patch)
        ax.text(
            x + w / 2,
            y + h / 2,
            label,
            ha="center",
            va="center",
            weight="bold",
            fontsize=fontsize,
            linespacing=1.15,
            zorder=4,
        )

    def container(x: float, y: float, w: float, h: float, label: str, fill: str) -> None:
        patch = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="square,pad=0",
            facecolor=fill, edgecolor=COLORS["ink"], linewidth=1.1,
            zorder=0,
        )
        ax.add_patch(patch)
        ax.text(
            x + 0.012,
            y + h - 0.027,
            label,
            ha="left",
            va="center",
            weight="bold",
            fontsize=6.6,
            zorder=4,
        )

    def route(
        points: list[tuple[float, float]],
        *,
        color: str | None = None,
        dashed: bool = False,
    ) -> None:
        line_color = color or COLORS["ink"]
        if len(points) > 2:
            ax.plot(
                [point[0] for point in points[:-1]],
                [point[1] for point in points[:-1]],
                color=line_color,
                lw=1.0,
                ls=(0, (4, 2)) if dashed else "-",
                solid_joinstyle="miter",
                solid_capstyle="butt",
                zorder=2,
            )
        ax.annotate(
            "",
            xy=points[-1],
            xytext=points[-2],
            arrowprops={
                "arrowstyle": "-|>",
                "color": line_color,
                "lw": 1.0,
                "linestyle": (0, (4, 2)) if dashed else "-",
                "shrinkA": 0,
                "shrinkB": 0,
                "mutation_scale": 8,
            },
            zorder=2,
        )

    ax.text(
        0.5,
        1.045,
        "Online SCF control with separate intervention analysis",
        ha="center",
        va="top",
        weight="bold",
        fontsize=7.4,
    )

    container(0.01, 0.39, 0.24, 0.55, "Proposal Sources", "white")
    container(0.26, 0.39, 0.45, 0.55, "Response Certificate", COLORS["light_blue"])
    container(0.72, 0.39, 0.27, 0.55, "State Controller", "white")

    box(0.045, 0.79, 0.17, 0.09, "Current State\n$\mathbf{H_t,\ D_t,\ M_t}$", "white")
    for y, label in ((0.68, "Low Cutoff"), (0.60, "Subspace Reuse"), (0.52, "Coarse k")):
        box(0.050, y, 0.16, 0.055, label, COLORS["light_orange"], edge=COLORS["orange"])
    box(0.026, 0.415, 0.10, 0.072, "Snapshot\n$\mathbf{S_t}$", "white")
    box(0.143, 0.415, 0.088, 0.072, "Candidate\n$\mathbf{D_c}$", COLORS["light_blue"], edge=COLORS["blue"])

    ax.plot([0.045, 0.03, 0.03], [0.835, 0.835, 0.548], color=COLORS["ink"], lw=0.95, zorder=2)
    for y in (0.7075, 0.6275, 0.5475):
        route([(0.03, y), (0.05, y)])
        route([(0.21, y), (0.225, y)])
    ax.plot([0.225, 0.225], [0.5475, 0.7075], color=COLORS["ink"], lw=0.95, zorder=2)
    route([(0.13, 0.79), (0.13, 0.487)])
    route([(0.225, 0.52), (0.225, 0.50), (0.187, 0.487)], color=COLORS["blue"])

    box(
        0.285,
        0.43,
        0.40,
        0.115,
        "Mermin Gradient\n"
        + r"$\mathbf{g=H_t+\tau\log[D_c(I-D_c)^{-1}]}$",
        "white",
        edge=COLORS["blue"],
        fontsize=5.8,
    )
    box(
        0.285,
        0.575,
        0.40,
        0.095,
        "Trace Projection\n$\mathbf{r=\Pi_0 g}$",
        "white",
        edge=COLORS["blue"],
        fontsize=6.0,
    )
    box(
        0.285,
        0.70,
        0.40,
        0.17,
        "Certified Gates\n"
        + r"$\mathbf{|\sum_k Tr D_{c,k}-N_e|\leq\epsilon_N}$"
        + "\n"
        + r"$\mathbf{\Delta_D\leq r/(4\tau),\quad \Delta_\Omega\leq r^2/(8\tau)}$",
        "white",
        edge=COLORS["orange"],
        fontsize=5.7,
    )
    route([(0.231, 0.451), (0.285, 0.451)], color=COLORS["blue"])
    route([(0.485, 0.545), (0.485, 0.575)])
    route([(0.485, 0.670), (0.485, 0.700)])

    box(0.755, 0.77, 0.14, 0.10, "Admissible?\nAll Pass", COLORS["light_orange"], edge=COLORS["orange"])
    box(0.915, 0.79, 0.058, 0.065, "Error\nStop", "white", edge=COLORS["red"], fontsize=5.6)
    box(0.745, 0.61, 0.105, 0.10, "Accept\n$\mathbf{Commit\ D_c}$", "white", edge=COLORS["green"])
    box(0.865, 0.61, 0.105, 0.10, "Reject / Abstain\nReference Map", "white", fontsize=5.5)
    box(0.79, 0.43, 0.14, 0.10, "Next State\n$\mathbf{S_{t+1}}$", COLORS["light_blue"], edge=COLORS["blue"])
    route([(0.685, 0.79), (0.755, 0.79)])
    route([(0.895, 0.82), (0.915, 0.82)], color=COLORS["red"])
    route([(0.80, 0.77), (0.80, 0.71)], color=COLORS["green"])
    route([(0.85, 0.77), (0.92, 0.71)])
    route([(0.7975, 0.61), (0.7975, 0.56), (0.86, 0.56), (0.86, 0.53)], color=COLORS["green"])
    route([(0.9175, 0.61), (0.9175, 0.56), (0.86, 0.56), (0.86, 0.53)])
    route([(0.93, 0.48), (0.998, 0.48), (0.998, 0.965), (0.13, 0.965), (0.13, 0.88)])
    ax.text(0.56, 0.97, "Next Iteration", ha="center", va="bottom", weight="bold", fontsize=6.2)

    container(0.01, 0.015, 0.98, 0.29, "Separate Intervention Analysis", "white")
    ax.text(0.035, 0.23, "Causal Response Decomposition", ha="left", va="center", weight="bold", fontsize=6.1)
    ax.text(0.55, 0.23, "Rank-One Intervention", ha="left", va="center", weight="bold", fontsize=6.1)
    box(0.045, 0.075, 0.085, 0.09, "Jacobian\n$\mathbf{J}$", COLORS["light_blue"], edge=COLORS["blue"])
    box(0.155, 0.075, 0.095, 0.09, "Response\n$\mathbf{A}$", COLORS["light_blue"], edge=COLORS["blue"])
    box(0.29, 0.075, 0.19, 0.09, "SCF Map\n$\mathbf{M_B=I-B(I-J)}$", "white", fontsize=5.8)
    route([(0.13, 0.12), (0.155, 0.12)])
    route([(0.25, 0.12), (0.29, 0.12)])
    ax.text(0.385, 0.185, "$\mathbf{Mixer\ B}$", ha="center", va="center", weight="bold", fontsize=6.0)
    route([(0.385, 0.175), (0.385, 0.165)])
    box(0.55, 0.075, 0.17, 0.09, "Perturbation\n$\mathbf{\delta J=u\,v^T}$", COLORS["light_orange"], edge=COLORS["orange"])
    box(0.82, 0.135, 0.14, 0.065, "Mode Shift", "white", edge=COLORS["blue"])
    box(0.82, 0.045, 0.14, 0.065, "Anchor Fixed", "white", edge=COLORS["green"])
    ax.plot([0.72, 0.76], [0.12, 0.12], color=COLORS["ink"], lw=0.95, zorder=2)
    route([(0.76, 0.12), (0.76, 0.1675), (0.82, 0.1675)], color=COLORS["blue"])
    route([(0.76, 0.12), (0.76, 0.0775), (0.82, 0.0775)], color=COLORS["green"])
    route([(0.47, 0.305), (0.47, 0.39)], color=COLORS["blue"], dashed=True)
    route([(0.74, 0.305), (0.74, 0.39)], color=COLORS["blue"], dashed=True)
    ax.text(0.455, 0.345, "Response Geometry", ha="right", va="center", weight="bold", fontsize=5.7)
    ax.text(0.725, 0.345, "Mixer Separation", ha="right", va="center", weight="bold", fontsize=5.7)


def render_fig1(inputs: dict[str, Any] | None, *, schematic_only: bool) -> Any:
    _, plt, FancyBboxPatch = plotting_modules()
    if schematic_only:
        fig, ax = plt.subplots(figsize=(7.2, 3.55), constrained_layout=True)
        method_schematic(ax, FancyBboxPatch)
        return fig

    assert inputs is not None
    trials = inputs["fig1_response_trials.csv"]
    shells = inputs["fig1_shell_weights.csv"]
    fig = plt.figure(figsize=(7.2, 6.85), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=(3.0, 1, 1))
    schematic = fig.add_subplot(grid[0, :])
    method_schematic(schematic, FancyBboxPatch)
    schematic.text(-0.03, 1.02, "a", transform=schematic.transAxes, weight="bold", fontsize=9)

    scopes = ["random_matrix", "material_window", "global_k"]
    titles = ["Random-matrix verification", "Spectral crossings and metals", "Shared chemical potential"]
    axes = [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1]), fig.add_subplot(grid[2, 0])]
    palette = [COLORS["blue"], COLORS["green"], COLORS["orange"], COLORS["purple"]]
    for panel_index, (scope, title, ax) in enumerate(zip(scopes, titles, axes), start=1):
        grouped: dict[str, list[float]] = defaultdict(list)
        for row in trials:
            if row["scope"] == scope:
                grouped[row["system"]].append(float(row["distance_ratio"]))
        if not grouped:
            raise ReproductionError(f"fig1_response_trials.csv has no rows for scope={scope}")
        for color, (system, values) in zip(palette, sorted(grouped.items())):
            x, y = ecdf(values)
            ax.step(x, y, where="post", color=color, lw=1.4, label=f"{system} (n={len(values)})")
        ax.axvline(1.0, color=COLORS["red"], lw=1.0, ls="--")
        ax.set(xlabel="achieved / bound ratio", ylabel="ECDF", title=title, xlim=(0, 1.04), ylim=(0, 1.02))
        ax.legend(loc="upper left")
        ax.text(-0.14, 1.03, chr(ord("a") + panel_index), transform=ax.transAxes, weight="bold", fontsize=9)

    shell_ax = fig.add_subplot(grid[2, 1])
    ordered = sorted(shells, key=lambda row: float(row["g2_shell"]))
    shell_ax.step(
        [float(row["g2_shell"]) for row in ordered],
        [float(row["cumulative_weight"]) for row in ordered],
        where="post", color=COLORS["blue"], lw=1.5,
    )
    shell_ax.axhline(0.90, color=COLORS["gray"], ls="--", lw=0.9)
    shell_ax.axhline(0.99, color=COLORS["gray"], ls=":", lw=0.9)
    shell_ax.set(xlabel=r"$|G|^2$ shell", ylabel="cumulative response weight", title="Response-shell concentration", ylim=(0, 1.02))
    shell_ax.text(-0.14, 1.03, "e", transform=shell_ax.transAxes, weight="bold", fontsize=9)
    for ax in (schematic, *axes, shell_ax):
        bold_axis_text(ax)
        ax.title.set_fontweight("bold")
        for text_artist in ax.texts:
            text_artist.set_fontweight("bold")
        legend = ax.get_legend()
        if legend is not None:
            for label in legend.get_texts():
                label.set_fontweight("bold")
    return fig


def render_fig2(inputs: dict[str, Any]) -> Any:
    _, plt, _ = plotting_modules()
    rows = inputs["fig2_interventions.csv"]
    fig, axes = plt.subplots(3, 1, figsize=(3.39, 5.83), constrained_layout=True)
    quantities = [
        ("update_error", "rank-one update error", True),
        ("intrinsic_margin", r"fixed-point margin $m_{\rm fp}$", False),
        ("spectral_radius", r"$\rho(M_B)$", False),
    ]
    palette = [COLORS["blue"], COLORS["green"], COLORS["orange"], COLORS["purple"], COLORS["red"]]
    for panel, (key, ylabel, log_scale) in enumerate(quantities):
        grouped: dict[str, list[tuple[float, float]]] = defaultdict(list)
        for row in rows:
            value = row[key]
            if value is None:
                continue
            group = str(row["mixer"] or row["system"])
            grouped[group].append((float(row["coupling_fraction"]), float(value)))
        if not grouped:
            raise ReproductionError(f"fig2_interventions.csv has no values for {key}")
        for color, (group, values) in zip(palette, sorted(grouped.items())):
            values.sort()
            axes[panel].plot(
                [item[0] for item in values], [item[1] for item in values],
                marker="o", ms=3, lw=1.2, color=color, label=group,
            )
        if log_scale:
            axes[panel].set_yscale("log")
        axes[panel].set(xlabel=r"coupling fraction $|\lambda|/\lambda_c$", ylabel=ylabel)
        legend = axes[panel].legend(loc="best", fontsize=6)
        for label in legend.get_texts():
            label.set_fontweight("bold")
        close_axes(axes[panel])
        bold_axis_text(axes[panel])
        axes[panel].text(
            -0.11,
            1.02,
            chr(ord("a") + panel),
            transform=axes[panel].transAxes,
            weight="bold",
            fontsize=9,
        )
    return fig


def render_fig3(inputs: dict[str, Any]) -> Any:
    _, plt, _ = plotting_modules()
    summary = inputs["fig3_worklaw_summary.json"]
    scaling = inputs["fig3_local_scaling.csv"]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), constrained_layout=True)

    models = summary["models"]
    y = list(range(len(models)))[::-1]
    effects = [100.0 * float(model["improvement_fraction"]) for model in models]
    low = [100.0 * float(model["ci_low"]) for model in models]
    high = [100.0 * float(model["ci_high"]) for model in models]
    axes[0].errorbar(
        effects, y,
        xerr=[[effect - lo for effect, lo in zip(effects, low)], [hi - effect for effect, hi in zip(effects, high)]],
        fmt="o", color=COLORS["blue"], ecolor=COLORS["gray"], capsize=2,
    )
    axes[0].axvline(0, color=COLORS["gray"], lw=0.8)
    criteria = sorted({100.0 * float(model["criterion"]) for model in models})
    for criterion in criteria:
        axes[0].axvline(criterion, color=COLORS["red"], lw=1.0, ls="--")
        axes[0].text(
            criterion,
            0.96,
            f"{criterion:+g}% gate",
            transform=axes[0].get_xaxis_transform(),
            ha="right",
            va="top",
            fontsize=6.3,
            fontweight="bold",
        )
    axes[0].set(
        yticks=y,
        yticklabels=[model["name"] for model in models],
        xlabel="held-out error reduction (%)",
    )
    for yi, model, effect in zip(y, models, effects):
        criterion = 100.0 * float(model["criterion"])
        verdict = "PASS" if effect >= criterion else "FAIL"
        axes[0].annotate(
            f"{effect:+.2f}%",
            xy=(effect, yi),
            xytext=(0, 7),
            textcoords="offset points",
            ha="center",
            va="bottom",
            color=COLORS["green"] if verdict == "PASS" else COLORS["red"],
            fontsize=6.6,
            fontweight="bold",
        )
    overall_verdict = "PASS" if all(
        effect >= 100.0 * float(model["criterion"])
        for model, effect in zip(models, effects)
    ) else "FAIL"
    axes[0].text(
        0.96,
        0.82,
        overall_verdict,
        transform=axes[0].transAxes,
        ha="right",
        va="top",
        color=COLORS["green"] if overall_verdict == "PASS" else COLORS["red"],
        fontsize=7.4,
        fontweight="bold",
    )

    categories = summary["categories"]
    category_y = list(range(len(categories)))[::-1]
    for yi, item in zip(category_y, categories):
        if item["effect_fraction"] is None:
            axes[1].scatter(
                0,
                yi,
                marker="D",
                s=22,
                facecolor="white",
                edgecolor=COLORS["gray"],
                linewidth=1.0,
                zorder=3,
            )
            axes[1].annotate(
                "insufficient data",
                xy=(0, yi),
                xytext=(5, 0),
                textcoords="offset points",
                va="center",
                fontsize=6.2,
                fontweight="bold",
            )
            continue
        effect = 100.0 * float(item["effect_fraction"])
        color = COLORS["red"] if effect < -5.0 else COLORS["blue"]
        axes[1].plot([0, effect], [yi, yi], color=color, lw=1.2)
        axes[1].scatter(effect, yi, s=25, color=color, zorder=3)
        axes[1].annotate(
            f"{effect:+.1f}",
            xy=(effect, yi),
            xytext=(5 if effect <= 35 else -5, 4 if effect < -5 else 0),
            textcoords="offset points",
            ha="left" if effect <= 35 else "right",
            va="center",
            fontsize=6.2,
            fontweight="bold",
        )
    axes[1].axvline(0, color=COLORS["gray"], lw=0.8)
    axes[1].axvline(-5.0, color=COLORS["orange"], lw=1.0, ls="--")
    axes[1].set(
        yticks=category_y,
        yticklabels=[item["name"] for item in categories],
        xlabel="category error reduction (%)",
    )

    x = [float(row["adiabaticity_squared"]) for row in scaling]
    y_scale = [float(row["projector_distance"]) for row in scaling]
    if any(value <= 0 for value in x + y_scale):
        raise ReproductionError("fig3_local_scaling.csv requires positive values for log axes")
    log_x = [math.log(value) for value in x]
    log_y = [math.log(value) for value in y_scale]
    x_mean = sum(log_x) / len(log_x)
    y_mean = sum(log_y) / len(log_y)
    denominator = sum((value - x_mean) ** 2 for value in log_x)
    if denominator == 0:
        raise ReproductionError("fig3_local_scaling.csv requires varying adiabaticity values")
    slope = sum(
        (x_value - x_mean) * (y_value - y_mean)
        for x_value, y_value in zip(log_x, log_y)
    ) / denominator
    intercept = y_mean - slope * x_mean
    fit_x = [min(x), max(x)]
    fit_y = [math.exp(intercept) * value ** slope for value in fit_x]
    mean_ratio = sum(y_value / x_value for x_value, y_value in zip(x, y_scale)) / len(x)
    axes[2].scatter(x, y_scale, s=10, color=COLORS["blue"], linewidths=0)
    axes[2].plot(fit_x, fit_y, color=COLORS["gray"], lw=1.1)
    axes[2].set(
        xscale="log",
        yscale="log",
        xlabel=r"coupling-weighted adiabaticity $A_k$",
        ylabel=r"occupied-projector distance $d_P$",
    )
    axes[2].text(
        0.03,
        0.95,
        f"slope {slope:.6f}\nmean $d_P/A_k$ {mean_ratio:.6f}\n$n={len(x)}$ positive pairs",
        transform=axes[2].transAxes,
        va="top",
        fontsize=6.8,
        fontweight="bold",
        linespacing=1.15,
    )

    for panel, ax in enumerate(axes):
        close_axes(ax)
        bold_axis_text(ax)
        ax.text(
            -0.12,
            1.04,
            chr(ord("a") + panel),
            transform=ax.transAxes,
            weight="bold",
            fontsize=9,
        )
    return fig


def render_fig5(inputs: dict[str, Any]) -> Any:
    _, plt, _ = plotting_modules()
    eos = inputs["fig5_eos.csv"]
    bands = inputs["fig5_bands.csv"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), constrained_layout=True)
    palette = [COLORS["blue"], COLORS["gray"], COLORS["green"], COLORS["orange"]]

    eos_groups: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in eos:
        eos_groups[str(row["method"])].append((float(row["volume_a3_per_atom"]), float(row["energy_ev_per_atom"])))
    for color, (method, values) in zip(palette, sorted(eos_groups.items())):
        values.sort()
        minimum = min(value[1] for value in values)
        axes[0].plot(
            [value[0] for value in values], [1000.0 * (value[1] - minimum) for value in values],
            marker="o", ms=3, lw=1.2, color=color, label=method,
        )
    axes[0].set(xlabel=r"volume ($\AA^3$/atom)", ylabel="relative energy (meV/atom)", title="Equation of state")
    axes[0].legend(loc="best")

    band_groups: dict[tuple[str, int], list[tuple[float, float]]] = defaultdict(list)
    method_order: list[str] = []
    for row in bands:
        method = str(row["method"])
        if method not in method_order:
            method_order.append(method)
        band_groups[(method, int(row["band_index"]))].append((float(row["k_distance"]), float(row["energy_ev"])))
    for method_index, method in enumerate(method_order):
        first = True
        for (row_method, _), values in sorted(band_groups.items()):
            if row_method != method:
                continue
            values.sort()
            axes[1].plot(
                [value[0] for value in values], [value[1] for value in values],
                color=palette[method_index % len(palette)], lw=1.0,
                label=method if first else None,
            )
            first = False
    axes[1].axhline(0, color=COLORS["gray"], lw=0.7, ls="--")
    axes[1].set(xlabel="wave-vector path", ylabel="energy (eV)", title="Band structure")
    axes[1].legend(loc="best")
    for panel, ax in enumerate(axes):
        ax.text(-0.14, 1.04, chr(ord("a") + panel), transform=ax.transAxes, weight="bold", fontsize=9)
    return fig


RENDERERS = {
    "fig1_screen": render_fig1,
    "fig2_causal": render_fig2,
    "fig3_worklaw": render_fig3,
    "fig5_anchor": render_fig5,
}


def main() -> int:
    args = parse_args()
    root = repository_root()
    try:
        manifest = load_manifest(args.manifest or root / "manifest.json")
        figures = manifest["figures"]
        if args.list:
            for figure in figures:
                print(f"{figure['id']}\t{figure['display']}")
            return 0
        requested = set(args.only or [])
        known = {figure["id"] for figure in figures}
        unknown = sorted(requested - known)
        if unknown:
            raise ReproductionError(f"unknown figure identifier(s): {', '.join(unknown)}")
        selected = [figure for figure in figures if not requested or figure["id"] in requested]
        if args.schematic_only and {figure["id"] for figure in selected} != {"fig1_screen"}:
            raise ReproductionError("--schematic-only requires --only fig1_screen")
        data_root = args.data_root or root / manifest["default_data_root"]
        output_dir = args.output_dir or root / "outputs" / "figures"
        for figure in selected:
            if figure["id"] == "fig1_screen" and args.schematic_only:
                fig = render_fig1(None, schematic_only=True)
            else:
                inputs = load_figure_inputs(figure, data_root)
                renderer = RENDERERS[figure["id"]]
                fig = renderer(inputs) if figure["id"] != "fig1_screen" else renderer(inputs, schematic_only=False)
            save_figure(fig, output_dir, figure["id"])
            import matplotlib.pyplot as plt
            plt.close(fig)
    except ReproductionError as exc:
        print(f"INCOMPLETE: {exc}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Build the five figures and three material tables from this version's inputs."""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "build")
    args = parser.parse_args()
    output = args.output.resolve()
    os.environ["IDMATE_FIGURE_OUTPUT"] = str(output / "figures")
    sys.path.insert(0, str(HERE / "figures"))
    import common
    import fig1
    import fig2_scope
    import fig3_physics
    import fig_supplement
    from style import plt, save_figure, use_style

    subprocess.run([sys.executable, str(HERE / "verify.py")], check=True)
    fig1.main()
    fig2_scope.main()
    plt.close("all")
    physics = common.load_physics()
    fig, axes = fig3_physics.draw(physics)
    fig3_physics.check_visuals(fig, axes, physics)
    save_figure(fig, common.OUTPUT / "fig5_anchor")
    plt.close(fig)
    use_style()
    plt.rcParams.update({"xtick.labelsize": 7.4, "ytick.labelsize": 7.4,
                         "legend.fontsize": 7.4, "path.simplify": True,
                         "svg.hashsalt": "idmate-supplement-2026-09-05"})
    supplemental = common.load_supplement()
    for name, builder in (("fig2_causal", fig_supplement.fig_s1),
                          ("fig3_worklaw", fig_supplement.fig_s2)):
        fig, series = builder(supplemental)
        fig_supplement.figure_checks(fig, series)
        save_figure(fig, common.OUTPUT / name)
        plt.close(fig)

    spec = importlib.util.spec_from_file_location("material_tables", HERE / "materials" / "reproduce.py")
    material = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(material)
    material.main(output / "material_tables")
    subprocess.run([sys.executable, str(HERE / "numerical_clarifications" / "reproduce.py"), "--check"], check=True)
    report = {"figures": 5, "material_tables": 3,
              "DFT_calculations_repeated": False, "stored_work_regression_interval_recomputed": False,
              "numerical_clarifications_checked": True,
              "material_separation_row_bootstrap_recomputed": True,
              "silicon_band_mae_eV": float(physics["bands"]["metrics"]["mae_eV"]),
              "silicon_band_max_error_eV": float(physics["bands"]["metrics"]["max_error_eV"]),
              "projector_log_slope": float(supplemental["slope"])}
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

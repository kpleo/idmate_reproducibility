# Data dictionary

All numerical files in this document are expected inputs and are not included
in the repository. Column requirements are machine-readable in
`manifest.json`. Missing values are permitted only for fields marked
`nullable`.

## Figure inputs

### `fig1_response_trials.csv`

One row per finite-temperature bound evaluation.

- `scope`: `random_matrix`, `material_window`, or `global_k`.
- `system`: physical system or ensemble name.
- `trial_id`: identifier unique within `scope` and `system`.
- `temperature_ha`: electronic temperature in Hartree.
- `distance_ratio`: density distance divided by its response bound.
- `gap_ratio`: free-energy gap divided by its response bound.
- `violation`: logical result of the declared numerical criterion.

### `fig1_shell_weights.csv`

Reciprocal-shell decomposition of the response derivative. `weight` is the
fraction assigned to one shell and `cumulative_weight` is its cumulative sum
in ascending `g2_shell` order. The same file supplies SI table
`tab:shell-sequence`.

### `fig2_interventions.csv`

One row per response intervention, optionally repeated for each mixer.
`coupling_fraction` is expressed relative to the critical coupling;
`update_error` tests the rank-one response update; `intrinsic_margin` describes
the electronic fixed-point margin; and `spectral_radius` describes the fixed-
mixer iteration matrix. Mixer-specific fields may be empty on rows that carry
only electronic quantities.

### `fig3_worklaw_summary.json`

Grouped prediction results, category effects, and experiment counts. The
required nested structure is defined by
`schemas/fig3_worklaw_summary.schema.json`. Improvement and interval values are
fractions, not percentages.

### `fig3_local_scaling.csv`

Positive per-step pairs used for the local log-log scaling test.
`adiabaticity_squared` is the squared perturbative predictor and
`projector_distance` is the occupied-projector distance.

### `fig5_eos.csv`

Equation-of-state points. Energies are in eV per atom. The plotting script
subtracts the minimum energy separately for each `method`; absolute energy
zeros are not compared.

### `fig5_bands.csv`

Long-form band data. Each row identifies a `method`, `k_distance`,
`band_index`, and energy in eV relative to the stated electronic reference.
`k_label` is populated only at symmetry points.

## Table inputs

Table CSV files live under `tables/` relative to the selected data root. Each
file corresponds to one block in `manifest.json`:

| Manuscript label | Expected file or files |
|---|---|
| `tab:three-zone-decisions` | `tables/main_controller_decisions.csv`, `tables/main_controller_terminal.csv` |
| `tab:mixing-baselines` | `tables/main_mixing_baselines.csv` |
| `tab:physical-convergence` | `tables/main_physical_convergence.csv` |
| `tab:mg-factorial-rungs` | `tables/si_mg_factorial_rungs.csv` |
| `tab:k2-nested-rungs` | `tables/si_k2_nested_rungs.csv` |
| `tab:numerical-parameters` | `tables/si_numerical_parameters.csv` |
| `tab:hgh` | `tables/si_hgh_parameters.csv` |
| `tab:shell-sequence` | `fig1_shell_weights.csv` |
| `tab:work-category` | `tables/si_work_category.csv` |
| `tab:h2-mixer-radii` | `tables/si_h2_mixer_radii.csv` |
| `tab:silicon-fixed-b` | `tables/si_silicon_mixer_radii.csv` |
| `tab:three-system-controller` | `tables/si_three_system_controller.csv` |
| `tab:slack-verification` | `tables/si_screen_slack.csv` |
| `tab:reference-decisions` | `tables/si_reference_decisions.csv` |
| `tab:accepted-reference-rows` | `tables/si_accepted_reference_rows.csv` |
| `tab:production-accepted` | `tables/si_production_accepted.csv` |
| `tab:production-cost` | `tables/si_production_cost.csv` |
| `tab:mixing-comparison` | `tables/si_mixing_comparison.csv` |

The scripts retain full input precision for validation and computation. Plotted
results are formatted to three decimal places, with scientific notation or a
fourth decimal used only when required by scale; table values are reproduced as
supplied. They do not reconstruct missing values or substitute values from
manuscript prose.

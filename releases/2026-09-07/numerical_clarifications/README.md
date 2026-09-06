# IDMate: numerical clarifications

This packet reconstructs four sets of numerical comparisons from the IDMate
study: historical window-bound comparisons, material-stratified candidate
separation, silicon response-shell profiles, and the H2 subspace-reuse example.
It contains processed numerical observations and independent analysis code.
Python 3.9 or later is sufficient; no third-party libraries or network access
are required.

## Reproduce

Run these commands from this directory:

```sh
python3 -B reproduce.py
python3 -B reproduce.py --check
python3 -B -m unittest -v test_reproduce.py
```

The first command writes `results.json` and `results.md`. The second verifies
input checksums and compares a fresh calculation with both saved results.
Counts and text must agree exactly; derived floating-point values allow a
relative difference of at most 2e-13 for differences between math libraries.
The tests cover the numerical results, input validation, ties, stratification,
rounding, and execution independently of the working directory. All output is
written beside the script. Numerical inputs retain their available precision;
the Markdown tables use three decimal places, including scientific-notation
mantissas. Statistics are calculated before rounding.

## Numerical inputs

| File | Contents |
| --- | --- |
| `historical_candidates.csv` | Ten accepted candidates: raw bound, window and full distances, and recorded trace deviation. |
| `material_candidates.csv` | The 58 accepted or rejected candidates used for the separation analysis. |
| `excluded_candidates.csv` | The eight abstentions excluded from that binary comparison. |
| `shell_profiles.csv` | All printed silicon shell shares: 37, 142, and 142 rows in three configurations. |
| `shell_probe_scattering.csv` | Nine input-shell groups from the exploratory silicon calculation. |
| `shell_configurations.json` | Numerical settings and reported shell counts for each configuration. |
| `h2_reuse.csv` | Sixty consecutive H2 reuse attempts, with spectral margins and individual decision conditions. |
| `source_coverage.json` | Original-record checksums, extraction checks, numerical definitions, and coverage limits. |
| `input_manifest.json` | Checksums of the supplied numerical inputs. |

Candidate identifiers distinguish separate single-proposal calculations from
the three calculations containing multiple proposals. `observation_order`
preserves the order used by the original bootstrap calculation; it is not a
ranking by distance. Proposal families are `low_cutoff`, `target_subspace`,
and `coarse_k`. The target-subspace proposals use a target-converged reference
subspace and are not prospective predictions. The ten historical rows overlap
the multi-proposal subset of the 58-row comparison and are not additional
independent samples.

## Definitions and interpretation

### Window bound and trace deviation

`raw_window_bound` is the historical bound without the trace correction.
`window_distance` measures the candidate's distance from the canonical
compressed current-Hamiltonian reference; `full_distance` uses the full current
reference map. Both use the k-weighted Frobenius norm of the real embedding of
the per-spin density matrix. `trace_deviation` is the absolute global embedded
trace mismatch. These quantities are dimensionless. The real embedding doubles
the per-spin trace mathematically; this is distinct from physical twofold spin
occupation.

The script reports both `window_distance > raw_window_bound` and
`window_distance > 0.05`. These are different tests. It also compares each
positive bound excess with that same row's recorded trace deviation. This last
comparison is empirical: it neither identifies the cause of the excess nor
reconstructs the corrected bound used in subsequent calculations. The original
window occupation margins and complete density matrices are not supplied.

### Material-stratified separation

The binary comparison contains 40 accepted and 18 rejected candidates from
silicon, graphene, and fcc aluminium. All eight abstentions are graphene
observations and are listed separately. The quantity called AUC here is

```text
P(full_distance_reject > full_distance_accept) + 0.5 P(tie).
```

It measures separation by full-space distance between decision groups, not the
predictive accuracy of a fitted classifier or a continuous-bound regression.
Every accepted/rejected pair contributes once. The report separates pairs
within the same material from pairs across materials and reports all nine
material-pair cells. A pair-weighted within-material AUC and an equal-material
average are different summaries and are labelled separately.

Deletion values exclude one material and recompute the statistic on the other
two; there is no training or prediction on a held-out material. Every rejected
candidate belongs to the coarse-k family, whereas accepted candidates belong
to the other two families. The decision and proposal-family effects therefore
cannot be isolated in these observations. The material-specific associations
differ, and the pooled result does not establish transferable separation.

The 95% percentile interval uses 10,000 independent resamples of the accepted
and rejected candidate lists with seed 20260829. It reproduces the original
candidate-level calculation. It does not resample materials or SCF trajectories
as clusters and is not an uncertainty interval for a population of materials.
Within-material associations vary rather than vanishing uniformly.

### Silicon response shells

All three configurations use diamond silicon with a conventional lattice
constant of 10.26 bohr. The shell key is the integer S in
`|G| = (2 pi / a) sqrt(S)`. The table contains both S and the originally printed
radius. A shell share is its output derivative-energy contribution divided by
the energy across all shells and all measured derivative columns. The
cumulative share includes all shells up to S, including unprinted tiny shares.

The exploratory calculation has 128 derivative columns at 6 Ha, a Gamma-only
mesh, and a 13 x 13 x 13 density grid. The selection and frozen calculations
have respectively 136 and 168 columns at 25 Ha, a 4 x 4 x 4 mesh, and a
33 x 33 x 33 grid. The shell projection measures the same type of quantity,
but the numerical settings and probe ensembles differ. These profiles are not
the separate 16-direction adaptive-response calculation.

Only shells with reported energy share greater than 1e-12 were printed. There
are 99 total shells in the exploratory calculation and 715 in each of the
other two. The supplied rows reproduce the reported 90% and 99% crossings,
the selected shell shares, and the exploratory decay fit. They do not reproduce
full matrix closure checks. The fit uses shares greater than 1e-8 and exact
radii derived from S and a; the original shares themselves were printed with
seven significant digits. Self-shell fractions are means over the columns of
each input-shell group, not cumulative output-shell shares.

### H2 spectral conditions

The H2 data are a single fixed-grid model trajectory, not a multi-material
validation. The retained spectral margin is

```text
margin = anchor_guard_energy - anchor_retained_edge_energy
         - 2 potential_change_linf.
```

Strict positivity is the retained-state separation condition obtained from
Weyl's eigenvalue perturbation inequality. The occupied-state margin is tested
separately. The full provisional reuse test additionally requires convergence
of the candidate residuals, positive candidate separations, a density bound
within tolerance, and orthogonality within tolerance. The script reconstructs
the composite decision from these recorded conditions and checks it against
the recorded decision.

The number of composite passes is not a count of locality or perturbation-law
agreement. The supplied first-order parameter `adiabatic_number_squared` and
normalized occupied-projector distance provide a separate comparison:
`A = sqrt(adiabatic_number_squared)` and `abs(d_P / A - 1)`. They are not used
by the reuse decision. The numerical summary labels the mean relative errors
with denominator A and with denominator d_P separately; the historical mean
discrepancy uses d_P. All decisions remain provisional floating-point tests,
not interval-enclosed spectral certificates. Their association with late SCF
iteration is visible in the supplied rows.

## Coverage

This is a self-contained reconstruction of the four numerical analyses above.
It does not regenerate the underlying SCF calculations, proposal density
matrices, finite-difference response matrices, all manuscript figures or
tables, VASP comparisons, or timing measurements. It contains neither the DFT
implementation nor licensed pseudopotentials. Checksums identify the original
records used in preparing the numerical extracts; those original files are
not required to run this packet and are not included. Exact historical source
versions cannot be inferred from a numerical match alone.

The reported bootstrap is descriptive at the candidate level. The shell
comparisons distinguish numerical configurations. The trace comparison does
not convert the old raw bound into a subsequently corrected one. These scope
distinctions are part of the numerical conclusions, not missing computations.

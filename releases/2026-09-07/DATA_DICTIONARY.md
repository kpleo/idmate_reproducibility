# Numerical inputs

All paths below are relative to this release. JSON is UTF-8 with finite
numerical values or explicit `null`. CSV uses a header row and comma delimiter.
Stored numerical precision is not a statement of physical accuracy.

## Distance distributions

`data/e0.csv`, `e3_graphene.csv`, `e3_al.csv`, `gk_graphene.csv`, and
`gk_al.csv` contain 10,000, 164, 1,614, 310 and 400 rows, respectively. Their
columns are `distance`, `bound`, `allowance` and `ratio`. The first three are
nonnegative, dimensionless density-matrix distances or their comparison
tolerances. `ratio` is the saved value of `distance/bound`; a blank means the
original calculation did not assign a ratio at its numerical resolution, not
a missing trial or a zero ratio. The corresponding bound can be positive;
the original missing-value mask is preserved without inferring a new threshold.
There are 114 such rows
in the matrix sample and two in the weighted graphene sample.

The `e0` rows are seeded finite matrices, `e3` rows use individual material
blocks, and `gk` rows use positive normalized k-point weights and a common
chemical potential. `allowance` is the recorded numerical comparison
tolerance, not an addition to the exact theorem. Row order is preserved.
The release recomputes the ratios and retains every defined ratio in the ECDF.

`data/shells.csv` contains the 142 printed non-negligible reciprocal shells
from the frozen 168-column silicon calculation. `shell` is squared reciprocal
radius in units `(2 pi/a)^2`, `g_bohr_inv` has inverse-bohr units, and `share`
and `cumulative` are derivative-energy fractions. `retained` identifies the
shells used by the truncated stencil. These printed shells are not the full
715-shell derivative matrix; no full-matrix reconstruction is implied.

## Historical window comparison

`data/screen_scope.csv` contains ten historical accepted candidates. It records
material, iteration, candidate family, uncorrected `window_bound`, measured
`window_error`, `full_error`, and `window_over_bound`. Distances use the recorded
real-embedded weighted norm. The `online_low_cutoff` family is distinct from
the target-converged subspace comparison. These are not the newer three
historical-subspace candidates.

## Silicon energy and band comparison

`data/silicon_physics.json` contains two five-point energy-volume curves and
their stored third-order Birch-Murnaghan parameters. Volume is in cubic
angstroms and energy is in eV per two-atom cell; the plot subtracts each fit's
energy minimum and divides by two atoms. Fit coefficients are retained, not
refitted in this version. The residuals and the plotted curves are recomputed.

Band energies are individually referenced to each method's valence-band
maximum. Five selected bands are supplied at 31 candidate points and 78
reference points. The main display uses the 39 reference points on L-Gamma-X;
the comparison retains all 78 reference samples for interpolation. The
`candidate_grid_ang_inv` and `reference_grid_ang_inv` arrays define the
physical reciprocal distances used for that interpolation. `cx` and `rx` are
segment-fractional display coordinates, not the metric coordinates. The mean
and maximum absolute errors are recomputed over all 155 paired values.

`figures/structure_cells.json` supplies real-space cells in bohr and fractional
atomic positions for the crystal motifs. Visual repeats identify the material;
they do not indicate the number of atoms in a calculation.

## Supplemental numerical displays

`data/supplement.json` includes rank-one intervention fractions and errors,
intrinsic margins, dimensionless mixer spectral radii, and stored regression
estimates. Couplings retain their original units; plotted coupling axes divide
by the recorded critical or maximum coupling. Fractional improvements and
interval endpoints are converted to percentages only for display. Category
identifiers C1-C8 retain the supplemental definitions. Parent identifiers are
anonymous integers used only to retain the number of lineages.

The 60 positive `ak` and `dp` values describe one H2 trajectory, with
dimensionless adiabaticity and occupied-projector distance. The displayed
log-log slope and intercept are refitted from these pairs. The regression
improvement interval is a stored hierarchical-bootstrap result: its endpoints
are displayed, not statistically reconstructed by this plotting package.

## Material tables

`materials/materials.json` contains three scalar comparisons, calculation
settings and candidate descriptions. Distances are defined in
`materials/README.md`. `null` window values denote abstention. The six
reference-map count means two maps per material. The separate
`numerical_crosscheck.json` records the comparison count and tolerances of an
independent reconstruction; it does not distribute the underlying matrices.

## Numerical clarifications

`numerical_clarifications/README.md` documents its eight numerical inputs,
their quantities, statistical units, original-record checksums and limits.
Its `results.json` and `results.md` are frozen expected scientific summaries,
included so `reproduce.py --check` can independently recalculate and compare
every result. They are not full calculation outputs. The 58 binary-decision
observations include the ten historical candidates and are not independent
additional observations. Eight abstentions are retained in a separate input.

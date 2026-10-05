# Numerical inputs, release 2026-10-05

Paths are relative to this release. JSON is UTF-8 with finite numbers or explicit `null`; CSV has a header row
and comma delimiters. Distances are density-matrix distances in the embedded weighted metric of the manuscript,
d(A,B) = [2 sum_k w_k ||A_k - B_k||_F^2]^(1/2) for per-spin matrices (dimensionless); energies are in hartree.
Stored precision is not a statement of physical accuracy.

## Certificate verification (Fig. 1c)

`data/e0.csv`, `e3_graphene.csv`, `e3_al.csv`, `gk_graphene.csv`, `gk_al.csv` (10 000, 164, 1 614, 310 and 400
rows) are identical to the files of release 2026-09-07: columns `distance`, `bound`, `allowance`, `ratio`. A blank
`ratio` means the original calculation assigned no ratio at its numerical resolution (114 rows of `e0`, two of
`gk_graphene`); 12 372 ratios are defined.

## Block matrices (Figs. 1a,b and S1)

`data/blocks_{si,graphene,al}.json`: for the first k point of each saved Hamiltonian, `m` and `n` (retained and
total dimensions), `tau_ha`, `mu_ha`, `ritz_values_ha` (the m Ritz values), `H_abs_ha` (n x n magnitudes of H in
the basis [Ritz vectors of P | eigenvectors of QHQ], hartree) and `dD_abs` (n x n magnitudes of D_c - D* in the
same basis), rounded to six significant digits.

## In-loop SCF candidates (Fig. 2a,c)

`data/material_candidates.csv`: the 58 screened candidates. Columns: `observation_order`, `candidate_id`, `run_id`
(`<material>_multiple` for the three multiple-insertion reference configurations, `<material>_single_<family>_<it>`
for single-insertion runs), `material`, `iteration` (insertion iteration), `proposal_family` (`low_cutoff`,
`target_subspace` = oracle subspace, `coarse_k`), `decision` (`accept`, `reject`), `raw_window_bound` (uncorrected
bound B_u), `window_distance` (d_W), `full_distance` (d_full).
`data/excluded_candidates.csv`: the eight abstentions (no window distance), same columns without the window ones.

## In-loop rerun with saved Hamiltonians (Fig. 2b, Sec. IV D, Tables S12 and S13)

`data/inloop_results.json`: one record per presented candidate (66). Identification: `leg` (harness run name),
`run_id`, `material`, `iteration`, `kind` (`low_cutoff_continuation`, `converged_subspace_reuse`,
`coarse_k_mesh_upsample`), `decision`, `abstain_reason`, `dump` (stored record name). Logged quantities:
`window_bound`, `window_bands`, `d_window`, `d_full_recorded`, `mu_c_recorded`, `mu_exact_recorded`. Setup: `tau`
(Ha), `N` (electrons), `nk`, `n` and `m` (plane waves and candidate orbitals per k point). Checks:
`d_full_reproduced` (recomputed from the stored frames), `ref_vs_exact_frame` (distance between D* from full
diagonalization of the stored H_k and the stored exact frame), `herm_max`; for Rayleigh-Ritz candidates (`ritz` =
true) also `orthonormality`, `ritz_offdiag_rel` (relative off-diagonal norm of P^dag H P), `ritz_vs_recorded_eigs`,
`candidate_reconstruction`, `identity_residual` (residual of the exact block identity), `mu_c_recomputed`.
Decomposition (Ritz candidates; all embedded unless noted): `d_full`; exact blocks `blk_P`, `blk_Q`, `blk_off`
(delta_P, delta_Q, delta_E); zeroth order `eP`, `eQ`; first order `d_dk` (d_E^(1)); `pyth` (leading-order sum);
second-order estimates `eP2`, `eQ2`, `pyth2`; coupling norm `c_w` (||QHP||_w, unembedded, Ha) and the gap-free
coupling bound `bound_cpl` (B_E); `rigorous` (bound with exact e_P, e_Q) and `rigorous_computable` (with the
computable estimates); indicator pieces `eP_est`, `eQ_est`, `d_cheap` (c_E hat), `C_hat`, `indicator` (I_1),
`indicator2` (I_1 plus the second-order near-Fermi estimate), `eP2_est`; `d_cpl_exact` (||D_0* - D*||); chemical
potentials `mu_c`, `mu0`, `mu_star` (Ha); `gapQ` (lowest complement level minus mu_0, Ha); `topP` (highest Ritz
level minus mu_0, Ha). Coarse-k records carry the identification, logged and check fields only.

`data/rerun_provenance.json`: code commit, patch and binary checksums, toolchain, test-suite result, hosts
(cluster A: Intel Xeon E5-2680 v4; cluster B: AMD EPYC 7H12), one record per run (material, insertion scheme,
family, iteration, cluster, wall time in s, exit code), and the comparison of the rerun ledgers with the original
records: for the single-insertion runs the number of arms compared and bitwise identical; for each
multiple-insertion configuration, per arm (`reference_scoring`, `integrated`), decision and iteration agreement,
whether the density trajectory is bitwise identical, the largest relative differences of the logged distances,
and the terminal free energies; and the size and format of the stored records (not distributed).

## Saved-Hamiltonian decomposition (Fig. 3, Tables II, S1, S2)

`data/decomposition_results.json`: per material (`si`, `graphene`, `al`): `tau`, `N`, `nk`, `w_min`,
`extra_sweep` (every buffer size b; field `extra`) and `tau_sweep` (15 temperatures at b = 2; field
`duplicate_of_buffer_sweep` marks the two temperatures that coincide with the buffer sweep). Row fields as in
`inloop_results.json`, plus `eps_rho` (normalized real-space density error), `eP_est_crude`, `mu_hat`.
`data/table2_decomposition.json`: the three saved candidates of Table II (`B_W`, `d_W`, the blocks, `dk`, `pyth`,
`d_full`, `cheap` = I_1, `rigorous`, `eps_rho`, and the saved-candidate cross-checks `d_full_saved`,
`eps_rho_saved`).

## Number ledgers

`data/stats_ledger.json` and `data/inloop_ledger.json`: every statistic quoted from the two decomposition studies,
as computed by `analysis/stats_ledger.py` and `analysis/inloop_ledger.py`.

## Physical comparisons and auxiliary studies (Figs. S2-S4)

`data/silicon_physics.json` (equation-of-state fits and points for IDMate and VASP, band paths and their
differences) and `data/supplement.json` (rank-one intervention errors, intrinsic margins, mixer spectral radii,
the stored work-regression verdict and the 60 projector-distance pairs) are identical to the inputs of the
earlier release where they overlap.

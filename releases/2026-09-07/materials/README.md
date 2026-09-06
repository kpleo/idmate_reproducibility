# Historical-subspace material comparison

This addition provides the scalar results and electronic-structure settings
for the three single-step comparisons in the IDMate manuscript. It concerns
diamond silicon, graphene, and face-centred cubic aluminium.

Run `python3 reproduce.py` to reconstruct the main-table rows, the two
Supplemental Material tables, and a numerical CSV from `materials.json`.
Python 3.6 or later suffices; no external package,
network access, or solver installation is needed for this table reconstruction.
Run `python3 -m unittest -v test_reproduce` for the reconstruction tests.
`numerical_crosscheck.json` summarizes an independent NumPy diagonalization
and density reconstruction: 1,135 numerical comparisons across all three
materials, with no failed comparisons. The largest embedded difference between
complete reference density matrices is 2.742e-13. Comparison tolerances are
1e-10 absolute and 1e-9 relative; the original screening decisions are unchanged.

## Quantities

- `window_bound`: the residual distance term plus the weighted particle-number
  correction, in the same real-embedded norm as the measured window distance.
- `window_distance`: distance to the canonical state of the compressed
  Hamiltonian at its prescribed window charge.
- `full_distance`: `sqrt(2 sum_k w_k ||Q_k-Qref_k||_F^2)`, with complex per-spin
  density matrices and positive normalized k-point weights.
- `density_distance`: `sqrt(Omega integral(delta rho)^2)/N_e`, a dimensionless
  real-space density difference.
- `next_mixed_density_distance`: the difference between two linear mixed
  inputs. Before separate charge normalization it equals `alpha*density_distance`.
- `delta_F_full_ha`: signed spin-degenerate fixed-Hamiltonian free-energy
  difference in hartree, not a nonlinear Kohn-Sham total-energy difference.
- `delta_F_window_ha`: the analogous target-window-charge difference. It is
  not assigned a residual-only energy bound when the candidate trace differs.
- `maps`: reference-map evaluations, not individual eigensolver calls.

All prescribed points are retained. Missing window quantities remain missing,
not zero. The three materials are individual numerical examples, not a
statistical sample establishing generalization to other materials.

## Calculation and scope

Each candidate uses retained orbitals from the first reference map, starting
from a neutral uniform density. Rayleigh-Ritz projection uses the Hamiltonian
at the next mixed input. It does not use the complete solution of that current
Hamiltonian. The candidate and screening decision precede reference scoring.

The settings include cell vectors, fractional atomic positions, the analytic
HGH LDA potential family, k points, grid dimensions, cutoffs, temperatures,
and mixing factors. Potential coefficient tables are not distributed in this
addition. The candidate keeps bands through the last historical
occupation above 1e-12, followed by exactly two additional bands. Its omitted
orthogonal complement is zero. No target-converged subspace is used.

The numerical solver and full matrix records are maintained separately. This
addition reconstructs the reported table; it does not independently regenerate
the underlying DFT maps without that solver. There is no nonlinear convergence
or runtime acceleration measurement in this three-point experiment.

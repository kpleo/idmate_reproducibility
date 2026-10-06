# IDMate figures, tables and numerical checks, release 2026-10-05

This release accompanies the revised manuscript "Gap-free error bounds for compressed Kohn-Sham states at
finite temperature". It reconstructs the three main figures, the four Supplemental figures and the
Supplemental tables of the decomposition and of the in-loop rerun from compact derived records, and it
recomputes the numbers quoted in the text. The manuscript itself is not included.

## Run

Python 3.10 or later with `requirements.txt` (NumPy, Matplotlib). The figures are typeset in Liberation Sans
(text and math); the build stops if that font is not installed. From this directory:

```sh
python3 build.py                                   # figures, tables, tests and input check
python3 -m unittest discover -s tests -v           # numbers quoted in the manuscript
python3 tables/make_tables.py                      # Supplemental table rows, compared with the typeset rows
python3 verify.py                                  # file inventory and checksums
```

Outputs are written under `build/` and are not distributed.

## Coverage

| Display | Inputs (data/) | Operation |
|---|---|---|
| Fig. 1 | `blocks_si.json`, `e0.csv`, `e3_*.csv`, `gk_*.csv` | Block magnitudes of H and of the full-state error; distributions of 12 372 distance/bound ratios |
| Fig. 2 | `material_candidates.csv`, `excluded_candidates.csv`, `inloop_results.json` | Window and full-state distances of the 58 screened and 8 abstained in-loop candidates; full-state distance versus the indicator for the 48 Rayleigh-Ritz candidates of the rerun; distance versus insertion iteration |
| Fig. 3 | `decomposition_results.json`, `table2_decomposition.json` | Buffer and temperature sweeps of the decomposition on the three saved Hamiltonians (296 candidates) |
| Fig. S1 | `blocks_{si,graphene,al}.json` | Full-state error blocks in the three materials |
| Fig. S2 | `silicon_physics.json` | Silicon equation-of-state and band comparisons |
| Figs. S3, S4 | `supplement.json` | Interventions and mixer response; work-counter regression |
| Tables S1, S2 | `decomposition_results.json` | Buffer and temperature sweeps |
| Tables S12, S13 | `inloop_results.json` | In-loop rerun: summary by family and the 66 presented candidates |

`DATA_DICTIONARY.md` defines every input. `analysis/` holds the scripts that produced the derived records;
see `analysis/README.md` for what they need.

## In-loop rerun with saved Hamiltonians

The SCF tests of the manuscript were rerun on 2026-10-05 with the code version that produced the original
records and a test-only hook that stores, at every insertion, the screen's Hamiltonian, the candidate orbitals
and occupations, and the exact reference frame. `data/rerun_provenance.json` records the hosts, the binary and
patch checksums, the per-run wall times, and the comparison with the original records: the 48 single-insertion
runs reproduce every logged decision, distance and density-trajectory fingerprint bit for bit in both arms; the
three multiple-insertion configurations, originally run on arm64, give on x86-64 the same decisions and iteration
counts, full-state distances that agree to 2e-11 (relative) in the reference-scoring arm analysed here and to
1.2e-8 in the integrated arm, window bounds and distances at the rounding floor that differ by up to 14%, and
terminal free energies that agree to 4.4e-15 Ha. `inloop_results.json` holds, for each of the 66
presented candidates, the logged screen and reference distances, the distance recomputed from the stored
records, the consistency checks, and every term of the decomposition and of the indicator.

## Scope

This is reconstruction of displayed and quoted results, not regeneration of the electronic-structure
calculations. The saved Hamiltonian and orbital records of the decomposition (about 27 MB) and the stored
in-loop records (about 0.27 GB) are wave-function-level data and remain outside this repository, as do the
solver, the potential files and all large calculation outputs; they are available from the authors. The
scripts in `analysis/` document how the derived records were computed from them.

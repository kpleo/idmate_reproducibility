# Analysis scripts

These scripts produced the derived records in `../data`. They are distributed to document the computation;
they do not run on the derived records alone.

| Script | Produces | Needs |
|---|---|---|
| `fullstate.py` | estimators shared by all analyses (D*, exact blocks, first- and second-order terms, I_1, bounds) | NumPy |
| `run_all.py` | `decomposition_results.json` (buffer and temperature sweeps) | the three saved Hamiltonian records (`IDMATE_H2_RUNS`) |
| `export_heatmap.py` | block matrices of Figs. 1 and S1 | the saved Hamiltonian records |
| `stats_ledger.py`, `make_si_tables.py` | number ledger and Tables S1, S2 of the saved-Hamiltonian study | `decomposition_results.json` |
| `inloop.py` | `inloop_results.json` | the stored in-loop records (`IDMATE_INLOOP_RUNS`) |
| `inloop_ledger.py`, `make_si_inloop_table.py` | number ledger and tables of the in-loop rerun | `inloop_results.json` |
| `verify_gr16.py`, `verify_secondorder_blocks.py`, `indicator_mu_check.py` | checks quoted in the Supplemental Material | the saved Hamiltonian records |

The records (about 27 MB for the saved Hamiltonians, 0.27 GB for the in-loop records) are available from the
authors; point `IDMATE_H2_RUNS` and `IDMATE_INLOOP_RUNS` at them (defaults: `./h2_records` and
`../inloop_rerun_2026-10-05/runs`). Stored in-loop matrices are little-endian complex128 arrays in column-major order with JSON metadata,
one directory per presented candidate. The scripts write next to themselves (the ledgers and tables of this
release were copied from those outputs; `../tables/make_tables.py` regenerates the tables from `../data`).

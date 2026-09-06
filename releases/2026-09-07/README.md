# IDMate numerical figures and material comparisons

This version reconstructs the three main figures and two supplemental figures
from compact, derived numerical records. It also reconstructs the three tables
for the historical-subspace comparison of silicon, graphene and aluminium.
The source files contain numerical inputs and plotting code, not the manuscript.

The earlier top-level interface is retained for historical compatibility. Its
four-figure numbering and missing-input status do not describe this version.

## Run

Python 3.10 or later is required. Install `requirements.txt` in a separate
environment, then run from this directory:

```sh
python3 build.py
python3 -m unittest discover -s tests -v
python3 -m unittest discover -s materials -v
python3 -m unittest discover -s numerical_clarifications -v
python3 verify.py
```

Generated figures, tables and numerical checks are written under `build/` and
are not distributed. The build has no network or electronic-structure software
dependency. Rendering requires Arial and Times New Roman installed locally;
the script stops with a specific font error when either is missing. Silent
font substitutions are not used because their different text widths can alter
the compact schematic layout. No font files are distributed. The numerical
tests and material-table reconstruction do not require those fonts.

## Coverage

| Display | Input | Operation |
|---|---|---|
| Main Fig. 1 | Five distance-series CSV files and shell weights | Empirical cumulative distributions and finite-difference shell profile; method schematic |
| Main Fig. 2 | Ten historical candidates | Window and full-state distances; uncorrected ratios `d_W/B_0` |
| Main Fig. 3 | Silicon energy-volume and band records | Stored-fit evaluation, energy alignment, interpolation, and crystal motifs |
| Supplemental Fig. S1 | Rank-one interventions and mixer radii | Coupling, intrinsic-margin and mixer-response comparison |
| Supplemental Fig. S2 | Stored regression results and 60 projector-distance pairs | Reported intervals and category estimates; local scaling regression |
| Three material tables | `materials/materials.json` | Main comparison, calculation settings and density/free-energy diagnostics |
| Numerical clarifications | `numerical_clarifications/` | Ten-row trace comparison, 58-row stratified separation, shell profiles and 60 H2 reuse attempts |

Every numerical point required for these five figures is included. This is
display reconstruction, not regeneration of the underlying calculations.
Intervals in Supplemental Fig. S2 are stored results; this version does not
refit its hierarchical regression or repeat its bootstrap. The three-material
crosscheck summary reports a separate calculation and does not contain the
matrices required to rerun that calculation. The historical fit and counter
results must not be interpreted as an acceleration measurement.

Other main-text and supplemental tables are not all covered by this version.
The solver, full density matrices, wave functions, large simulation outputs,
and potential files remain outside this repository. The analytic potential
family is identified in the material metadata, without its coefficient tables.

`DATA_DICTIONARY.md` defines the numerical inputs and units.
The numerical-clarification subdirectory includes an independent standard-
library implementation and frozen expected summaries. It recomputes the
material-pair decomposition, the recorded 10,000-draw row-resampling interval,
the shell thresholds and provisional reuse counts. This row-resampling interval
is distinct from the stored hierarchical work-regression interval in Fig. S2.
`manifest.json` records the included file checksums and scope. Precision in the
inputs is retained for reconstruction; displayed values follow the figure and
table formatting. File checksums establish input identity, not byte-identical
PDF output across font or software versions.

## Scientific interpretation

The finite-temperature distance statements concern the specified retained
window, not the full nonlinear self-consistent state. The historical ten-row
comparison uses the uncorrected bound. Numerical allowances in the distribution
records are separate from an exact mathematical inequality. The later material
comparison includes the weighted particle-number correction and preserves the
graphene abstention as missing window quantities, never as zeros.

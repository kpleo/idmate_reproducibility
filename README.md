# IDMate manuscript reproducibility

Repository name: `idmate_reproducibility`.

This repository contains the lightweight scripts and metadata needed to
reconstruct the figures and tables accompanying the manuscript
"IDMate: Finite-temperature response bounds and reference-map fallback for
self-consistent-field proposal screening." It does not distribute numerical
records or an electronic-structure implementation.

## Repository scope

The repository contains documentation, Python scripts, environment
descriptions, machine-readable mappings, and empty input/output directory
placeholders. It excludes:

- numerical input records and simulation outputs;
- source code for the electronic-structure solver;
- pseudopotential files, including PAW, POTCAR, PSP8, UPF, and HGH datasets;
- executable files, fitted parameters, scheduler records, and
  machine-specific configuration;
- proprietary reference-code inputs and outputs.

The numerical records required by the scripts are therefore marked
`not_distributed`. Running a complete build without supplying those records
must fail with an explicit list of missing files; it is not reported as a
successful reproduction.

## Display map

| Manuscript item | Script | Required derived inputs | Expected output |
|---|---|---|---|
| Fig. 1 | `scripts/generate_figures.py --only fig1_screen` | `fig1_response_trials.csv`, `fig1_shell_weights.csv` | `fig1_screen.pdf`, `fig1_screen.svg` |
| Fig. 2 | `scripts/generate_figures.py --only fig2_causal` | `fig2_interventions.csv` | `fig2_causal.pdf`, `fig2_causal.svg` |
| Fig. 3 | `scripts/generate_figures.py --only fig3_worklaw` | `fig3_worklaw_summary.json`, `fig3_local_scaling.csv` | `fig3_worklaw.pdf`, `fig3_worklaw.svg` |
| Fig. 4 in the current four-figure layout | `scripts/generate_figures.py --only fig5_anchor` | `fig5_eos.csv`, `fig5_bands.csv` | `fig5_anchor.pdf`, `fig5_anchor.svg` |
| Main-text and SI tables | `scripts/generate_tables.py` | the 18 table inputs listed in `manifest.json` | one LaTeX fragment per table |

Figure numbering follows the current manuscript filenames: the fourth main
display retains the source identifier `fig5_anchor` because earlier low-data
displays were converted to tables. `config/display_map.json` provides the
machine-readable manuscript-label mapping.

## Input data contract

Place author-supplied, derived numerical records under a directory outside
this Git checkout, or under the ignored `data/derived/` directory. Input paths
are relative to the directory passed with `--data-root`.

CSV column names, types, units, and required fields are defined in
`manifest.json` and summarized in `docs/DATA_DICTIONARY.md`. The nested JSON
record for Fig. 3 is specified by
`schemas/fig3_worklaw_summary.schema.json`. These schemas describe derived
quantities only; they do not prescribe or expose the solver that produced
them.

## Minimal environment

- Python 3.11 or newer
- Matplotlib 3.9.2

Using `venv`:

```bash
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

Alternatively:

```bash
conda env create -f environment.yml
conda activate idmate-paper-reproduction
```

## Reproduction commands

Inspect the expected inputs without treating absent data as an error:

```bash
python scripts/validate_inputs.py --data-root data/derived --allow-missing
```

Render the data-independent method schematic used in Fig. 1a:

```bash
python scripts/generate_figures.py \
  --only fig1_screen --schematic-only --output-dir build/figures
```

After the derived records have been supplied, build all four figures and all
mapped main-text and SI tables:

```bash
python scripts/build_all.py \
  --data-root data/derived \
  --output-root build
```

Build one display or one table:

```bash
python scripts/generate_figures.py \
  --data-root data/derived \
  --output-dir build/figures \
  --only fig5_anchor

python scripts/generate_tables.py \
  --data-root data/derived \
  --output-dir build/tables \
  --only tab:physical-convergence
```

Verify that the repository contains no local paths, credentials, binary files,
large files, restricted numerical formats, or accidentally included numerical
records:

```bash
python scripts/verify_repository.py
```

## Interpretation of build status

- `PASS`: every required input exists, conforms to the declared columns and
  types, and the requested output was written.
- `INCOMPLETE`: one or more non-distributed input files are absent. This is
  the expected state of the repository as distributed.
- `ERROR`: an input is malformed or an output could not be generated.

See `docs/REPRODUCIBILITY_SCOPE.md` for the boundary between numerical
reproduction and regeneration from an electronic-structure calculation.

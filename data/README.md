# Numerical input placement

The repository does not distribute numerical records. Author-supplied derived
records may be placed in `data/derived/` for local use; the directory contents
are ignored by Git.

The expected filenames and columns are defined in `../manifest.json`. The same
files may instead remain outside the checkout and be selected with a relative
location such as `--data-root ../derived-records`. Paths recorded in scripts
and metadata remain relative to that root.

Only derived quantities used directly in manuscript figures and tables belong
in this interface. Wave functions, charge-density grids, pseudopotentials,
simulation checkpoints, executable files, fitted records, and complete
solver outputs are outside its scope.

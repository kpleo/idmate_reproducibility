# Reproducibility scope

The [2026-09-07 release](../releases/2026-09-07/README.md) includes its own
numerical inputs, coverage statement and reconstruction commands. The sections
below describe the earlier top-level interface, whose numerical inputs remain
separate. They do not limit or extend the versioned release's stated coverage.

## What this repository can reproduce

When the declared derived numerical inputs are supplied, the repository can:

1. validate filenames, headers, basic types, finite numerical values, and the
   nested Fig. 3 summary structure;
2. reconstruct the four main-text displays;
3. reconstruct LaTeX fragments for every table currently labelled in the main
   manuscript and SI;
4. verify that the repository contains no numerical records, local paths,
   credentials, restricted file types, binary files, or oversized artifacts.

## What this repository cannot reproduce by itself

The repository does not regenerate electronic-structure trajectories from
atomic structures. That stage requires an electronic-structure implementation,
licensed or separately distributed pseudopotentials, and the numerical
settings described in the manuscript. Those components are not distributed
here.

Accordingly, two distinct claims must not be conflated:

- **display reconstruction**: deterministic transformation of supplied derived
  records into manuscript figures and tables;
- **calculation regeneration**: repetition of the underlying self-consistent-
  field calculations and response experiments.

This repository supports the first claim and documents the interface needed
for a future calculation-regeneration stage.

## Data release boundary

The intended future data deposit is a separate, citable archive of derived
records matching the schemas in this repository. Until that deposit exists,
the manifest records `data_distribution: not_distributed`, and a full build is
expected to stop with an `INCOMPLETE` message.

No checksum is supplied for unavailable data. Checksums can be added when a
data deposit is available and the deposited files can be verified.

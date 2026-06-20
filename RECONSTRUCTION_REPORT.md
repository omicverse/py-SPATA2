# SPATA2 Python Reconstruction Report

Source reference: <https://themilolab.github.io/SPATA2/index.html>

SPATA2 v3 is a broad R toolkit for spatial expression analysis. The pkgdown
reference lists object initiation, platform readers, molecular variables,
image handling, spatial segmentation, spatial measures, spatial annotations,
spatial trajectories, plotting layers, and spatial gradient screening.

## Implemented first slice

The first Python slice focuses on stable infrastructure that maps cleanly to
AnnData and can later be used in `omicverse.space`:

- coordinate extraction from `adata.obsm["spatial"]` or `adata.obs`;
- expression and metadata variable extraction;
- joined coordinate-variable tables;
- convex tissue outline generation;
- spatial outlier detection and filtering;
- explicit pixel/physical-unit distance conversion.

## Namespace audit

The SPATA2 3.1.4 `NAMESPACE` currently exports 751 symbols. This repository now
ships a strict symbol-level audit:

- `NAMESPACE_PARITY.md`
- `references/spata2_namespace_parity.csv`

Current strict R-name coverage is 16 / 751 exports (2.1%). This is not a full
port yet. It is an auditable starting point for a rewrite.

## Benchmark

The benchmark fixture contains 6,402 observations x 100 genes and two isolated
artifact-like coordinates.

| Operation | Seconds |
|---|---:|
| `getCoordsDf` | 0.000418 |
| `getFeatureDf` | 0.001804 |
| `identifyTissueOutline` | 0.001991 |
| `identifySpatialOutliers` | 0.010231 |
| `removeSpatialOutliers` | 0.005265 |

Benchmark outputs are stored in `benchmarks/results/`.

## Deferred

The following areas are intentionally not reconstructed in this first pass:

- R S4 object model compatibility;
- Shiny/interactive workflows;
- image containers and pixel-content detection;
- spatial annotations and trajectories;
- spatial gradient screening statistics;
- ggplot-style plotting layers.

Those areas should be reconstructed only after API boundaries are clear and
after fixture data or parity outputs are available.

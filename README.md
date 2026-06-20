# py-SPATA2

AnnData-native Python reconstruction of the SPATA2 spatial transcriptomics
toolkit.

This is a rewrite project, not an `rpy2` bridge. The goal is namespace-level and
behavior-level parity with SPATA2 while using AnnData as the Python object model.

Current status is **early reconstruction**. The repository now tracks strict
SPATA2 3.1.4 NAMESPACE parity so incomplete areas are visible instead of hidden.

## First API slice

```python
import spata2py as spata

adata = spata.create_spata2_fixture()

coords = spata.get_coords_df(adata)
table = spata.join_with_variables(adata, ["region", "GeneA", "GeneB"])

outliers = spata.identify_spatial_outliers(adata, min_neighbors=2)
filtered = spata.remove_spatial_outliers(adata)

outline = spata.identify_tissue_outline(adata)
```

SPATA2-compatible R-style aliases are also exposed for the implemented subset:

```python
coords = spata.getCoordsDf(adata)
outline = spata.identifyTissueOutline(adata)
outliers = spata.identifySpatialOutliers(adata, radius=1.6, min_neighbors=2)
filtered = spata.removeSpatialOutliers(adata)
```

The current implemented slice covers the SPATA2 object foundation rather than
the full R toolkit:

| Python API | Purpose |
|---|---|
| `get_coords_df` | Extract AnnData spatial coordinates as a tidy table |
| `get_coords_mtr`, `get_coords_range`, `get_coords_center` | Coordinate matrix and summaries |
| `extract_variables` | Fetch observation metadata and molecular variables |
| `join_with_variables` | Combine coordinates with metadata/expression |
| `identify_tissue_outline` | Build a convex tissue hull from coordinates |
| `get_tissue_outline_df`, `contains_tissue_outline` | Stored tissue-outline accessors |
| `identify_spatial_outliers` | Flag isolated spatial observations |
| `remove_spatial_outliers`, `contains_spatial_outliers` | Filter and inspect outlier flags |
| `pixels_to_unit`, `unit_to_pixels` | Convert pixel distances with explicit scale |
| `create_spatial_data` | Lightweight `SpatialData` coordinate container |
| `is_outlier` | Robust MAD-based numeric outlier helper |

## Namespace Parity

Strict parity is generated from the SPATA2 3.1.4 `NAMESPACE` file.

Current audit:

| Metric | Value |
|---|---:|
| R exports audited | 751 |
| Matching Python R-style symbols implemented | 16 |
| Strict namespace coverage | 2.1% |

Implemented R-compatible symbols:

`SpatialData`, `as_micrometer`, `as_pixel`, `containsSpatialOutliers`,
`containsTissueOutline`, `createSpatialData`, `getCoordsCenter`, `getCoordsDf`,
`getCoordsMtr`, `getCoordsRange`, `getFeatureDf`, `getTissueOutlineDf`,
`identifySpatialOutliers`, `identifyTissueOutline`, `is_outlier`,
`removeSpatialOutliers`.

Full audit:

- [`NAMESPACE_PARITY.md`](NAMESPACE_PARITY.md)
- [`references/spata2_namespace_parity.csv`](references/spata2_namespace_parity.csv)

Regenerate it with:

```bash
python scripts/audit_namespace.py \
  --namespace /path/to/SPATA2/NAMESPACE \
  --csv-out references/spata2_namespace_parity.csv \
  --md-out NAMESPACE_PARITY.md
```

## Benchmark

Synthetic spatial fixture: 6,402 observations x 100 genes, including two
artifact-like isolated spots.

| Operation | Seconds | Result shape |
|---|---:|---|
| `getCoordsDf` | 0.000418 | `(6402, 3)` |
| `getFeatureDf` | 0.001804 | `(6402, 6)` |
| `identifyTissueOutline` | 0.001991 | `(5, 2)` |
| `identifySpatialOutliers` | 0.010231 | `(6402,)` |
| `removeSpatialOutliers` | 0.005265 | `(6400, 100)` |

Benchmark artifacts:

- [`benchmarks/benchmark_spata2py.py`](benchmarks/benchmark_spata2py.py)
- [`benchmarks/results/spata2py_benchmark.md`](benchmarks/results/spata2py_benchmark.md)
- [`benchmarks/results/spata2py_benchmark.json`](benchmarks/results/spata2py_benchmark.json)

Regenerate with:

```bash
python benchmarks/benchmark_spata2py.py
```

## Install

From GitHub:

```bash
pip install git+https://github.com/omicverse/py-SPATA2.git
```

From a local checkout:

```bash
pip install -e ".[dev]"
```

## Reconstruction Notes

SPATA2 is a broad R/S4 framework covering object initiation, molecular
variables, images, segmentation, spatial measures, spatial annotations, spatial
trajectories, plotting, and gradient screening. This Python reconstruction
starts with the low-level pieces that can be validated against AnnData and then
expands module by module:

- AnnData is the Python object model.
- Coordinates are read from `adata.obsm["spatial"]` by default.
- Expression variables come from `adata.X`, `adata.layers[...]`, or
  `adata.raw.X`.
- Metadata variables come from `adata.obs`.
- Tissue outline and outlier detection are deterministic numerical routines
  based on `scipy`.

## Parity Roadmap

The next reconstruction passes should follow SPATA2's actual export surface, not
similar-looking convenience functions:

1. Object model: `SPATA2`, `MolecularAssay`, `SpatialData`, image containers.
2. Readers/initiation: Visium, VisiumHD, Xenium, MERFISH, Slide-seq helpers.
3. Get/set/add/contains accessors with matching semantics.
4. Spatial preprocessing: coordinate transforms, tissue outlines, outliers.
5. Spatial annotations and trajectories.
6. Spatial gradient screening and spatial measures.
7. Plotting parity using matplotlib/plotly equivalents where appropriate.

## Testing

```bash
pytest
```

Current local gate:

```text
12 passed
```

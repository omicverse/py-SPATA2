# py-SPATA2

AnnData-native Python reconstruction of selected stable utilities from the
SPATA2 spatial transcriptomics toolkit.

This repository is intentionally scoped to Python APIs that map cleanly onto
AnnData and can be tested without an R runtime. It does not vendor SPATA2, use
`rpy2`, or require an installed R session.

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

The first slice covers the SPATA2 object foundation rather than the full R
toolkit:

| Python API | Purpose |
|---|---|
| `get_coords_df` | Extract AnnData spatial coordinates as a tidy table |
| `extract_variables` | Fetch observation metadata and molecular variables |
| `join_with_variables` | Combine coordinates with metadata/expression |
| `identify_tissue_outline` | Build a convex tissue hull from coordinates |
| `identify_spatial_outliers` | Flag isolated spatial observations |
| `remove_spatial_outliers` | Filter flagged observations |
| `pixels_to_unit`, `unit_to_pixels` | Convert pixel distances with explicit scale |

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
starts with the low-level pieces that are stable across platforms and useful in
OmicVerse:

- AnnData is the Python object model.
- Coordinates are read from `adata.obsm["spatial"]` by default.
- Expression variables come from `adata.X`, `adata.layers[...]`, or
  `adata.raw.X`.
- Metadata variables come from `adata.obs`.
- Tissue outline and outlier detection are deterministic numerical routines
  based on `scipy`.

## Parity Roadmap

The following SPATA2 areas are intentionally deferred until their Python API
shape is clear:

- platform-specific raw readers beyond what Scanpy/OmicVerse already provide;
- histology image containers and pixel-content segmentation;
- interactive spatial annotations and trajectories;
- spatial gradient screening models;
- ggplot-style visualization layers.

## Testing

```bash
pytest
```


import numpy as np
from anndata import AnnData

from spata2py import identify_spatial_outliers, remove_spatial_outliers


def test_regular_lattice_retains_tissue_and_removes_isolated_control():
    x, y = np.meshgrid(np.arange(10), np.arange(10))
    coords = np.column_stack([x.ravel(), y.ravel()]).astype(float) * 0.13
    coords += np.random.default_rng(7).normal(0, 1e-14, coords.shape)
    coords = np.vstack([coords, [50, 50]])
    adata = AnnData(np.ones((len(coords), 2)))
    adata.obsm['spatial'] = coords
    flags = identify_spatial_outliers(adata)
    assert not flags.iloc[:-1].any()
    assert flags.iloc[-1]
    assert remove_spatial_outliers(adata).n_obs == 100


def test_mad_does_not_treat_roundoff_as_spatial_isolation():
    from spata2py._core import _robust_threshold
    distances = 0.13 + np.linspace(-1e-14, 1e-14, 100)
    assert _robust_threshold(distances, 3.5) == distances.max()

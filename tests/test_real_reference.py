from pathlib import Path
import numpy as np
from anndata import AnnData
from scipy import sparse
from spata2py import identify_spatial_outliers


def test_real_dlpfc_and_controls_match_r_dbscan():
    f = np.load(Path(__file__).parents[1]/'examples/data/dlpfc_coords.npz')
    a = AnnData(sparse.csr_matrix((len(f['coords']),1)))
    a.obsm['spatial'] = f['coords']
    flags = identify_spatial_outliers(a).to_numpy()
    np.testing.assert_array_equal(flags,f['expected'].astype(bool))
    assert flags[-5:].all()

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData


def create_spata2_fixture() -> AnnData:
    """Create a small spatial AnnData fixture with one isolated observation."""

    coords = np.array(
        [
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.5, 0.5],
            [12.0, 12.0],
        ],
        dtype=float,
    )
    x = np.array(
        [
            [1.0, 0.0, 5.0],
            [2.0, 1.0, 4.0],
            [3.0, 1.0, 3.0],
            [4.0, 2.0, 2.0],
            [5.0, 3.0, 1.0],
            [8.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    obs = pd.DataFrame(
        {
            "region": ["core", "core", "edge", "edge", "core", "artifact"],
            "total_counts": x.sum(axis=1),
        },
        index=[f"spot_{idx}" for idx in range(coords.shape[0])],
    )
    var = pd.DataFrame(index=["GeneA", "GeneB", "GeneC"])
    adata = AnnData(X=x, obs=obs, var=var)
    adata.obsm["spatial"] = coords
    return adata

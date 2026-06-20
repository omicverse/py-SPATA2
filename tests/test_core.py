import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import spata2py as spata


def test_get_coords_df_reads_obsm_spatial():
    adata = spata.create_spata2_fixture()

    coords = spata.get_coords_df(adata, include_obs=["region"])

    assert list(coords.columns) == ["barcode", "x", "y", "region"]
    assert coords.loc["spot_2", "x"] == 1.0
    assert coords.loc["spot_2", "barcode"] == "spot_2"


def test_extract_variables_from_obs_and_expression():
    adata = spata.create_spata2_fixture()

    values = spata.extract_variables(adata, ["region", "GeneA", "GeneC"])

    assert values.loc["spot_0", "region"] == "core"
    assert values.loc["spot_3", "GeneA"] == 4.0
    assert values.loc["spot_4", "GeneC"] == 1.0


def test_extract_variables_supports_sparse_layers():
    adata = spata.create_spata2_fixture()
    adata.layers["scaled"] = sparse.csr_matrix(adata.X * 2.0)

    values = spata.extract_variables(adata, "GeneB", layer="scaled")

    assert values["GeneB"].tolist() == [0.0, 2.0, 2.0, 4.0, 6.0, 0.0]


def test_join_with_variables_combines_coords_and_values():
    adata = spata.create_spata2_fixture()

    table = spata.join_with_variables(adata, ["region", "GeneB"])

    assert list(table.columns) == ["barcode", "x", "y", "region", "GeneB"]
    assert table.loc["spot_1", "GeneB"] == 1.0


def test_identify_tissue_outline_writes_hull_to_uns():
    adata = spata.create_spata2_fixture()

    outline = spata.identify_tissue_outline(adata)

    assert {"x", "y"} == set(outline.columns)
    assert len(outline) >= 3
    assert "spata2py_tissue_outline" in adata.uns


def test_identify_spatial_outliers_with_radius():
    adata = spata.create_spata2_fixture()

    outliers = spata.identify_spatial_outliers(adata, radius=1.6, min_neighbors=2)

    assert outliers.loc["spot_5"]
    assert not outliers.loc["spot_0"]
    assert adata.obs["spata2py_spatial_outlier"].dtype == bool


def test_remove_spatial_outliers_returns_filtered_copy():
    adata = spata.create_spata2_fixture()
    spata.identify_spatial_outliers(adata, radius=1.6, min_neighbors=2)

    filtered = spata.remove_spatial_outliers(adata)

    assert filtered.n_obs == 5
    assert "spot_5" not in set(filtered.obs_names)
    assert adata.n_obs == 6


def test_unit_conversion_roundtrip():
    units = spata.pixels_to_unit(np.array([0.0, 10.0, 25.0]), pixels_per_unit=5.0)
    pixels = spata.unit_to_pixels(units, pixels_per_unit=5.0)

    np.testing.assert_allclose(units, [0.0, 2.0, 5.0])
    np.testing.assert_allclose(pixels, [0.0, 10.0, 25.0])


def test_missing_variable_raises_keyerror():
    adata = spata.create_spata2_fixture()

    with pytest.raises(KeyError):
        spata.extract_variables(adata, "MissingGene")

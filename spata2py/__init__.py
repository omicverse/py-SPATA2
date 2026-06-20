"""AnnData-native reconstruction of selected SPATA2 spatial utilities."""

from ._core import (
    extract_variables,
    get_coords_df,
    identify_spatial_outliers,
    identify_tissue_outline,
    join_with_variables,
    pixels_to_unit,
    remove_spatial_outliers,
    unit_to_pixels,
)
from .fixtures import create_spata2_fixture

spata_get_coords = get_coords_df
spata_identify_spatial_outliers = identify_spatial_outliers
spata_identify_tissue_outline = identify_tissue_outline

__all__ = [
    "get_coords_df",
    "extract_variables",
    "join_with_variables",
    "identify_tissue_outline",
    "identify_spatial_outliers",
    "remove_spatial_outliers",
    "pixels_to_unit",
    "unit_to_pixels",
    "spata_get_coords",
    "spata_identify_spatial_outliers",
    "spata_identify_tissue_outline",
    "create_spata2_fixture",
]

__version__ = "0.1.0"

"""
Coordinate extraction for MPAS plots.
"""

import numpy as np
from .common import (
    _rad2deg,
    _lon_wrap_to_center,
    _normalize_lon_window,
    _infer_center_lon,
    _infer_extent_from_points,
)


__all__ = ["extract_cell_latlon"]


def extract_cell_latlon(
    ds,
    *,
    da=None,
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    auto_extent=True,
    extent_margin=0.05,
    global_plot=False,
    center_lon=None,
    input_radians=True,
    lon_name_cell="lonCell",
    lat_name_cell="latCell",
    lon_name_vertex="lonVertex",
    lat_name_vertex="latVertex",
):
    """
    Extract and wrap MPAS cell and vertex lon/lat.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset <xarray.Dataset>`
        MPAS dataset containing cell and vertex coordinate variables.
    da : :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.ndarray <numpy.ndarray>`, optional
        Cell-centered data. If provided, finite values are used for extent
        inference and ``center_lon`` inference.
    lon_min, lon_max, lat_min, lat_max : :py:class:`float <float>`, optional
        User-defined plotting window. Supports cross-dateline windows, e.g.
        ``lon_min=90, lon_max=-60`` means 90E to 300E.
    auto_extent : :py:class:`bool <bool>`, default: True
        If True and no explicit extent is supplied, infer the extent from valid
        cell centers.
    extent_margin : :py:class:`float <float>`, default: 0.05
        Fractional margin added to the automatically inferred extent.
    global_plot : :py:class:`bool <bool>`, default: False
        If True, return a global extent centered on ``center_lon``.
    center_lon : :py:class:`float <float>`, optional
        Longitude wrapping center. If None, inferred automatically.
    input_radians : :py:class:`bool <bool>`, default: True
        If True, MPAS longitude and latitude variables are interpreted as
        radians and converted to degrees.
    lon_name_cell, lat_name_cell, lon_name_vertex, lat_name_vertex : :py:class:`str <str>`
        MPAS coordinate variable names.

    Returns
    -------
    dict
        Dictionary with wrapped cell and vertex coordinates, the inferred
        longitude center, and plotting extent values.
    """

    # -----------------------------
    # Read raw lon/lat
    # -----------------------------
    lon_cell_raw = np.asarray(ds[lon_name_cell].values)
    lat_cell = np.asarray(ds[lat_name_cell].values)

    lon_vertex_raw = np.asarray(ds[lon_name_vertex].values)
    lat_vertex = np.asarray(ds[lat_name_vertex].values)

    if input_radians:
        lon_cell_raw = _rad2deg(lon_cell_raw)
        lat_cell = _rad2deg(lat_cell)
        lon_vertex_raw = _rad2deg(lon_vertex_raw)
        lat_vertex = _rad2deg(lat_vertex)

    # -----------------------------
    # Optional valid-data mask
    # -----------------------------
    if da is not None:
        if hasattr(da, "values"):
            values = np.asarray(da.values)
        else:
            values = np.asarray(da)

        if values.ndim != 1:
            raise ValueError(
                f"`da` must be 1D on nCells after selection, got shape {values.shape}."
            )

        if values.size != lon_cell_raw.size:
            raise ValueError(
                f"`da` length does not match nCells: "
                f"len(da)={values.size}, nCells={lon_cell_raw.size}"
            )

        valid = np.isfinite(values)
    else:
        valid = np.isfinite(lon_cell_raw) & np.isfinite(lat_cell)

    # -----------------------------
    # Infer center longitude
    # -----------------------------
    if center_lon is None:
        center_lon = _infer_center_lon(
            lon_min=lon_min,
            lon_max=lon_max,
            lon_cell_raw=lon_cell_raw,
            valid=valid,
        )

    # -----------------------------
    # Wrap lon around center_lon
    # -----------------------------
    lon_cell = _lon_wrap_to_center(lon_cell_raw, center_lon)
    lon_vertex = _lon_wrap_to_center(lon_vertex_raw, center_lon)

    # -----------------------------
    # Decide plotting extent
    # -----------------------------
    user_has_extent = (
        lon_min is not None
        and lon_max is not None
        and lat_min is not None
        and lat_max is not None
    )

    if global_plot:
        lon_min_plot = center_lon - 180.0
        lon_max_plot = center_lon + 180.0
        lat_min_plot = -90.0
        lat_max_plot = 90.0

    elif user_has_extent:
        lon_min_plot, lon_max_plot = _normalize_lon_window(lon_min, lon_max)
        lat_min_plot = float(lat_min)
        lat_max_plot = float(lat_max)

    elif auto_extent:
        lon_min_plot, lon_max_plot, lat_min_plot, lat_max_plot = (
            _infer_extent_from_points(
                lon_cell[valid],
                lat_cell[valid],
                margin=extent_margin,
            )
        )

    else:
        lon_min_plot = None
        lon_max_plot = None
        lat_min_plot = None
        lat_max_plot = None

    return {
        "lon_cell": lon_cell,
        "lat_cell": lat_cell,
        "lon_vertex": lon_vertex,
        "lat_vertex": lat_vertex,
        "center_lon": center_lon,
        "lon_min_plot": lon_min_plot,
        "lon_max_plot": lon_max_plot,
        "lat_min_plot": lat_min_plot,
        "lat_max_plot": lat_max_plot,
    }

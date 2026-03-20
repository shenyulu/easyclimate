"""
Custom mask methods
"""

import numpy as np
import xarray as xr
from typing import Literal

__all__ = ["mask_custom_rectangular_box"]


def mask_custom_rectangular_box(
    da: xr.DataArray,
    lon1: float,
    lon2: float,
    lat1: float,
    lat2: float,
    angle: float = 0.0,
    center: Literal["center", "lowerleft"] = "center",
    use_coslat: bool = True,
    lat_dim: str = "lat",
    lon_dim: str = "lon",
):
    """
    Generate a mask for a custom geographical rectangular box with optional rotation.

    Parameters
    ----------
    da : :py:class:`xarray.DataArray <xarray:xarray.DataArray>`
        Input data array containing latitude and longitude coordinates.
        The mask is generated on the horizontal grid defined by `lat_dim`
        and `lon_dim`.
    lon1, lon2: :py:class:`float <float>`.
        Rectangular box longitude point. The applicable value should be between
        -180 :math:`^\circ` and 360 :math:`^\circ`. `lon1` and `lon2` are used
        to define the unrotated box extent, and do not strictly require the size
        relationship between them.
    lat1, lat2: :py:class:`float <float>`.
        Rectangular box latitude point. The applicable value should be between
        -90 :math:`^\circ` and 90 :math:`^\circ`. `lat1` and `lat2` are used to
        define the unrotated box extent, and do not strictly require the size
        relationship between them.
    angle: :py:class:`float <float>`, default: `0.0`.
        Counterclockwise rotation angle in degrees.
    center: :py:class:`str <str>`, default: `"center"`.
        Rotation pivot of the rectangular box.

        - `"center"`: rotate around the box center.
        - `"lowerleft"`: rotate around the lower-left corner of the box.
    use_coslat: :py:class:`bool <bool>`, default: `True`.
        Whether to apply local geographical metric correction before rotation.
        If `True`, longitude is scaled by :math:`\cos(lat0)` using the selected
        pivot latitude so that longitude and latitude have a locally comparable
        distance scale.
    lat_dim: :py:class:`str <str>`, default: `"lat"`.
        Name of the latitude coordinate in `da`.
    lon_dim: :py:class:`str <str>`, default: `"lon"`.
        Name of the longitude coordinate in `da`.

    Returns
    -------
    mask : :py:class:`xarray.DataArray <xarray:xarray.DataArray>`
        Boolean mask on the same horizontal grid as `da`, where values are
        `True` inside the rotated rectangular box and `False` outside.
    """

    # Normalize order
    x1, x2 = (lon1, lon2) if lon1 <= lon2 else (lon2, lon1)
    y1, y2 = (lat1, lat2) if lat1 <= lat2 else (lat2, lat1)

    lon = da[lon_dim]
    lat = da[lat_dim]

    # Construct 2D grid
    Lon, Lat = xr.broadcast(lon, lat)
    Lon = Lon.transpose(lat_dim, lon_dim)
    Lat = Lat.transpose(lat_dim, lon_dim)

    # Rotation center
    if center.lower() in ("center", "c"):
        px = (x1 + x2) / 2
        py = (y1 + y2) / 2
    elif center.lower() in ("lowerleft", "ll"):
        px = x1
        py = y1
    else:
        raise ValueError('center must be "center" or "lowerleft"')

    theta = np.deg2rad(angle)

    # Latitude-longitude scale correction
    if use_coslat:
        lat0 = py
        s = np.cos(np.deg2rad(lat0))
        s = max(s, 1e-8)
    else:
        s = 1.0

    # ===== Key step: reverse rotation of grid points =====
    X = Lon.values * s
    Y = Lat.values
    PX = px * s
    PY = py

    # Reverse rotation (-theta)
    Xr = np.cos(theta) * (X - PX) + np.sin(theta) * (Y - PY) + PX
    Yr = -np.sin(theta) * (X - PX) + np.cos(theta) * (Y - PY) + PY

    # Restore longitude scale
    Xr = Xr / s

    # Determine if inside axis-aligned rectangle
    mask_np = (Xr >= x1) & (Xr <= x2) & (Yr >= y1) & (Yr <= y2)

    mask = xr.DataArray(
        mask_np,
        coords=Lon.coords,
        dims=Lon.dims,
        name="mask",
    )

    return mask

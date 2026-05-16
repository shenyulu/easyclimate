"""
Divergence and Vorticity
"""

import xarray as xr
import numpy as np

from typing import Literal

from ..core.utility import (
    compare_multi_dataarray_coordinate,
)
from ..backend import (
    dvrfidf_ncl,
    ddvfidf_ncl,
    ddvfidf_rs,
    dvrfidf_rs,
    ddvfidf_batch_rs,
    dvrfidf_batch_rs,
)

__all__ = [
    "calc_divergence_rs",
    "calc_vorticity_rs",
    "calc_divergence_ncl",
    "calc_vorticity_ncl",
    "calc_divergence",
    "calc_vorticity",
]


def calc_divergence_rs(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    R: float = 6.37122e6,
    cyclic_boundary_setting: Literal["nan", "cyclic", "cyclic+diff", "diff"] = "nan",
    method: Literal["rust_batch", "rust_raw"] = "rust_batch",
) -> xr.DataArray:
    """
    Calculate the horizontal divergence term using the Rust backend.

    This function wraps the Rust finite-difference implementation for spherical
    wind fields and supports both batched and raw execution modes.

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name.
    R: :py:class:`float <float>`, default: `6371220.0`.
        Radius of the Earth in meters.
    cyclic_boundary_setting: {"nan", "cyclic", "cyclic+diff", "diff"}, default: `nan`.
        A scalar integer equal to the boundary condition option:

        - ``nan``: Boundary points are set to the missing value.
        - ``cyclic``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic point.) The upper and lower boundaries will be set to missing.
        - ``cyclic+diff``: Boundary points are estimated using one-sided difference schemes normal to the boundary.
        - ``diff``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic points.) The upper and lower boundaries are estimated using a one-sided difference scheme normal to the boundary.
    method: {`rust_batch`, `rust_raw`}, default: `rust_batch`.
        Rust execution mode. ``rust_batch`` uses the batched backend function,
        while ``rust_raw`` loops over non-core dimensions and calls the raw backend.

    Returns
    -------
    The horizontal divergence term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    xmsg = np.nan
    compare_multi_dataarray_coordinate([u_data, v_data])

    if lon_dim not in u_data.dims or lat_dim not in u_data.dims:
        raise ValueError(
            f"u_data and v_data must have dimensions {lon_dim} and {lat_dim}."
        )

    match cyclic_boundary_setting:
        case "nan":
            iopt = 0
        case "cyclic":
            iopt = 1
        case "cyclic+diff":
            iopt = 2
        case "diff":
            iopt = 3

    transpose_order_u = [d for d in u_data.dims if d not in [lat_dim, lon_dim]] + [
        lat_dim,
        lon_dim,
    ]
    u_trans = u_data.transpose(*transpose_order_u)
    v_trans = v_data.transpose(*transpose_order_u)

    def _core_batch(u_vals, v_vals, lon_vals, lat_vals, xmsg, R):
        orig_lat_vals = lat_vals.copy()
        flip_lat = orig_lat_vals[0] > orig_lat_vals[-1]
        if flip_lat:
            u_vals = np.flip(u_vals, axis=-2)
            v_vals = np.flip(v_vals, axis=-2)
            lat_vals = np.flip(lat_vals)

        glat_vals = np.asarray(lat_vals, dtype=np.float64, order="C")
        glon_vals = np.asarray(lon_vals, dtype=np.float64, order="C")

        # Ensure C-contiguous + float64 to avoid implicit copy/stride pitfalls in PyTorch 3/NumPy.
        u_batch = np.asarray(u_vals, dtype=np.float64, order="C")
        v_batch = np.asarray(v_vals, dtype=np.float64, order="C")

        xmsg_rs = float(xmsg) if not np.isnan(xmsg) else np.nan

        dv_out, ier = ddvfidf_batch_rs(
            u_batch, v_batch, glat_vals, glon_vals, int(iopt), xmsg_rs, float(R)
        )
        if ier != 0:
            raise ValueError(f"Rust backend error in ddvfidf_batch: ier={ier}")

        dv_out = np.asarray(dv_out)

        if flip_lat:
            dv_out = np.flip(dv_out, axis=-2)

        return dv_out

    def _core(u_vals, v_vals, lon_vals, lat_vals, xmsg, R):
        # u_vals/v_vals: (..., nlat, mlon)
        orig_lat_vals = lat_vals.copy()
        flip_lat = orig_lat_vals[0] > orig_lat_vals[-1]
        if flip_lat:
            u_vals = np.flip(u_vals, axis=-2)
            v_vals = np.flip(v_vals, axis=-2)
            lat_vals = np.flip(lat_vals)

        # Rust side expectations (nlat, mlon), no swapaxes required here
        glat_vals = lat_vals.astype(np.float64)
        glon_vals = lon_vals.astype(np.float64)

        # dtype/contiguous (avoids implicit pyo3/numpy copy overhead and strange strides)
        u_batch = np.asarray(u_vals, dtype=np.float64, order="C")
        v_batch = np.asarray(v_vals, dtype=np.float64, order="C")

        batch_shape = u_batch.shape[:-2]
        nlat_out = u_batch.shape[-2]
        mlon_out = u_batch.shape[-1]
        dv_out = np.empty((*batch_shape, nlat_out, mlon_out), dtype=np.float64)

        # Rust now supports `xmsg=np.nan` => missing measurements are replaced with NaN; therefore, sentinels are no longer used.
        xmsg_rs = float(xmsg) if not np.isnan(xmsg) else np.nan

        for idx in np.ndindex(*batch_shape):
            u_slice = u_batch[idx]
            v_slice = v_batch[idx]

            dv_slice, ier = ddvfidf_rs(
                u_slice, v_slice, glat_vals, glon_vals, int(iopt), xmsg_rs, R
            )
            if ier != 0:
                raise ValueError(f"Rust backend error in ddvfidf: ier={ier}")

            dv_out[idx] = np.asarray(dv_slice)

        if flip_lat:
            dv_out = np.flip(dv_out, axis=-2)

        return dv_out

    if method == "rust_batch":
        div = xr.apply_ufunc(
            _core_batch,
            u_trans,
            v_trans,
            u_trans[lon_dim],
            u_trans[lat_dim],
            xmsg,
            R,
            input_core_dims=[
                (lat_dim, lon_dim),
                (lat_dim, lon_dim),
                (lon_dim,),
                (lat_dim,),
                [],
                [],
            ],
            output_core_dims=[(lat_dim, lon_dim)],
            output_dtypes=[np.float64],
            keep_attrs=True,
            dask="allowed",
            vectorize=False,
        )
        div = div.transpose(*u_data.dims)

    elif method == "rust_raw":
        div = xr.apply_ufunc(
            _core,
            u_trans,
            v_trans,
            u_trans[lon_dim],
            u_trans[lat_dim],
            xmsg,
            R,
            input_core_dims=[
                (lat_dim, lon_dim),
                (lat_dim, lon_dim),
                (lon_dim,),
                (lat_dim,),
                [],
                [],
            ],
            output_core_dims=[(lat_dim, lon_dim)],
            output_dtypes=[np.float64],
            keep_attrs=True,
            dask="allowed",
            vectorize=False,
        )
        div = div.transpose(*u_data.dims)

    div.name = "divergence"
    div.attrs["long_name"] = "divergence"
    div.attrs["units"] = "s^-1"
    return div


def calc_vorticity_rs(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    R: float = 6.37122e6,
    cyclic_boundary_setting: Literal["nan", "cyclic", "cyclic+diff", "diff"] = "nan",
    method: Literal["rust_batch", "rust_raw"] = "rust_batch",
) -> xr.DataArray:
    """
    Calculate the horizontal relative vorticity term using the Rust backend.

    This function wraps the Rust finite-difference implementation for spherical
    wind fields and supports both batched and raw execution modes.

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name.
    R: :py:class:`float <float>`, default: `6371220.0`.
        Radius of the Earth in meters.
    cyclic_boundary_setting: {"nan", "cyclic", "cyclic+diff", "diff"}, default: `nan`.
        A scalar integer equal to the boundary condition option:

        - ``nan``: Boundary points are set to the missing value.
        - ``cyclic``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic point.) The upper and lower boundaries will be set to missing.
        - ``cyclic+diff``: Boundary points are estimated using one-sided difference schemes normal to the boundary.
        - ``diff``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic points.) The upper and lower boundaries are estimated using a one-sided difference scheme normal to the boundary.
    method: {`rust_batch`, `rust_raw`}, default: `rust_batch`.
        Rust execution mode. ``rust_batch`` uses the batched backend function,
        while ``rust_raw`` loops over non-core dimensions and calls the raw backend.

    Returns
    -------
    The horizontal relative vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    xmsg = np.nan
    compare_multi_dataarray_coordinate([u_data, v_data])

    match cyclic_boundary_setting:
        case "nan":
            iopt = 0
        case "cyclic":
            iopt = 1
        case "cyclic+diff":
            iopt = 2
        case "diff":
            iopt = 3

    transpose_order_u = [d for d in u_data.dims if d not in [lat_dim, lon_dim]] + [
        lat_dim,
        lon_dim,
    ]
    u_trans = u_data.transpose(*transpose_order_u)
    v_trans = v_data.transpose(*transpose_order_u)

    def _core_batch(u_vals, v_vals, lon_vals, lat_vals, xmsg, R):
        orig_lat_vals = lat_vals.copy()
        flip_lat = orig_lat_vals[0] > orig_lat_vals[-1]
        if flip_lat:
            u_vals = np.flip(u_vals, axis=-2)
            v_vals = np.flip(v_vals, axis=-2)
            lat_vals = np.flip(lat_vals)

        glat_vals = np.asarray(lat_vals, dtype=np.float64, order="C")
        glon_vals = np.asarray(lon_vals, dtype=np.float64, order="C")

        u_batch = np.asarray(u_vals, dtype=np.float64, order="C")
        v_batch = np.asarray(v_vals, dtype=np.float64, order="C")

        xmsg_rs = float(xmsg) if not np.isnan(xmsg) else np.nan

        rv_out, ier = dvrfidf_batch_rs(
            u_batch, v_batch, glat_vals, glon_vals, int(iopt), xmsg_rs, float(R)
        )
        if ier != 0:
            raise ValueError(f"Rust backend error in dvrfidf_batch: ier={ier}")

        rv_out = np.asarray(rv_out)

        if flip_lat:
            rv_out = np.flip(rv_out, axis=-2)

        return rv_out

    def _core(u_vals, v_vals, lon_vals, lat_vals, xmsg, R):
        orig_lat_vals = lat_vals.copy()
        flip_lat = orig_lat_vals[0] > orig_lat_vals[-1]
        if flip_lat:
            u_vals = np.flip(u_vals, axis=-2)
            v_vals = np.flip(v_vals, axis=-2)
            lat_vals = np.flip(lat_vals)

        glat_vals = lat_vals.astype(np.float64)
        glon_vals = lon_vals.astype(np.float64)

        u_batch = np.asarray(u_vals, dtype=np.float64, order="C")
        v_batch = np.asarray(v_vals, dtype=np.float64, order="C")

        batch_shape = u_batch.shape[:-2]
        nlat_out = u_batch.shape[-2]
        mlon_out = u_batch.shape[-1]
        rv_out = np.empty((*batch_shape, nlat_out, mlon_out), dtype=np.float64)

        xmsg_rs = float(xmsg) if not np.isnan(xmsg) else np.nan

        for idx in np.ndindex(*batch_shape):
            u_slice = u_batch[idx]
            v_slice = v_batch[idx]

            rv_slice, ier = dvrfidf_rs(
                u_slice, v_slice, glat_vals, glon_vals, int(iopt), xmsg_rs, R
            )
            if ier != 0:
                raise ValueError(f"Rust backend error in dvrfidf: ier={ier}")

            rv_out[idx] = np.asarray(rv_slice)

        if flip_lat:
            rv_out = np.flip(rv_out, axis=-2)

        return rv_out

    if method == "rust_batch":
        rv = xr.apply_ufunc(
            _core_batch,
            u_trans,
            v_trans,
            u_trans[lon_dim],
            u_trans[lat_dim],
            xmsg,
            R,
            input_core_dims=[
                (lat_dim, lon_dim),
                (lat_dim, lon_dim),
                (lon_dim,),
                (lat_dim,),
                [],
                [],
            ],
            output_core_dims=[(lat_dim, lon_dim)],
            output_dtypes=[np.float64],
            keep_attrs=True,
            dask="allowed",
            vectorize=False,
        )
        rv = rv.transpose(*u_data.dims)

    elif method == "rust_raw":
        rv = xr.apply_ufunc(
            _core,
            u_trans,
            v_trans,
            u_trans[lon_dim],
            u_trans[lat_dim],
            xmsg,
            R,
            input_core_dims=[
                (lat_dim, lon_dim),
                (lat_dim, lon_dim),
                (lon_dim,),
                (lat_dim,),
                [],
                [],
            ],
            output_core_dims=[(lat_dim, lon_dim)],
            output_dtypes=[np.float64],
            keep_attrs=True,
            dask="allowed",
            vectorize=False,
        )
        rv = rv.transpose(*u_data.dims)

    rv.name = "relative_vorticity"
    rv.attrs["long_name"] = "relative_vorticity"
    rv.attrs["units"] = "s^-1"
    return rv


def calc_divergence_ncl(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    cyclic_boundary_setting: Literal["nan", "cyclic", "cyclic+diff", "diff"] = "nan",
) -> xr.DataArray:
    """
    Calculate the horizontal divergence term using the NCL-compatible backend.

    This function wraps the NCL-style finite-difference implementation and keeps
    the input dimension order unchanged in the returned result.

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name.
    cyclic_boundary_setting: {"nan", "cyclic", "cyclic+diff", "diff"}, default: `nan`.
        A scalar integer equal to the boundary condition option:

        - ``nan``: Boundary points are set to the missing value.
        - ``cyclic``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic point.) The upper and lower boundaries will be set to missing.
        - ``cyclic+diff``: Boundary points are estimated using one-sided difference schemes normal to the boundary.
        - ``diff``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic points.) The upper and lower boundaries are estimated using a one-sided difference scheme normal to the boundary.

    Returns
    -------
    The horizontal divergence term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    compare_multi_dataarray_coordinate([u_data, v_data])

    if lon_dim not in u_data.dims or lat_dim not in u_data.dims:
        raise ValueError(
            f"u_data and v_data must have dimensions {lon_dim} and {lat_dim}."
        )

    missing = None
    xmsg = missing if missing is not None else u_data.attrs.get("_FillValue", np.nan)

    match cyclic_boundary_setting:
        case "nan":
            iopt = 0
        case "cyclic":
            iopt = 1
        case "cyclic+diff":
            iopt = 2
        case "diff":
            iopt = 3

    # Transpose to ensure lat_dim is second-to-last, lon_dim is last
    transpose_order_u = [d for d in u_data.dims if d not in [lat_dim, lon_dim]] + [
        lat_dim,
        lon_dim,
    ]
    u_trans = u_data.transpose(*transpose_order_u)
    v_trans = v_data.transpose(*transpose_order_u)  # Assume same dims

    def _core(u_vals, v_vals, lon_vals, lat_vals, xmsg):
        # Core slices: u_vals.shape = (..., nlat, mlon)
        orig_lat_vals = lat_vals.copy()
        flip_lat = orig_lat_vals[0] > orig_lat_vals[-1]
        if flip_lat:
            u_vals = np.flip(u_vals, axis=-2)
            v_vals = np.flip(v_vals, axis=-2)
            lat_vals = np.flip(lat_vals)

        # Swap to (..., mlon, nlat) for Fortran
        u_vals_t = np.swapaxes(u_vals, -2, -1)
        v_vals_t = np.swapaxes(v_vals, -2, -1)
        glat_vals = lat_vals.astype(np.float64)
        glon_vals = lon_vals.astype(np.float64)

        # Handle missing values
        if np.isnan(xmsg):
            sentinel = 1e20
            u_fort_batch = np.nan_to_num(u_vals_t, nan=sentinel).astype(np.float64)
            v_fort_batch = np.nan_to_num(v_vals_t, nan=sentinel).astype(np.float64)
            xmsg_fort = float(sentinel)
        else:
            u_fort_batch = u_vals_t.astype(np.float64)
            v_fort_batch = v_vals_t.astype(np.float64)
            xmsg_fort = float(xmsg)

        # Batch processing
        batch_shape = u_fort_batch.shape[:-2]
        nlat_out = u_fort_batch.shape[-1]
        mlon_out = u_fort_batch.shape[-2]
        dv_out = np.empty((*batch_shape, nlat_out, mlon_out), dtype=np.float64)

        for idx in np.ndindex(*batch_shape):
            u_slice = u_fort_batch[idx]
            v_slice = v_fort_batch[idx]

            dv_slice_fort, ier = ddvfidf_ncl(
                u_slice, v_slice, glat_vals, glon_vals, xmsg_fort, iopt
            )
            if ier != 0:
                raise ValueError(f"easyclimate-backend error in ddvfidf: ier={ier}")

            # Restore missing values if using sentinel
            if np.isnan(xmsg):
                dv_slice_fort = np.where(
                    dv_slice_fort == xmsg_fort, np.nan, dv_slice_fort
                )

            # Swap back to (nlat, mlon)
            dv_slice = np.swapaxes(dv_slice_fort, -2, -1)
            dv_out[idx] = dv_slice

        # Flip back if original was decreasing
        if flip_lat:
            dv_out = np.flip(dv_out, axis=-2)

        return dv_out

    div = xr.apply_ufunc(
        _core,
        u_trans,
        v_trans,
        u_trans[lon_dim],
        u_trans[lat_dim],
        xmsg,  # Pass as scalar, broadcasted
        input_core_dims=[
            (lat_dim, lon_dim),
            (lat_dim, lon_dim),
            (lon_dim,),
            (lat_dim,),
            [],
        ],
        output_core_dims=[(lat_dim, lon_dim)],
        output_dtypes=[np.float64],
        keep_attrs=True,
        dask="allowed",  # Allow dask on non-core dimensions
        vectorize=False,
    )

    # Transpose back to original dimension order
    div = div.transpose(*u_data.dims)

    div.name = "divergence"
    div.attrs["long_name"] = "divergence"
    div.attrs["units"] = "s^-1"
    return div


def calc_vorticity_ncl(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    cyclic_boundary_setting: Literal["nan", "cyclic", "cyclic+diff", "diff"] = "nan",
) -> xr.DataArray:
    """
    Calculate the horizontal relative vorticity term using the NCL-compatible backend.

    This function wraps the NCL-style finite-difference implementation and keeps
    the input dimension order unchanged in the returned result.

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name.
    cyclic_boundary_setting: {"nan", "cyclic", "cyclic+diff", "diff"}, default: `nan`.
        A scalar integer equal to the boundary condition option:

        - ``nan``: Boundary points are set to the missing value.
        - ``cyclic``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic point.) The upper and lower boundaries will be set to missing.
        - ``cyclic+diff``: Boundary points are estimated using one-sided difference schemes normal to the boundary.
        - ``diff``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic points.) The upper and lower boundaries are estimated using a one-sided difference scheme normal to the boundary.

    Returns
    -------
    The horizontal relative vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).
    """
    compare_multi_dataarray_coordinate([u_data, v_data])

    missing = None
    xmsg = missing if missing is not None else u_data.attrs.get("_FillValue", np.nan)

    match cyclic_boundary_setting:
        case "nan":
            iopt = 0
        case "cyclic":
            iopt = 1
        case "cyclic+diff":
            iopt = 2
        case "diff":
            iopt = 3

    # Transpose to ensure lat_dim is second-to-last, lon_dim is last
    transpose_order_u = [d for d in u_data.dims if d not in [lat_dim, lon_dim]] + [
        lat_dim,
        lon_dim,
    ]
    u_trans = u_data.transpose(*transpose_order_u)
    v_trans = v_data.transpose(*transpose_order_u)  # Assume same dims

    def _core(u_vals, v_vals, lon_vals, lat_vals, xmsg):
        # Core slices: u_vals.shape = (..., nlat, mlon)
        orig_lat_vals = lat_vals.copy()
        flip_lat = orig_lat_vals[0] > orig_lat_vals[-1]
        if flip_lat:
            u_vals = np.flip(u_vals, axis=-2)
            v_vals = np.flip(v_vals, axis=-2)
            lat_vals = np.flip(lat_vals)

        # Swap to (..., mlon, nlat) for Fortran
        u_vals_t = np.swapaxes(u_vals, -2, -1)
        v_vals_t = np.swapaxes(v_vals, -2, -1)
        glat_vals = lat_vals.astype(np.float64)
        glon_vals = lon_vals.astype(np.float64)

        # Handle missing values
        if np.isnan(xmsg):
            sentinel = 1e20
            u_fort_batch = np.nan_to_num(u_vals_t, nan=sentinel).astype(np.float64)
            v_fort_batch = np.nan_to_num(v_vals_t, nan=sentinel).astype(np.float64)
            xmsg_fort = float(sentinel)
        else:
            u_fort_batch = u_vals_t.astype(np.float64)
            v_fort_batch = v_vals_t.astype(np.float64)
            xmsg_fort = float(xmsg)

        # Batch processing
        batch_shape = u_fort_batch.shape[:-2]
        nlat_out = u_fort_batch.shape[-1]
        mlon_out = u_fort_batch.shape[-2]
        rv_out = np.empty((*batch_shape, nlat_out, mlon_out), dtype=np.float64)

        for idx in np.ndindex(*batch_shape):
            u_slice = u_fort_batch[idx]
            v_slice = v_fort_batch[idx]

            rv_slice_fort, ier = dvrfidf_ncl(
                u_slice, v_slice, glat_vals, glon_vals, xmsg_fort, iopt
            )
            if ier != 0:
                raise ValueError(f"Fortran error in dvrfidf: ier={ier}")

            # Restore missing values if using sentinel
            if np.isnan(xmsg):
                rv_slice_fort = np.where(
                    rv_slice_fort == xmsg_fort, np.nan, rv_slice_fort
                )

            # Swap back to (nlat, mlon)
            rv_slice = np.swapaxes(rv_slice_fort, -2, -1)
            rv_out[idx] = rv_slice

        # Flip back if original was decreasing
        if flip_lat:
            rv_out = np.flip(rv_out, axis=-2)

        return rv_out

    rv = xr.apply_ufunc(
        _core,
        u_trans,
        v_trans,
        u_trans[lon_dim],
        u_trans[lat_dim],
        xmsg,  # Pass as scalar, broadcasted
        input_core_dims=[
            (lat_dim, lon_dim),
            (lat_dim, lon_dim),
            (lon_dim,),
            (lat_dim,),
            [],
        ],
        output_core_dims=[(lat_dim, lon_dim)],
        output_dtypes=[np.float64],
        keep_attrs=True,
        dask="allowed",  # Allow dask on non-core dimensions
        vectorize=False,
    )

    # Transpose back to original dimension order
    rv = rv.transpose(*u_data.dims)

    rv.name = "relative_vorticity"
    rv.attrs["long_name"] = "relative_vorticity"
    rv.attrs["units"] = "s^-1"
    return rv


def calc_divergence(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    R: float = 6.37122e6,
    spherical_coord=True,
) -> xr.DataArray:
    """
    Calculate the horizontal divergence term.

    rectangular coordinates

    .. math::
        \\mathrm{D} = \\frac{\\partial u}{\\partial x} + \\frac{\\partial v}{\\partial y}

    Spherical coordinates

    .. math::
        \\mathrm{D} = \\frac{\\partial u}{\\partial x} + \\frac{\\partial v}{\\partial y} - \\frac{v}{R} \\tan \\varphi

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    R: :py:class:`float <float>`, default: `6371200.0`.
        Radius of the Earth.
    spherical_coord: :py:class:`bool<bool>`, default: `True`.
        Whether or not to compute the horizontal Laplace term in spherical coordinates.

    Returns
    -------
    The horizontal divergence term. (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2dv_cfd.shtml
        - Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1. p113-114

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    from .diff import calc_dx_gradient, calc_dy_gradient
    from .utility import transfer_deg2rad

    compare_multi_dataarray_coordinate([u_data, v_data])

    lat_array = u_data[lat_dim]
    dudx = calc_dx_gradient(u_data, lon_dim=lon_dim, lat_dim=lat_dim, R=R)
    dvdy = calc_dy_gradient(v_data, lat_dim=lat_dim, R=R)

    if spherical_coord == True:
        term3 = v_data / R * np.tan(transfer_deg2rad(lat_array))
        div = dudx + dvdy - term3
    elif spherical_coord == False:
        div = dudx + dvdy

    div.name = "divergence"
    div.attrs["long_name"] = "divergence"
    div.attrs["units"] = "s^-1"
    return div


def calc_vorticity(
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    R: float = 6.37122e6,
    spherical_coord: bool = True,
) -> xr.DataArray:
    """
    Calculate the horizontal relative vorticity term.

    rectangular coordinates

    .. math::
        \\zeta = \\frac{\\partial v}{\\partial x} - \\frac{\\partial u}{\\partial y}

    Spherical coordinates

    .. math::
        \\zeta = \\frac{\\partial v}{\\partial x} - \\frac{\\partial u}{\\partial y} + \\frac{u}{R} \\tan \\varphi

    Parameters
    ----------
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    R: :py:class:`float <float>`, default: `6370000`.
        Radius of the Earth.
    spherical_coord: :py:class:`bool<bool>`, default: `True`.
        Whether or not to compute the horizontal Laplace term in spherical coordinates.

    Returns
    -------
    The horizontal relative vorticity term. (:py:class:`xarray.DataArray<xarray.DataArray>`).

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2vr_cfd.shtml
        - Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1. p113-114

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    from .diff import calc_dx_gradient, calc_dy_gradient
    from .utility import transfer_deg2rad

    compare_multi_dataarray_coordinate([u_data, v_data])

    dvdx = calc_dx_gradient(v_data, lon_dim=lon_dim, lat_dim=lat_dim, R=R)
    dudy = calc_dy_gradient(u_data, lat_dim=lat_dim, R=R)

    if spherical_coord == True:
        term3 = u_data / R * np.tan(transfer_deg2rad(u_data[lat_dim]))
        vor = dvdx - dudy + term3
    elif spherical_coord == False:
        vor = dvdx - dudy

    vor.name = "relative_vorticity"
    vor.attrs["long_name"] = "relative_vorticity"
    vor.attrs["units"] = "s^-1"
    return vor

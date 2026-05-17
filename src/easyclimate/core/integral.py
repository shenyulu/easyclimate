"""
Vertical integration using beta factors.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import xarray as xr

from .units import transfer_data_multiple_units, transfer_units_coeff

__all__ = [
    "calc_top2surface_integral",
    "calc_top2surface_average",
    "calc_top2surface_integral_rs",
    "calc_top2surface_average_rs",
    "calc_delta_pressure",
    "calc_p_integral",
]


def _prepare_top2surface_inputs(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
) -> tuple[xr.DataArray, xr.DataArray, np.ndarray]:
    vertical_dim_base = transfer_units_coeff(vertical_dim_units, "Pa")
    data_input = data_input.assign_coords(
        {vertical_dim: data_input[vertical_dim] * vertical_dim_base}
    )

    surface_pressure_data = transfer_data_multiple_units(
        surface_pressure_data, surface_pressure_data_units, "Pa"
    )

    level = np.asarray(data_input.coords[vertical_dim].values, dtype=np.float64)
    if not np.all(np.diff(level) < 0):
        raise ValueError(
            "Vertical levels must be in decreasing order after conversion to Pa."
        )

    return data_input, surface_pressure_data, level


def _get_vibeta_core(level: np.ndarray):
    from ..backend import dvibeta_ncl

    ptop = np.float64(level[-1])
    nlev = np.int32(len(level))
    xmsg = np.float64(1e30)
    linlog = np.int32(2)
    level64 = level.astype(np.float64)

    def vibeta_core(x_profile, psfc):
        if np.isnan(psfc) or psfc <= ptop:
            return np.nan

        x_profile = np.asarray(x_profile, dtype=np.float64)
        x_profile = np.where(np.isnan(x_profile), xmsg, x_profile)

        xsfc = x_profile[0] if abs(x_profile[0] - xmsg) <= 1e10 else 0.0
        pbot = np.float64(psfc)
        plvcrt = ptop

        vint, ier = dvibeta_ncl(
            level64,
            x_profile,
            xmsg,
            linlog,
            np.float64(psfc),
            np.float64(xsfc),
            pbot,
            ptop,
            plvcrt,
            nlev,
        )

        if ier == 0:
            return vint
        return np.nan

    return vibeta_core


def _apply_vibeta_integral(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    vibeta_core,
) -> xr.DataArray:
    return xr.apply_ufunc(
        vibeta_core,
        data_input,
        surface_pressure_data,
        input_core_dims=[(vertical_dim,), ()],
        output_core_dims=[()],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.float64],
        join="outer",
        dataset_fill_value=np.nan,
        keep_attrs=False,
    )


def _build_top2surface_output(
    integral: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    data_input: xr.DataArray,
    units: str,
    integral_type: str,
    normalize: bool,
    mass_weighted: bool,
    method: str | None = None,
) -> xr.DataArray:
    data_name = data_input.name or "data"
    if method is None:
        long_name = f"Vertical integral of {data_name} from top to surface"
    else:
        long_name = (
            f"Vertical integral of {data_name} from top to surface using {method}"
        )

    output = xr.DataArray(
        integral,
        dims=surface_pressure_data.dims,
        coords=surface_pressure_data.coords,
        attrs={
            "long_name": long_name,
            "units": units,
            "integral_type": integral_type,
            "normalize": str(normalize),
            "mass_weighted": str(mass_weighted),
            "method": method or "",
        },
        name=data_input.name,
    )
    return output


def _get_mass_weighted_units(input_units: str) -> str:
    if input_units in {"m s-1", "m s**-1", "m s^-1"}:
        return "kg m-1 s-1"
    return f"({input_units} Pa)/(m s-2)" if input_units else "Pa/(m s-2)"


def _import_rust_vibeta_functions():
    try:
        from ..backend import dvibeta_rs, dvibeta_batch, dvibeta_batch_sum_norm
    except ImportError as exc:
        raise ImportError(
            "Rust vibeta backend is unavailable. Please install `easyclimate_rust`."
        ) from exc

    return dvibeta_rs, dvibeta_batch, dvibeta_batch_sum_norm


def _get_rust_profile_block(level: np.ndarray):
    dvibeta_rs, _, _ = _import_rust_vibeta_functions()

    level_f64 = np.ascontiguousarray(level, dtype=np.float64)
    ptop = float(level_f64[-1])
    xmsg = np.float64(1e30)
    linlog = np.int32(2)

    def vibeta_block(x, psfc):
        x = np.asarray(x, dtype=np.float64)
        psfc = np.asarray(psfc, dtype=np.float64)

        out = np.full(psfc.shape, np.nan, dtype=np.float64)
        valid = np.isfinite(psfc) & (psfc > ptop)
        if not np.any(valid):
            return out

        nlev = x.shape[-1]
        x2 = x.reshape(-1, nlev)
        ps2 = psfc.reshape(-1)
        out2 = out.reshape(-1)
        valid2 = valid.reshape(-1)

        xwork = np.array(x2, copy=True, dtype=np.float64)
        bad = ~np.isfinite(xwork)
        if bad.any():
            xwork[bad] = xmsg

        xsfc_all = xwork[:, 0].copy()
        xsfc_all[np.abs(xsfc_all - xmsg) <= 1e10] = 0.0

        for i in range(xwork.shape[0]):
            if not valid2[i]:
                continue

            ps = float(ps2[i])
            vint, ier = dvibeta_rs(
                level_f64,
                np.ascontiguousarray(xwork[i]),
                xmsg,
                linlog,
                np.float64(ps),
                np.float64(xsfc_all[i]),
                np.float64(ps),
                np.float64(ptop),
                np.float64(ptop),
            )

            if ier == 0:
                out2[i] = vint

        return out

    return vibeta_block


def _calc_top2surface_rust(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    method: Literal["rust", "rust-block"],
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    _, dvibeta_batch, dvibeta_batch_sum_norm = _import_rust_vibeta_functions()

    data_input, surface_pressure_data, level = _prepare_top2surface_inputs(
        data_input=data_input,
        surface_pressure_data=surface_pressure_data,
        vertical_dim=vertical_dim,
        surface_pressure_data_units=surface_pressure_data_units,
        vertical_dim_units=vertical_dim_units,
    )

    if method == "rust":
        vibeta_block = _get_rust_profile_block(level)
        integral = xr.apply_ufunc(
            vibeta_block,
            data_input,
            surface_pressure_data,
            input_core_dims=[(vertical_dim,), ()],
            output_core_dims=[()],
            vectorize=False,
            dask="parallelized",
            output_dtypes=[np.float64],
            join="outer",
            keep_attrs=False,
        )
        layer_thickness = xr.apply_ufunc(
            vibeta_block,
            xr.ones_like(data_input, dtype=np.float64),
            surface_pressure_data,
            input_core_dims=[(vertical_dim,), ()],
            output_core_dims=[()],
            vectorize=False,
            dask="parallelized",
            output_dtypes=[np.float64],
            join="outer",
            keep_attrs=False,
        )
        return integral, layer_thickness, surface_pressure_data

    if method == "rust-block":
        level_f64 = np.ascontiguousarray(level, dtype=np.float64)
        ptop = np.float64(level_f64[-1])
        xmsg = np.float64(1e30)
        linlog = np.int32(2)

        def vibeta_block(x, psfc):
            x = np.ascontiguousarray(x, dtype=np.float64)
            psfc = np.ascontiguousarray(psfc, dtype=np.float64)
            vint, _ = dvibeta_batch(level_f64, x, psfc, xmsg, linlog, ptop)
            return vint

        def vibeta_block_sum_norm(x, psfc):
            x = np.ascontiguousarray(x, dtype=np.float64)
            psfc = np.ascontiguousarray(psfc, dtype=np.float64)
            vsum, vnorm, _ = dvibeta_batch_sum_norm(
                level_f64, x, psfc, xmsg, linlog, ptop
            )
            return vsum, vnorm

        integral = xr.apply_ufunc(
            vibeta_block,
            data_input,
            surface_pressure_data,
            input_core_dims=[(vertical_dim,), ()],
            output_core_dims=[()],
            vectorize=False,
            dask="parallelized",
            output_dtypes=[np.float64],
            join="outer",
            keep_attrs=False,
        )
        layer_integral, layer_thickness = xr.apply_ufunc(
            vibeta_block_sum_norm,
            data_input,
            surface_pressure_data,
            input_core_dims=[(vertical_dim,), ()],
            output_core_dims=[(), ()],
            vectorize=False,
            dask="parallelized",
            output_dtypes=[np.float64, np.float64],
            join="outer",
            keep_attrs=False,
        )
        return integral, layer_thickness, surface_pressure_data

    raise ValueError("The parameter `method` should be `rust` or `rust-block`.")


def _calc_top2surface_vibeta_ncl(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    data_input, surface_pressure_data, level = _prepare_top2surface_inputs(
        data_input=data_input,
        surface_pressure_data=surface_pressure_data,
        vertical_dim=vertical_dim,
        surface_pressure_data_units=surface_pressure_data_units,
        vertical_dim_units=vertical_dim_units,
    )
    vibeta_core = _get_vibeta_core(level)
    integral = _apply_vibeta_integral(
        data_input=data_input,
        surface_pressure_data=surface_pressure_data,
        vertical_dim=vertical_dim,
        vibeta_core=vibeta_core,
    )
    layer_thickness = _apply_vibeta_integral(
        data_input=xr.ones_like(data_input, dtype=np.float64),
        surface_pressure_data=surface_pressure_data,
        vertical_dim=vertical_dim,
        vibeta_core=vibeta_core,
    )
    return integral, layer_thickness, surface_pressure_data


def calc_top2surface_integral(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    normalize: bool = False,
    mass_weighted: bool = False,
    gravity: float = 9.80665,
) -> xr.DataArray:
    """
    Calculate the vertical integral in the p-coordinate system from the ground
    to the zenith along the barometric pressure direction with Fortran vibeta backends.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data to be calculated.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`
        Surface level pressure. Must be local surface pressure, not MSLP.
    vertical_dim: :py:class:`str <str>`
        Vertical coordinate dimension name.
    surface_pressure_data_units: {"hPa", "Pa", "mbar"}
        Unit of `surface_pressure_data`.
    vertical_dim_units: {"hPa", "Pa", "mbar"}
        Unit of vertical pressure coordinate.
    mass_weighted: :py:class:`bool <bool>`, default False
        If True, convert pressure integral :math:`\\int x\\,dp` to mass-weighted
        integral :math:`(1/g)\\int x\\,dp`.

        .. note::

            For moisture flux quantities such as :math:`q\\cdot u` or :math:`q\\cdot v`, set
            `mass_weighted=True` to obtain the commonly used vertically integrated
            moisture flux with units :math:`\\mathrm{kg \, m^{-1} \, s^{-1}}`.

    gravity: :py:class:`float <float>`, default 9.80665
        Gravitational acceleration used for mass weighting, in :math:`\\mathrm{m \, s^{-2}}.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        Vertical pressure integral, or mass-weighted integral depending on the
        selected options.
    """
    if gravity <= 0:
        raise ValueError("`gravity` must be positive.")
    if normalize and mass_weighted:
        raise ValueError("`normalize` and `mass_weighted` cannot both be True.")
    if normalize:
        return calc_top2surface_average(
            data_input=data_input,
            surface_pressure_data=surface_pressure_data,
            vertical_dim=vertical_dim,
            surface_pressure_data_units=surface_pressure_data_units,
            vertical_dim_units=vertical_dim_units,
        )

    integral, _, surface_pressure_data = _calc_top2surface_vibeta_ncl(
        data_input=data_input,
        surface_pressure_data=surface_pressure_data,
        vertical_dim=vertical_dim,
        surface_pressure_data_units=surface_pressure_data_units,
        vertical_dim_units=vertical_dim_units,
    )

    integral_type = "integral summation"
    input_units = data_input.attrs.get("units", "").strip()
    units = f"{input_units} Pa".strip()

    if mass_weighted:
        integral = integral / gravity
        integral_type = "mass-weighted integral"
        if input_units in {"m s-1", "m s**-1", "m s^-1"}:
            units = "kg m-1 s-1"
        else:
            units = f"({input_units} Pa)/(m s-2)" if input_units else "Pa/(m s-2)"

    return _build_top2surface_output(
        integral=integral,
        surface_pressure_data=surface_pressure_data,
        data_input=data_input,
        units=units,
        integral_type=integral_type,
        normalize=False,
        mass_weighted=mass_weighted,
        method="ncl",
    )


def calc_top2surface_average(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
) -> xr.DataArray:
    """
    Calculate the pressure-thickness-weighted layer average from the top to surface
    with Fortran vibeta backends.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data to be calculated.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`
        Surface level pressure. Must be local surface pressure, not MSLP.
    vertical_dim: :py:class:`str <str>`
        Vertical coordinate dimension name.
    surface_pressure_data_units: {"hPa", "Pa", "mbar"}
        Unit of `surface_pressure_data`.
    vertical_dim_units: {"hPa", "Pa", "mbar"}
        Unit of vertical pressure coordinate.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        The pressure-thickness-weighted layer average from the top to surface.
    """
    integral, layer_thickness, surface_pressure_data = _calc_top2surface_vibeta_ncl(
        data_input=data_input,
        surface_pressure_data=surface_pressure_data,
        vertical_dim=vertical_dim,
        surface_pressure_data_units=surface_pressure_data_units,
        vertical_dim_units=vertical_dim_units,
    )
    average = integral / layer_thickness
    units = data_input.attrs.get("units", "")

    return _build_top2surface_output(
        integral=average,
        surface_pressure_data=surface_pressure_data,
        data_input=data_input,
        units=units,
        integral_type="integral average",
        normalize=True,
        mass_weighted=False,
        method="ncl",
    )


def calc_top2surface_integral_rs(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    normalize: bool = False,
    mass_weighted: bool = False,
    gravity: float = 9.80665,
    method: Literal["rust", "rust-block"] = "rust-block",
) -> xr.DataArray:
    """
    Calculate the vertical integral in the p-coordinate system from the ground
    to the zenith along the barometric pressure direction with Rust vibeta backends.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data to be calculated.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`
        Surface level pressure. Must be local surface pressure, not MSLP.
    vertical_dim: :py:class:`str <str>`
        Vertical coordinate dimension name.
    surface_pressure_data_units: {"hPa", "Pa", "mbar"}
        Unit of `surface_pressure_data`.
    vertical_dim_units: {"hPa", "Pa", "mbar"}
        Unit of vertical pressure coordinate.
    mass_weighted: :py:class:`bool <bool>`, default False
        If True, convert pressure integral :math:`\\int x\\,dp` to mass-weighted
        integral :math:`(1/g)\\int x\\,dp`.

        .. note::

            For moisture flux quantities such as :math:`q\\cdot u` or :math:`q\\cdot v`, set
            `mass_weighted=True` to obtain the commonly used vertically integrated
            moisture flux with units :math:`\\mathrm{kg \, m^{-1} \, s^{-1}}`.

    gravity: :py:class:`float <float>`, default 9.80665
        Gravitational acceleration used for mass weighting, in :math:`\\mathrm{m \, s^{-2}}.
    method: {"rust", "rust-block"}
        Rust backend engine.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        Vertical pressure integral, or mass-weighted integral depending on the
        selected options.
    """
    if gravity <= 0:
        raise ValueError("`gravity` must be positive.")
    if normalize and mass_weighted:
        raise ValueError("`normalize` and `mass_weighted` cannot both be True.")
    if normalize:
        return calc_top2surface_average_rs(
            data_input=data_input,
            surface_pressure_data=surface_pressure_data,
            vertical_dim=vertical_dim,
            surface_pressure_data_units=surface_pressure_data_units,
            vertical_dim_units=vertical_dim_units,
            method=method,
        )

    integral, _, surface_pressure_data = _calc_top2surface_rust(
        data_input=data_input,
        surface_pressure_data=surface_pressure_data,
        vertical_dim=vertical_dim,
        surface_pressure_data_units=surface_pressure_data_units,
        vertical_dim_units=vertical_dim_units,
        method=method,
    )

    input_units = data_input.attrs.get("units", "").strip()
    units = f"{input_units} Pa".strip()
    integral_type = "integral summation"

    if mass_weighted:
        integral = integral / gravity
        units = _get_mass_weighted_units(input_units)
        integral_type = "mass-weighted integral"

    return _build_top2surface_output(
        integral=integral,
        surface_pressure_data=surface_pressure_data,
        data_input=data_input,
        units=units,
        integral_type=integral_type,
        normalize=False,
        mass_weighted=mass_weighted,
        method=method,
    )


def calc_top2surface_average_rs(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    method: Literal["rust", "rust-block"] = "rust-block",
) -> xr.DataArray:
    """
    Calculate the pressure-thickness-weighted layer average from the top to surface
    with Rust vibeta backends.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
        The spatio-temporal data to be calculated.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`
        Surface level pressure. Must be local surface pressure, not MSLP.
    vertical_dim: :py:class:`str <str>`
        Vertical coordinate dimension name.
    surface_pressure_data_units: {"hPa", "Pa", "mbar"}
        Unit of `surface_pressure_data`.
    vertical_dim_units: {"hPa", "Pa", "mbar"}
        Unit of vertical pressure coordinate.
    method: {"rust", "rust-block"}
        Rust backend engine.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        The pressure-thickness-weighted layer average from the top to surface.
    """
    integral, layer_thickness, surface_pressure_data = _calc_top2surface_rust(
        data_input=data_input,
        surface_pressure_data=surface_pressure_data,
        vertical_dim=vertical_dim,
        surface_pressure_data_units=surface_pressure_data_units,
        vertical_dim_units=vertical_dim_units,
        method=method,
    )
    average = integral / layer_thickness

    return _build_top2surface_output(
        integral=average,
        surface_pressure_data=surface_pressure_data,
        data_input=data_input,
        units=data_input.attrs.get("units", ""),
        integral_type="integral average",
        normalize=True,
        mass_weighted=False,
        method=method,
    )


def calc_delta_pressure(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
) -> xr.DataArray:
    """
    Calculates the pressure layer thickness (delta pressure) of a constant
    pressure level coordinate system.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The spatio-temporal data to be calculated.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Mean surface sea level pressure.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    surface_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.

    Returns
    -------
    The pressure layer thickness (delta pressure) of a constant pressure level coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso::
        - :py:func:`geocat.comp.meteorology.delta_pressure <geocat-comp:geocat.comp.meteorology.delta_pressure>`
        - `dpres_plevel - NCL <https://www.ncl.ucar.edu/Document/Functions/Built-in/dpres_plevel.shtml>`__

    Examples
    --------
    The results in :py:func:`geocat.comp.meteorology.delta_pressure <geocat-comp:geocat.comp.meteorology.delta_pressure>`:

    >>> from geocat.comp.meteorology import delta_pressure
    >>> dp = delta_pressure(
    ...     pressure_lev= np.array([1000.,925.,850.,700.,600.,500., 400.,300.,250.,200.,150.,100., 70.,50.,30.,20.,10.]),
    ...     surface_pressure = np.array([1013]),
    ... )
    >>> print(dp)
    [[ 50.5  75.  112.5 125.  100.  100.  100.   75.   50.   50.   50.   40.
       25.   20.   15.   10.    5. ]]

    For comparison, the results in :py:func:`easyclimate.calc_delta_pressure <calc_delta_pressure>`:

    >>> temp_sample = xr.DataArray(
    ...     np.array([[292.,285.,283.,277.,270.,260., 250.,235.,225.,215.,207.,207., 213.,220.,225.,228.,230.]]),
    ...     dims = ("lat", "plev"),
    ...     coords = {"plev": np.array([1000.,925.,850.,700.,600.,500., 400.,300.,250.,200.,150.,100., 70.,50.,30.,20.,10.]),
    ...             "lat": np.array([0])}
    ... )
    >>> dp = ecl.calc_delta_pressure(
    ...     data_input = temp_sample,
    ...     surface_pressure_data = xr.DataArray([1013], dims = "lat"),
    ...     vertical_dim = "plev",
    ...     surface_pressure_data_units = "Pa",
    ...     vertical_dim_units = "Pa",
    ... ).transpose("lat", "plev")
    >>> print(dp)
    <xarray.DataArray 'plev' (lat: 1, plev: 17)> Size: 136B
    array([[ 50.5,  75. , 112.5, 125. , 100. , 100. , 100. ,  75. ,  50. ,
            50. ,  50. ,  40. ,  25. ,  20. ,  15. ,  10. ,   5. ]])
    Coordinates:
    * lat      (lat) int64 8B 0
    * plev     (plev) float64 136B 1e+03 925.0 850.0 700.0 ... 50.0 30.0 20.0 10.0
    """
    from .diff import calc_gradient

    # `vertical_dim` data is forced into descending order
    data_input = data_input.sortby(vertical_dim, ascending=False)

    # Change the vertical coordinate unit to `Pa`
    vertical_dim_base = transfer_units_coeff(vertical_dim_units, "Pa")
    data_input = data_input.assign_coords(
        {vertical_dim: data_input[vertical_dim] * vertical_dim_base}
    )

    # Change the surface pressure data unit to `Pa`
    surface_pressure_data_Pa = transfer_data_multiple_units(
        surface_pressure_data, surface_pressure_data_units, "Pa"
    )

    pressure_lev = xr.broadcast(data_input[vertical_dim], data_input)[0]
    pressure_top = pressure_lev.min(dim=vertical_dim)
    surface_pressure = surface_pressure_data_Pa

    delta_pressure = np.abs(calc_gradient(pressure_lev, dim=vertical_dim))
    # top level
    delta_pressure[{vertical_dim: -1}] = (
        pressure_lev[{vertical_dim: -1}] + pressure_lev[{vertical_dim: -2}]
    ) / 2 - pressure_top
    # bottom level
    delta_pressure[{vertical_dim: 0}] = (
        surface_pressure
        - (pressure_lev[{vertical_dim: 0}] + pressure_lev[{vertical_dim: 1}]) / 2
    )

    # Add attributes to the output DataArray
    delta_pressure.attrs = {
        "units": "Pa",
        "long_name": "pressure layer thickness",
        "description": "Thickness of pressure layers in the vertical coordinate system",
        "standard_name": "pressure_thickness",
    }
    delta_pressure[vertical_dim].attrs = {"units": "Pa"}
    return delta_pressure


def calc_p_integral(
    data_input: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    normalize: bool = False,
) -> xr.DataArray:
    """
    Calculate the vertical integral along the barometric pressure direction in the p-coordinate system.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The spatio-temporal data to be calculated.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    normalize: :py:class:`bool<bool>`, default: `True`.
        Whether or not the integral results are averaged over the entire layer.

    Returns
    -------
    The vertical integral along the barometric pressure direction in the p-coordinate system (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`).

    .. attention::
        This method ignores the effect of topography, so it applies to altitudes **above 900hPa** and is **NOT applicable to the Tibetan Plateau region**.
        For a fully accurate vertical integration, please use the :py:func:`calc_top2surface_integral <calc_top2surface_integral>` function to calculate,
        but the speed of the calculation is slightly slowed down.
    """
    # Change the vertical coordinate unit to `Pa`
    vertical_dim_base = transfer_units_coeff(vertical_dim_units, "Pa")
    data_input = data_input.assign_coords(
        {vertical_dim: data_input[vertical_dim] * vertical_dim_base}
    )

    part1 = data_input.sortby(vertical_dim).integrate(
        coord=vertical_dim
    )  # unit: [data_input] *kg *s^-2

    if normalize == False:
        fint = part1

    elif normalize == True:
        part2 = (
            xr.ones_like(data_input).sortby(vertical_dim).integrate(coord=vertical_dim)
        )
        fint = part1 / part2  # unit: [data_input] *kg *s^-2

    else:
        raise ValueError("The parameter `normalize` should be `True` or `False`.")

    # Create output with coords from surface_pressure_data
    output = xr.DataArray(
        fint,
        dims=fint.dims,
        coords=fint.coords,
        attrs={
            "long_name": f"Vertical integral of {data_input.name} from top to surface",
            "units": f'{data_input.attrs.get("units", "")} Pa',
        },
        name=data_input.name,
    )
    return output

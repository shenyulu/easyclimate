"""
Water Flux
"""

from __future__ import annotations

from .units import (
    transfer_data_multiple_units,
)
import xarray as xr
from typing import Literal

__all__ = [
    "calc_horizontal_water_flux",
    "calc_vertical_water_flux",
    "calc_water_flux_top2surface_integral",
    "calc_divergence_watervaporflux",
    "calc_divergence_watervaporflux_top2surface_integral",
]


def calc_horizontal_water_flux(
    specific_humidity_data: xr.DataArray,
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    g: float = 9.8,
) -> xr.Dataset:
    """
    Calculate horizontal water vapor flux at each vertical level.

    .. math::
        \\frac{1}{g} q \\mathbf{V} = \\frac{1}{g} (u q\\ \\mathbf{i} + vq\\ \\mathbf{j})

    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.

    Returns
    -------
    The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`).

    - :math:`qu`: zonal water vapor flux.
    - :math:`qv`: meridional water vapor flux.

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    water_flux = xr.Dataset(
        data_vars={
            "qu": specific_humidity_data * u_data / g,
            "qv": specific_humidity_data * v_data / g,
        }
    )
    return water_flux


def calc_vertical_water_flux(
    specific_humidity_data: xr.DataArray, omega_data: xr.DataArray, g: float = 9.8
) -> xr.DataArray:
    """
    Calculate vertical water vapor flux.

    .. math::
        -\\omega \\frac{q}{g}

    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    omega_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The vertical velocity data (:math:`\\frac{\\mathrm{d} p}{\\mathrm{d} t}`).
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.

    Returns
    -------
    The vertical water flux. (:py:class:`xarray.DataArray <xarray.DataArray>`).

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    water_flux = -omega_data * specific_humidity_data / g
    return water_flux


def calc_water_flux_top2surface_integral(
    specific_humidity_data: xr.DataArray,
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    specific_humidity_data_units: Literal["kg/kg", "g/kg", "g/g"],
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    method: Literal["ncl", "rust", "rust-block"] = "rust-block",
    g: float = 9.8,
) -> xr.DataArray:
    """
    Calculate the water vapor flux across the vertical level.

    .. math::

        \\frac{1}{g} \\int_0^{p_s} (q\\mathbf{v}),dp

    Parameters
    ----------
    specific_humidity: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Mean surface sea level pressure.
    surface_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    specific_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    method: {"ncl", "rust", "rust-block"}, default: `rust-block`.
        Vertical integration backend.
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.

    Returns
    -------
    The water vapor flux. (:py:class:`xarray.Dataset<xarray.Dataset>`, :math:`\\mathrm{kg \\cdot m^-1 \\cdot s^-1 }`).

    - :math:`qu`: zonal water vapor flux.
    - :math:`qv`: meridional water vapor flux.

    .. seealso::
        :py:func:`calc_top2surface_integral <calc_top2surface_integral>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    from .integral import calc_top2surface_integral, calc_top2surface_integral_rs

    specific_humidity_kg_kg = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_data_units, "kg/kg"
    )

    # Calculate the single-layer water flux
    water_flux_single_layer = calc_horizontal_water_flux(
        specific_humidity_kg_kg, u_data, v_data, g=g
    )  # 1/g *(q \mathbf{v})
    qu = water_flux_single_layer["qu"]
    qv = water_flux_single_layer["qv"]

    if method == "ncl":
        # Calculate the integral of `qu` over the whole atmosphere
        qu_top2surface_integral = calc_top2surface_integral(
            data_input=qu,
            vertical_dim=vertical_dim,
            vertical_dim_units=vertical_dim_units,
            surface_pressure_data=surface_pressure_data,
            surface_pressure_data_units=surface_pressure_data_units,
            mass_weighted=False,
            gravity=g,
        )

        # Calculate the integral of `qv` over the whole atmosphere
        qv_top2surface_integral = calc_top2surface_integral(
            data_input=qv,
            vertical_dim=vertical_dim,
            vertical_dim_units=vertical_dim_units,
            surface_pressure_data=surface_pressure_data,
            surface_pressure_data_units=surface_pressure_data_units,
            mass_weighted=False,
            gravity=g,
        )
    elif method in {"rust", "rust-block"}:
        # Calculate the integral of `qu` over the whole atmosphere
        qu_top2surface_integral = calc_top2surface_integral_rs(
            data_input=qu,
            vertical_dim=vertical_dim,
            vertical_dim_units=vertical_dim_units,
            surface_pressure_data=surface_pressure_data,
            surface_pressure_data_units=surface_pressure_data_units,
            mass_weighted=False,
            gravity=g,
            method=method,
        )

        # Calculate the integral of `qv` over the whole atmosphere
        qv_top2surface_integral = calc_top2surface_integral_rs(
            data_input=qv,
            vertical_dim=vertical_dim,
            vertical_dim_units=vertical_dim_units,
            surface_pressure_data=surface_pressure_data,
            surface_pressure_data_units=surface_pressure_data_units,
            mass_weighted=False,
            gravity=g,
            method=method,
        )
    else:
        raise ValueError("method should be `ncl`, `rust`, or `rust-block`.")

    quv_top2surface_integral = xr.Dataset(
        data_vars={"qu": qu_top2surface_integral, "qv": qv_top2surface_integral}
    )

    quv_top2surface_integral.attrs = dict()
    quv_top2surface_integral.attrs["long_name"] = (
        "Vertical integral of water vapour flux (1/g \int q U dp)"
    )
    quv_top2surface_integral.attrs["units"] = "kg m**-1 s**-1"

    return quv_top2surface_integral


def calc_divergence_watervaporflux(
    specific_humidity_data: xr.DataArray,
    qu_data: xr.DataArray,
    qv_data: xr.DataArray,
    specific_humidity_data_units: Literal["kg/kg", "g/kg", "g/g"],
    cyclic_boundary_setting: Literal["nan", "cyclic", "cyclic+diff", "diff"] = "nan",
    method: Literal["ncl", "rust-batch", "rust-raw"] = "rust-batch",
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    R: float = 6.37122e6,
) -> xr.DataArray:
    """
    Calculate water vapor flux divergence at each vertical level.

    .. math::
        \\nabla \\left( \\frac{1}{g} q \\mathbf{V} \\right) = \\frac{1}{g} \\nabla \\cdot \\left( q \\mathbf{V} \\right)


    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    qu_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal horizontal water flux (i.e., :math:`qu`).
    qv_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional horizontal water flux (i.e., :math:`qv`).
    specific_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
    cyclic_boundary_setting: {"nan", "cyclic", "cyclic+diff", "diff"}, default: `nan`.
        A scalar integer equal to the boundary condition option:

        - ``nan``: Boundary points are set to the missing value.
        - ``cyclic``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic point.) The upper and lower boundaries will be set to missing.
        - ``cyclic+diff``: Boundary points are estimated using one-sided difference schemes normal to the boundary.
        - ``diff``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic points.) The upper and lower boundaries are estimated using a one-sided difference scheme normal to the boundary.

    method: {"ncl", "rust-batch", "rust-raw"}, default: `rust-batch`.
        The method to calculate horizontal divergence term. Optional values are ``ncl``, ``rust-batch``, or ``rust-raw``.
    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    R: :py:class:`float <float>`, default: `6.37122e6`.
        Radius of the Earth.

        .. note::
            The parameter is applicable only when ``method = "rust-batch"`` or ``method = "rust-raw"``.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        The water vapor flux divergence (i.e., :math:`q\\mathbf{V}`).

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    from .rvdv import calc_divergence_ncl, calc_divergence_rs

    specific_humidity_data_kgkg = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_data_units, "kg/kg"
    )

    flux_u = specific_humidity_data_kgkg * qu_data
    flux_v = specific_humidity_data_kgkg * qv_data

    if method in ["ncl"]:
        div_watervaporflux = calc_divergence_ncl(
            u_data=flux_u,
            v_data=flux_v,
            lon_dim=lon_dim,
            lat_dim=lat_dim,
            cyclic_boundary_setting=cyclic_boundary_setting,
        )
    elif method in ["rust-batch", "rust-raw"]:
        div_watervaporflux = calc_divergence_rs(
            u_data=flux_u,
            v_data=flux_v,
            lon_dim=lon_dim,
            lat_dim=lat_dim,
            R=R,
            cyclic_boundary_setting=cyclic_boundary_setting,
            method=method,
        )
    else:
        raise ValueError(
            f"Unsupported method={method!r}. Expected one of 'ncl', 'rust-batch', or 'rust-raw'."
        )

    return div_watervaporflux


def calc_divergence_watervaporflux_top2surface_integral(
    specific_humidity_data: xr.DataArray,
    u_data: xr.DataArray,
    v_data: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    vertical_dim: str,
    specific_humidity_data_units: Literal["kg/kg", "g/kg", "g/g"],
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    cyclic_boundary_setting: Literal["nan", "cyclic", "cyclic+diff", "diff"] = "nan",
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    integral_method: Literal["ncl", "rust", "rust-block"] = "rust-block",
    div_method: Literal["raw", "ncl", "rust-batch", "rust-raw"] = "rust-batch",
    g: float = 9.8,
    R: float = 6.37122e6,
) -> xr.DataArray:
    """
    Calculate water vapor flux divergence across the vertical level.

    .. math::

        \\nabla \\cdot \\frac{1}{g} \\int_0^{p_s} (q\\mathbf{v}),dp

    Parameters
    ----------
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    u_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The zonal wind data.
    v_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The meridional wind data.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Mean surface sea level pressure.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    specific_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
    surface_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    cyclic_boundary_setting: {"nan", "cyclic", "cyclic+diff", "diff"}, default: `nan`.
        A scalar integer equal to the boundary condition option:

        - ``nan``: Boundary points are set to the missing value.
        - ``cyclic``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic point.) The upper and lower boundaries will be set to missing.
        - ``cyclic+diff``: Boundary points are estimated using one-sided difference schemes normal to the boundary.
        - ``diff``: The u and v arrays are cyclic in longitude. (The arrays should **NOT** include the cyclic points.) The upper and lower boundaries are estimated using a one-sided difference scheme normal to the boundary.

        .. note::
            The parameter is applicable only when ``div_method = ncl`` or ``div_method = rust-batch`` or ``div_method = rust-raw``.

    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    integral_method: {"ncl", "rust", "rust-block"}, default: `rust-block`.
        The vertical integration backend.
    div_method: {"raw", "ncl", "rust-batch", "rust-raw"}, default: `rust-batch`.
        The method to calculate horizontal divergence term. Optional values are ``ncl``, ``rust-batch``, or ``rust-raw``.
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.
    R: :py:class:`float <float>`, default: `6.37122e6`.
        Radius of the Earth.

        .. note::
            The parameter is applicable only when ``div_method = "rust-batch"`` or ``div_method = "rust-raw"``.

    Returns
    -------
    The water vapor flux divergence. (:py:class:`xarray.DataArray<xarray.DataArray>`, :math:`\\mathrm{kg \\cdot m^-2 \\cdot s^-1 }`).

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_geographic_finite_difference.py
    """
    from .rvdv import calc_divergence, calc_divergence_ncl, calc_divergence_rs

    specific_humidity_kg_kg = transfer_data_multiple_units(
        specific_humidity_data, specific_humidity_data_units, "kg/kg"
    )

    # Calculation water vapor flux of the whole layer
    quv = calc_water_flux_top2surface_integral(
        specific_humidity_data=specific_humidity_kg_kg,
        u_data=u_data,
        v_data=v_data,
        surface_pressure_data=surface_pressure_data,
        surface_pressure_data_units=surface_pressure_data_units,
        specific_humidity_data_units="kg/kg",
        vertical_dim=vertical_dim,
        vertical_dim_units=vertical_dim_units,
        method=integral_method,
        g=g,
    )  # 1/g \int q \mathbf{V} dp

    # Calculation of water vapor flux divergence
    if div_method == "raw":
        div_quv = calc_divergence(
            quv["qu"],
            quv["qv"],
            lon_dim=lon_dim,
            lat_dim=lat_dim,
            R=R,
        )
    elif div_method in ["ncl"]:
        div_quv = calc_divergence_ncl(
            quv["qu"],
            quv["qv"],
            lon_dim=lon_dim,
            lat_dim=lat_dim,
            cyclic_boundary_setting=cyclic_boundary_setting,
        )
    elif div_method in ["rust-batch", "rust-raw"]:
        div_quv = calc_divergence_rs(
            quv["qu"],
            quv["qv"],
            lon_dim=lon_dim,
            lat_dim=lat_dim,
            R=R,
            cyclic_boundary_setting=cyclic_boundary_setting,
            method=div_method,
        )
    else:
        raise ValueError(
            f"Unsupported div_method={div_method!r}. Expected one of 'ncl', 'rust-batch', or 'rust-raw'."
        )

    div_quv.attrs = dict()
    div_quv.name = "wvdiv"
    div_quv.attrs["long_name"] = "Divergence of vertical integral of water vapour flux"
    div_quv.attrs["units"] = "kg m**-2 s**-1"

    result = xr.Dataset(data_vars={"qu": quv["qu"], "qv": quv["qv"], "wvdiv": div_quv})

    return result

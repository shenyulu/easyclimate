"""
Interpolates Community Atmosphere Model (CAM) or Community Earth System Model (CESM) hybrid coordinates to pressure coordinates

.. seealso::

    - https://ncar.github.io/CAM/doc/build/html/
    - https://www.cesm.ucar.edu/models/cam
    - https://www.cesm.ucar.edu/
"""

from __future__ import annotations

import numpy as np
import xarray as xr
from ..backend import _vinth2p_dp, _vinth2p_ecmwf, _vintp2p_ecmwf
from ..core.utility import (
    compare_multi_dataarray_coordinate,
    transfer_data_multiple_units,
)
from typing import Literal, Optional

__all__ = ["interp_vinth2p_dp", "interp_vinth2p_ecmwf", "interp_vintp2p_ecmwf"]


def interp_vinth2p_dp(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    hybrid_A_coefficients: xr.DataArray,
    hybrid_B_coefficients: xr.DataArray,
    vertical_output_level: list[int | float],
    vertical_input_dim: str,
    vertical_output_dim: str,
    vertical_output_dim_units: str,
    interp_method: Literal["linear", "log", "loglog"] = "linear",
    extrapolation: bool = False,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
    p0_hPa=1000.0,
) -> xr.DataArray:
    """
    Interpolate atmospheric data from Community Atmosphere Model (CAM) hybrid sigma-pressure coordinates to pressure levels.

    This function performs vertical interpolation of atmospheric data (typically temperature)
    from hybrid sigma-pressure coordinates to specified pressure levels using the NCL's
    ``vinth2p`` algorithm implemented in Fortran.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray <xarray.DataArray>`.
        Input 3D temperature field on hybrid levels with dimensions, e.g., (time, lev, lat, lon).
    surface_pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`.
        Surface pressure field with dimensions, e.g., (time, lat, lon).
    surface_pressure_data_units : Literal["hPa", "Pa", "mbar"]
        Units of the surface pressure data.
    hybrid_A_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`.
        Hybrid A coefficients (pressure term) for the model levels.
    hybrid_B_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`.
        Hybrid B coefficients (sigma term) for the model levels.
    vertical_output_level : list[:py:class:`int <in>`t | :py:class:`float <float>`]
        List of target pressure levels for interpolation.
    vertical_input_dim : :py:class:`str <str>`.
        Name of the vertical dimension in the input data.
    vertical_output_dim : :py:class:`str <str>`.
        Name to use for the vertical dimension in the output data.
    vertical_output_dim_units : :py:class:`str <str>`.
        Units for the output pressure levels (must be convertible to hPa).
    interp_method : Literal["linear", "log", "loglog"], optional
        Interpolation method:

        - ``"linear"``: Linear interpolation
        - ``"log"``: Logarithmic interpolation
        - ``"loglog"``: Log-log interpolation

        Default is ``"linear"``.
    extrapolation : :py:class:`bool <bool>`., optional
        Whether to extrapolate below the lowest model level when needed.
        Default is ``False``.
    lon_dim : :py:class:`str <str>`., optional
        Name of the longitude dimension. Default is ``"lon"``.
    lat_dim : :py:class:`str <str>`., optional
        Name of the latitude dimension. Default is ``"lat"``.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name. Default is ``"time"``.
    p0_hPa : :py:class:`float <float>`., optional
        Reference pressure in hPa for hybrid level calculation. Default is ``1000.0`` hPa.

    Returns
    -------
    :py:class:`xarray.DataArray <xarray.DataArray>`
        Interpolated data on pressure levels with dimensions, e.g., (time, plev, lat, lon),
        where plev corresponds to vertical_output_level.

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p.shtml

    Examples
    --------
    >>> # Interpolate temperature to pressure levels
    >>> interp_vinth2p_dp(
    ...     data_input=temp_data,
    ...     surface_pressure_data=psfc_data,
    ...     surface_pressure_data_units="Pa",
    ...     hybrid_A_coefficients=hyam,
    ...     hybrid_B_coefficients=hybm,
    ...     vertical_output_level=[1000, 850, 700, 500, 300],
    ...     vertical_input_dim="lev",
    ...     vertical_output_dim="plev",
    ...     vertical_output_dim_units="hPa",
    ...     interp_method="log"
    ... )
    """
    if interp_method == "linear":
        intyp = 1
    elif interp_method == "log":
        intyp = 2
    elif interp_method == "loglog":
        intyp = 3
    else:
        raise ValueError("interp_method must be `linear`, `log`, or `loglog`.")

    if extrapolation == True:
        kxtrp = 1
    elif extrapolation == False:
        kxtrp = 0
    else:
        raise ValueError(
            "extrapolation must be `False` (no extrapolation) or `True` (extrapolate)"
        )

    surface_pressure_data = transfer_data_multiple_units(
        surface_pressure_data, surface_pressure_data_units, "Pa"
    )
    vertical_output_level = xr.DataArray(
        vertical_output_level,
        dims=[vertical_output_dim],
        coords={vertical_output_dim: vertical_output_level},
    )
    vertical_output_level = transfer_data_multiple_units(
        vertical_output_level, vertical_output_dim_units, "hPa"
    )
    vertical_output_level.attrs["standard_name"] = "air_pressure"
    vertical_output_level.attrs["long_name"] = "Level"
    vertical_output_level.attrs["positive"] = "down"
    vertical_output_level.attrs["axis"] = "Z"

    # Compute plevi using mean surface pressure for simplicity
    psfc_mean = surface_pressure_data.mean().values
    plevi = hybrid_A_coefficients * p0_hPa + hybrid_B_coefficients * (psfc_mean * 0.01)
    plevi = xr.DataArray(plevi.data, dims=vertical_input_dim)

    # Validate input shapes
    compare_multi_dataarray_coordinate(
        [hybrid_A_coefficients, hybrid_B_coefficients, plevi],
        time_dim=time_dim,
        exclude_dims=[vertical_input_dim],
    )
    compare_multi_dataarray_coordinate(
        [data_input, surface_pressure_data],
        time_dim=time_dim,
        exclude_dims=[vertical_input_dim],
    )

    def vinth2p(
        dati,
        hbcofa,
        hbcofb,
        p0,
        plevi,
        plevo,
        psfc,
        intyp=1,
        ilev=2,
        kxtrp=0,
        spvl=-99999.0,
    ):
        """
        Interpolate data from hybrid coordinate levels to pressure levels.

        This function interpolates 3D data from hybrid coordinates to specified pressure levels
        using hybrid coefficients and surface pressure. It supports linear, logarithmic, or
        log-log interpolation, with optional extrapolation below the lowest hybrid level.
        Inputs can be either NumPy arrays objects, and the output will
        preserve xarray metadata if provided.

        Parameters
        ----------
        dati : ndarray
            Input 3D array of shape (nlon, nlat, nlevi) containing data on hybrid surfaces.
            The vertical dimension is ordered from top to bottom.
        hbcofa : ndarray
            1D array of shape (nlevip1,) containing the 'A' (pressure) hybrid coefficients.
        hbcofb : ndarray
            1D array of shape (nlevip1,) containing the 'B' (sigma) hybrid coefficients.
        p0 : float
            Reference pressure in millibars (mb) for computing hybrid pressure levels.
        plevi : ndarray
            1D array of shape (nlevip1,) containing pressure values of hybrid surfaces in mb.
        plevo : ndarray
            1D array of shape (nlevo,) containing output pressure levels in mb.
            Ordered from low to high pressure (top to bottom).
        psfc : ndarray
            2D array of shape (nlon, nlat) containing surface pressure in Pascals.
        intyp : int, optional
            Interpolation type (default: 1).
            - 1: Linear interpolation
            - 2: Logarithmic interpolation
            - 3: Log-log interpolation
        ilev : int, optional
            Level type flag (default: 2).
            - 1: Data is on level interfaces (half levels)
            - 2: Data is on level midpoints (full levels)
        kxtrp : int, optional
            Extrapolation flag (default: 0).
            - 0: Do not extrapolate; use spvl for values below lowest hybrid level
            - 1: Extrapolate data below lowest hybrid level
        spvl : float, optional
            Special value to use when extrapolation is disabled (default: -99999.).

        Returns
        -------
        dato : ndarray or xarray.DataArray
            3D array of shape (nlon, nlat, nlevo) containing interpolated data on pressure levels.

        Notes
        -----
        - The input array `dati` must be 3D with dimensions (longitude, latitude, level).
        - The hybrid pressure at a level is computed as: P(k) = hbcofa(k)*p0 + hbcofb(k)*psfc*0.01
        where psfc is converted from Pascals to millibars.
        """
        dati = np.asarray(dati, dtype=np.float64)
        hbcofa = np.asarray(hbcofa, dtype=np.float64)
        hbcofb = np.asarray(hbcofb, dtype=np.float64)
        p0 = np.float64(p0)
        plevi = np.asarray(plevi, dtype=np.float64)
        plevo = np.asarray(plevo, dtype=np.float64)
        psfc = np.asarray(psfc, dtype=np.float64)
        intyp = int(intyp)
        ilev = int(ilev)
        kxtrp = int(kxtrp)
        spvl = np.float64(spvl)

        imax, nlat, nlevi = dati.shape
        nlevip1 = hbcofa.shape[0]
        nlevo = plevo.shape[0]

        # Call the Fortran subroutine and get dato directly
        dato = _vinth2p_dp(
            dati,
            hbcofa,
            hbcofb,
            p0,
            plevi,
            plevo,
            intyp,
            ilev,
            psfc,
            spvl,
            kxtrp,
            imax=imax,
            nlat=nlat,
            nlevi=nlevi,
            nlevip1=nlevip1,
            nlevo=nlevo,
        )
        return dato

    # Define wrapper function for apply_ufunc
    def vinth2p_wrapper(
        datai_slice, psfc_slice, hyam, hybm, p0, plevi, plevo, intyp, ilev, kxtrp
    ):
        """Wrapper function to handle single time slice"""
        # Transpose dimensions to (lon, lat, lev) for Fortran
        datai_slice = np.transpose(datai_slice, (2, 1, 0))  # lon, lat, lev
        psfc_slice = np.transpose(psfc_slice, (1, 0))  # lon, lat

        # Call the Fortran function
        result = vinth2p(
            datai_slice,
            hyam,
            hybm,
            p0,
            plevi,
            plevo,
            psfc_slice,
            intyp=intyp,
            ilev=ilev,
            kxtrp=kxtrp,
        )

        # Transpose back to (lev, lat, lon)
        return np.transpose(result, (2, 1, 0))  # plev, lat, lon

    # Apply vinth2p to all time slices using apply_ufunc
    result = xr.apply_ufunc(
        vinth2p_wrapper,
        data_input,  # Input data (time, lev, lat, lon)
        surface_pressure_data,  # Surface pressure (time, lat, lon)
        input_core_dims=[
            [vertical_input_dim, lat_dim, lon_dim],
            [lat_dim, lon_dim],
        ],  # Core dimensions to reduce
        output_core_dims=[
            [vertical_output_dim, lat_dim, lon_dim]
        ],  # Output core dimensions
        exclude_dims=set((vertical_input_dim,)),  # Remove lev from input
        vectorize=True,  # Handle multiple time steps
        dask="parallelized",  # Enable parallel processing
        output_dtypes=[np.float64],  # Output data type
        kwargs={
            "hyam": hybrid_A_coefficients.values,
            "hybm": hybrid_B_coefficients.values,
            "p0": p0_hPa,
            "plevi": plevi.values,
            "plevo": vertical_output_level.values,
            "intyp": intyp,
            "ilev": 2,
            "kxtrp": kxtrp,
        },
        dask_gufunc_kwargs={
            "allow_rechunk": True,
            "output_sizes": {vertical_output_dim: len(vertical_output_level)},
        },
    )

    result = result.where(result != -99999.0, np.nan)
    # Assign coordinates to the output
    result = result.assign_coords({vertical_output_dim: vertical_output_level})
    return result


def interp_vinth2p_ecmwf(
    data_input: xr.DataArray,
    surface_pressure_data: xr.DataArray,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    hybrid_A_coefficients: xr.DataArray,
    hybrid_B_coefficients: xr.DataArray,
    vertical_output_level: list[int | float],
    vertical_input_dim: str,
    vertical_output_dim: str,
    vertical_output_dim_units: str,
    variable_flag: Literal["T", "Z", "other"],
    temperature_bottom_data: Optional[xr.DataArray] = None,
    surface_geopotential_data: Optional[xr.DataArray] = None,
    interp_method: Literal["linear", "log", "loglog"] = "linear",
    extrapolation: bool = True,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
    p0_hPa: float = 1000.0,
) -> xr.DataArray:
    """
    Interpolate atmospheric data from Community Atmosphere Model (CAM) hybrid sigma-pressure
    coordinates to pressure levels using ECMWF extrapolation methods.

    This function performs vertical interpolation of atmospheric data from hybrid sigma-pressure
    coordinates to specified pressure levels using the NCL's ``vinth2p_ecmwf`` algorithm implemented
    in Fortran. It supports ECMWF-specific extrapolation for temperature ('T') and geopotential
    height ('Z') below the lowest hybrid level.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray <xarray.DataArray>`
        Input 3D field (e.g., temperature or geopotential) on hybrid levels with dimensions,
        e.g., (time, lev, lat, lon).
    surface_pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`
        Surface pressure field with dimensions, e.g., (time, lat, lon).
    surface_pressure_data_units : ``Literal["hPa", "Pa", "mbar"]``
        Units of the surface pressure data.
    hybrid_A_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`
        Hybrid A coefficients (pressure term) for the model levels.
    hybrid_B_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`
        Hybrid B coefficients (sigma term) for the model levels.
    vertical_output_level : list[:py:class:`int <int>` | :py:class:`float <float>`]
        List of target pressure levels for interpolation (in specified units).
    vertical_input_dim : :py:class:`str <str>`
        Name of the vertical dimension in the input data.
    vertical_output_dim : :py:class:`str <str>`
        Name to use for the vertical dimension in the output data.
    vertical_output_dim_units : :py:class:`str <str>`
        Units for the output pressure levels (must be convertible to hPa).
    variable_flag : ``Literal["T", "Z", "other"]``
        Indicates the type of variable being interpolated:
        - "T": Temperature (uses ECMWF extrapolation if enabled)
        - "Z": Geopotential height (uses ECMWF extrapolation if enabled)
        - "other": Any other variable (uses lowest level value for extrapolation)
    temperature_bottom_data : Optional[:py:class:`xarray.DataArray <xarray.DataArray>`], optional
        Temperature at the lowest model level (required for 'Z' extrapolation).
        Dimensions, e.g., (time, lat, lon). Default is None.
    surface_geopotential_data : Optional[:py:class:`xarray.DataArray <xarray.DataArray>`], optional
        Surface geopotential (required for 'T' or 'Z' extrapolation).
        Dimensions, e.g., (time, lat, lon). Default is None.
    interp_method : ``Literal["linear", "log", "loglog"]``, optional
        Interpolation method:

        - ``"linear"``: Linear interpolation
        - ``"log"``: Logarithmic interpolation
        - ``"loglog"``: Log-log interpolation

        Default is ``"linear"``.
    extrapolation : :py:class:`bool <bool>`, optional
        Whether to extrapolate below the lowest model level when needed.
        Default is False.
    lon_dim : :py:class:`str <str>`, optional
        Name of the longitude dimension. Default is "lon".
    lat_dim : :py:class:`str <str>`, optional
        Name of the latitude dimension. Default is "lat".
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name. Default is "time".
    p0_hPa : :py:class:`float <float>`, optional
        Reference pressure in **hPa** for hybrid level calculation. Default is ``1000.0 hPa``.

    Returns
    -------
    :py:class:`xarray.DataArray <xarray.DataArray>`
        Interpolated data on pressure levels with dimensions, e.g., (time, plev, lat, lon),
        where plev corresponds to ``vertical_output_level``.

    Notes
    -----
    - The hybrid level pressure is calculated as: :math:`P = A \\cdot p_0 + B \\cdot \\mathrm{psfc}`
    - Output pressure levels are converted to hPa internally for calculations.
    - Missing values are converted to NaN in the output.
    - ECMWF extrapolation is applied only when ``extrapolation=True`` and variable_flag is 'T' or 'Z'.
    - For 'T' or 'Z', tbot and phis must be provided when ``extrapolation=True``.

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p_ecmwf.shtml

    Examples
    --------
    >>> # Interpolate temperature to pressure levels with ECMWF extrapolation
    >>> interp_vinth2p_ecmwf(
    ...     data_input=temp_data,
    ...     surface_pressure_data=psfc_data,
    ...     surface_pressure_data_units="Pa",
    ...     hybrid_A_coefficients=hyam,
    ...     hybrid_B_coefficients=hybm,
    ...     vertical_output_level=[1000, 850, 700, 500, 300],
    ...     vertical_input_dim="lev",
    ...     vertical_output_dim="plev",
    ...     vertical_output_dim_units="hPa",
    ...     variable_flag="T",
    ...     temperature_bottom_data=tbot_data,
    ...     surface_geopotential_data=phis_data,
    ...     interp_method="log",
    ...     extrapolation=True
    ... )
    """
    if interp_method == "linear":
        intyp = 1
    elif interp_method == "log":
        intyp = 2
    elif interp_method == "loglog":
        intyp = 3
    else:
        raise ValueError("interp_method must be 'linear', 'log', or 'loglog'.")

    if extrapolation:
        kxtrp = 1
    else:
        kxtrp = 0

    if variable_flag == "T":
        varflg = 1
    elif variable_flag == "Z":
        varflg = -1
    elif variable_flag == "other":
        varflg = 0
    else:
        raise ValueError("variable_flag must be 'T', 'Z', or 'other'.")

    if extrapolation and (variable_flag == "T" or variable_flag == "Z"):
        if temperature_bottom_data is None or surface_geopotential_data is None:
            raise ValueError(
                "`temperature_bottom_data` and `surface_geopotential_data` must be provided when extrapolation=True for 'T' or 'Z'."
            )

    # Convert units
    surface_pressure_data = transfer_data_multiple_units(
        surface_pressure_data, surface_pressure_data_units, "Pa"
    )
    vertical_output_level = xr.DataArray(
        vertical_output_level,
        dims=[vertical_output_dim],
        coords={vertical_output_dim: vertical_output_level},
    )
    vertical_output_level = transfer_data_multiple_units(
        vertical_output_level, vertical_output_dim_units, "hPa"
    )
    vertical_output_level.attrs["standard_name"] = "air_pressure"
    vertical_output_level.attrs["long_name"] = "Level"
    vertical_output_level.attrs["positive"] = "down"
    vertical_output_level.attrs["axis"] = "Z"

    # Compute plevi using mean surface pressure
    psfc_mean = surface_pressure_data.mean().values
    plevi = hybrid_A_coefficients * p0_hPa + hybrid_B_coefficients * (psfc_mean * 0.01)
    plevi = xr.DataArray(plevi.data, dims=vertical_input_dim)

    # Validate input shapes
    compare_multi_dataarray_coordinate(
        [hybrid_A_coefficients, hybrid_B_coefficients, plevi],
        time_dim=time_dim,
        exclude_dims=[vertical_input_dim],
    )
    compare_multi_dataarray_coordinate(
        [
            data_input,
            surface_pressure_data,
            temperature_bottom_data,
        ],
        time_dim=time_dim,
        exclude_dims=[vertical_input_dim],
    )
    compare_multi_dataarray_coordinate(
        [
            surface_pressure_data,
            surface_geopotential_data,
        ],
        time_dim=time_dim,
        exclude_dims=[vertical_input_dim, time_dim],
    )

    def vinth2pecmwf(
        dati,
        hbcofa,
        hbcofb,
        p0,
        plevi,
        plevo,
        psfc,
        varflg,
        tbot,
        phis,
        intyp=1,
        ilev=2,
        kxtrp=0,
        spvl=-99999.0,
    ):
        """
        Interpolate data from hybrid coordinate levels to pressure levels with ECMWF extrapolation.

        This function interpolates 3D data from hybrid coordinates to specified pressure levels
        using hybrid coefficients and surface pressure. It supports linear, logarithmic, or
        log-log interpolation, with optional ECMWF extrapolation below the lowest hybrid level
        for temperature ('T') or geopotential height ('Z').

        Parameters
        ----------
        dati : ndarray
            Input 3D array of shape (nlon, nlat, nlevi) containing data on hybrid surfaces.
        hbcofa : ndarray
            1D array of shape (nlevip1,) containing the 'A' (pressure) hybrid coefficients.
        hbcofb : ndarray
            1D array of shape (nlevip1,) containing the 'B' (sigma) hybrid coefficients.
        p0 : float
            Reference pressure in millibars (mb) for computing hybrid pressure levels.
        plevi : ndarray
            1D array of shape (nlevip1,) containing pressure values of hybrid surfaces in mb.
        plevo : ndarray
            1D array of shape (nlevo,) containing output pressure levels in mb.
        psfc : ndarray
            2D array of shape (nlon, nlat) containing surface pressure in Pascals.
        varflg : int
            Flag indicating the variable type:
            - 1: Temperature ('T')
            - -1: Geopotential height ('Z')
            - 0: Other variables
        tbot : ndarray
            2D array of shape (nlon, nlat) containing temperature at the lowest model level.
        phis : ndarray
            2D array of shape (nlon, nlat) containing surface geopotential.
        intyp : int, optional
            Interpolation type (default: 1).
            - 1: Linear interpolation
            - 2: Logarithmic interpolation
            - 3: Log-log interpolation
        ilev : int, optional
            Level type flag (default: 2).
            - 1: Data is on level interfaces (half levels)
            - 2: Data is on level midpoints (full levels)
        kxtrp : int, optional
            Extrapolation flag (default: 0).
            - 0: Do not extrapolate; use spvl for values below lowest hybrid level
            - 1: Extrapolate data using ECMWF formulation
        spvl : float, optional
            Special value to use when extrapolation is disabled (default: -99999.).

        Returns
        -------
        dato : ndarray or xarray.DataArray
            3D array of shape (nlon, nlat, nlevo) containing interpolated data on pressure levels.
        """
        imax, nlat, nlevi = dati.shape
        nlevip1 = hbcofa.shape[0]
        nlevo = plevo.shape[0]

        dati = np.asarray(dati, dtype=np.float64)
        hbcofa = np.asarray(hbcofa, dtype=np.float64)
        hbcofb = np.asarray(hbcofb, dtype=np.float64)
        p0 = np.float64(p0)
        plevi = np.asarray(plevi, dtype=np.float64)
        plevo = np.asarray(plevo, dtype=np.float64)
        psfc = np.asarray(psfc, dtype=np.float64)
        varflg = np.int32(varflg)
        tbot = np.asarray(tbot, dtype=np.float64)
        phis = np.asarray(phis, dtype=np.float64)
        intyp = np.int32(intyp)
        ilev = np.int32(ilev)
        kxtrp = np.int32(kxtrp)
        spvl = np.float64(spvl)

        imax, nlat, nlevi = dati.shape
        nlevip1 = hbcofa.shape[0]
        nlevo = plevo.shape[0]

        # Call the Fortran subroutine
        dato = _vinth2p_ecmwf(
            dati,
            hbcofa,
            hbcofb,
            p0,
            plevi,
            plevo,
            intyp,
            ilev,
            psfc,
            spvl,
            kxtrp,
            imax=imax,
            nlat=nlat,
            nlevi=nlevi,
            nlevip1=nlevip1,
            nlevo=nlevo,
            varflg=varflg,
            tbot=tbot,
            phis=phis,
        )
        return dato

    # Define wrapper function for apply_ufunc
    def vinth2pecmwf_wrapper(
        datai_slice,
        psfc_slice,
        tbot_slice,
        phis_slice,
        hyam,
        hybm,
        p0,
        plevi,
        plevo,
        varflg,
        intyp,
        ilev,
        kxtrp,
    ):
        """Wrapper function to handle single time slice"""
        # Transpose dimensions to (lon, lat, lev) for Fortran
        datai_slice = np.transpose(datai_slice, (2, 1, 0))  # lon, lat, lev
        psfc_slice = np.transpose(psfc_slice, (1, 0))  # lon, lat
        tbot_slice = np.transpose(tbot_slice, (1, 0))
        phis_slice = np.transpose(phis_slice, (1, 0))

        # Call the Fortran function
        result = vinth2pecmwf(
            datai_slice,
            hyam,
            hybm,
            p0,
            plevi,
            plevo,
            psfc_slice,
            varflg,
            tbot_slice,
            phis_slice,
            intyp=intyp,
            ilev=ilev,
            kxtrp=kxtrp,
        )

        # Transpose back to (lev, lat, lon)
        return np.transpose(result, (2, 1, 0))  # plev, lat, lon

    # Prepare inputs for apply_ufunc
    input_core_dims = [
        [vertical_input_dim, lat_dim, lon_dim],
        [lat_dim, lon_dim],
        [lat_dim, lon_dim],
        [lat_dim, lon_dim],
    ]

    # First two elements are straightforward
    input_arrays = [data_input.transpose(lon_dim, lat_dim, vertical_input_dim, ...)]

    surface_list_item = surface_pressure_data.transpose(lon_dim, lat_dim, ...)
    input_arrays.append(surface_list_item)

    # Handle temperature_bottom_data
    if temperature_bottom_data is None:
        input_arrays.append(surface_list_item)
    else:
        input_arrays.append(temperature_bottom_data.transpose(lon_dim, lat_dim, ...))

    # Handle surface_geopotential_data
    if surface_geopotential_data is None:
        input_arrays.append(surface_list_item)
    else:
        input_arrays.append(surface_geopotential_data.transpose(lon_dim, lat_dim, ...))

    # Apply vinth2pecmwf to all time slices using apply_ufunc
    result = xr.apply_ufunc(
        vinth2pecmwf_wrapper,
        *input_arrays,
        input_core_dims=input_core_dims,
        output_core_dims=[[vertical_output_dim, lat_dim, lon_dim]],
        exclude_dims=set((vertical_input_dim,)),
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.float64],
        # output_sizes={vertical_output_dim: len(vertical_output_level)},
        kwargs={
            "hyam": hybrid_A_coefficients.values,
            "hybm": hybrid_B_coefficients.values,
            "p0": p0_hPa,
            "plevi": plevi.values,
            "plevo": vertical_output_level.values,
            "varflg": varflg,
            "intyp": intyp,
            "ilev": 2,
            "kxtrp": kxtrp,
        },
        dask_gufunc_kwargs={
            "allow_rechunk": True,
            "output_sizes": {vertical_output_dim: len(vertical_output_level)},
        },
    )

    result = result.where(result != -99999.0, np.nan)
    # Assign coordinates to the output
    result = result.assign_coords(
        {vertical_output_dim: vertical_output_level}
    ).drop_vars(vertical_input_dim)
    return result


def interp_vintp2p_ecmwf(
    data_input: xr.DataArray,
    pressure_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    surface_pressure_data: xr.DataArray,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    vertical_output_level: list[int | float],
    vertical_input_dim: str,
    vertical_output_dim: str,
    vertical_output_dim_units: str,
    variable_flag: Literal["T", "Z", "other"],
    temperature_bottom_data: Optional[xr.DataArray] = None,
    surface_geopotential_data: Optional[xr.DataArray] = None,
    interp_method: Literal["linear", "log", "loglog"] = "linear",
    extrapolation: bool = False,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    time_dim: str = "time",
) -> xr.DataArray:
    """
    Interpolates data at multidimensional pressure levels to constant pressure coordinates and uses an ECMWF formulation to extrapolate values below ground.

    This function performs vertical interpolation of atmospheric data from input pressure levels
    to specified output pressure levels using the NCL's `vintp2p_ecmwf` algorithm implemented
    in Fortran. It supports ECMWF-specific extrapolation for temperature ('T') and geopotential
    height ('Z') below the lowest pressure level.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray <xarray.DataArray>`
        Input 3D field (e.g., temperature or geopotential) on pressure levels with dimensions,
        e.g., (time, lev, lat, lon).
    pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`
        3D pressure field corresponding to the input data levels, with dimensions,
        e.g., (time, lev, lat, lon).
    pressure_data_units : Literal["hPa", "Pa", "mbar"]
        Units of the pressure data.
    surface_pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`
        Surface pressure field with dimensions, e.g., (time, lat, lon).
    surface_pressure_data_units : Literal["hPa", "Pa", "mbar"]
        Units of the surface pressure data.
    vertical_output_level : list[:py:class:`int <int>` | :py:class:`float <float>`]
        List of target pressure levels for interpolation (in specified units).
    vertical_input_dim : :py:class:`str <str>`
        Name of the vertical dimension in the input data.
    vertical_output_dim : :py:class:`str <str>`
        Name to use for the vertical dimension in the output data.
    vertical_output_dim_units : :py:class:`str <str>`
        Units for the output pressure levels (must be convertible to hPa).
    variable_flag : ``Literal["T", "Z", "other"]``
        Indicates the type of variable being interpolated:
        - "T": Temperature (uses ECMWF extrapolation if enabled)
        - "Z": Geopotential height (uses ECMWF extrapolation if enabled)
        - "other": Any other variable (uses lowest level value for extrapolation)
    temperature_bottom_data : Optional[:py:class:`xarray.DataArray <xarray.DataArray>`], optional
        Temperature at the lowest model level (required for 'Z' extrapolation).
        Dimensions, e.g., (time, lat, lon). Default is None.
    surface_geopotential_data : Optional[:py:class:`xarray.DataArray <xarray.DataArray>`], optional
        Surface geopotential (required for 'T' or 'Z' extrapolation).
        Dimensions, e.g., (time, lat, lon). Default is None.
    interp_method : ``Literal["linear", "log", "loglog"]``, optional
        Interpolation method:

        - "linear": Linear interpolation
        - "log": Logarithmic interpolation
        - "loglog": Log-log interpolation

        Default is "linear".
    extrapolation : :py:class:`bool <bool>`, optional
        Whether to extrapolate below the lowest pressure level when needed.
        Default is False.
    lon_dim : :py:class:`str <str>`, optional
        Name of the longitude dimension. Default is "lon".
    lat_dim : :py:class:`str <str>`, optional
        Name of the latitude dimension. Default is "lat".
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name. Default is "time".

    Returns
    -------
    :py:class:`xarray.DataArray <xarray.DataArray>`
        Interpolated data on pressure levels with dimensions, e.g., (time, plev, lat, lon),
        where plev corresponds to vertical_output_level.

    Notes
    -----
    - The input pressure levels are provided directly via ``pressure_data``.
    - Output pressure levels and surface pressure are converted to hPa internally for calculations.
    - Missing values are converted to NaN in the output.
    - ECMWF extrapolation is applied only when ``extrapolation=True`` and variable_flag is 'T' or 'Z'.
    - For 'T' or 'Z', ``temperature_bottom_data`` and ``surface_geopotential_data`` must be provided when ``extrapolation=True``.

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Built-in/vintp2p_ecmwf.shtml

    Examples
    --------
    >>> # Interpolate temperature to pressure levels with ECMWF extrapolation
    >>> interp_vintp2p_ecmwf(
    ...     data=temp_data,
    ...     pressure_data=pres_data,
    ...     pressure_data_units="Pa",
    ...     surface_pressure_data=psfc_data,
    ...     surface_pressure_data_units="Pa",
    ...     vertical_output_level=[1000, 850, 700, 500, 300],
    ...     vertical_input_dim="lev",
    ...     vertical_output_dim="plev",
    ...     vertical_output_dim_units="hPa",
    ...     variable_flag="T",
    ...     temperature_bottom_data=tbot_data,
    ...     surface_geopotential_data=phis_data,
    ...     interp_method="log",
    ...     extrapolation=True
    ... )
    """
    # Map interpolation method to integer
    if interp_method == "linear":
        intyp = 1
    elif interp_method == "log":
        intyp = 2
    elif interp_method == "loglog":
        intyp = 3
    else:
        raise ValueError("interp_method must be 'linear', 'log', or 'loglog'.")

    # Map extrapolation flag
    if extrapolation:
        kxtrp = 1
    else:
        kxtrp = 0

    # Map variable flag
    if variable_flag == "T":
        varflg = 1
    elif variable_flag == "Z":
        varflg = -1
    elif variable_flag == "other":
        varflg = 0
    else:
        raise ValueError("variable_flag must be 'T', 'Z', or 'other'.")

    # Validate tbot and phis for extrapolation
    if extrapolation and (variable_flag == "T" or variable_flag == "Z"):
        if temperature_bottom_data is None or surface_geopotential_data is None:
            raise ValueError(
                "`temperature_bottom_data` and `surface_geopotential_data` must be provided when extrapolation=True for 'T' or 'Z'."
            )

    # Convert units
    pressure_data = transfer_data_multiple_units(
        pressure_data, pressure_data_units, "hPa"
    )
    surface_pressure_data = transfer_data_multiple_units(
        surface_pressure_data, surface_pressure_data_units, "Pa"
    )
    vertical_output_level = xr.DataArray(
        vertical_output_level,
        dims=[vertical_output_dim],
        coords={vertical_output_dim: vertical_output_level},
    )
    vertical_output_level = transfer_data_multiple_units(
        vertical_output_level, vertical_output_dim_units, "hPa"
    )
    vertical_output_level.attrs["standard_name"] = "air_pressure"
    vertical_output_level.attrs["long_name"] = "Level"
    vertical_output_level.attrs["positive"] = "down"
    vertical_output_level.attrs["axis"] = "Z"

    # Compute plevi using mean pressure for simplicity
    presi_mean = pressure_data.mean(
        dim=[dim for dim in pressure_data.dims if dim != vertical_input_dim]
    )
    plevi = presi_mean.values
    plevi = xr.DataArray(plevi, dims=vertical_input_dim)

    # Validate input shapes
    compare_multi_dataarray_coordinate(
        [data_input, pressure_data, surface_pressure_data],
        time_dim=time_dim,
        exclude_dims=[vertical_input_dim],
    )
    if extrapolation and (variable_flag == "T" or variable_flag == "Z"):
        compare_multi_dataarray_coordinate(
            [
                data_input,
                pressure_data,
                surface_pressure_data,
                temperature_bottom_data,
                surface_geopotential_data,
            ],
            time_dim=time_dim,
            exclude_dims=[vertical_input_dim],
        )

    def vintp2pecmwf(
        dati,
        presi,
        plevi,
        plevo,
        psfc,
        varflg,
        tbot,
        phis,
        intyp=1,
        kxtrp=0,
        spvl=-99999.0,
    ):
        """
        Interpolate data from pressure levels to specified pressure levels with ECMWF extrapolation.

        Parameters
        ----------
        dati : ndarray
            Input 3D array of shape (nlon, nlat, nlevi) containing data on pressure levels.
        presi : ndarray
            3D array of shape (nlon, nlat, nlevi) containing pressure values of input levels in hPa.
        plevi : ndarray
            1D array of shape (nlevip1,) containing pressure values for a vertical column in hPa.
        plevo : ndarray
            1D array of shape (nlevo,) containing output pressure levels in hPa.
        psfc : ndarray
            2D array of shape (nlon, nlat) containing surface pressure in Pascals.
        varflg : int
            Flag indicating the variable type:
            - 1: Temperature ('T')
            - -1: Geopotential height ('Z')
            - 0: Other variables
        tbot : ndarray
            2D array of shape (nlon, nlat) containing temperature at the lowest level.
        phis : ndarray
            2D array of shape (nlon, nlat) containing surface geopotential.
        intyp : int, optional
            Interpolation type (default: 1).
            - 1: Linear interpolation
            - 2: Logarithmic interpolation
            - 3: Log-log interpolation
        kxtrp : int, optional
            Extrapolation flag (default: 0).
            - 0: Do not extrapolate; use spvl for values below lowest level
            - 1: Extrapolate data using ECMWF formulation
        spvl : float, optional
            Special value to use when extrapolation is disabled (default: -99999.).

        Returns
        -------
        dato : ndarray
            3D array of shape (nlon, nlat, nlevo) containing interpolated data on pressure levels.
        """
        dati = np.asarray(dati, dtype=np.float64)
        presi = np.asarray(presi, dtype=np.float64)
        plevi = np.asarray(plevi, dtype=np.float64)
        plevo = np.asarray(plevo, dtype=np.float64)
        psfc = np.asarray(psfc, dtype=np.float64)
        varflg = np.int32(varflg)
        tbot = (
            np.asarray(tbot, dtype=np.float64)
            if tbot is not None
            else np.zeros_like(psfc, dtype=np.float64)
        )
        phis = (
            np.asarray(phis, dtype=np.float64)
            if phis is not None
            else np.zeros_like(psfc, dtype=np.float64)
        )
        intyp = np.int32(intyp)
        kxtrp = np.int32(kxtrp)
        spvl = np.float64(spvl)

        imax, nlat, nlevi = dati.shape
        nlevip1 = plevi.shape[0]
        nlevo = plevo.shape[0]

        # Call the Fortran subroutine
        dato = _vintp2p_ecmwf(
            dati,
            presi=presi,
            plevi=plevi,
            plevo=plevo,
            intyp=intyp,
            psfc=psfc,
            spvl=spvl,
            kxtrp=kxtrp,
            imax=imax,
            nlat=nlat,
            nlevi=nlevi,
            nlevip1=nlevip1,
            nlevo=nlevo,
            varflg=varflg,
            tbot=tbot,
            phis=phis,
        )
        return dato

    # Define wrapper function for apply_ufunc
    def vintp2pecmwf_wrapper(
        datai_slice,
        presi_slice,
        psfc_slice,
        plevi,
        plevo,
        varflg,
        tbot_slice,
        phis_slice,
        intyp,
        kxtrp,
    ):
        """Wrapper function to handle single time slice"""
        # Transpose dimensions to (lon, lat, lev) for Fortran
        datai_slice = np.transpose(datai_slice, (2, 1, 0))  # lon, lat, lev
        presi_slice = np.transpose(presi_slice, (2, 1, 0))  # lon, lat, lev
        psfc_slice = np.transpose(psfc_slice, (1, 0))  # lon, lat
        tbot_slice = (
            np.transpose(tbot_slice, (1, 0)) if tbot_slice is not None else None
        )
        phis_slice = (
            np.transpose(phis_slice, (1, 0)) if phis_slice is not None else None
        )

        # Call the Fortran function
        result = vintp2pecmwf(
            datai_slice,
            presi_slice,
            plevi,
            plevo,
            psfc_slice,
            varflg,
            tbot_slice,
            phis_slice,
            intyp=intyp,
            kxtrp=kxtrp,
        )

        # Transpose back to (lev, lat, lon)
        return np.transpose(result, (2, 1, 0))  # plev, lat, lon

    # Prepare inputs for apply_ufunc
    input_core_dims = [
        [vertical_input_dim, lat_dim, lon_dim],
        [vertical_input_dim, lat_dim, lon_dim],
        [lat_dim, lon_dim],
    ]
    input_arrays = [data_input, pressure_data, surface_pressure_data]
    kwargs = {
        "plevi": plevi.values,
        "plevo": vertical_output_level.values,
        "varflg": varflg,
        "intyp": intyp,
        "kxtrp": kxtrp,
    }
    if extrapolation and (variable_flag == "T" or variable_flag == "Z"):
        input_core_dims.append([lat_dim, lon_dim])
        input_core_dims.append([lat_dim, lon_dim])
        input_arrays.append(temperature_bottom_data)
        input_arrays.append(surface_geopotential_data)
        kwargs["tbot_slice"] = temperature_bottom_data.values
        kwargs["phis_slice"] = surface_geopotential_data.values
    else:
        kwargs["tbot_slice"] = np.zeros_like(surface_pressure_data.values[0, ...])
        kwargs["phis_slice"] = np.zeros_like(surface_pressure_data.values[0, ...])

    # Apply vintp2pecmwf to all time slices using apply_ufunc
    result = xr.apply_ufunc(
        vintp2pecmwf_wrapper,
        *input_arrays,
        input_core_dims=input_core_dims,
        output_core_dims=[[vertical_output_dim, lat_dim, lon_dim]],
        exclude_dims=set((vertical_input_dim,)),
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.float64],
        kwargs=kwargs,
        dask_gufunc_kwargs={
            "allow_rechunk": True,
            "output_sizes": {vertical_output_dim: len(vertical_output_level)},
        },
    )

    result = result.where(result != -99999.0, np.nan)
    # Assign coordinates to the output
    result = result.assign_coords({vertical_output_dim: vertical_output_level})
    return result

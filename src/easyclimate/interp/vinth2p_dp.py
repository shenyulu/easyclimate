"""
Interpolates CAM hybrid coordinates to pressure coordinates
"""

from __future__ import annotations

import numpy as np
import xarray as xr
from easyclimate_backend.vinth2p_dp._vinth2p_dp import vinth2p as _vinth2p_dp
from ..core.utility import (
    compare_multi_dataarray_coordinate,
    transfer_data_multiple_units,
)
from typing import Literal

__all__ = ["interp_vinth2p_dp"]


def interp_vinth2p_dp(
    temperature_data: xr.DataArray,
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
    lon_dim="lon",
    lat_dim="lat",
    p0_hPa=1000.0,
) -> xr.DataArray:
    """
    Interpolate atmospheric data from Community Atmosphere Model (CAM) hybrid sigma-pressure coordinates to pressure levels.

    This function performs vertical interpolation of atmospheric data (typically temperature)
    from hybrid sigma-pressure coordinates to specified pressure levels using the NCL's
    ``vinth2p`` algorithm implemented in Fortran.

    Parameters
    ----------
    temperature_data : :py:class:`xarray.DataArray <xarray.DataArray>`.
        Input 3D temperature field on hybrid levels with dimensions, e.g., (time, lev, lat, lon).
    surface_pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`.
        Surface pressure field with dimensions, e.g., (time, lat, lon).
    surface_pressure_data_units : Literal["hPa", "Pa", "mbar"]
        Units of the surface pressure data.
    hybrid_A_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`.
        Hybrid A coefficients (pressure term) for the model levels.
    hybrid_B_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`.
        Hybrid B coefficients (sigma term) for the model levels.
    vertical_output_level : list[int | float]
        List of target pressure levels for interpolation.
    vertical_input_dim : :py:class:`str <str>`.
        Name of the vertical dimension in the input data.
    vertical_output_dim : :py:class:`str <str>`.
        Name to use for the vertical dimension in the output data.
    vertical_output_dim_units : :py:class:`str <str>`.
        Units for the output pressure levels (must be convertible to hPa).
    interp_method : Literal["linear", "log", "loglog"], optional
        Interpolation method:

        - "linear": Linear interpolation
        - "log": Logarithmic interpolation
        - "loglog": Log-log interpolation

        Default is "linear".
    extrapolation : :py:class:`bool <bool>`., optional
        Whether to extrapolate below the lowest model level when needed.
        Default is ``False``.
    lon_dim : :py:class:`str <str>`., optional
        Name of the longitude dimension. Default is ``"lon"``.
    lat_dim : :py:class:`str <str>`., optional
        Name of the latitude dimension. Default is ``"lat"``.
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
    ...     temperature_data=temp_data,
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
        exclude_dims=[vertical_input_dim],
    )
    compare_multi_dataarray_coordinate(
        [temperature_data, surface_pressure_data], exclude_dims=[vertical_input_dim]
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
        temperature_data,  # Input data (time, lev, lat, lon)
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
    )

    result = result.where(result != -99999.0, np.nan)
    # Assign coordinates to the output
    result = result.assign_coords({vertical_output_dim: vertical_output_level})
    return result

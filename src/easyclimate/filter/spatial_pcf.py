"""
2D spatial parabolic cylinder function
"""

from __future__ import annotations
import numpy as np
import xarray as xr
from scipy import signal
from ..core.utility import validate_dataarrays
from ..core.stat import calc_detrend_spatial

__all__ = ["filter_2D_spatial_parabolic_cylinder_function"]


def fft2_along_dims(da, dims, time_dim="time", lon_dim="lon"):
    """Apply 2D FFT along specified dimensions in xarray DataArray."""
    # Determines the order of the original dimensions
    original_dims = da.dims

    # Check whether the input dimensions are valid
    for dim in dims:
        if dim not in original_dims:
            raise ValueError(
                f"Dimension {dim} not found in input dimensions {original_dims}"
            )

    # Gets the position of the dimension in the original array (axes parameter)
    axes = [original_dims.index(dim) for dim in dims]

    # Apply np.fft.fft2 using apply_ufunc
    result = xr.apply_ufunc(
        np.fft.fft2,
        da,
        kwargs={"axes": tuple(axes)},
        input_core_dims=[list(original_dims)],
        output_core_dims=[list(original_dims)],
        keep_attrs=True,
        dask="parallelized",
        output_dtypes=[complex],
    )

    # Frequency and wavenumber parameters
    n_time = da[time_dim].shape[0]
    n_lon = da[lon_dim].shape[0]
    frequency = np.fft.fftfreq(n_time)
    wavenumber = np.fft.fftfreq(n_lon) * n_lon

    result = result.rename({lon_dim: "k", time_dim: "freq"})
    result = result.assign_coords({"k": wavenumber, "freq": frequency})
    return result


def ifft2_along_dims(
    da, dims, original_data, freq_dim="freq", k_dim="k", lon_dim="lon", time_dim="time"
):
    """Apply 2D invert FFT along specified dimensions in xarray DataArray."""
    # Determines the order of the original dimensions
    original_dims = da.dims

    # Check whether the input dimensions are valid
    for dim in dims:
        if dim not in original_dims:
            raise ValueError(
                f"Dimension {dim} not found in input dimensions {original_dims}"
            )

    # Gets the position of the dimension in the original array (axes parameter)
    axes = [original_dims.index(dim) for dim in dims]

    # Apply np.fft.ifft2 using apply_ufunc
    result = xr.apply_ufunc(
        np.fft.ifft2,
        da,
        kwargs={"axes": tuple(axes)},
        input_core_dims=[list(original_dims)],
        output_core_dims=[list(original_dims)],
        keep_attrs=True,
        dask="parallelized",
        output_dtypes=[complex],
    )

    result = np.real(result)
    result = result.rename({k_dim: "lon", freq_dim: "time"})
    result = result.assign_coords(
        {"lon": original_data[lon_dim], "time": original_data[time_dim]}
    )
    new_dims = ("wave_type",) + original_data.dims
    result = result.transpose(*new_dims)
    return result


def filter_2D_spatial_parabolic_cylinder_function(
    zonal_wind_speed_data: xr.DataArray,
    meridional_wind_speed_data: xr.DataArray,
    z_data: xr.DataArray,
    period_min=3.0,
    period_max=30.0,
    wavenumber_min=2,
    wavenumber_max=20,
    trapping_scale_deg=6.0,
    lon_dim="lon",
    lat_dim="lat",
    time_dim: str = "time",
):
    """
    Perform space-time spectral analysis to extract equatorial wave components.

    This function filters atmospheric fields to isolate equatorial wave modes, including Kelvin waves,
    westward-moving mixed Rossby-gravity (WMRG) waves, and equatorial Rossby waves of the first and
    second kind (R1 and R2), within user-specified period and wavenumber ranges. It processes input
    data using a combination of detrending, windowing, Fourier transforms, and projection onto
    parabolic cylinder functions to decompose the fields into these wave types.

    Parameters
    ----------

    zonal_wind_speed_data: :class:`xarray.DataArray <xarray.DataArray>`
        The zonal (east-west) wind speed component.
    meridional_wind_speed_data: :class:`xarray.DataArray <xarray.DataArray>`
        The meridional (north-south) wind speed component.
    z_data: :class:`xarray.DataArray <xarray.DataArray>`
        The geopotential height.
    period_min: :class:`float`, optional (default=3.0)
        The minimum period (in days) of the waves to be extracted.
    period_max: :class:`float`, optional (default=30.0)
        The maximum period (in days) of the waves to be extracted.
    wavenumber_min: :class:`int`, optional (default=2)
        The minimum zonal wavenumber of the waves to be extracted.
    wavenumber_max: :class:`int`, optional (default=20)
        The maximum zonal wavenumber of the waves to be extracted.
    trapping_scale_deg: :class:`float`, optional (default=6.0)
        The meridional trapping scale in degrees, defining the width of the equatorial waveguide
        for the parabolic cylinder functions.
    lon_dim: :class:`str`, optional (default="lon")
        The name of the longitude dimension in the input data arrays.
    lat_dim: :class:`str`, optional (default="lat")
        The name of the latitude dimension in the input data arrays.
    time_dim: :class:`str`, optional (default="time")
        The name of the time dimension in the input data arrays.

    Returns
    -------

    :class:`xarray.Dataset <xarray.Dataset>`

    A dataset containing the filtered fields for each wave type:
        - ``'u'``: Zonal wind component (m/s)
        - ``'v'``: Meridional wind component (m/s)
        - ``'z'``: Geopotential height (m)
    Each variable includes an additional ``'wave_type'`` dimension with values:
    ``['kelvin', 'wmrg', 'r1', 'r2']``.

    Mathematical Explanation
    ------------------------

    The function implements a space-time spectral analysis method to decompose atmospheric fields into
    equatorial wave modes, based on the normal modes of the shallow water equations on an equatorial
    beta-plane. The meridional structure of these waves is represented by parabolic cylinder functions,
    which arise as solutions to the quantum harmonic oscillator problem adapted to the equatorial
    waveguide.

    Key Steps
    ~~~~~~~~~

    1. Detrending and Windowing:
        - The input fields are spatially detrended to remove large-scale trends.
        - A Tukey window (alpha=0.1) is applied along the time dimension to minimize spectral leakage.

    2. Variable Transformation:
        Following Yang et al. (2003), new variables :math:`q` and :math:`r` are defined:

        .. math::
            q = \\frac{g}{c_e} z + u, \\quad r = \\frac{g}{c_e} z - u

        where:
            - :math:`g = 9.8 \\text{m/s}^2` is gravitational acceleration.
            - :math:`c_e = 2 \\beta L^2` is a characteristic speed.
            - :math:`\\beta = 2.3 \\times 10^{-11}, \\text{m}^{-1}\\text{s}^{-1}` is the beta parameter.
            - :math:`L = \\frac{2\\pi R_e \\cdot \\text{trapping_scale_deg}}{360}` is the trapping scale in meters (:math:`R_e = 6.371 \\times 10^6, \\text{m}` is Earth's radius).
            - :math:`z` is geopotential height, and :math:`u` is zonal wind speed.

    3. Fourier Transform:
        A 2D Fast Fourier Transform (FFT) is applied along the time and longitude dimensions to convert
        the fields into frequency-wavenumber space.

    4. Projection onto Parabolic Cylinder Functions:
        - The spectral coefficients are projected onto parabolic cylinder functions :math:`D_n(y)`, where
            :math:`y = \\text{latitude} / \\text{trapping_scale_deg}` is the scaled latitude.
        - The first four modes (:math:`n = 0, 1, 2, 3`) are computed with normalization factors:

        .. math::

            D_0(y) = e^{-y^2/4}, \\quad D_1(y) = y e^{-y^2/4}, \\quad D_2(y) = (y^2 - 1) e^{-y^2/4}, \\quad D_3(y) = y (y^2 - 3) e^{-y^2/4}

        - These functions define the meridional structure of the waves.

    5. Wave Type Identification and Filtering:
        - Kelvin Wave: Uses :math:`n = 0` mode, eastward propagating within specified frequency (:math:`1/\\text{period_max}` to :math:`1/\\text{period_min}`) and wavenumber ranges.
        - WMRG Wave: Combines :math:`n = 1` for :math:`q` and :math:`n = 0` for meridional wind :math:`v`, westward propagating.
        - R1 Wave: Uses :math:`n = 2` for :math:`q`, :math:`n = 0` for :math:`r`, and :math:`n = 1` for :math:`v`, westward propagating.
        - R2 Wave: Uses :math:`n = 3` for :math:`q`, :math:`n = 1` for :math:`r`, and :math:`n = 2` for :math:`v`, westward propagating.
        - Frequencies and wavenumbers are selected based on ``period_min``, ``period_max``, ``wavenumber_min``, and ``wavenumber_max``.

    6. Reconstruction:
        - The filtered spectral coefficients are recombined with the parabolic cylinder functions and transformed back to physical space using an inverse 2D FFT.
        - The physical fields are reconstructed as:
            - :math:`u = (q + r) / 2`.
            - :math:`v` directly from its coefficients.
            - :math:`z = (q - r) \\cdot c_e / (2g)`.

    References
    ----------
    - Gill, A.E. (1980), Some simple solutions for heat-induced tropical circulation. Q.J.R. Meteorol. Soc., 106: 447-462. https://doi.org/10.1002/qj.49710644905
    - Li, X.f., Cho, HR. Development and propagation of equatorial waves. Adv. Atmos. Sci. 14, 323–338 (1997). https://doi.org/10.1007/s00376-997-0053-6
    - Yang, G.-Y., Hoskins, B., & Slingo, J. (2003). Convectively coupled equatorial waves: A new methodology for identifying wave structures in observational data. *Journal of the Atmospheric Sciences*, 60(14), 1637-1654.
    - Knippertz, P., Gehne, M., Kiladis, G.N., Kikuchi, K., Rasheeda Satheesh, A., Roundy, P.E., et al. (2022) The intricacies of identifying equatorial waves. Quarterly Journal of the Royal Meteorological Society, 148(747), 2814–2852. Available from: https://doi.org/10.1002/qj.4338

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_spatial_pcf.py
    """
    # Constants
    GRAVITY = 9.8
    EARTH_RADIUS = 6.371e6
    BETA = 2.3e-11

    # Wave parameters
    trapping_scale_meters = 2.0 * np.pi * EARTH_RADIUS * trapping_scale_deg / 360.0
    ce = 2.0 * trapping_scale_meters**2 * BETA
    gravity_over_ce = GRAVITY / ce
    ce_over_gravity = ce / GRAVITY

    validate_dataarrays([zonal_wind_speed_data, meridional_wind_speed_data, z_data])

    # Detrend data
    zonal_wind_speed_data_detrend = calc_detrend_spatial(
        zonal_wind_speed_data, time_dim=time_dim
    )
    meridional_wind_speed_data_detrend = calc_detrend_spatial(
        meridional_wind_speed_data, time_dim=time_dim
    )
    z_data_detrend = calc_detrend_spatial(z_data, time_dim=time_dim)

    # Array shape
    time_shape = zonal_wind_speed_data_detrend.shape[0]
    time_data = zonal_wind_speed_data_detrend[time_dim].data
    latitude_data = zonal_wind_speed_data_detrend[lat_dim].data

    # Tukey window
    tukey_window_raw = signal.windows.tukey(time_shape, alpha=0.1)
    tukey_window = xr.DataArray(
        tukey_window_raw, dims=time_dim, coords={time_dim: time_data}
    )

    # Apply window
    zonal_wind_windowed = zonal_wind_speed_data_detrend * tukey_window
    meridional_windowed = meridional_wind_speed_data_detrend * tukey_window
    z_windowed = z_data_detrend * tukey_window

    # Transform variables
    # see (2.2) in (Yang et all., 2003)
    q = z_windowed * gravity_over_ce + zonal_wind_windowed
    r = z_windowed * gravity_over_ce - zonal_wind_windowed

    q_fft = fft2_along_dims(q, dims=(time_dim, lon_dim))
    v_fft = fft2_along_dims(meridional_windowed, dims=(time_dim, lon_dim))
    r_fft = fft2_along_dims(r, dims=(time_dim, lon_dim))
    q_fft_reordered = q_fft.transpose("freq", "k", ...)
    v_fft_reordered = v_fft.transpose("freq", "k", ...)
    r_fft_reordered = r_fft.transpose("freq", "k", ...)

    freq = q_fft["freq"]
    n_time = q[time_dim].shape[0]
    n_lon = q[lon_dim].shape[0]
    n_lat = q[lat_dim].shape[0]

    f_min_idx = np.where(freq >= 1.0 / period_max)[0][0]
    f_max_idx = (np.where(freq > 1.0 / period_min)[0][0]) - 1

    # Index
    f_pos_start = f_min_idx
    f_pos_end = f_max_idx + 1
    f_neg_start = n_time - f_max_idx
    f_neg_end = n_time - f_min_idx + 1

    k_pos_start = wavenumber_min
    k_pos_end = wavenumber_max + 1
    k_neg_start = n_lon - wavenumber_max
    k_neg_end = n_lon - wavenumber_min + 1

    # Parabolic cylinder functions
    sqrt_2pi = np.sqrt(2.0 * np.pi)
    normalization_factors = np.array(
        [sqrt_2pi, sqrt_2pi, 2.0 * sqrt_2pi, 6.0 * sqrt_2pi]
    )
    normalization_factors_shape = normalization_factors.size

    y = latitude_data / trapping_scale_deg

    cylinder_functions = np.zeros((normalization_factors_shape, n_lat))
    exp_term = np.exp(-(y**2) / 4.0)

    cylinder_functions[0, :] = exp_term
    cylinder_functions[1, :] = y * exp_term
    cylinder_functions[2, :] = (y**2 - 1.0) * exp_term
    cylinder_functions[3, :] = y * (y**2 - 3.0) * exp_term

    cylinder_functions = xr.DataArray(
        cylinder_functions,
        dims=("mode", "lat"),
        coords={
            "mode": normalization_factors,
            lat_dim: zonal_wind_speed_data[lat_dim].data,
        },
    )

    d_lat = np.abs(latitude_data[1] - latitude_data[0]) * np.pi / 180.0

    #
    dims = ("mode",) + q_fft_reordered.dims
    shape = (normalization_factors_shape,) + q_fft_reordered.shape

    kelvin_q_fft = xr.zeros_like(q_fft_reordered, dtype="complex")
    q_wave_components = xr.DataArray(
        data=np.zeros(shape, dtype="complex"),
        dims=dims,
        coords={"mode": normalization_factors, **q_fft_reordered.coords},
    )
    v_wave_components = xr.DataArray(
        data=np.zeros(shape, dtype="complex"),
        dims=dims,
        coords={"mode": normalization_factors, **q_fft_reordered.coords},
    )
    r_wave_components = xr.DataArray(
        data=np.zeros(shape, dtype="complex"),
        dims=dims,
        coords={"mode": normalization_factors, **q_fft_reordered.coords},
    )

    # Calculate wave components
    for mode in range(normalization_factors_shape):
        if mode == 0:
            # Kelvin wave
            part_neg = (
                q_fft_reordered.isel(
                    freq=slice(f_neg_start, f_neg_end), k=slice(k_pos_start, k_pos_end)
                )
                * cylinder_functions.isel(mode=mode)
                * d_lat
            )
            part_neg = part_neg.sum(dim=lat_dim)
            part_pos = (
                q_fft_reordered.isel(
                    freq=slice(f_pos_start, f_pos_end), k=slice(k_neg_start, k_neg_end)
                )
                * cylinder_functions.isel(mode=mode)
                * d_lat
            )
            part_pos = part_pos.sum(dim=lat_dim)
            part_below = cylinder_functions["mode"].isel(mode=mode) / trapping_scale_deg
            kelvin_q_fft[
                {
                    "freq": slice(f_neg_start, f_neg_end),
                    "k": slice(k_pos_start, k_pos_end),
                }
            ] = (
                part_neg / part_below
            )
            kelvin_q_fft[
                {
                    "freq": slice(f_pos_start, f_pos_end),
                    "k": slice(k_neg_start, k_neg_end),
                }
            ] = (
                part_pos / part_below
            )

        # Others
        part_neg = (
            q_fft_reordered.isel(
                freq=slice(f_neg_start, f_neg_end), k=slice(k_neg_start, k_neg_end)
            )
            * cylinder_functions.isel(mode=mode)
            * d_lat
        )
        part_neg = part_neg.sum(dim=lat_dim)
        part_pos = (
            q_fft_reordered.isel(
                freq=slice(f_pos_start, f_pos_end), k=slice(k_pos_start, k_pos_end)
            )
            * cylinder_functions.isel(mode=mode)
            * d_lat
        )
        part_pos = part_pos.sum(dim=lat_dim)
        part_below = cylinder_functions["mode"].isel(mode=mode) / trapping_scale_deg
        q_wave_components[
            {
                "mode": mode,
                "freq": slice(f_neg_start, f_neg_end),
                "k": slice(k_neg_start, k_neg_end),
            }
        ] = (
            part_neg / part_below
        )
        q_wave_components[
            {
                "mode": mode,
                "freq": slice(f_pos_start, f_pos_end),
                "k": slice(k_pos_start, k_pos_end),
            }
        ] = (
            part_pos / part_below
        )

        part_neg = (
            v_fft_reordered.isel(
                freq=slice(f_neg_start, f_neg_end), k=slice(k_neg_start, k_neg_end)
            )
            * cylinder_functions.isel(mode=mode)
            * d_lat
        )
        part_neg = part_neg.sum(dim=lat_dim)
        part_pos = (
            v_fft_reordered.isel(
                freq=slice(f_pos_start, f_pos_end), k=slice(k_pos_start, k_pos_end)
            )
            * cylinder_functions.isel(mode=mode)
            * d_lat
        )
        part_pos = part_pos.sum(dim=lat_dim)
        part_below = cylinder_functions["mode"].isel(mode=mode) / trapping_scale_deg
        v_wave_components[
            {
                "mode": mode,
                "freq": slice(f_neg_start, f_neg_end),
                "k": slice(k_neg_start, k_neg_end),
            }
        ] = (
            part_neg / part_below
        )
        v_wave_components[
            {
                "mode": mode,
                "freq": slice(f_pos_start, f_pos_end),
                "k": slice(k_pos_start, k_pos_end),
            }
        ] = (
            part_pos / part_below
        )

        part_neg = (
            r_fft_reordered.isel(
                freq=slice(f_neg_start, f_neg_end), k=slice(k_neg_start, k_neg_end)
            )
            * cylinder_functions.isel(mode=mode)
            * d_lat
        )
        part_neg = part_neg.sum(dim=lat_dim)
        part_pos = (
            r_fft_reordered.isel(
                freq=slice(f_pos_start, f_pos_end), k=slice(k_pos_start, k_pos_end)
            )
            * cylinder_functions.isel(mode=mode)
            * d_lat
        )
        part_pos = part_pos.sum(dim=lat_dim)
        part_below = cylinder_functions["mode"].isel(mode=mode) / trapping_scale_deg
        r_wave_components[
            {
                "mode": mode,
                "freq": slice(f_neg_start, f_neg_end),
                "k": slice(k_neg_start, k_neg_end),
            }
        ] = (
            part_neg / part_below
        )
        r_wave_components[
            {
                "mode": mode,
                "freq": slice(f_pos_start, f_pos_end),
                "k": slice(k_pos_start, k_pos_end),
            }
        ] = (
            part_pos / part_below
        )

    # Reconstruct physical space fields for each wave type
    dims = ("wave_type",) + q_fft_reordered.dims
    shape = (normalization_factors_shape,) + q_fft_reordered.shape
    wave_types = np.array(["kelvin", "wmrg", "r1", "r2"])

    u_wave_fields = xr.DataArray(
        data=np.zeros(shape, dtype="complex"),
        dims=dims,
        coords={"wave_type": wave_types, **q_fft_reordered.coords},
    )
    v_wave_fields = xr.DataArray(
        data=np.zeros(shape, dtype="complex"),
        dims=dims,
        coords={"wave_type": wave_types, **q_fft_reordered.coords},
    )
    z_wave_fields = xr.DataArray(
        data=np.zeros(shape, dtype="complex"),
        dims=dims,
        coords={"wave_type": wave_types, **q_fft_reordered.coords},
    )

    for wave_idx, wave_type in enumerate(wave_types):
        if wave_type == "kelvin":
            u_wave_fields[{"wave_type": wave_idx}] = (
                0.5 * kelvin_q_fft * cylinder_functions.isel(mode=0)
            )
            z_wave_fields[{"wave_type": wave_idx}] = (
                0.5 * kelvin_q_fft * cylinder_functions.isel(mode=0) * ce_over_gravity
            )

        elif wave_type == "wmrg":
            u_wave_fields[{"wave_type": wave_idx}] = (
                0.5 * q_wave_components.isel(mode=1) * cylinder_functions.isel(mode=1)
            )
            v_wave_fields[{"wave_type": wave_idx}] = (
                0.5 * v_wave_components.isel(mode=0) * cylinder_functions.isel(mode=0)
            )
            z_wave_fields[{"wave_type": wave_idx}] = (
                0.5
                * q_wave_components.isel(mode=1)
                * cylinder_functions.isel(mode=1)
                * ce_over_gravity
            )

        elif wave_type == "r1":
            u_wave_fields[{"wave_type": wave_idx}] = 0.5 * (
                q_wave_components.isel(mode=2) * cylinder_functions.isel(mode=2)
                - r_wave_components.isel(mode=0) * cylinder_functions.isel(mode=0)
            )
            v_wave_fields[{"wave_type": wave_idx}] = (
                0.5 * v_wave_components.isel(mode=1) * cylinder_functions.isel(mode=1)
            )
            z_wave_fields[{"wave_type": wave_idx}] = (
                0.5
                * (
                    q_wave_components.isel(mode=2) * cylinder_functions.isel(mode=2)
                    + r_wave_components.isel(mode=0) * cylinder_functions.isel(mode=0)
                )
                * ce_over_gravity
            )

        elif wave_type == "r2":
            u_wave_fields[{"wave_type": wave_idx}] = 0.5 * (
                q_wave_components.isel(mode=3) * cylinder_functions.isel(mode=3)
                - r_wave_components.isel(mode=1) * cylinder_functions.isel(mode=1)
            )
            v_wave_fields[{"wave_type": wave_idx}] = (
                0.5 * v_wave_components.isel(mode=2) * cylinder_functions.isel(mode=2)
            )
            z_wave_fields[{"wave_type": wave_idx}] = (
                0.5
                * (
                    q_wave_components.isel(mode=3) * cylinder_functions.isel(mode=3)
                    + r_wave_components.isel(mode=1) * cylinder_functions.isel(mode=1)
                )
                * ce_over_gravity
            )

    # Inverse FFT to get physical fields
    u_wave = ifft2_along_dims(
        u_wave_fields, dims=("freq", "k"), original_data=zonal_wind_speed_data
    )
    v_wave = ifft2_along_dims(
        v_wave_fields, dims=("freq", "k"), original_data=meridional_wind_speed_data
    )
    z_wave = ifft2_along_dims(z_wave_fields, dims=("freq", "k"), original_data=z_data)

    results = xr.Dataset()
    u_wave.attrs["long_name"] = "zonal wind"
    u_wave.attrs["units"] = "m/s"
    results["u"] = u_wave

    v_wave.attrs["long_name"] = "zonal wind"
    v_wave.attrs["units"] = "m/s"
    results["v"] = v_wave

    z_wave.attrs["long_name"] = "geopotential height"
    z_wave.attrs["units"] = "m"
    results["z"] = z_wave
    return results

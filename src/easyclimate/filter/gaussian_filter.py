"""
Gaussian filter
"""

from __future__ import annotations
import numpy as np
import xarray as xr
from scipy.ndimage import gaussian_filter1d


__all__ = ["calc_gaussian_filter"]


def calc_gaussian_filter(
    data: xr.DataArray,
    window_length: float,
    sigma: float = None,
    dim: str = "time",
    keep_attrs: bool = False,
) -> xr.DataArray:
    """
    Apply a Gaussian filter to data along a specified dimension.

    Parameters
    ----------
    da : :py:class:`xarray.DataArray<xarray.DataArray>`
        Input data array
    window_length: :py:class:`int <int>`, optional
        The window width.
    sigma : :py:class:`float <float>`, optional
        Standard deviation for Gaussian kernel. If None, calculated as :math:`\\mathrm{window_length} / \\sqrt{8 log2}`.
    dim : :py:class:`str <str>`, optional
        Dimension along which to filter (default: 'time')
    keep_attrs : :py:class:`bool <bool>`, optional
        Whether to preserve attributes (default: False)

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`
        Smoothed data with the same dimensions as input.

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_ao_index.py
    """
    # Calculate sigma if not provided
    if sigma is None:
        # np.sqrt(8*np.log(2)) \\approx 2.3548
        sigma = window_length / np.sqrt(
            8 * np.log(2)
        )  # Convert FWHM (window; Full Width at Half Maximum) to sigma (math parameter)

    # Apply Gaussian filter along the specified dimension
    smoothed_data = gaussian_filter1d(
        data.values, sigma=sigma, axis=data.get_axis_num(dim)
    )

    # Create new DataArray with smoothed values
    result = xr.DataArray(
        smoothed_data,
        dims=data.dims,
        coords=data.coords,
        attrs=data.attrs if keep_attrs else None,
    )

    return result

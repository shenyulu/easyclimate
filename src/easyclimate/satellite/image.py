"""
Satellite Map
"""

from __future__ import annotations
import numpy as np
from scipy.interpolate import interp1d
import xarray as xr

__all__ = [
    "get_stretched_rgb_data",
]


def get_stretched_rgb_data(
    data_input: xr.DataArray, r_band: str, g_band: str, b_band: str
) -> xr.DataArray:
    """
    Extract and process RGB bands from an xarray DataArray to create a stretched RGB composite.

    This function takes three spectral bands from a multispectral dataset, applies contrast stretching
    with a predefined interpolation curve, and combines them into an RGB composite DataArray.

    Parameters:
    -----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`.
        Input DataArray containing multiple spectral bands (must include red, green, and blue channels).
    r_band : :py:class:`str <str>`.
        Name of the band to use for the Red channel.
    g_band : :py:class:`str <str>`.
        Name of the band to use for the Green channel.
    b_band : :py:class:`str <str>`.
        Name of the band to use for the Blue channel.

    Returns:
    --------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
        A 3-band DataArray (RGB) with stretched values (0-255) and 'band' dimension
        containing ``'r'``, ``'g'``, ``'b'`` coordinates

    Notes:
    ------
    The stretching function applies a piecewise linear interpolation with the following breakpoints:

    - Input: ``[0, 30, 60, 120, 190, 255]``
    - Output: ``[0, 110, 160, 210, 240, 255]``

    Values are first normalized to 0-255 range before stretching.
    """

    def stretch(data: np.ndarray | xr.DataArray):
        """
        Apply contrast stretching to input data (numpy array or xarray DataArray).

        The stretching includes:
        1. Normalization to 0-255 range
        2. Piecewise linear interpolation with predefined breakpoints
        3. Conversion to uint8

        Parameters
        -----------
        data : Union[np.ndarray, xr.DataArray]
            Input data to be stretched

        Returns
        --------
        Union[np.ndarray, xr.DataArray]
            Stretched data in the same format as input
        """
        # Extract numpy array if input is DataArray
        if isinstance(data, xr.DataArray):
            data_values = data.values
        else:
            data_values = data

        # Normalize to 0-255 range (assuming input is 0-1)
        data_norm = (data_values - 0) / (1 - 0) * 255

        # Define stretching curve breakpoints
        x = [0, 30, 60, 120, 190, 255]  # Input values
        y = [0, 110, 160, 210, 240, 255]  # Output values
        interp = interp1d(x, y, bounds_error=False, fill_value=255)

        # Apply interpolation and convert to 8-bit unsigned integer
        stretched_data = interp(data_norm).astype(np.uint8)

        # Return DataArray if input was DataArray, otherwise return numpy array
        if isinstance(data, xr.DataArray):
            return xr.DataArray(
                stretched_data, dims=data.dims, coords=data.coords, attrs=data.attrs
            )
        else:
            return stretched_data

    # Extract the three bands from input DataArray
    R = data_input[r_band]
    G = data_input[g_band]
    B = data_input[b_band]

    # Apply stretching to each band
    R_stretched = stretch(R)
    G_stretched = stretch(G)
    B_stretched = stretch(B)

    # Add band dimension and coordinate to each channel
    R_stretched = R_stretched.assign_coords({"band": "r"}).expand_dims("band")
    G_stretched = G_stretched.assign_coords({"band": "g"}).expand_dims("band")
    B_stretched = B_stretched.assign_coords({"band": "b"}).expand_dims("band")

    # Combine all bands into single RGB DataArray
    RGB_stretched = xr.concat([R_stretched, G_stretched, B_stretched], dim="band")
    RGB_stretched.name = "rgb"
    return RGB_stretched

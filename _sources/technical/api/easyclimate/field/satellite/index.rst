easyclimate.field.satellite
===========================

.. py:module:: easyclimate.field.satellite


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/field/satellite/image/index


Functions
---------

.. autoapisummary::

   easyclimate.field.satellite.get_stretched_rgb_data


Package Contents
----------------

.. py:function:: get_stretched_rgb_data(data_input: xarray.DataArray, r_band: str, g_band: str, b_band: str) -> xarray.DataArray

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



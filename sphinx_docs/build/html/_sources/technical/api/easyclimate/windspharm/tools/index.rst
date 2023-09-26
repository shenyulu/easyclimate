:py:mod:`easyclimate.windspharm.tools`
======================================

.. py:module:: easyclimate.windspharm.tools

.. autoapi-nested-parse::

   Tools for managing data for use with `~windspharm.standard.VectorWind`
   (or indeed `spharm.Spharmt`).



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.windspharm.tools.__order_dims
   easyclimate.windspharm.tools.__reshape
   easyclimate.windspharm.tools.prep_data
   easyclimate.windspharm.tools.recover_data
   easyclimate.windspharm.tools.get_recovery
   easyclimate.windspharm.tools.reverse_latdim
   easyclimate.windspharm.tools.order_latdim



Attributes
~~~~~~~~~~

.. autoapisummary::

   easyclimate.windspharm.tools.__recover_docstring_template


.. py:function:: __order_dims(d, inorder)


.. py:function:: __reshape(d)


.. py:function:: prep_data(data, dimorder)

   Prepare data for input to `~windspharm.standard.VectorWind` (or to
   `spharm.Spharmt` method calls).

   Returns a dictionary of intermediate information that can be passed
   to `recover_data` or `get_recovery` to recover the original shape
   and order of the data.

   **Arguments:**

   *data*
       Data array. The array must be at least 2D.

   *dimorder*
       String specifying the order of dimensions in the data array. The
       characters 'x' and 'y' represent longitude and latitude
       respectively. Any other characters can be used to represent
       other dimensions.

   **Returns:**

   *pdata*
       *data* reshaped/reordered to (latitude, longitude, other).

   *info*
       A dictionary of information required to recover *data*.

   **See also:**

   `recover_data`, `get_recovery`.

   **Examples:**

   Prepare an array with dimensions (12, 17, 73, 144) where the
   dimensions are (time, level, latitude, longitude)::

       pdata, info = prep_data(data, 'tzyx')

   Prepare an array with dimensions (144, 16, 73, 21) where the first
   dimension is longitude and the third dimension is latitude. The
   characters used to represent the other dimensions are arbitrary::

       pdata, info = prep_data(data, 'xayb')



.. py:function:: recover_data(pdata, info)

   Recover the shape and dimension order of an array output from
   `~windspharm.standard.VectorWind` methods (or from `spharm.Spharmt`
   methods).

   This function performs the opposite of `prep_data`.

   For recovering the shape of multiple variables, see `get_recovery`.

   **Arguments:**

   *pdata*
       Data array with either 2 or 3 dimensions. The first two
       dimensions are latitude and longitude respectively.

   *info*
       Information dictionary output from `prep_data`.

   **Returns:**

   *data*
       The data reshaped/reordered.

   **See also:**

   `prep_data`, `get_recovery`.

   **Example:**

   Recover the original input shape and dimension order of an array
   processed with `prep_data` or an output of
   `~windspharm.standard.VectorWind` or `sparm.Spharmt` method calls on
   such data::

       data = recover_data(pdata, info)



.. py:data:: __recover_docstring_template
   :value: Multiline-String

    .. raw:: html

        <details><summary>Show Value</summary>

    .. code-block:: python

        """Shape/dimension recovery.
        
        Recovers variable shape/dimension according to:
        
        {!s}
        
        Returns a `list` of variables.
        
        """

    .. raw:: html

        </details>

   

.. py:function:: get_recovery(info)

   Return a function that can be used to recover the shape and
   dimension order of multiple arrays output from
   `~windspharm.standard.VectorWind` methods (or from `spharm.Spharmt`
   methods) according to a single dictionary of recovery information.

   **Argument:**

   *info*
       Information dictionary output from `prep_data`.

   **Returns:**

   *recover*
       A function used to recover arrays.

   **See also:**

   `recover_data`, `prep_data`.

   **Example:**

   Generate a function to recover the original input shape and
   dimension order of arrays processed with `prep_data` and outputs of
   `~windspharm.standard.VectorWind` method calls on this data::

       u, info = prep_data(u, 'tzyx')
       v, info = prep_data(v, 'tzyx')
       w = VectorWind(u, v)
       sf, vp = w.sfvp()
       recover = get_recovery(info)
       u, v, sf, vp = recover(u, v, sf, vp)



.. py:function:: reverse_latdim(u, v, axis=0)

   Reverse the order of the latitude dimension of zonal and meridional
   wind components.

   **Arguments:**

   *u*, *v*
       Zonal and meridional wind components respectively.

   **Optional argument:**

   *axis*
       Index of the latitude dimension. This dimension will be reversed
       in the input arrays. Defaults to 0 (the first dimension).

   **Returns:**

   *ur*, *vr*
       Zonal and meridional wind components with the latitude dimensions
       reversed. These are always copies of the input.

   **See also:**

   `order_latdim`.

   **Examples:**

   Reverse the dimension corresponding to latitude when it is the first
   dimension of the inputs::

       u, v = reverse_latdim(u, v)

   Reverse the dimension corresponding to latitude when it is the third
   dimension of the inputs::

       u, v = reverse_latdim(u, v, axis=2)



.. py:function:: order_latdim(latdim, u, v, axis=0)

   Ensure the latitude dimension is north-to-south.

   Returns copies of the latitude dimension and wind components
   with the latitude dimension going from north to south. If the
   latitude dimension is already in this order then the output will
   just be copies of the input.

   **Arguments:**

   *latdim*
       Array of latitude values.

   *u*, *v*
       Zonal and meridional wind components respectively.

   **Keyword argument:**

   *axis*
       Index of the latitude dimension in the zonal and meridional wind
       components. Defaults to 0 (the first dimension).

   **Returns:**

   *latdimr*
       Possibly reversed *latdim*, always a copy of *latdim*.

   *ur*, *vr*
       Possibly reversed *u* and *v* respectively. Always copies of *u*
       and *v* respectively.

   **See also:**

   `reverse_latdim`.

   **Examples:**

   Order the latitude dimension when latitude is the first dimension of
   the wind components::

       latdim, u, v = order_latdim(latdim, u, v)

   Order the latitude dimension when latitude is the third dimension of
   the wind components::

       latdim, u, v = order_latdim(latdim, u, v, axis=2)




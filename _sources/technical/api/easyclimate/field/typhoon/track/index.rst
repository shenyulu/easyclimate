easyclimate.field.typhoon.track
===============================

.. py:module:: easyclimate.field.typhoon.track

.. autoapi-nested-parse::

   Track cyclone center



Functions
---------

.. autoapisummary::

   easyclimate.field.typhoon.track.track_cyclone_center_msl_only


Module Contents
---------------

.. py:function:: track_cyclone_center_msl_only(msl_data: xarray.DataArray, sample_point: Tuple[float, float], index_value: List[int] | List[float] = [0], lon_dim: str = 'lon', lat_dim: str = 'lat') -> pandas.DataFrame

   Tracks the center of a cyclone using biquadratic interpolation on mean sea level pressure (MSL) data.

   This function identifies local minima in the MSL data using a minimum filter with cyclic boundary conditions
   in longitude, selects the minimum closest to the provided sample point in geodesic distance, and applies
   biquadratic interpolation to estimate the precise location of the cyclone center. If interpolation fails,
   the grid-based minimum is used.

   This is a simple approach used by the author in `pytrack <https://github.com/tenomoto/pytrack>`__ to identify local minima in sea-level pressure.

   Local minima are identified using :py:func:`scipy.ndimage.minimum_filter() <scipy.ndimage.minimum_filter>`.
   To ensure consistency with the subsequent quadratic interpolation, we search for minima within a :math:`3 \cdot 3` grid.
   Since the region is cropped, the ``mode`` parameter is set to ``nearest``,
   which extends the boundary values for both latitude and longitude dimensions.

   Given an estimated position, e.g., :math:`\lambda = 140^\circ, \phi = 20^\circ` in ``sample_point``,
   we calculate the great-circle distance :math:`d = a\alpha` (:math:`a` is the Earth's radius) to the identified local minima using the formula:

   .. math::

       \cos \alpha = \sin\theta_0 \sin\theta + \cos\theta_0 \cos\theta \cos(\lambda - \lambda_0).

   Since the relative magnitude remains unchanged when comparing the central angle :math:`\alpha` on a unit sphere, the Earth's radius is omitted.

   The location of the local minima closest to the given longitude and latitude is identified.
   This point and its eight neighboring points, totaling nine points, are stored in a one-dimensional array.
   The center is indexed as 0, and the points are stored counterclockwise starting from the bottom-left corner.

   Let :math:`f` be a quadratic function of :math:`x` and :math:`y`:

   .. math::

       f(x, y) = c_0 + c_1x + c_2y + c_3xy + c_4x^2 + c_5y^2 + c_6x^2y + c_7xy^2 + c_8x^2y^2.

   A necessary condition for the quadratic function to have an extremum is that its gradient is zero.
   At grid points, assume :math:`f(x_0, y_0)` is a local minimum.
   The extremum of :math:`f(x, y)` may not necessarily lie on a grid point.
   Suppose the extremum of :math:`f(x, y)` is at :math:`x_0 + \Delta x, y_0 + \Delta y`. Expanding around :math:`x_0, y_0` using a Taylor series gives:

   .. math::

       f(x, y) = f(x_0, y_0) + f_x\Delta x + f_y\Delta y + \frac{1}{2}f_{xx}(\Delta x)^2 + \frac{1}{2}f_{yy}(\Delta y)^2 + f_{xy}\Delta x\Delta y.


   Define:

   .. math::

       \begin{align}
       \mathbf{b} &= -\begin{bmatrix} f_x \ f_y \end{bmatrix}, \\
       \mathbf{A} &= \begin{bmatrix} f_{xx} & f_{xy} \\ f_{xy} & f_{yy} \end{bmatrix}, \\
       \mathbf{x} &= \begin{bmatrix} \Delta x \\ \Delta y \end{bmatrix}.
       \end{align}

   Then:

   .. math::
       \begin{align}
       f(x, y) &= f(x_0, y_0) - \mathbf{b}^T\mathbf{x} + \frac{1}{2}\mathbf{x}^T\mathbf{A}\mathbf{x}, \\
       \nabla f &= \mathbf{A}\mathbf{x} - \mathbf{b} = 0, \\
       \mathbf{x} &= \mathbf{A}^{-1}\mathbf{b}.
       \end{align}

   When :math:`d \equiv f_{xx}f_{yy} - f_{xy}^2 \ne 0`:

   .. math::

       \mathbf{A}^{-1} = \frac{1}{d} \begin{pmatrix} f_{yy} & -f_{xy} \\ -f_{xy} & f_{xx} \end{pmatrix},

   allowing :math:`\Delta x` and :math:`\Delta y` to be determined, thus locating the extremum.

   Parameters
   ----------
   msl_data : :py:class:`xarray.DataArray<xarray.DataArray>`
       Input mean sea level pressure data with latitude and longitude dimensions.
   sample_point : :py:class:`tuple[float, float] <tuple>`
       Initial guess for the cyclone center as (longitude, latitude) in degrees.
   index_value : :py:class:`List[int] | List[float] <list>`, optional
       Index value(s) for the output DataFrame, default is [0].
   lon_dim : :py:class:`str <str>`, optional
       Name of the longitude dimension, default is 'lon'.
   lat_dim : :py:class:`str <str>`, optional
       Name of the latitude dimension, default is 'lat'.

   Returns
   -------
   :py:class:`pandas.DataFrame<pandas.DataFrame>`
       DataFrame containing the longitude, latitude, and minimum MSL pressure of the cyclone center.

   .. seealso::

       - https://github.com/tenomoto/pytrack
       - 台風: https://www.dpac.dpri.kyoto-u.ac.jp/enomoto/pymetds/Typhoon.html

   Example
   -------
   >>> import xarray as xr
   >>> import numpy as np
   >>> slp = xr.DataArray(np.random.rand(20, 30), dims=['lat', 'lon'],
   ...                    coords={'lat': np.linspace(-10, 10, 20), 'lon': np.linspace(100, 130, 30)})
   >>> result = track_cyclone_center_msl_only(slp, (110, 0), index_value = [0])

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_tc_track_axis.py



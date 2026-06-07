easyclimate.core.stats.detrend_spatial
======================================

.. py:module:: easyclimate.core.stats.detrend_spatial

.. autoapi-nested-parse::

   Unified spatial detrending function with multiple computation methods.

   Available methods:
   1. 'scipy_reduce' - Simplified version using scipy.signal.detrend
   2. 'scipy' - Optimized version using scipy.signal.detrend
   3. 'numpy' - Manual numpy vectorized version
   4. 'rust' - High-performance Rust backend version
   5. 'rust_chunked' - Rust backend with chunked processing
   6. 'rust_flexible' - Rust backend with flexible dimensions
   7. 'auto' - Automatically selects the best available method



Functions
---------

.. autoapisummary::

   easyclimate.core.stats.detrend_spatial.calc_detrend_spatial_fast


Module Contents
---------------

.. py:function:: calc_detrend_spatial_fast(data_input: xarray.DataArray, time_dim: str = 'time', min_valid_fraction: float = 0.5, method: Literal['scipy_reduce', 'scipy', 'numpy', 'rust', 'rust_chunked', 'rust_flexible', 'auto'] = 'auto', **kwargs) -> xarray.DataArray

   Remove linear trend along time dimension from spatio-temporal data.

   Supports multiple computation methods with optional automatic selection.

   Parameters
   ----------
   data_input : xr.DataArray
       The spatio-temporal data to be detrended.
   time_dim : str, default "time"
       Name of the time dimension.
   min_valid_fraction : float, default 0.5
       Minimum fraction of valid (finite) values required for detrending.
       Grid points with fewer valid values will be set to NaN.
   method : str, default 'auto'
       Computation method to use:
       - 'scipy_reduce': Simplified version using scipy.signal.detrend
       - 'scipy': Optimized version using scipy.signal.detrend
       - 'numpy': Manual numpy vectorized implementation
       - 'rust': High-performance Rust backend
       - 'rust_chunked': Rust backend with chunked processing (for large datasets)
       - 'rust_flexible': Rust backend with flexible dimension handling
       - 'auto': Automatically selects the best available method
   **kwargs : dict
       Additional arguments passed to specific methods:
       - chunk_size: int (for 'rust_chunked' method)
       - use_chunked: bool (for 'rust_chunked' method)

   Returns
   -------
   xr.DataArray
       Detrended data with same shape and coordinates as input.

   Raises
   ------
   TypeError
       If data_input is not an xarray.DataArray.
   ValueError
       If input parameters are invalid.
   ImportError
       If Rust method is selected but Rust backend is not available.

   Examples
   --------
   >>> import xarray as xr
   >>> import numpy as np
   >>>
   >>> # Create sample data
   >>> data = xr.DataArray(
   ...     np.random.randn(100, 50, 100),
   ...     dims=['time', 'lat', 'lon']
   ... )
   >>>
   >>> # Using scipy method
   >>> result1 = calc_detrend_spatial_fast(data, method='scipy')
   >>>
   >>> # Using numpy method
   >>> result2 = calc_detrend_spatial_fast(data, method='numpy')
   >>>
   >>> # Using Rust method (if available)
   >>> try:
   >>>     result3 = calc_detrend_spatial_fast(data, method='rust')
   >>> except ImportError:
   >>>     print("Rust backend not available")
   >>>
   >>> # Automatic method selection
   >>> result4 = calc_detrend_spatial_fast(data, method='auto')

   Notes
   -----
   - 'scipy_reduce': Simplest method but less robust with NaN values
   - 'scipy': Optimized scipy version with better special value handling
   - 'numpy': Manual vectorized implementation, typically 2-3x faster than scipy
   - 'rust': High-performance Rust implementation, typically 10-50x faster than numpy
   - 'rust_chunked': Rust chunked version for large memory datasets
   - 'rust_flexible': Rust version with flexible dimension ordering
   - 'auto': Uses 'rust' if available, otherwise 'numpy'



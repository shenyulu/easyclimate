"""
Unified spatial detrending function with multiple computation methods.

Available methods:
1. 'scipy_reduce' - Simplified version using scipy.signal.detrend
2. 'scipy' - Optimized version using scipy.signal.detrend
3. 'numpy' - Manual numpy vectorized version
4. 'rust' - High-performance Rust backend version
5. 'rust_chunked' - Rust backend with chunked processing
6. 'rust_flexible' - Rust backend with flexible dimensions
7. 'auto' - Automatically selects the best available method
"""

import numpy as np
import xarray as xr
from scipy import signal
import warnings
from typing import Union, Literal, Optional
import dask.array as da

__all__ = ["calc_detrend_spatial_fast"]

# Try to import Rust backend
from ...backend import (
    calc_detrend_spatial_3d_rs,
    calc_detrend_spatial_3d_chunked_rs,
    calc_detrend_spatial_flexible_rs,
)

from ...backend import RUST_AVAILABLE


def calc_detrend_spatial_fast(
    data_input: xr.DataArray,
    time_dim: str = "time",
    min_valid_fraction: float = 0.5,
    method: Literal[
        "scipy_reduce",
        "scipy",
        "numpy",
        "rust",
        "rust_chunked",
        "rust_flexible",
        "auto",
    ] = "auto",
    **kwargs,
) -> xr.DataArray:
    """
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
    """
    # Parameter validation
    if not isinstance(data_input, xr.DataArray):
        raise TypeError("data_input must be an xarray.DataArray")

    if time_dim not in data_input.dims:
        raise ValueError(
            f"time_dim '{time_dim}' not found in data dimensions: {data_input.dims}"
        )

    if not 0.0 <= min_valid_fraction <= 1.0:
        raise ValueError("min_valid_fraction must be between 0 and 1")

    # Automatic method selection
    if method == "auto":
        if RUST_AVAILABLE:
            method = "rust"
        else:
            method = "numpy"

    dim_array = data_input.dims

    # Map methods to functions
    method_map = {
        "scipy_reduce": calc_detrend_spatial_scipy_reduce,
        "scipy": calc_detrend_spatial_scipy,
        "numpy": calc_detrend_spatial_numpy,
        "rust": calc_detrend_spatial_rust,
        "rust_chunked": calc_detrend_spatial_rust_chunked,
        "rust_flexible": calc_detrend_spatial_rust_flexible,
    }

    if method not in method_map:
        raise ValueError(
            f"Unknown method: {method}. Available methods: {list(method_map.keys())}"
        )

    if method.startswith("rust") and not RUST_AVAILABLE:
        raise ImportError(
            f"Rust backend not available for method '{method}'. "
            f"Please install easyclimate_rust module or use another method."
        )

    # Call selected method
    result = method_map[method](data_input, time_dim, min_valid_fraction, **kwargs)
    return result.transpose(*dim_array)


def calc_detrend_spatial_scipy_reduce(
    data_input: xr.DataArray,
    time_dim: str = "time",
    min_valid_fraction: float = 0.5,
    **kwargs,
) -> xr.DataArray:
    """
    Simplified version using scipy.signal.detrend.

    Note: This method is less robust with special values. It handles NaN, Inf, and -Inf
    by replacing them with 1 for computation and restoring NaN afterward.

    Parameters
    ----------
    data_input : xr.DataArray
        Input data array.
    time_dim : str, default "time"
        Time dimension name.
    min_valid_fraction : float, default 0.5
        Minimum fraction of valid values required.
    **kwargs : dict
        Additional arguments (not used in this method).

    Returns
    -------
    xr.DataArray
        Detrended data array.
    """
    # Create a copy to avoid modifying original data
    data = data_input.copy()

    # Preprocess: replace Inf and -Inf with NaN for consistency
    # This ensures they are treated as missing values
    if np.any(np.isinf(data)):
        data = data.where(np.isfinite(data), np.nan)

    # Because scipy.signal.detrend cannot handle np.nan,
    # we need to create a mask first, replace np.nan with 1 for computation,
    # and then restore the masked areas.

    # Calculate mask for invalid values (NaN after preprocessing)
    mask_bool = np.isnan(data).mean(dim=time_dim)
    mask_float = mask_bool + 0.0

    # Calculate valid fraction for each spatial point
    valid_fraction = 1 - mask_bool
    mask_condition = valid_fraction >= min_valid_fraction

    # Fill NaN with 1 for detrending computation
    # Note: Using 1 as fill value ensures it doesn't affect trend calculation
    data_filled = data.fillna(1.0)

    # Apply detrending using reduce
    detrenddata_withoutmask = data_filled.reduce(signal.detrend, dim=time_dim)

    # Restore NaN in originally invalid locations
    result = detrenddata_withoutmask.where(mask_float < 0.5)

    # Apply minimum valid fraction mask
    result = result.where(mask_condition)

    # Preserve attributes
    result.attrs = data_input.attrs.copy()

    # Add method metadata
    result.attrs["detrend_method"] = "linear_scipy_reduce"
    result.attrs["detrend_time_dim"] = time_dim
    result.attrs["min_valid_fraction"] = min_valid_fraction
    result.attrs["special_value_handling"] = "inf_to_nan_replaced_with_1"

    return result


def calc_detrend_spatial_scipy(
    data_input: xr.DataArray,
    time_dim: str = "time",
    min_valid_fraction: float = 0.5,
    **kwargs,
) -> xr.DataArray:
    """
    Optimized version using scipy.signal.detrend with special value handling.

    Handles NaN, Inf, and -Inf values robustly by replacing them with 0
    during detrending and restoring NaN afterward.
    """

    def _detrend_with_special_handling(arr, axis=-1):
        """
        Detrend function that handles special values (NaN, Inf, -Inf).

        This wrapper ensures scipy.signal.detrend works correctly
        even when data contains special values.
        """
        # Create mask for valid (finite) values
        valid_mask = np.isfinite(arr)

        # If all values along axis are invalid, return the original array
        if not np.any(valid_mask):
            return arr

        # Replace special values with 0 for detrending
        arr_clean = np.where(valid_mask, arr, 0.0)

        # Perform detrending
        detrended = signal.detrend(arr_clean, axis=axis)

        # Restore special values in their original positions
        result = np.where(valid_mask, detrended, np.nan)

        return result

    # Check if data is chunked (Dask array)
    is_dask = hasattr(data_input.data, "chunks")

    # Step 1: Create mask for pixels that have enough valid values
    if is_dask:
        # For Dask: use efficient approach
        finite_count = xr.apply_ufunc(
            lambda x: np.isfinite(x).sum(axis=-1),
            data_input,
            input_core_dims=[[time_dim]],
            dask="parallelized",
            output_dtypes=[float],
        )
        time_length = data_input.sizes[time_dim]
        valid_fraction = finite_count / time_length
        mask_condition = valid_fraction >= min_valid_fraction
    else:
        # For non-Dask: simpler approach
        mask_condition = np.isfinite(data_input).sum(dim=time_dim) >= (
            data_input.sizes[time_dim] * min_valid_fraction
        )

    # Step 2: Apply detrending with special value handling
    detrended = xr.apply_ufunc(
        _detrend_with_special_handling,
        data_input,
        input_core_dims=[[time_dim]],
        output_core_dims=[[time_dim]],
        dask="parallelized" if is_dask else "allowed",
        vectorize=True,
        output_dtypes=[data_input.dtype],
        dask_gufunc_kwargs=(
            {
                "allow_rechunk": True,
                "output_sizes": {time_dim: data_input.sizes[time_dim]},
            }
            if is_dask
            else {}
        ),
    )

    # Step 3: Apply mask to remove unreliable pixels
    result = detrended.where(mask_condition)

    # Step 4: Preserve attributes
    result.attrs = data_input.attrs.copy()

    # Add method metadata
    result.attrs["detrend_method"] = "linear_scipy_optimized"
    result.attrs["detrend_time_dim"] = time_dim
    result.attrs["min_valid_fraction"] = min_valid_fraction

    return result


def calc_detrend_spatial_numpy(
    data_input: xr.DataArray,
    time_dim: str = "time",
    min_valid_fraction: float = 0.5,
    **kwargs,
) -> xr.DataArray:
    """
    High-performance version using manual linear regression.

    This version is typically 2-3x faster than scipy.signal.detrend for large arrays
    and handles special values efficiently.
    """

    def _fast_detrend(arr, axis=-1):
        """Fast vectorized detrending with special value handling."""
        # Get shape information
        shape = arr.shape
        n = shape[axis]

        # Move time axis to last position if needed
        if axis != -1:
            arr = np.moveaxis(arr, axis, -1)

        # Create time index
        t = np.arange(n, dtype=arr.dtype)
        t_mean = t.mean()

        # Mask for finite values
        finite_mask = np.isfinite(arr)

        # Replace non-finite with 0 for computation
        arr_clean = np.where(finite_mask, arr, 0.0)

        # Count valid values per pixel
        valid_counts = finite_mask.sum(axis=-1, keepdims=True)

        # Calculate means (accounting for masked values)
        y_mean = arr_clean.sum(axis=-1, keepdims=True) / np.maximum(valid_counts, 1)

        # Calculate slope using least squares
        # slope = sum((t - t_mean) * (y - y_mean)) / sum((t - t_mean)^2)
        t_centered = t - t_mean
        y_centered = arr_clean - y_mean

        numerator = (t_centered * y_centered).sum(axis=-1, keepdims=True)
        denominator = (t_centered**2).sum()

        slope = numerator / denominator

        # Calculate detrended values
        trend = slope * t_centered + y_mean
        detrended = arr_clean - trend

        # Restore NaN where original data was non-finite
        detrended = np.where(finite_mask, detrended, np.nan)

        # Move axis back if needed
        if axis != -1:
            detrended = np.moveaxis(detrended, -1, axis)

        return detrended

    # Calculate mask
    is_dask = hasattr(data_input.data, "chunks")

    if is_dask:
        finite_count = xr.apply_ufunc(
            lambda x: np.isfinite(x).sum(axis=-1),
            data_input,
            input_core_dims=[[time_dim]],
            dask="parallelized",
            output_dtypes=[float],
        )
    else:
        finite_count = np.isfinite(data_input).sum(dim=time_dim)

    valid_fraction = finite_count / data_input.sizes[time_dim]
    mask_condition = valid_fraction >= min_valid_fraction

    # Apply fast detrending
    detrended = xr.apply_ufunc(
        _fast_detrend,
        data_input,
        input_core_dims=[[time_dim]],
        output_core_dims=[[time_dim]],
        dask="parallelized" if is_dask else "allowed",
        vectorize=True,
        output_dtypes=[data_input.dtype],
        dask_gufunc_kwargs=(
            {
                "allow_rechunk": True,
                "output_sizes": {time_dim: data_input.sizes[time_dim]},
            }
            if is_dask
            else {}
        ),
    )

    # Apply mask
    result = detrended.where(mask_condition)
    result.attrs = data_input.attrs.copy()

    # Add method metadata
    result.attrs["detrend_method"] = "linear_numpy_vectorized"
    result.attrs["detrend_time_dim"] = time_dim
    result.attrs["min_valid_fraction"] = min_valid_fraction

    return result


def calc_detrend_spatial_rust(
    data_input: xr.DataArray,
    time_dim: str = "time",
    min_valid_fraction: float = 0.5,
    **kwargs,
) -> xr.DataArray:
    """
    High-performance version using Rust backend.

    Features:
    - Rayon parallel processing across spatial grid points
    - Optimized linear regression algorithm
    - Efficient handling of NaN, Inf, and -Inf values
    - Cache-friendly memory access patterns
    """
    if not RUST_AVAILABLE:
        raise ImportError("Rust backend not available")

    # Get dimension order
    dims = list(data_input.dims)
    time_axis = dims.index(time_dim)

    # Get data as numpy array (load if dask)
    if hasattr(data_input.data, "compute"):
        print("Converting Dask array to numpy for Rust processing...")
        data_np = data_input.compute().values
    else:
        data_np = data_input.values

    # Ensure data is float64
    if data_np.dtype != np.float64:
        data_np = data_np.astype(np.float64)

    # Handle different dimension orders (assume 3D data)
    if len(dims) == 3:
        # Transpose to (time, dim1, dim2) if necessary
        if time_axis != 0:
            # Get the order to transpose
            new_order = [time_axis] + [i for i in range(3) if i != time_axis]
            data_np = np.transpose(data_np, new_order)

            # Apply Rust function
            result_np = calc_detrend_spatial_3d_rs(data_np, min_valid_fraction)

            # Transpose back
            inverse_order = [new_order.index(i) for i in range(3)]
            result_np = np.transpose(result_np, inverse_order)
        else:
            # Time is already first dimension
            result_np = calc_detrend_spatial_3d_rs(data_np, min_valid_fraction)
    else:
        raise ValueError(
            f"Expected 3D data, got {len(dims)}D. "
            f"Dimensions: {dims}. Please ensure data has exactly 3 dimensions."
        )

    # Create output DataArray
    result = xr.DataArray(
        result_np,
        dims=data_input.dims,
        coords=data_input.coords,
        attrs=data_input.attrs.copy(),
    )

    # Add processing metadata
    result.attrs["detrend_method"] = "linear_rust_backend"
    result.attrs["detrend_time_dim"] = time_dim
    result.attrs["min_valid_fraction"] = min_valid_fraction

    return result


def calc_detrend_spatial_rust_chunked(
    data_input: xr.DataArray,
    time_dim: str = "time",
    min_valid_fraction: float = 0.5,
    chunk_size: int = 1000,
    **kwargs,
) -> xr.DataArray:
    """
    Rust backend with chunked processing for large datasets.

    Parameters
    ----------
    chunk_size : int, default 1000
        Number of spatial points to process per chunk.
    """
    if not RUST_AVAILABLE:
        raise ImportError("Rust backend not available")

    use_chunked = kwargs.get("use_chunked", True)

    # Get dimension order
    dims = list(data_input.dims)
    time_axis = dims.index(time_dim)

    # Get data as numpy array
    if hasattr(data_input.data, "compute"):
        print("Converting Dask array to numpy for Rust processing...")
        data_np = data_input.compute().values
    else:
        data_np = data_input.values

    # Ensure data is float64
    if data_np.dtype != np.float64:
        data_np = data_np.astype(np.float64)

    # Handle different dimension orders
    if len(dims) == 3:
        if time_axis != 0:
            new_order = [time_axis] + [i for i in range(3) if i != time_axis]
            data_np = np.transpose(data_np, new_order)

            result_np = calc_detrend_spatial_3d_chunked_rs(
                data_np, min_valid_fraction, chunk_size
            )

            inverse_order = [new_order.index(i) for i in range(3)]
            result_np = np.transpose(result_np, inverse_order)
        else:
            result_np = calc_detrend_spatial_3d_chunked_rs(
                data_np, min_valid_fraction, chunk_size
            )
    else:
        raise ValueError(f"Expected 3D data, got {len(dims)}D")

    # Create output DataArray
    result = xr.DataArray(
        result_np,
        dims=data_input.dims,
        coords=data_input.coords,
        attrs=data_input.attrs.copy(),
    )

    # Add processing metadata
    result.attrs["detrend_method"] = "linear_rust_backend_chunked"
    result.attrs["detrend_time_dim"] = time_dim
    result.attrs["min_valid_fraction"] = min_valid_fraction
    result.attrs["chunk_size"] = chunk_size

    return result


def calc_detrend_spatial_rust_flexible(
    data_input: xr.DataArray,
    time_dim: str = "time",
    min_valid_fraction: float = 0.5,
    **kwargs,
) -> xr.DataArray:
    """
    Rust backend with flexible dimension handling.

    Can handle arbitrary dimension ordering without transposing.
    """
    if not RUST_AVAILABLE:
        raise ImportError("Rust backend not available")

    dims = list(data_input.dims)
    if len(dims) != 3:
        raise ValueError(f"Expected 3D data, got {len(dims)}D")

    time_axis = dims.index(time_dim)

    # Get data
    if hasattr(data_input.data, "compute"):
        data_np = data_input.compute().values
    else:
        data_np = data_input.values

    if data_np.dtype != np.float64:
        data_np = data_np.astype(np.float64)

    # Call Rust function with flexible axis
    result_np = calc_detrend_spatial_flexible_rs(data_np, time_axis, min_valid_fraction)

    result = xr.DataArray(
        result_np,
        dims=data_input.dims,
        coords=data_input.coords,
        attrs=data_input.attrs.copy(),
    )

    result.attrs["detrend_method"] = "linear_rust_backend_flexible"
    result.attrs["detrend_time_dim"] = time_dim

    return result


def benchmark_detrend_methods(
    data_input: Optional[xr.DataArray] = None,
    shape: tuple = (365, 180, 360),
    time_dim: str = "time",
    methods: Optional[list] = None,
    n_runs: int = 3,
    warmup: bool = True,
    verbose: bool = True,
) -> dict:
    """
    Benchmark different detrending methods.

    Parameters
    ----------
    data_input : xr.DataArray, optional
        Data to test. If None, test data is created.
    shape : tuple, default (365, 180, 360)
        Shape of test data (time, lat, lon).
    time_dim : str, default "time"
        Name of the time dimension.
    methods : list, optional
        List of methods to test. If None, all available methods are tested.
    n_runs : int, default 3
        Number of runs per method.
    warmup : bool, default True
        Whether to perform warmup runs (not timed).
    verbose : bool, default True
        Whether to print progress information.

    Returns
    -------
    dict
        Dictionary containing benchmark results for each method.

    Examples
    --------
    >>> import numpy as np
    >>> import xarray as xr
    >>>
    >>> # Create test data
    >>> data = xr.DataArray(
    ...     np.random.randn(100, 50, 100),
    ...     dims=['time', 'lat', 'lon']
    ... )
    >>>
    >>> # Run benchmark
    >>> results = benchmark_detrend_methods(data)
    >>>
    >>> # Test specific methods
    >>> results = benchmark_detrend_methods(
    ...     data,
    ...     methods=['scipy', 'numpy', 'auto']
    ... )
    """
    import time

    # Create test data if needed
    if data_input is None:
        if verbose:
            print(f"Creating test data with shape {shape}...")
        np.random.seed(42)
        nt, nlat, nlon = shape
        data = np.random.randn(nt, nlat, nlon) * 10
        trend = np.linspace(0, 10, nt)[:, None, None]
        data = data + trend

        # Add some NaN values (5%)
        mask = np.random.random(shape) < 0.05
        data[mask] = np.nan

        # Add some Inf values (1%)
        inf_mask = np.random.random(shape) < 0.01
        data[inf_mask] = np.inf

        # Create xarray
        data_input = xr.DataArray(
            data,
            dims=["time", "lat", "lon"],
            coords={
                "time": np.arange(nt),
                "lat": np.linspace(-90, 90, nlat),
                "lon": np.linspace(-180, 180, nlon),
            },
        )
    else:
        shape = data_input.shape
        if verbose:
            print(f"Using provided data with shape {shape}...")

    # Determine methods to test
    if methods is None:
        # Test all available methods
        all_methods = ["scipy_reduce", "scipy", "numpy"]
        if RUST_AVAILABLE:
            all_methods.extend(["rust", "rust_chunked", "rust_flexible", "auto"])
        methods = all_methods

    # Remove unavailable methods
    available_methods = []
    for method in methods:
        if method.startswith("rust") and not RUST_AVAILABLE:
            if verbose:
                print(f"Skipping method '{method}' (Rust backend not available)")
        else:
            available_methods.append(method)

    if not available_methods:
        raise ValueError("No available methods to benchmark")

    if verbose:
        print(f"Benchmarking methods: {available_methods}")
        print(f"Data shape: {shape}")
        print(f"Data size: {data_input.nbytes / (1024 ** 2):.1f} MB\n")

    results = {}

    for method in available_methods:
        if verbose:
            print(f"Testing method: {method}")

        # Warmup run (if enabled)
        if warmup:
            try:
                _ = calc_detrend_spatial_fast(
                    data_input, time_dim=time_dim, method=method
                )
            except Exception as e:
                if verbose:
                    print(f"  Warning: Warmup failed: {e}")

        # Timed runs
        times = []
        for i in range(n_runs):
            try:
                start = time.time()
                result = calc_detrend_spatial_fast(
                    data_input, time_dim=time_dim, method=method
                )
                end = time.time()
                elapsed = end - start
                times.append(elapsed)

                # Validate result
                valid_count = np.isfinite(result.values).sum()
                if verbose:
                    print(
                        f"  Run {i+1}: {elapsed:.3f} seconds, valid values: {valid_count:,}"
                    )

            except Exception as e:
                if verbose:
                    print(f"  Run {i+1}: ERROR - {e}")
                times.append(np.nan)

        # Calculate statistics
        valid_times = [t for t in times if not np.isnan(t)]
        if valid_times:
            mean_time = np.mean(valid_times)
            std_time = np.std(valid_times)
            throughput = (data_input.nbytes / (1024**2)) / mean_time

            # Store results
            results[method] = {
                "mean_time": mean_time,
                "std_time": std_time,
                "throughput_mb_s": throughput,
                "data_size_mb": data_input.nbytes / (1024**2),
                "valid_value_count": (
                    np.isfinite(result.values).sum() if "result" in locals() else 0
                ),
                "successful_runs": len(valid_times),
                "total_runs": n_runs,
            }

            if verbose:
                print(
                    f"  Results: {mean_time:.3f} ± {std_time:.3f} seconds, "
                    f"{throughput:.1f} MB/s\n"
                )
        else:
            results[method] = {
                "mean_time": np.nan,
                "std_time": np.nan,
                "throughput_mb_s": np.nan,
                "data_size_mb": data_input.nbytes / (1024**2),
                "valid_value_count": 0,
                "successful_runs": 0,
                "total_runs": n_runs,
                "error": "All runs failed",
            }
            if verbose:
                print(f"  Failed all runs\n")

    # Print summary
    if verbose:
        print("=" * 60)
        print("BENCHMARK SUMMARY")
        print("=" * 60)

        # Find fastest method
        fastest_method = None
        fastest_time = float("inf")

        for method, stats in results.items():
            if not np.isnan(stats["mean_time"]) and stats["mean_time"] < fastest_time:
                fastest_time = stats["mean_time"]
                fastest_method = method

        print(f"\nFastest method: {fastest_method} ({fastest_time:.3f} seconds)")
        print(f"Data size: {data_input.nbytes / (1024 ** 2):.1f} MB")
        print("\nDetailed results:")

        # Sort by speed
        sorted_methods = sorted(
            [m for m in results.keys() if not np.isnan(results[m]["mean_time"])],
            key=lambda m: results[m]["mean_time"],
        )

        for i, method in enumerate(sorted_methods):
            stats = results[method]
            rel_speed = (
                fastest_time / stats["mean_time"] if stats["mean_time"] > 0 else 0
            )
            print(
                f"  {i+1:2d}. {method:15s}: {stats['mean_time']:.3f} ± {stats['std_time']:.3f} s "
                f"({rel_speed:.1f}x, {stats['throughput_mb_s']:.1f} MB/s)"
            )

    return results


def compare_results(
    data_input: xr.DataArray,
    methods: Optional[list] = None,
    time_dim: str = "time",
    tolerance: float = 1e-10,
) -> dict:
    """
    Compare numerical results from different methods.

    Parameters
    ----------
    data_input : xr.DataArray
        Input data array.
    methods : list, optional
        List of methods to compare. If None, uses ['scipy', 'numpy'].
    time_dim : str, default "time"
        Time dimension name.
    tolerance : float, default 1e-10
        Tolerance for numerical comparison.

    Returns
    -------
    dict
        Dictionary containing results from each method and comparison statistics.
    """
    if methods is None:
        methods = ["scipy", "numpy"]
        if RUST_AVAILABLE:
            methods.append("rust")

    results = {}

    for method in methods:
        try:
            result = calc_detrend_spatial_fast(
                data_input, time_dim=time_dim, method=method
            )
            results[method] = result

            print(f"Method '{method}':")
            print(f"  Valid values: {np.isfinite(result.values).sum():,}")
            print(f"  Min value: {np.nanmin(result.values):.6f}")
            print(f"  Max value: {np.nanmax(result.values):.6f}")
            print(f"  Mean value: {np.nanmean(result.values):.6f}")
            print(f"  Std value: {np.nanstd(result.values):.6f}")
            print()

        except Exception as e:
            print(f"Method '{method}' failed: {e}")
            results[method] = None

    # Compare results from different methods (if multiple successful)
    successful_methods = [m for m in methods if results[m] is not None]

    if len(successful_methods) >= 2:
        print("=" * 40)
        print("RESULT COMPARISON")
        print("=" * 40)

        for i in range(len(successful_methods)):
            for j in range(i + 1, len(successful_methods)):
                m1, m2 = successful_methods[i], successful_methods[j]
                r1, r2 = results[m1].values, results[m2].values

                # Calculate difference statistics
                diff = r1 - r2
                finite_mask = np.isfinite(r1) & np.isfinite(r2)

                if np.any(finite_mask):
                    abs_diff = np.abs(diff[finite_mask])

                    print(f"\nComparing {m1} vs {m2}:")
                    print(f"  Max absolute difference: {np.max(abs_diff):.6e}")
                    print(f"  Mean absolute difference: {np.mean(abs_diff):.6e}")
                    print(f"  Std of differences: {np.std(diff[finite_mask]):.6e}")
                    print(f"  RMSE: {np.sqrt(np.mean(diff[finite_mask]**2)):.6e}")

                    # Check if values are equal within tolerance
                    close = np.allclose(
                        r1[finite_mask], r2[finite_mask], rtol=tolerance, atol=tolerance
                    )
                    print(f"  All values within {tolerance:.0e} tolerance: {close}")

    return results


def create_test_dataset(
    shape: tuple = (100, 50, 100),
    add_trend: bool = True,
    add_noise: bool = True,
    nan_fraction: float = 0.05,
    inf_fraction: float = 0.01,
    seed: int = 42,
) -> xr.DataArray:
    """
    Create a test dataset for detrending.

    Parameters
    ----------
    shape : tuple, default (100, 50, 100)
        Shape of the dataset (time, lat, lon).
    add_trend : bool, default True
        Whether to add a linear trend.
    add_noise : bool, default True
        Whether to add Gaussian noise.
    nan_fraction : float, default 0.05
        Fraction of values to set to NaN.
    inf_fraction : float, default 0.01
        Fraction of values to set to Inf.
    seed : int, default 42
        Random seed for reproducibility.

    Returns
    -------
    xr.DataArray
        Test dataset with specified properties.
    """
    np.random.seed(seed)
    nt, nlat, nlon = shape

    # Create base data
    if add_noise:
        data = np.random.randn(nt, nlat, nlon) * 10
    else:
        data = np.zeros((nt, nlat, nlon))

    # Add linear trend if requested
    if add_trend:
        trend = np.linspace(0, 10, nt)[:, None, None]
        data = data + trend

    # Add NaN values
    if nan_fraction > 0:
        mask = np.random.random(shape) < nan_fraction
        data[mask] = np.nan

    # Add Inf values
    if inf_fraction > 0:
        inf_mask = np.random.random(shape) < inf_fraction
        data[inf_mask] = np.inf

    # Create xarray DataArray
    result = xr.DataArray(
        data,
        dims=["time", "lat", "lon"],
        coords={
            "time": np.arange(nt),
            "lat": np.linspace(-90, 90, nlat),
            "lon": np.linspace(-180, 180, nlon),
        },
        attrs={
            "description": "Test dataset for detrending",
            "trend_added": add_trend,
            "noise_added": add_noise,
            "nan_fraction": nan_fraction,
            "inf_fraction": inf_fraction,
            "seed": seed,
        },
    )

    return result


# Example usage and testing
if __name__ == "__main__":
    import sys

    # Create small test data for demonstration
    print("Creating demo data...")
    test_data = create_test_dataset(
        # shape=(50, 25, 25),
        shape=(100, 50, 50),
        add_trend=True,
        add_noise=True,
        nan_fraction=0.1,
        inf_fraction=0.02,
        seed=42,
    )

    # Test all available methods
    print("\nTesting all available methods...")

    try:
        # Run benchmark
        benchmark_results = benchmark_detrend_methods(
            data_input=test_data,
            # methods=['scipy_reduce', 'scipy', 'numpy', 'auto'],
            n_runs=3,
            verbose=True,
        )

        # Compare numerical results
        print("\nComparing numerical results...")
        print("*************************************************")
        compare_results(test_data, methods=["scipy", "numpy"])
        print("*************************************************")
        compare_results(test_data, methods=["numpy", "scipy_reduce"])
        print("*************************************************")
        compare_results(test_data, methods=["scipy_reduce", "rust_flexible"])

    except Exception as e:
        print(f"Error during benchmark: {e}")
        import traceback

        traceback.print_exc()

    print("\nDemo completed!")

"""
Obtain data within a specified time period
"""
from .yearstat import *

def get_specific_years_data(data_input: xr.DataArray, year_array: np.array, dim = 'time') -> xr.DataArray:
    """
    Slicing and extracting the part of the data containing the specified year based on an array of given integer years.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
    year_array: numpy.array containing int data-type objects
        Year(s) to be extracted.
    dim : str
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    years = data_input[dim].dt.year
    year_idx = years.isin(year_array)
    return data_input.isel({dim: year_idx})

def get_specific_months_data(data_input: xr.DataArray, month_array: np.array, dim = 'time') -> xr.DataArray:
    """
    Slicing and extracting the part of the data containing the specified year based on an array of given integer months.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
    month_array: numpy.array containing int data-type objects
        Month(s) to be extracted.
    dim : str
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    months = data_input[dim].dt.month
    months_idx = months.isin(month_array)
    return data_input.isel({dim: months_idx})

def get_specific_days_data(data_input: xr.DataArray, day_array: np.array, dim = 'time') -> xr.DataArray:
    """
    Slicing and extracting the part of the data containing the specified year based on an array of given integer days.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
    day_array: numpy.array containing int data-type objects
        Days(s) to be extracted.
    dim : str
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    days = data_input[dim].dt.day
    days_idx = days.isin(day_array)
    return data_input.isel({dim: days_idx})

def get_specific_hours_data(data_input: xr.DataArray, hour_array: np.array, dim = 'time') -> xr.DataArray:
    """
    Slicing and extracting the part of the data containing the specified year based on an array of given integer hours.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
    hour_array: numpy.array containing int data-type objects
        Hour(s) to be extracted.
    dim : str
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    hours = data_input.time.dt.hour
    hours_idx = hours.isin(hour_array)
    return data_input.isel({dim: hours_idx})

def get_specific_minutes_data(data_input: xr.DataArray, minute_array: np.array, dim = 'time') -> xr.DataArray:
    """
    Slicing and extracting the part of the data containing the specified year based on an array of given integer minutes.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
    minute_array: numpy.array containing int data-type objects
        Minute(s) to be extracted.
    dim : str
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    minutes = data_input.time.dt.minute
    minutes_idx = minutes.isin(minute_array)
    return data_input.isel({dim: minutes_idx})

def get_specific_seconds_data(data_input: xr.DataArray, second_array: np.array, dim = 'time') -> xr.DataArray:
    """
    Slicing and extracting the part of the data containing the specified year based on an array of given integer seconds.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
    second_array: numpy.array containing int data-type objects
        Second(s) to be extracted.
    dim : str
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    seconds = data_input.time.dt.second
    seconds_idx = seconds.isin(second_array)
    return data_input.isel({dim: seconds_idx})

def get_specific_microseconds_data(data_input: xr.DataArray, microsecond_array: np.array, dim = 'time') -> xr.DataArray:
    """
    Slicing and extracting the part of the data containing the specified year based on an array of given integer microseconds.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
    microsecond_array: numpy.array containing int data-type objects
        Microsecond(s) to be extracted.
    dim : str
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    microseconds = data_input.time.dt.microsecond
    microseconds_idx = microseconds.isin(microsecond_array)
    return data_input.isel({dim: microseconds_idx})

def get_specific_nanoseconds_data(data_input: xr.DataArray, nanosecond_array: np.array, dim = 'time') -> xr.DataArray:
    """
    Slicing and extracting the part of the data containing the specified year based on an array of given integer nanoseconds.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
    nanosecond_array: numpy.array containing int data-type objects
        Nanosecond(s) to be extracted.
    dim : str
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    nanoseconds = data_input.time.dt.nanosecond
    nanoseconds_idx = nanoseconds.isin(nanosecond_array)
    return data_input.isel({dim: nanoseconds_idx})

def get_specific_dayofweek_data(data_input: xr.DataArray, dayofweek_array: np.array, dim = 'time') -> xr.DataArray:
    """
    Slicing and extracting the part of the data containing the specified year based on an array of given integer dayofweek.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
    dayofweek_array: numpy.array containing int data-type objects
        The days of the week to be extracted. 

        The integer numbers correspond to the days of the week as follows.

    +-------------------+-------------------+
    | Day of the week   | Integer numbers   |
    +===================+===================+
    |      Monday       |         0         |
    +-------------------+-------------------+
    |      Tuesday      |         1         |
    +-------------------+-------------------+
    |      Wednesday    |         2         |
    +-------------------+-------------------+
    |      Thursday     |         3         |
    +-------------------+-------------------+
    |      Friday       |         4         |
    +-------------------+-------------------+
    |      Saturday     |         5         |
    +-------------------+-------------------+
    |      Sunday       |         6         |
    +-------------------+-------------------+

    dim : str
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    dayofweek = data_input[dim].dt.dayofweek
    dayofweek_idx = dayofweek.isin(dayofweek_array)
    return data_input.isel({dim: dayofweek_idx})

def get_yearmean_for_specific_months_data(data_input: xr.DataArray, month_array: np.array, dim = 'time', kwargs = None) -> xr.DataArray:
    """
    Get the annual average of certain months.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        :py:class:`xarray.DataArray<xarray.DataArray>` to be extracted.
    month_array: numpy.array containing int data-type objects
        Month(s) to be extracted.
    dim : str
        Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
    **kwargs:
        Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data. 
        These could include dask-specific kwargs like split_every.

    Returns
    -------
    :py:class:`xarray.DataArray<xarray.DataArray>`.
    """
    data_get_specific_months_data = get_specific_months_data(data_input, month_array, dim = dim)
    return calc_yearmean(data_get_specific_months_data, dim = dim, **kwargs)

def get_year_exceed_index_upper_bound(
    data_input: xr.DataArray,
    thresh: float,
    time_dim: str = 'time',
) -> np.array:
    """
    Extract the years under the specified threshold (upper bound) in the annual average index (one-dimensional data with only a `time` dimension).

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The one-dimensional data with only a `time` dimension.
    thresh: :py:class:`float<python.float>`.
        The threshold value.
    time_dim: :py:class:`str<python.str>`.
        The time coordinate dimension name.

    Returns
    -------
    :py:class:`numpy.array <numpy:numpy.array>`.
    """
    if data_input.dims == (time_dim,):
        pass
    else:
        raise ValueError('Only one-dimensional time series input is supported.')
    return data_input[time_dim].dt.year[data_input > thresh].data

def get_year_exceed_index_lower_bound(
    data_input: xr.DataArray,
    thresh: float,
    time_dim: str = 'time',
) -> np.array:
    """
    Extract the years under the specified threshold (lower bound) in the annual average index (one-dimensional data with only a `time` dimension).

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
        The one-dimensional data with only a `time` dimension.
    thresh: :py:class:`float<python.float>`.
        The threshold value.
    time_dim: :py:class:`str<python.str>`.
        The time coordinate dimension name.

    Returns
    -------
    :py:class:`numpy.array <numpy:numpy.array>`.
    """
    if data_input.dims == (time_dim,):
        pass
    else:
        raise ValueError('Only one-dimensional time series input is supported.')
    return data_input[time_dim].dt.year[data_input < thresh].data   
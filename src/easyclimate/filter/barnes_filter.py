"""
https://github.com/LinOuyang/pybarnes
"""

from warnings import warn
import xarray as xr
import numpy as np
from metpy.calc import lat_lon_grid_deltas

def field_grids(data, grids):
    """
    For each point of the grids, its nearby region are chosen.

    Parameters
    ----------
    data : array
        An N-dimensional array.
    grids : int or 2-size tuple
        The total number of grid points of the x and y directions.

    Returns
    -------
    D : ndarray
        A ndarray with extra two dimensions. The last two dimensions are
        each points's nearby region.
    """
    radius_grid = (np.array(grids) - 1)// 2
    ny, nx = data.shape[-2:]
    D = np.zeros(data.shape + tuple(grids[::-1]))
    d = np.lib.stride_tricks.sliding_window_view(data, grids[::-1], axis=(-2, -1))
    D[..., radius_grid[1]:-radius_grid[1], radius_grid[0]:-radius_grid[0], :, :] = d
    for i in np.arange(nx)[radius_grid[0]:-radius_grid[0]]:
        idx = slice(i - radius_grid[0], i + radius_grid[0] + 1)
        D[..., :radius_grid[1], i, :, :] = data[..., None, :grids[1], idx]
        D[..., -radius_grid[1]:, i, :, :] = data[..., None, -grids[1]:, idx]

    for i in np.arange(ny)[radius_grid[1]:-radius_grid[1]]:
        idx = slice(i - radius_grid[1], i + radius_grid[1] + 1)
        D[..., i, :radius_grid[0], :, :] = data[..., None, idx, :grids[0]]
        D[..., i, -radius_grid[0]:, :, :] = data[..., None, idx, -grids[0]:]
    D[..., :radius_grid[1], :radius_grid[0], :, :] = data[..., None, None, :grids[1], :grids[0]]
    D[..., :radius_grid[1], -radius_grid[0]:, :, :] = data[..., None, None, :grids[1], -grids[0]:]
    D[..., -radius_grid[1]:, :radius_grid[0], :, :] = data[..., None, None, -grids[1]:, :grids[0]]
    D[..., -radius_grid[1]:, -radius_grid[0]:, :, :] = data[..., None, None, -grids[1]:, -grids[0]:]
    return D


class BarnesFilter:

    """
    The Barnes method performs grid point interpolation by selecting appropriate
    filtering parameters *c* and *g* to filter out shortwave noise in the original field,
    making the analysis results stable and smooth. In addition, it can form a bandpass filter
    to separate various sub weather scales that affect weather processes according to actual needs,
    achieving the purpose of scale separation.

    Reference:
    DOI : https://doi.org/10.1175/1520-0493(1980)108<1108:AOTFSM>2.0.CO;2
    """


    def __init__(self, data_arr, lon=None, lat=None, radius_degree=10):
        """
        Initializing the data and caculate the distance.


        Parameters
        ----------
        data_arr : numpy.array or xarray.Dataset or xarray.DataArray
            An N-dimensional array which to be filtered.
        lon : array
            If the data_arr are numpy.array which has no longitude and latitude infomation,
            then the longitude infomation must be specified.
        lat : array
            If the data_arr are numpy.array which has no longitude and latitude infomation,
            then the latitude infomation must be specified.
        radius_degree : int or tuple
            The radius of each point when caculating the distance of each other.
            Units : degree.

        Returns
        -------
        out : object
            A Barnes filter object.
        """
        self.radius_degree = radius_degree
        if isinstance(data_arr, xr.Dataset):
            self.data_vars = list(data_arr.data_vars)
            self.dims = data_arr[self.data_vars[0]].dims
            self.coords = data_arr.coords
            lat_var, lon_var = self.dims[-2:]
            lon, lat = data_arr[lon_var].data, data_arr[lat_var].data
            self.data = np.stack([data_arr[v] for v in self.data_vars], axis=0)
        elif isinstance(data_arr, xr.DataArray):
            self.dims = data_arr.dims
            self.coords = data_arr.coords
            lat_var, lon_var = self.dims[-2:]
            lon, lat = data_arr[lon_var].data, data_arr[lat_var].data
            self.data = data_arr.data[None]
            self.data_vars = [data_arr.name]
        else:
            if lat is None:
                raise KeyError("The longitude and latitude of the data are missing.")
            self.data = data_arr.copy()
        self.kilometer_distance2 = self.__calculate_distance(lon, lat)

    def __convert_data(self, data):
        if hasattr(self, "dims"):
            nv = data.shape[0]
            ds = []
            for i in range(nv):
                da = xr.DataArray(data[i], coords=self.coords, dims=self.dims, name=self.data_vars[i])
                if nv < 2:
                    return da
                ds.append(da)
            return xr.merge(ds)
        else:
            return data

    def __calculate_distance(self, lon, lat):
        nx, ny = len(lon), len(lat)
        if not np.iterable(self.radius_degree):
            radius_degree = np.array([self.radius_degree]*2)
        if (lon.ndim < 2) & (lat.ndim < 2):
            lon, lat = np.meshgrid(lon, lat)
        mean_delta_lon = np.abs(np.mean(lon[:, :-1] - lon[:, 1:]))
        mean_delta_lat = np.abs(np.mean(lat[:-1] - lat[1:]))
        self.grids = 2 * (radius_degree//np.array([mean_delta_lon, mean_delta_lat])).astype(int) + 1
        radius_grid = (self.grids - 1)// 2

        dx, dy = lat_lon_grid_deltas(lon, lat)
        x = np.array(np.concatenate([np.zeros(ny,)[:, None], np.cumsum(dx, 1)], 1))/1000
        y = np.array(np.concatenate([np.zeros(nx,)[None], np.cumsum(dy, 0)], 0))/1000
        del dx, dy

        gridx = field_grids(x, self.grids) 
        gridy = field_grids(y, self.grids)

        dX = x[:, :, None, None] - gridx
        dY = y[:, :, None, None] - gridy
        del x, y, gridx, gridy

        kilometer_distance2 = dX ** 2 + dY ** 2
        return kilometer_distance2

    def __lowpass(self, g=0.3, c=150000):
        data = field_grids(self.data, self.grids)
        sum_func = np.sum
        if np.isnan(self.data).any():
            warn("The input data contains NANs which may result in unexpected result")
            sum_func = np.nansum
        weights1 = np.exp(-self.kilometer_distance2/(4*c))
        sum_weights1 = np.sum(weights1, (-1, -2), keepdims=True)
        normed_weights1 = weights1/sum_weights1
        del weights1, sum_weights1
        weights2 = np.exp(-self.kilometer_distance2/(4*g*c))
        sum_weights2 = np.sum(weights2, (-1, -2), keepdims=True)
        normed_weights2 = weights2/sum_weights2
        del weights2, sum_weights2
        revision1 = sum_func(data * normed_weights1, (-1, -2))
        diff = field_grids(self.data - revision1, self.grids)
        weighted_diff = sum_func(diff * normed_weights2, (-1, -2))
        revision2 = revision1 + weighted_diff
        return revision2

    def lowpass(self, g=0.3, c=150000):
        """
        Selecting different parameters *g* and *c*
        will result in different filtering characteristics.

        Reference:
        DOI : https://doi.org/10.1175/1520-0493(1980)108<1108:AOTFSM>2.0.CO;2

        Parameters
        ----------
        g : float, generally between (0, 1]
            Constant parameter.
        c : int
            Constant parameter. When *c* takes a larger value, the filter function converges
            at a larger wavelength, and the response function slowly approaches the maximum value,
            which means that high-frequency fluctuations have been filtered out.

        Returns
        -------
        data_vars : array
            Data field after filtering out high-frequency fluctuations

        """
        return self.__convert_data(self.__lowpass(g=g, c=c))

    def bandpass(self, g1=0.3, c1=30000, g2=0.3, c2=150000):
        """
        Select two different filtering schemes 1 and 2, and perform the filtering separately.
        And then perform the difference, that means *scheme1 - scheme2*.
        The mesoscale fluctuations are thus preserved.

        Parameters
        ----------
        g1 : float, generally between (0, 1]
            Constant parameter of scheme1.
        c1 : int
            Constant parameterof scheme1.
        g2 : float, generally between (0, 1]
            Constant parameter of scheme2.
        c2 : int
            Constant parameterof scheme2.

        Returns
        -------
        data_vars : array
            Mesoscale wave field filtered out from raw data
        """
        lowpass1 = self.__lowpass(g1, c1)
        lowpass2 = self.__lowpass(g2, c2)
        return self.__convert_data(lowpass1 - lowpass2)



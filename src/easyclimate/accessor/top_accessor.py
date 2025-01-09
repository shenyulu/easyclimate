import xarray as xr
from .plot_accessor import PlotAccessor


@xr.register_dataset_accessor("easyclimate")
@xr.register_dataset_accessor("ecl")
class EasyClimateAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._plot = PlotAccessor(xarray_obj)

    @property
    def plot(self):
        return self._plot

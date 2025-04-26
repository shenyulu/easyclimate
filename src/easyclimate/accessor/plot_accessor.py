import xarray as xr
import matplotlib.pyplot as plt
from typing import Hashable
from matplotlib.axes import Axes
from .. import plot as eclplot


class PlotAccessor:
    def __init__(self, ds: xr.Dataset):
        self.ds = ds

    def curved_quiver(self, **kwargs):
        return eclplot.curved_quiver(self.ds, **kwargs)

    def bar_plot_with_threshold(self, **kwargs):
        return eclplot.bar_plot_with_threshold(self.ds, **kwargs)

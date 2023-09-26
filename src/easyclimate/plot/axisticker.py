"""
ticker
"""
import cartopy.mpl.ticker as geoticker
import matplotlib.ticker as ticker

def set_lon_format_axis(ax, axis = 'x' , *args):
    """
    axis = ax.xaxis
    """
    if axis == 'x':
        axis = ax.xaxis
    elif axis == 'y':
        axis = ax.yaxis
    
    axis.set_major_formatter(geoticker.LongitudeFormatter(*args))

def set_lat_format_axis(ax, axis = 'y', *args):
    """
    axis = ax.yaxis
    """
    if axis == 'x':
        axis = ax.xaxis
    elif axis == 'y':
        axis = ax.yaxis
    
    axis.set_major_formatter(geoticker.LatitudeFormatter(*args))

def set_p_format_axis(ax, axis = 'y', axis_limits = (1000, 100), ticker_step = 100):
    """
    a
    """
    if axis == 'x':
        axis = ax.xaxis
        ax.set_xscale('log')
        ax.set_xlim(axis_limits)
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.xaxis.set_major_locator(ticker.MultipleLocator(ticker_step))
    elif axis == 'y':
        axis = ax.yaxis
        ax.set_yscale('log')
        ax.set_ylim(axis_limits)
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.yaxis.set_major_locator(ticker.MultipleLocator(ticker_step))

    
    
    
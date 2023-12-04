import numpy as np

def round_sf_np(
    x: np.array,
    significant_figure: int = 4
) -> np.array:
    """
    Take four significant figures
    """
    r=np.ceil(np.log(x)/np.log(10))
    f=significant_figure
    return np.round(x*(10**(f-r)),0)*(10**(r-f))
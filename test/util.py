import numpy as np
from pathlib import Path


def round_sf_np(x: np.array, significant_figure: int = 4) -> np.array:
    """
    Take four significant figures
    """
    r = np.ceil(np.log(x) / np.log(10))
    f = significant_figure
    return np.round(x * (10 ** (f - r)), 0) * (10 ** (r - f))


def round_sf_np_new(
    arr: np.array,
) -> np.array:
    """
    Take Two significant figures (NEW)
    """
    return np.array(["{:.2g}".format(num) for num in arr], dtype=float)


def assert_path_dir_exist(folder_path):
    if not Path(folder_path).exists():
        Path(folder_path).mkdir()
        print(r"The path {folder_path} is create!")
    else:
        pass

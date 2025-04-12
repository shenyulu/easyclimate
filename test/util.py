import numpy as np
from pathlib import Path
from typing import Union, List


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


def is_mostly_true(
    bool_array: Union[List[bool], np.ndarray], threshold: float = 0.8
) -> bool:
    """
    Checks if the proportion of True values in a boolean array exceeds the given threshold.
    Supports both Python lists and NumPy arrays.

    Args:
        bool_array: A list or NumPy array of boolean values
        threshold: The minimum proportion required (default: 0.8 or 80%)

    Returns:
        bool: True if proportion exceeds threshold, False otherwise

    Examples:
        >>> is_mostly_true([True, True, False])
        False
        >>> is_mostly_true(np.array([True, True, True, False]), threshold=0.7)
        True
        >>> is_mostly_true(np.array([True, False]))
        False
    """
    if len(bool_array) == 0:
        return False

    # Convert to NumPy array if it's a list
    arr = np.array(bool_array) if isinstance(bool_array, list) else bool_array

    true_count = np.sum(arr)
    proportion = true_count / len(arr)
    return proportion > threshold  # Strictly greater than

"""
pytest for core/utility.py

Part3
"""

import pytest

from packaging import version
from functools import wraps
import warnings
from unittest.mock import patch
from easyclimate.core.utility import check_deprecation_status, deprecated

# Test version for the module
__version__ = "1.2.3"


# Test cases for check_deprecation_status
@pytest.mark.parametrize(
    "current,removal,expected",
    [
        ("1.0.0", "2.0.0", False),  # current < removale
        ("2.0.0", "2.0.0", True),  # current == removale
        ("2.0.1", "2.0.0", True),  # current > removale
        ("1.9.9", "2.0.0", False),  # just below
        ("2.0.0", "1.9.9", True),  # just above
        ("1.0.0", "1.0.0", True),  # equal versions
        ("1.0.0-alpha", "1.0.0", False),  # pre-release
        ("1.0.0+local", "1.0.0", True),  # local version
    ],
)
def test_check_deprecation_status(current, removal, expected):
    assert check_deprecation_status(current, removal) == expected


# Test the deprecated decorator
def test_deprecated_decorator_warning():
    """Test that the decorator issues a warning when version < removal_version"""
    removal_version = "2.0.0"
    test_version = "1.0.0"

    @deprecated(
        version="1.0.0", removal_version=removal_version, replacement="new_func"
    )
    def test_func():
        return "original"

    with patch("easyclimate.core.utility.__version__", test_version):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")  # ensure all warnings are caught
            result = test_func()

            # Check the warning was issued
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "is deprecated since 1.0.0" in str(w[0].message)
            assert "will be removed in 2.0.0" in str(w[0].message)
            assert "Use `new_func` instead" in str(w[0].message)

            # Check the function still works
            assert result == "original"


def test_deprecated_decorator_error():
    """Test that the decorator raises an error when version >= removal_version"""
    removal_version = "1.0.0"
    test_version = "1.0.0"

    @deprecated(version="0.9.0", removal_version=removal_version)
    def test_func():
        return "original"

    with patch("easyclimate.core.utility.__version__", test_version):
        with pytest.raises(RuntimeError) as excinfo:
            test_func()

        assert "was removed in version 1.0.0" in str(excinfo.value)


def test_deprecated_decorator_no_replacement():
    """Test the decorator without replacement function"""
    removal_version = "2.0.0"
    test_version = "1.0.0"

    @deprecated(version="1.0.0", removal_version=removal_version)
    def test_func():
        return "original"

    with patch("easyclimate.core.utility.__version__", test_version):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            test_func()

            assert len(w) == 1
            assert "is deprecated since 1.0.0" in str(w[0].message)
            assert "will be removed in 2.0.0" in str(w[0].message)
            assert "Use `" not in str(w[0].message)  # No replacement mentioned


def test_deprecated_decorator_preserves_metadata():
    """Test that the decorator preserves the original function's metadata"""

    @deprecated(version="1.0.0", removal_version="2.0.0")
    def func_with_docstring():
        """Original docstring"""
        pass

    assert func_with_docstring.__name__ == "func_with_docstring"
    assert func_with_docstring.__doc__ == "Original docstring"

"""
pytest for core/spec/plot.py
"""

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr

from easyclimate.core.spec.plot import (
    plot_spectral_matrix,
    plot_spectral_power,
    plot_spectral_real_imag,
    spectral_coefficients_to_matrix,
    spectral_power_by_total_wavenumber,
)


@pytest.fixture
def packed_spectral_coefficients():
    """Create deterministic T4 packed triangular spectral coefficients."""
    ntrunc = 4
    nmode = (ntrunc + 1) * (ntrunc + 2) // 2
    values = np.arange(1, nmode + 1, dtype=np.float64) + 1j * np.arange(nmode)
    return xr.DataArray(
        values,
        dims=("spec_dim",),
        coords={"spec_dim": np.arange(nmode)},
        name="sample_spec",
    )


def test_spectral_coefficients_to_matrix_values(packed_spectral_coefficients):
    matrix = spectral_coefficients_to_matrix(packed_spectral_coefficients, value="real")

    assert matrix.dims == ("m", "n")
    assert matrix.shape == (5, 5)
    assert matrix.attrs["spectral_grid"] == "triangular"
    assert matrix.attrs["ntrunc"] == 4
    assert np.isnan(matrix.sel(m=1, n=0))
    assert matrix.sel(m=0, n=0) == pytest.approx(1.0)
    assert matrix.sel(m=0, n=4) == pytest.approx(5.0)
    assert matrix.sel(m=1, n=1) == pytest.approx(6.0)


def test_spectral_coefficients_to_matrix_quantities(packed_spectral_coefficients):
    coeff = packed_spectral_coefficients.values

    abs_matrix = spectral_coefficients_to_matrix(
        packed_spectral_coefficients, value="abs"
    )
    log_matrix = spectral_coefficients_to_matrix(
        packed_spectral_coefficients, value="logabs"
    )
    power_matrix = spectral_coefficients_to_matrix(
        packed_spectral_coefficients, value="power"
    )
    imag_matrix = spectral_coefficients_to_matrix(
        packed_spectral_coefficients, value="imag"
    )

    assert abs_matrix.sel(m=0, n=0) == pytest.approx(np.abs(coeff[0]))
    assert log_matrix.sel(m=0, n=0) == pytest.approx(np.log10(np.abs(coeff[0]) + 1e-30))
    assert power_matrix.sel(m=0, n=0) == pytest.approx(np.abs(coeff[0]) ** 2)
    assert imag_matrix.sel(m=0, n=0) == pytest.approx(coeff[0].imag)


def test_spectral_power_by_total_wavenumber(packed_spectral_coefficients):
    power = spectral_power_by_total_wavenumber(packed_spectral_coefficients)
    total_wavenumber = np.array([0, 1, 2, 3, 4, 1, 2, 3, 4, 2, 3, 4, 3, 4, 4])
    expected = np.bincount(
        total_wavenumber,
        weights=np.abs(packed_spectral_coefficients.values) ** 2,
        minlength=5,
    )

    assert power.dims == ("n",)
    np.testing.assert_allclose(power.values, expected)


def test_plot_spectral_power_draws_line(packed_spectral_coefficients):
    fig, ax = plot_spectral_power(packed_spectral_coefficients, logy=False, color="k")

    assert fig is ax.figure
    assert len(ax.lines) == 1
    assert ax.lines[0].get_color() == "k"
    assert ax.get_xlabel() == "Total wavenumber n"
    plt.close(fig)


def test_plot_spectral_matrix_draws_image(packed_spectral_coefficients):
    fig, ax = plt.subplots()
    returned_ax = plot_spectral_matrix(
        packed_spectral_coefficients, value="real", ax=ax, add_colorbar=False
    )

    assert returned_ax is ax
    assert len(ax.images) == 1
    assert ax.images[0].get_array().shape == (5, 5)
    plt.close(fig)


def test_plot_spectral_real_imag_draws_two_images(packed_spectral_coefficients):
    fig, axes = plot_spectral_real_imag(packed_spectral_coefficients)

    assert axes.shape == (2,)
    assert len(axes[0].images) == 1
    assert len(axes[1].images) == 1
    plt.close(fig)


def test_spectral_plot_input_validation(packed_spectral_coefficients):
    with pytest.raises(ValueError, match="non-mode dimension"):
        spectral_power_by_total_wavenumber(
            packed_spectral_coefficients.expand_dims(time=[0])
        )

    with pytest.raises(ValueError, match="must be 'abs'"):
        spectral_coefficients_to_matrix(packed_spectral_coefficients, value="bad")

    with pytest.raises(ValueError, match="triangular spectral truncation"):
        spectral_coefficients_to_matrix(
            xr.DataArray(np.arange(5), dims=("spec_dim",)), value="abs"
        )

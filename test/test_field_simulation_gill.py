"""Tests for the Gill-model Python API and plotting helpers."""

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
from matplotlib.quiver import Quiver
import numpy as np
import pytest
import xarray as xr
from xarray.plot.accessor import DataArrayPlotAccessor

from easyclimate.field.simulation.gill import api, plotting


@pytest.fixture
def raw_solver_result():
    """Return a small backend-shaped result with stored snapshots."""
    final = np.arange(9, dtype=float).reshape(3, 3)
    staggered = np.arange(6, dtype=float).reshape(3, 2)
    history = np.stack((final, final + 10.0))
    return {
        "case": "gill_sym",
        "case_label": "Gill: symmetric heating",
        "short_name": "Gill_sym",
        "t": 2.0,
        "nt": 4,
        "nx": 3,
        "ny": 2,
        "nz": 3,
        "lx": 6.0,
        "ly": 4.0,
        "lz": np.pi,
        "damping": 0.1,
        "kernel": "serial",
        "dt": 0.5,
        "dx": 2.0,
        "dy": 4.0 / 3.0,
        "dz": np.pi / 2.0,
        "x": [0.0, 2.0, 4.0],
        "y": [-1.0, 0.0, 1.0],
        "z": [0.0, np.pi / 2.0, np.pi],
        "q": final,
        "eta": final + 20.0,
        "u": final + 30.0,
        "v": final + 40.0,
        "v_staggered": staggered,
        "time": [1.0, 2.0],
        "eta_history": history,
        "u_history": history + 20.0,
        "v_history": history + 40.0,
    }


@pytest.fixture
def gill_result():
    """Return a labelled Gill result suitable for plotting tests."""
    x = np.array([-1.0, 0.0, 1.0, 2.0])
    y = np.array([-1.0, 0.0, 1.0])
    z = np.array([0.0, np.pi / 2.0, np.pi])
    u = np.arange(12, dtype=float).reshape(3, 4)
    v = u + 2.0
    return xr.Dataset(
        data_vars={
            "eta": (("y", "x"), u - 6.0),
            "u": (("y", "x"), u),
            "v": (("y", "x"), v),
        },
        coords={"x": x, "y": y, "z": z},
        attrs={"lz": np.pi, "case_label": "test case", "t": 1.25},
    )


def test_simulate_converts_backend_result_and_forwards_options(
    monkeypatch, raw_solver_result
):
    calls = []

    def fake_gill_simulate(**kwargs):
        calls.append(kwargs)
        return raw_solver_result

    monkeypatch.setattr(api, "gill_simulate", fake_gill_simulate)

    result = api.simulate(
        case="gill_asym",
        t=2.0,
        nt=4,
        nx=3,
        ny=2,
        nz=3,
        lx=6.0,
        ly=4.0,
        lz=np.pi,
        damping=0.1,
        store_every=2,
        kernel="serial",
    )

    assert calls == [
        {
            "case": "gill_asym",
            "t": 2.0,
            "nt": 4,
            "nx": 3,
            "ny": 2,
            "nz": 3,
            "lx": 6.0,
            "ly": 4.0,
            "lz": np.pi,
            "damping": 0.1,
            "store_every": 2,
            "kernel": "serial",
        }
    ]
    assert result.sizes == {"x": 3, "y": 3, "y_staggered": 2, "z": 3, "time": 2}
    assert result["eta"].dims == ("y", "x")
    assert result["v_staggered"].dims == ("y_staggered", "x")
    assert result["eta_history"].dims == ("time", "y", "x")
    np.testing.assert_array_equal(result["eta"].values, raw_solver_result["eta"].T)
    np.testing.assert_array_equal(
        result["u_history"].values, raw_solver_result["u_history"].transpose(0, 2, 1)
    )
    assert result.attrs["kernel"] == "serial"
    assert result.attrs["dt"] == pytest.approx(0.5)


def test_plot_snapshot_draws_on_supplied_axes(gill_result):
    fig, ax = plt.subplots()

    returned_ax = plotting.plot_snapshot(
        gill_result,
        ax=ax,
        nq=2,
        title="Custom title",
        contourf_kwargs={"levels": 3},
        quiver_kwargs={"color": "black"},
    )

    assert returned_ax is ax
    assert ax.get_title() == "Custom title"
    assert ax.get_xlabel() == "x"
    assert ax.get_ylabel() == "y"
    assert ax.get_aspect() == 1.0
    assert any(isinstance(child, Quiver) for child in ax.get_children())
    plt.close(fig)


def test_plot_cross_sections_reconstructs_expected_diagnostics(
    monkeypatch, gill_result
):
    captured = []

    def capture_contourf(self, *args, **kwargs):
        captured.append((self._da.copy(), kwargs))
        return None

    monkeypatch.setattr(DataArrayPlotAccessor, "contourf", capture_contourf)
    fig = plt.figure()

    returned_fig = plotting.plot_cross_sections(
        gill_result,
        fig=fig,
        contourf_uw_kwargs={"levels": 5},
        contourf_vw_kwargs={"cmap": "Blues"},
        contourf_u_kwargs={"vmax": 10.0},
    )

    u_ymean = gill_result["u"].mean(dim="y")
    v_xmean = gill_result["v"].mean(dim="x")
    u_xmean = gill_result["u"].mean(dim="x")
    z = gill_result["z"]
    expected_uw = (u_ymean * np.sin(z)).transpose("z", "x")
    expected_vw = (-v_xmean * np.sin(z)).transpose("z", "y")
    expected_u = (u_xmean * np.cos(z)).transpose("z", "y")

    assert returned_fig is fig
    assert len(captured) == 3
    xr.testing.assert_allclose(captured[0][0], expected_uw)
    xr.testing.assert_allclose(captured[1][0], expected_vw)
    xr.testing.assert_allclose(captured[2][0], expected_u)
    assert captured[0][1]["levels"] == 5
    assert captured[1][1]["cmap"] == "Blues"
    assert captured[2][1]["vmax"] == 10.0
    assert [ax.get_title().splitlines()[-1] for ax in fig.axes] == ["u-w", "v-w", "u"]
    assert [ax.get_yticklabels()[0].get_text() for ax in fig.axes] == ["0", "0", "0"]
    plt.close(fig)

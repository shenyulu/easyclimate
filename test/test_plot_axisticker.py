"""
pytest for plot.axisticker.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import matplotlib.pyplot as plt


def test_set_lon_format_axis():
    fig, ax = plt.subplots()

    ax.xaxis.set_ticks(np.linspace(0, 360, 10))
    ax.yaxis.set_ticks(np.linspace(-90, 90, 10))

    ecl.plot.set_lon_format_axis(ax, axis="x")
    result_data1 = str(ax.get_xticklabels())
    refer_data1 = "[Text(0.0, 0, '0°'), Text(40.0, 0, '40°E'), Text(80.0, 0, '80°E'), Text(120.0, 0, '120°E'), Text(160.0, 0, '160°E'), Text(200.0, 0, '160°W'), Text(240.0, 0, '120°W'), Text(280.0, 0, '80°W'), Text(320.0, 0, '40°W'), Text(360.0, 0, '0°')]"
    assert result_data1 == refer_data1

    ecl.plot.set_lon_format_axis(ax, axis="y")
    result_data2 = str(ax.get_yticklabels())
    refer_data2 = "[Text(0, -90.0, '90°W'), Text(0, -70.0, '70°W'), Text(0, -50.0, '50°W'), Text(0, -30.0, '30°W'), Text(0, -10.0, '10°W'), Text(0, 10.0, '10°E'), Text(0, 30.0, '30°E'), Text(0, 50.0, '50°E'), Text(0, 70.0, '70°E'), Text(0, 90.0, '90°E')]"
    assert result_data2 == refer_data2

    with pytest.raises(ValueError):
        ecl.plot.set_lon_format_axis(ax, axis="This is a sample test!")
        assert 1 == 1


def test_set_lat_format_axis():
    fig, ax = plt.subplots()

    ax.xaxis.set_ticks(np.linspace(0, 360, 10))
    ax.yaxis.set_ticks(np.linspace(-90, 90, 10))

    ecl.plot.set_lat_format_axis(ax, axis="x")
    result_data1 = str(ax.get_xticklabels())
    refer_data1 = "[Text(0.0, 0, '0°'), Text(40.0, 0, '40°N'), Text(80.0, 0, '80°N'), Text(120.0, 0, '120°N'), Text(160.0, 0, '160°N'), Text(200.0, 0, '200°N'), Text(240.0, 0, '240°N'), Text(280.0, 0, '280°N'), Text(320.0, 0, '320°N'), Text(360.0, 0, '360°N')]"
    assert result_data1 == refer_data1

    ecl.plot.set_lat_format_axis(ax, axis="y")
    result_data2 = str(ax.get_yticklabels())
    refer_data2 = "[Text(0, -90.0, '90°S'), Text(0, -70.0, '70°S'), Text(0, -50.0, '50°S'), Text(0, -30.0, '30°S'), Text(0, -10.0, '10°S'), Text(0, 10.0, '10°N'), Text(0, 30.0, '30°N'), Text(0, 50.0, '50°N'), Text(0, 70.0, '70°N'), Text(0, 90.0, '90°N')]"
    assert result_data2 == refer_data2

    with pytest.raises(ValueError):
        ecl.plot.set_lat_format_axis(ax, axis="This is a sample test!")
        assert 1 == 1


def test_set_p_format_axis():
    fig, ax = plt.subplots()
    ecl.plot.set_p_format_axis(ax, axis="y")
    result_data1 = str(ax.get_yticklabels())
    refer_data1 = "[Text(0, 0.0, '0'), Text(0, 100.0, '100'), Text(0, 200.0, '200'), Text(0, 300.0, '300'), Text(0, 400.0, '400'), Text(0, 500.0, '500'), Text(0, 600.0, '600'), Text(0, 700.0, '700'), Text(0, 800.0, '800'), Text(0, 900.0, '900'), Text(0, 1000.0, '1000'), Text(0, 1100.0, '1100')]"
    result_data2 = ax.get_yscale()
    refer_data2 = "log"
    assert result_data1 == refer_data1
    assert result_data2 == refer_data2

    fig, ax = plt.subplots()
    ecl.plot.set_p_format_axis(ax, axis="x")
    result_data3 = str(ax.get_xticklabels())
    refer_data3 = "[Text(0.0, 0, '0'), Text(100.0, 0, '100'), Text(200.0, 0, '200'), Text(300.0, 0, '300'), Text(400.0, 0, '400'), Text(500.0, 0, '500'), Text(600.0, 0, '600'), Text(700.0, 0, '700'), Text(800.0, 0, '800'), Text(900.0, 0, '900'), Text(1000.0, 0, '1000'), Text(1100.0, 0, '1100')]"
    result_data4 = ax.get_xscale()
    refer_data4 = "log"
    assert result_data3 == refer_data3
    assert result_data4 == refer_data4

    with pytest.raises(ValueError):
        ecl.plot.set_p_format_axis(ax, axis="This is a sample test!")
        assert 1 == 1

"""
pytest for plot.geolabel.py
"""

import matplotlib.pyplot as plt

import easyclimate as ecl


def test_add_geolatitude_label_default_kwargs():
    fig, ax = plt.subplots()

    text = ecl.plot.add_geolatitude_label(ax)

    assert text in ax.texts
    assert text.get_position() == (-0.07, 0.55)
    assert text.get_text() == ""
    assert text.get_va() == "bottom"
    assert text.get_ha() == "center"
    assert text.get_rotation() == 90.0
    assert text.get_rotation_mode() == "anchor"

    plt.close(fig)


def test_add_geolongitude_label_default_kwargs():
    fig, ax = plt.subplots()

    text = ecl.plot.add_geolongitude_label(ax)

    assert text in ax.texts
    assert text.get_position() == (0.5, -0.2)
    assert text.get_text() == ""
    assert text.get_va() == "bottom"
    assert text.get_ha() == "center"
    assert text.get_rotation() == 90.0
    assert text.get_rotation_mode() == "anchor"

    plt.close(fig)


def test_geolabel_kwargs_override_defaults():
    fig, ax = plt.subplots()

    text = ecl.plot.add_geolatitude_label(
        ax,
        x=-0.1,
        y=0.6,
        s="Latitude",
        va="top",
        ha="right",
        rotation=0,
        rotation_mode="default",
        color="red",
    )

    assert text.get_position() == (-0.1, 0.6)
    assert text.get_text() == "Latitude"
    assert text.get_va() == "top"
    assert text.get_ha() == "right"
    assert text.get_rotation() == 0.0
    assert text.get_rotation_mode() == "default"
    assert text.get_color() == "red"

    plt.close(fig)

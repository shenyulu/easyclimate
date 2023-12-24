.. _api:

List of functions and classes (API)
===================================

.. note::

    **All functions and classes should be accessed from the** :mod:`easyclimate`
    **top-level namespace.**

    Modules inside of the :mod:`easyclimate` package are meant mostly for internal
    organization. Please **avoid importing** directly from submodules since
    functions/classes may be moved around.

.. automodule:: easyclimate

Core
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.core.diagnosis
    easyclimate.core.diff
    easyclimate.core.eddy
    easyclimate.core.eof
    easyclimate.core.extract
    easyclimate.core.mk_test
    easyclimate.core.read
    easyclimate.core.stat
    easyclimate.core.tutorial
    easyclimate.core.utility
    easyclimate.core.variability
    easyclimate.core.yearstat

Filter
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.filter.butter_filter
    easyclimate.filter.smooth
    .. easyclimate.filter.barnes_filter

Interpolation
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.interp.point2mesh
    easyclimate.interp.mesh2mesh

Plot
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.plot.axisticker
    easyclimate.plot.projection
    easyclimate.plot.significance_plot
    easyclimate.plot.taylor_diagram
    easyclimate.plot.wind

Wavelet Transform
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.wavelet


Spherical Harmonics of Wind Fields
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.windspharm.top

Climate Field
----------------------------------------

Air–Sea Interaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
    :toctree: generated/

    easyclimate.field.air_sea_interaction.index_enso

Atmosphere
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
    :toctree: generated/

    easyclimate.field.atm.index_PNA

Land–Atmosphere Interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
    :toctree: generated/

    easyclimate.field.land_atm_interaction

Monsoon
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
    :toctree: generated/

    easyclimate.field.monsoon.index_npwi

Ocean
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
    :toctree: generated/

    easyclimate.field.ocean.mixlayer
    easyclimate.field.ocean.oceanic_front
    easyclimate.field.ocean.stability
    easyclimate.field.ocean.thermal
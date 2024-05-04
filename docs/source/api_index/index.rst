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
    easyclimate.core.windspharm

Filter
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.filter.butter_filter
    easyclimate.filter.smooth
    easyclimate.filter.wavelet

Interpolation
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.interp.point2mesh
    easyclimate.interp.mesh2mesh
    easyclimate.interp.interp1d_vertical_model2pressure

Plot
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.plot.axisticker
    easyclimate.plot.projection
    easyclimate.plot.significance_plot
    easyclimate.plot.taylor_diagram


Climate Field
----------------------------------------

Air–Sea Interaction
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.air_sea_interaction
    easyclimate.field.air_sea_interaction.index_enso
    easyclimate.field.air_sea_interaction.index_iod
    easyclimate.field.air_sea_interaction.index_iobm
    easyclimate.field.air_sea_interaction.index_amm
    easyclimate.field.air_sea_interaction.index_atlantic_nino

Teleconnections
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.teleconnection
    easyclimate.field.teleconnection.index_pna
    easyclimate.field.teleconnection.index_nao
    easyclimate.field.teleconnection.index_ea
    easyclimate.field.teleconnection.index_wa
    easyclimate.field.teleconnection.index_wp
    easyclimate.field.teleconnection.index_eu
    easyclimate.field.teleconnection.index_srp
    easyclimate.field.teleconnection.index_cgt

Land–Atmosphere Interactions
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.land

Monsoon
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.monsoon
    easyclimate.field.monsoon.index_npwi

Ocean
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.ocean
    easyclimate.field.ocean.mixlayer
    easyclimate.field.ocean.oceanic_front
    easyclimate.field.ocean.stability
    easyclimate.field.ocean.thermal

.. _api:

List of Functions and Classes (API)
===================================

.. note::

    **All functions and classes should be accessed from the** :mod:`easyclimate`
    **top-level namespace.**

    Modules inside of the :mod:`easyclimate` package are meant mostly for internal
    organization. Please **avoid importing** directly from submodules since
    functions/classes may be moved around.

.. automodule:: easyclimate

Core❤️
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.core.datanode
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
    easyclimate.core.spharm
    easyclimate.core.windspharm

Physics🗺️
----------------------------------------

Geophysical Parameters Related
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.physics.geo.coriolis
    easyclimate.physics.geo.linrood_latwgt

Dry Related
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.physics.dry

Moisture Related
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.physics.moisture.dewpoint
    easyclimate.physics.moisture.lapse_rate
    easyclimate.physics.moisture.mix
    easyclimate.physics.moisture.vapor_pressure
    easyclimate.physics.moisture.wet_bulb

Temperature Related
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.physics.temperature.equivalent_potential_temperature
    easyclimate.physics.temperature.potential_temperature
    easyclimate.physics.temperature.virtual_temperature

Condensation Related
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.physics.condensation.lcl

Convection Related
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.physics.convection.stability

Energy Related
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.physics.energy.enthalpy
    easyclimate.physics.energy.latent_heat_water

Others
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.physics.transfer

Filter🎹
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.filter.butter_filter
    easyclimate.filter.lanczos_filter
    easyclimate.filter.kf_filter
    easyclimate.filter.barnes_filter
    easyclimate.filter.gaussian_filter
    easyclimate.filter.smooth
    easyclimate.filter.wavelet
    easyclimate.filter.redfit
    easyclimate.filter.spatial_pcf
    easyclimate.filter.emd

Interpolation🔗
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.interp.barnes
    easyclimate.interp.mesh2mesh
    easyclimate.interp.mesh2point
    easyclimate.interp.interp1d_vertical_model2pressure
    easyclimate.interp.interp1d_vertical_pressure2altitude
    easyclimate.interp.vinth2p_dp

Plot🖊️
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.plot.axisticker
    easyclimate.plot.projection
    easyclimate.plot.significance_plot
    easyclimate.plot.taylor_diagram
    easyclimate.plot.quick_draw
    easyclimate.plot.curved_quiver_plot

WRF-python🌌
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.wrf.interface

Satellite🪐
----------------------------------------

.. autosummary::
    :toctree: generated/

    easyclimate.satellite.image

Meteorology Field🌂
----------------------------------------

Air–Sea Interaction
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.air_sea_interaction.index_enso
    easyclimate.field.air_sea_interaction.index_iod
    easyclimate.field.air_sea_interaction.index_iobm
    easyclimate.field.air_sea_interaction.index_amm
    easyclimate.field.air_sea_interaction.index_atlantic_nino

Teleconnections
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.teleconnection.index_pna
    easyclimate.field.teleconnection.index_nao
    easyclimate.field.teleconnection.index_ea
    easyclimate.field.teleconnection.index_wa
    easyclimate.field.teleconnection.index_wp
    easyclimate.field.teleconnection.index_eu
    easyclimate.field.teleconnection.index_srp
    easyclimate.field.teleconnection.index_cgt
    easyclimate.field.teleconnection.index_ao_nam

Land
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.land

Monsoon
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.monsoon.index_npwi

Ocean
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.ocean.mixlayer
    easyclimate.field.ocean.oceanic_front
    easyclimate.field.ocean.stability
    easyclimate.field.ocean.thermal

Typhoon
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.typhoon.potential_intensity

Equatorial Wave
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.equatorial_wave.mjo
    easyclimate.field.equatorial_wave.wk_spectra

Heat Stress
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.heat_stress.humanindexmod_2020

Boundary Layer
::::::::::::::::::::::::::::::::::::::::

.. autosummary::
    :toctree: generated/

    easyclimate.field.boundary_layer.aerobulk

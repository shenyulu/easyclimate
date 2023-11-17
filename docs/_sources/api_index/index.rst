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
----

.. autosummary::
    :toctree: generated/

    easyclimate.core.diagnosis
    easyclimate.core.diff
    easyclimate.core.eddy
    easyclimate.core.eof
    easyclimate.core.extract
    easyclimate.core.read
    easyclimate.core.stat
    easyclimate.core.utility
    easyclimate.core.variability
    easyclimate.core.yearstat

Filter
----

.. autosummary::
    :toctree: generated/

    .. easyclimate.filter.barnes_filter
    .. easyclimate.filter.filter

Interpolation
----

.. autosummary::
    :toctree: generated/

    easyclimate.interp.point2mesh
    easyclimate.interp.mesh2mesh

Ocean
----

.. autosummary::
    :toctree: generated/

    easyclimate.ocean.mixlayer
    easyclimate.ocean.stability
    easyclimate.ocean.thermal

Plot
----

.. autosummary::
    :toctree: generated/

    easyclimate.plot.axisticker
    easyclimate.plot.projection
    easyclimate.plot.significance_plot
    easyclimate.plot.taylor_diagram
    easyclimate.plot.wind

Spherical harmonics of wind fields
----

.. autosummary::
    :toctree: generated/

    easyclimate.windspharm.top

Index
----

.. autosummary::
    :toctree: generated/

    easyclimate.index.enso
    easyclimate.index.npwi
    easyclimate.index.oceanic_front
    easyclimate.index.pna
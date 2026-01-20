.. easyclimate documentation master file, created by
   sphinx-quickstart on Mon Mar 20 14:33:30 2023.

.. Easy Climate
.. =======================================

.. .. image:: _static/easyclimate_logo_mini.png

.. Easy Climate is a Python library for processing spatial data (bathymetry, geophysics surveys, etc)
.. and interpolating it on regular grids (i.e., gridding).

.. Our core interpolation methods are inspired by machine-learning.
.. As such, Verde implements an interface that is similar to the popular scikit-learn library.
.. We also provide other analysis methods that are often used in combination with gridding,
.. like trend removal, blocked/windowed operations, cross-validation, and more!

.. title:: Home

========
|banner|
========

.. |banner| image:: _static/easyclimate_logo_mini.png
    :alt: Easy-climate Documentation
    :align: middle

.. raw:: html

    <p class="lead centered front-page-callout">
        A line of code to analyze climate
        <br>
        Trying to make it easily without complicated writing code?
        <br>
        <strong>Easy Climate is just here to help you!</strong>
    </p>


You can directly install it via `pip` by using 🛒

.. code-block:: bash

    $ pip install easyclimate


----

.. grid:: 1 2 1 2
    :margin: 5 5 0 0
    :padding: 0 0 0 0
    :gutter: 4

    .. grid-item-card:: :octicon:`info` Try Online🤗
        :text-align: center
        :class-title: sd-fs-5
        :class-card: sd-p-3

        New to Easy Climate? Try!

        .. button-link:: https://mybinder.org/v2/gh/shenyulu/easyclimate/main?labpath=docs%2Fexample
            :click-parent:
            :color: primary
            :outline:
            :expand:

            Binder Online Engine :octicon:`rocket`

    .. grid-item-card:: :octicon:`comment-discussion` Need help?
        :text-align: center
        :class-title: sd-fs-5
        :class-card: sd-p-3

        Ask on our community channels.

        .. button-link:: https://github.com/shenyulu/easyclimate/discussions
            :click-parent:
            :color: primary
            :outline:
            :expand:

            Join the conversation :octicon:`link-external`

    .. grid-item-card:: :octicon:`file-badge` Reference documentation
        :text-align: center
        :class-title: sd-fs-5
        :class-card: sd-p-3

        A list of modules and functions.

        .. button-ref:: api
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

    .. grid-item-card:: :octicon:`bookmark` Using Easy Climate for research?
        :text-align: center
        :class-title: sd-fs-5
        :class-card: sd-p-3

        Citations help support our work!

        .. button-ref:: citenote
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

            Cite our repository

----

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Getting Started

    overview.rst
    install.rst
    ./auto_gallery/index

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Data Processing

   ./auto_gallery/plot_basic_statistical_analysis
   ./auto_gallery/plot_time_scale_average
   ./auto_gallery/plot_geographic_finite_difference
   ./auto_gallery/plot_interp
   ./auto_gallery/plot_wrf_tutorial

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Plotting

   ./auto_gallery/plot_formatting_coordinates
   ./auto_gallery/plot_taylor_diagram

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Meteorology Field

   ./static_docs/air_sea_interaction
   ./static_docs/teleconnections
   ./static_docs/ocean
   ./static_docs/monsoon
   ./static_docs/land
   ./static_docs/typhoons

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Reference documentation

    ./api_index/index.rst
    changes.md
    ./static_docs/cite
    Open Source Licenses <https://easyclimate-backend.readthedocs.io/en/latest/src/softlist.html>

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Community

    ./contributing.rst
    How to Contribute <https://github.com/shenyulu/easyclimate/blob/main/CONTRIBUTING.md>
    Source Code on GitHub <https://github.com/shenyulu/easyclimate>
    sponsor.rst


.. caution::
    🚨 This package is still undergoing rapid development. 🚨

    All of the API (functions/classes/interfaces) is subject to change.
    There may be non-backward compatible changes as we experiment with new design ideas and implement new features.
    This is not a finished product, use with caution.

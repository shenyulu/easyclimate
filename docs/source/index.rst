.. easyclimate documentation master file, created by
   sphinx-quickstart on Mon Mar 20 14:33:30 2023.

.. Easy climate
.. =======================================

.. .. image:: _static/easyclimate-logo.svg

.. Easy climate is a Python library for processing spatial data (bathymetry, geophysics surveys, etc)
.. and interpolating it on regular grids (i.e., gridding).

.. Our core interpolation methods are inspired by machine-learning.
.. As such, Verde implements an interface that is similar to the popular scikit-learn library.
.. We also provide other analysis methods that are often used in combination with gridding,
.. like trend removal, blocked/windowed operations, cross-validation, and more!

.. title:: Home

========
|banner|
========

.. |banner| image:: _static/easyclimate-logo.svg
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

Installation🛒
------------------------

You can do a direct install via `pip` by using:

.. code-block:: bash

    $ pip install easyclimate

Online experience💻
------------------------
Just click on the link below and wait for the online environment to be configured before using the Easyclimate package on Jupyter notebook online.

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/shenyulu/easyclimate/main?labpath=docs%2Fexample

How to cite📣
-----------------------

If you are using **Easy Climate** and would like to cite it in academic publication, we would certainly appreciate it. We recommend the following citations.
We provide a `Zenodo citation and DOI <https://zenodo.org/doi/10.5281/zenodo.10279567>`__ for this purpose:

An example BibTeX entry:

.. code:: BibTeX

    @misc{easyclimate_v2023_12_1,
        author = {Yulu Shen},
        title  = {easyclimate: v2023.12.1},
        month  = dec,
        year   = 2023,
        doi    = {10.5281/zenodo.10279567},
        url    = {https://doi.org/10.5281/zenodo.10279567}
        }


.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Getting Started

    overview.rst
    install.rst
    .. ./auto_gallery_output/index

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Data Processing

   ./auto_gallery_output/plot_basic_statistical_analysis
   ./auto_gallery_output/plot_time_scale_average
   ./auto_gallery_output/plot_geographic_finite_difference
   ./auto_gallery_output/plot_interp

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Plotting

   ./auto_gallery_output/plot_formatting_coordinates
   ./auto_gallery_output/plot_taylor_diagram

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Climate Field

   ./static_docs/air_sea_interaction
   ./static_docs/ocean
   ./static_docs/monsoon
   ./static_docs/land

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Reference documentation

    ./api_index/index.rst
    changes.md

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Community

    How to contribute <https://github.com/shenyulu/easyclimate/blob/main/CONTRIBUTING.md>
    Source code on GitHub <https://github.com/shenyulu/easyclimate>


Indices and tables🧭
------------------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. caution::
    🚨 This package is still undergoing rapid development. 🚨

    All of the API (functions/classes/interfaces) is subject to change.
    There may be non-backward compatible changes as we experiment with new design ideas and implement new features.
    This is not a finished product, use with caution.

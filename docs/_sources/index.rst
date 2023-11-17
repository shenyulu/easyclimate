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
        A line of code to finish the analysis of climatology
        <br>
        Trying to make it easily without complicated writing code?
        <br>
        <strong>Easy Climate is just here to help you!</strong>
    </p>

Installation
------------

You can do a direct install via `pip` by using:

.. code-block:: bash

    $ pip install easyclimate


How to cite
-----------
If you would like to cite `Easy-climate` you can do so using our `Zenodo deposit <https://zenodo.org/>`__.


.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Getting Started
    
    overview.rst
    install.rst

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Data Processing

   ./auto_gallery_output/plot_basic_statistical_analysis
   ./auto_gallery_output/plot_time_scale_average
   ./auto_gallery_output/plot_geographic_finite_difference   

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Plotting

   ./auto_gallery_output/plot_formatting_coordinates
   ./auto_gallery_output/plot_taylor_diagram

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Indices:

   ./static_docs/air_sea_interaction.md
   ./static_docs/monsoon.md

   
.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Community

    How to contribute <https://github.com/shenyulu/easyclimate/blob/main/CONTRIBUTING.md>
    Source code on GitHub <https://github.com/shenyulu/easyclimate>

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Reference documentation

    ./api_index/index.rst
    changes.rst
    

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. caution::
    🚨 This package is still undergoing rapid development. 🚨

    All of the API (functions/classes/interfaces) is subject to change. 
    There may be non-backward compatible changes as we experiment with new design ideas and implement new features. 
    This is not a finished product, use with caution.
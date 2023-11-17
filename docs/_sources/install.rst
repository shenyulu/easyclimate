.. _install:

Installation
====================

There are different ways to install **easyclimate**:

Pypi
------------

Using the `pip <https://pypi.org/project/pip/>`__ package manager:

.. code:: bash

    $ python -m pip install easyclimate

Development version
------------------------

You can use ``pip`` to install the latest **unreleased** version from
GitHub (**not recommended** in most situations):

.. code:: bash

    $ python -m pip install --upgrade git+https://github.com/shenyulu/easyclimate

.. note::

    The commands above should be executed in a terminal. On Windows, use the
    ``cmd.exe`` or the "Anaconda Prompt" app if you're using Anaconda.


Which Python?
-------------

You'll need **Python** :math:`\geq` **3.10**.


.. _dependencies:

Dependencies
------------

The required dependencies should be installed automatically when you install
Verde using ``conda`` or ``pip``.

Required:

* `numpy <http://www.numpy.org/>`__ (1.24.3 or later)
* `scipy <https://docs.scipy.org/doc/scipy/reference/>`__
* `pandas <http://pandas.pydata.org/>`__
* `xarray <http://xarray.pydata.org/>`__
* `dask <https://dask.org/>`__
* `geocat.viz <https://github.com/NCAR/geocat-viz>`__
* `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`__ (0.20 or later) for plotting maps
* `pymannkendall <https://github.com/mmhs013/pymannkendall>`__
* `xarray-regrid <https://github.com/EXCITED-CO2/xarray-regrid>`__
* `matplotlib <https://matplotlib.org/>`__
* `cartopy <https://scitools.org.uk/cartopy/>`__


Our examples use other packages as well which are not used within Verde itself.
If you wish to **run the examples in the documentation**, you will also have to
install:


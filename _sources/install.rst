.. _install:

Installation
====================================

There are different ways to install **easyclimate**:

Pypi
------------------------------------

Using the `pip <https://pypi.org/project/pip/>`__ package manager to obtain latest version of easyclimate:

.. code:: bash

    $ pip install easyclimate

Development version
------------------------------------

You can use ``pip`` to install the latest **unreleased** version from
GitHub (**not recommended** in most situations):

.. code:: bash

    $ pip install --upgrade git+https://github.com/shenyulu/easyclimate@dev

.. note::

    The commands above should be executed in a terminal. On Windows, use the
    ``cmd.exe`` or the "Anaconda Prompt" app if you're using Anaconda.


Which Python?
------------------------------------

You'll need **Python** :math:`\geq` **3.10**.


.. _dependencies:

Dependencies
------------------------------------

The required dependencies should be installed automatically when you install
easyclimate using ``conda`` or ``pip``.

Required:

* `numpy <http://www.numpy.org/>`__ (1.24.3 or later)
* `xarray <http://xarray.pydata.org/>`__ (0.17.0 or later)
* `geocat.viz <https://github.com/NCAR/geocat-viz>`__ (2023.10.0 or before)
* `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`__ (0.20 or later)
* `xeofs <https://github.com/nicrie/xeofs>`__ (2.2.2 or later)

* `matplotlib <https://matplotlib.org/>`__
* `pandas <http://pandas.pydata.org/>`__
* `fast-barnes-py <https://github.com/MeteoSwiss/fast-barnes-py>`__
* `python-oceans <https://github.com/pyoceans/python-oceans>`__
* `intel-fortran-rt <https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html>`__
* `scipy <https://docs.scipy.org/doc/scipy/reference/>`__ (1.8.0 or later)
* `statsmodels <https://github.com/statsmodels/statsmodels>`__
* `geocat-comp <https://github.com/NCAR/geocat-comp>`__
* `pyspharm-syl <https://github.com/shenyulu/pyspharm-syl>`__
* `windspharm-syl <https://github.com/shenyulu/windspharm-syl>`__
* `dask <https://dask.org/>`__

* `pymannkendall <https://github.com/mmhs013/pymannkendall>`__
* `flox <https://github.com/xarray-contrib/flox>`__
* `xarray-datatree <https://github.com/xarray-contrib/datatree>`__
* `xarray-regrid <https://github.com/EXCITED-CO2/xarray-regrid>`__
* `gsw-xarray <https://github.com/DocOtak/gsw-xarray>`__
* `pooch <https://github.com/fatiando/pooch>`__
* `tqdm <https://github.com/tqdm/tqdm>`__
* `zarr <https://github.com/zarr-developers/zarr-python>`__
* `metpy <https://github.com/Unidata/MetPy>`__

Our examples use other packages as well which are not used within easyclimate itself.
If you wish to **run the examples in the documentation**, you will also have to
install:

Acknowledgement
------------------------------------
- point2mesh: https://github.com/MeteoSwiss/fast-barnes-py by `MeteoSwiss <https://github.com/MeteoSwiss>`__.
- windspharm: https://github.com/ajdawson/windspharm by `ajdawson <https://github.com/ajdawson>`__.
- wavelets: https://github.com/ct6502/wavelets by `ct6502 <https://github.com/ct6502>`__ and `A Practical Guide to Wavelet Analysis <http://paos.colorado.edu/research/wavelets/>`__
- pybarnes: https://github.com/LinOuyang/pybarnes by `LinOuyang <https://github.com/LinOuyang>`__.

:py:mod:`easyclimate.core.read`
===============================

.. py:module:: easyclimate.core.read

.. autoapi-nested-parse::

   Functions for read data.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.core.read.open_muliti_dataset



.. py:function:: open_muliti_dataset(files: str, dim: str, **kwargs) -> xarray.Dataset

   Compare python library versions.

   .. attention::
       - Only for incoming version numbers without alphabetic characters.
       - Based on this method, the version number comparison should result in the following `"10.12.2.6.5">"10.12.2.6"`.

   Parameters
   ----------
   - ver1: Version number 1
   - ver2: Version number 2

   Returns
   -------
   :py:class:`int<python.int>`.

   .. note::
       If `ver1<ver2`, return `-1`; If `ver1=ver2`, return `0`; If `ver1>ver2`, return `1`.

   Examples
   --------

   .. code:: python
       >>> import easyclimate as ecl
       >>> result = assert_compared_version("10.12.2.6.5", "10.12.2.6")
       >>> print(result)
       1



easyclimate.field.equatorial_wave.mjo
=====================================

.. py:module:: easyclimate.field.equatorial_wave.mjo

.. autoapi-nested-parse::

   Madden Julian Oscillation (MJO)



Functions
---------

.. autoapisummary::

   easyclimate.field.equatorial_wave.mjo.draw_mjo_phase_space_basemap
   easyclimate.field.equatorial_wave.mjo.draw_mjo_phase_space


Module Contents
---------------

.. py:function:: draw_mjo_phase_space_basemap(ax=None, add_phase_text: bool = True, add_location_text: bool = True)

   Create a base map for Madden-Julian Oscillation (MJO) phase space visualization.

   This function generates the standardized background diagram for MJO analysis,
   including phase division lines (1-8), amplitude circles, and optional geographical
   annotations. The diagram uses RMM1/RMM2 indices as axes in a Cartesian coordinate
   system (-4 to 4 range).

   Parameters
   ----------
   ax : :py:class:`matplotlib.axes.Axes`, optional
       Target axes object for plotting. If None, uses current axes (gca()).
   add_phase_text : :py:class:`bool`, default: True
       Whether to annotate MJO phase numbers (1-8) on the diagram.
   add_location_text : :py:class:`bool`, default: True
       Whether to annotate geographical regions corresponding to MJO phases:
       - Pacific Ocean (top)
       - Indian Ocean (bottom)
       - Western Hemisphere/Africa (left)
       - Maritime Continent (right)

   Returns
   -------
   :py:class:`matplotlib.axes.Axes`
       The configured axes object with MJO phase space basemap elements.

   Notes
   -----
   The diagram includes:
   1. 8 radial lines dividing phases (45° intervals)
   2. Unit circle indicating amplitude threshold
   3. Axes labeled as RMM1 (x-axis) and RMM2 (y-axis)

   Examples
   --------
   >>> import matplotlib.pyplot as plt
   >>> fig, ax = plt.subplots()
   >>> ecl.field.equatorial_wave.draw_mjo_phase_space_basemap(ax=ax)
   >>> plt.show()

   .. seealso::

       https://www.ncl.ucar.edu/Applications/mjoclivar.shtml

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_mjo_phase.py


.. py:function:: draw_mjo_phase_space(mjo_data: xarray.Dataset, rmm1_dim: str = 'RMM1', rmm2_dim: str = 'RMM2', time_dim: str = 'time', ax=None, color='blue', start_text='START', add_start_text: bool = True)

   Visualize Madden-Julian Oscillation (MJO) phase space trajectory from RMM indices.

   Plots the temporal evolution of MJO phases as a parametric curve in RMM1-RMM2 space,
   with optional starting point marker. Combines scatter points (individual observations)
   with connecting lines to show temporal progression.

   Parameters
   ----------
   mjo_data : :py:class:`xarray.Dataset`
       Input dataset containing RMM indices. Must include:
       - RMM1 component (real-time multivariate MJO index 1)
       - RMM2 component (real-time multivariate MJO index 2)
       - Time coordinate
   rmm1_dim : :py:class:`str`, default: "RMM1"
       Variable name for RMM1 component in the dataset
   rmm2_dim : :py:class:`str`, default: "RMM2"
       Variable name for RMM2 component in the dataset
   time_dim : :py:class:`str`, default: "time"
       Time coordinate dimension name
   ax : :py:class:`matplotlib.axes.Axes`, optional
       Target axes for plotting. Creates new axes if None.
   color : :py:class:`str` or tuple, default: "blue"
       Color specification for trajectory and points
   start_text : :py:class:`str`, default: "START"
       Annotation text for trajectory starting point
   add_start_text : :py:class:`bool`, default: True
       Toggle for displaying start point annotation

   Returns
   -------
   :py:class:`matplotlib.axes.Axes`
       Configured axes object with MJO phase space trajectory

   Returns
   -------
   :py:class:`matplotlib.axes.Axes`
       Configured axes object with MJO phase space trajectory

   Examples
   --------
   >>> import xarray as xr
   >>> import matplotlib.pyplot as plt
   >>> # Load RMM indices dataset
   >>> mjo_ds = xr.open_dataset('http://iridl.ldeo.columbia.edu/SOURCES/.BoM/.MJO/.RMM/dods', decode_times=False)
   >>> T = mjo_ds.T.values
   >>> mjo_ds['T'] = pd.date_range("1974-06-01", periods=len(T))
   >>> mjo_ds = ecl.utility.get_compress_xarraydata(mjo_ds)
   >>> fig, ax = plt.subplots(figsize=(8,8))
   >>> ecl.field.equatorial_wave.draw_mjo_phase_space_basemap()
   >>> ecl.field.equatorial_wave.draw_mjo_phase_space(
   ...:    mjo_data = mjo_data.sel(time = slice('2024-12-01', '2024-12-31')),
   ...:    rmm1_dim = "RMM1",
   ...:    rmm2_dim = "RMM2",
   ...:    time_dim = "time"rmm_data, ax=ax, color="red")
   >>> plt.title("MJO Phase Space Trajectory")
   >>> plt.show()

   .. note::
       1. Recommended to use with :py:func:`draw_mjo_phase_space_basemap` for basic map.
       2. Temporal resolution affects trajectory smoothness (daily data recommended).

   .. seealso::

       https://www.ncl.ucar.edu/Applications/mjoclivar.shtml

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_mjo_phase.py



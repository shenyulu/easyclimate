easyclimate.plot.line
=====================

.. py:module:: easyclimate.plot.line

.. autoapi-nested-parse::

   Bar plot for xarray dataset



Functions
---------

.. autoapisummary::

   easyclimate.plot.line.line_plot_with_threshold


Module Contents
---------------

.. py:function:: line_plot_with_threshold(da: xarray.DataArray, threshold: float = 0, pos_color: str = 'red', neg_color: str = 'blue', ax=None, line_plot: bool = True, fill_pos_plot: bool = True, fill_neg_plot: bool = True, line_kwargs=None, fill_kwargs=None) -> tuple

   Plot a line chart with proper shading at threshold crossings.

   Parameters:
   -----------
   da : :py:class:`xarray.DataArray <xarray.DataArray>`
       1-dimensional data array
   threshold : :py:class:`float <float>`, optional
       Color separation threshold (default: 0)
   pos_color : :py:class:`str <str>`, optional
       Color for values ≥ threshold (default: 'red')
   neg_color : :py:class:`str <str>`, optional
       Color for values < threshold (default: 'blue')
   ax : matplotlib axes, optional
       Axes to plot on (default: current axes)
   line_kwargs : :py:class:`dict <dict>`, optional
       Arguments for plt.plot
   fill_kwargs : :py:class:`dict <dict>`, optional
       Arguments for plt.fill_between

   Returns:
   --------
   tuple
       (line plot, fill objects)

   .. seealso::

       :py:func:`matplotlib.lines.Line2D <matplotlib.lines.Line2D>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_basic_statistical_analysis.py
       ./dynamic_docs/plot_corr_reg.py



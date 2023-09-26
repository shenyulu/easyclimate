:py:mod:`easyclimate.plot.significance_plot`
============================================

.. py:module:: easyclimate.plot.significance_plot

.. autoapi-nested-parse::

   Functions for draw significant area.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.plot.significance_plot.draw_significant_area
   easyclimate.plot.significance_plot.get_significance_point



.. py:function:: draw_significant_area(p_value_data, threshold=0.05, lon_dim='lon', lat_dim='lat', ax=None, hatch_colors='k', point_density='...', reverse_level_plot=False)

   绘制显著性区域（contourf hatch 方法）
   data: 包含 p 值的 DataArray
   ax: 绘制的 axes
   hatch_colors: 更改 hatch 图案颜色，类似于 hatch_colors = ['maroon', 'red', 'darkorange']
   point_density: 点密度或者点的绘制类型，可选值有 '.', '..', '...'


.. py:function:: get_significance_point(p_value_data, threshold=0.05, lon_dim='lon', lat_dim='lat')

       
       



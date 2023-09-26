:py:mod:`easyclimate.core.mk_test`
==================================

.. py:module:: easyclimate.core.mk_test

.. autoapi-nested-parse::

   Mann-Kendall Test 



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.core.mk_test.original_test
   easyclimate.core.mk_test.hamed_rao_modification_test
   easyclimate.core.mk_test.yue_wang_modification_test
   easyclimate.core.mk_test.pre_whitening_modification_test
   easyclimate.core.mk_test.trend_free_pre_whitening_modification_test
   easyclimate.core.mk_test.seasonal_test
   easyclimate.core.mk_test.regional_test
   easyclimate.core.mk_test.correlated_seasonal_test
   easyclimate.core.mk_test.partial_test
   easyclimate.core.mk_test.sens_slope
   easyclimate.core.mk_test.seasonal_sens_slope



.. py:function:: original_test(data_input, dim, alpha=0.05, returns_type='dataset_returns')


.. py:function:: hamed_rao_modification_test(data_input, dim, alpha=0.05, lag=None, returns_type='dataset_returns')


.. py:function:: yue_wang_modification_test(data_input, dim, alpha=0.05, lag=None, returns_type='dataset_returns')


.. py:function:: pre_whitening_modification_test(data_input, dim, alpha=0.05, returns_type='dataset_returns')


.. py:function:: trend_free_pre_whitening_modification_test(data_input, dim, alpha=0.05, returns_type='dataset_returns')


.. py:function:: seasonal_test(data_input, dim, alpha=0.05, period=12, returns_type='dataset_returns')


.. py:function:: regional_test(data_input, dim, alpha=0.05, returns_type='dataset_returns')


.. py:function:: correlated_seasonal_test(data_input, dim, alpha=0.05, period=12, returns_type='dataset_returns')


.. py:function:: partial_test(data_input, dim, alpha=0.05, returns_type='dataset_returns')


.. py:function:: sens_slope(data_input, dim, alpha=0.05, returns_type='dataset_returns')


.. py:function:: seasonal_sens_slope(data_input, dim, alpha=0.05, period=12, returns_type='dataset_returns')



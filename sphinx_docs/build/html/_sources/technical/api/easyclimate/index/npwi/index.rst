:py:mod:`easyclimate.index.npwi`
================================

.. py:module:: easyclimate.index.npwi


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.index.npwi.calc_index_NPWI
   easyclimate.index.npwi.monsoon_onsetdate_cal
   easyclimate.index.npwi.monsoon_detreatdate_cal
   easyclimate.index.npwi.monsoonregion



.. py:function:: calc_index_NPWI(precipitable_water_data)

   https://journals.ametsoc.org/view/journals/clim/17/11/1520-0442_2004_017_2241_gumoar_2.0.co_2.xml


.. py:function:: monsoon_onsetdate_cal(NPWI, thresh=0.618, rollingday=3, n=7)

   NPWI: 三维数组


.. py:function:: monsoon_detreatdate_cal(data, monsoon_onset, thresh=0.618, rollingday=3, n=7)

   NPWI: 三维数组


.. py:function:: monsoonregion(PW)



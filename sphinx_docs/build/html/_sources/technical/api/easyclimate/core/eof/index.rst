:py:mod:`easyclimate.core.eof`
==============================

.. py:module:: easyclimate.core.eof

.. autoapi-nested-parse::

   EOF MCA

   https://xeofs.readthedocs.io/en/latest/



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.core.eof.get_EOF_model
   easyclimate.core.eof.calc_EOF_analysis
   easyclimate.core.eof.get_EOF_projection
   easyclimate.core.eof.save_EOF_model
   easyclimate.core.eof.load_EOF_model
   easyclimate.core.eof.get_MCA_model
   easyclimate.core.eof.calc_MCA_analysis
   easyclimate.core.eof.get_MCA_projection
   easyclimate.core.eof.save_MCA_model
   easyclimate.core.eof.load_MCA_model



.. py:function:: get_EOF_model(data_input: xarray.DataArray, time_dim: str = 'time', n_modes=10, standardize=False, use_coslat=False, use_weights=False, weights=None, solver='auto', **solver_kwargs)

       
       


.. py:function:: calc_EOF_analysis(model: xeofs.models.eof.EOF, mini_summary=False)


.. py:function:: get_EOF_projection(model: xeofs.models.eof.EOF, data: xarray.DataArray)


.. py:function:: save_EOF_model(model: xeofs.models.eof.EOF, path: save_EOF_model.str)

       
       


.. py:function:: load_EOF_model(path: str)

       
       


.. py:function:: get_MCA_model(data_input1: xarray.DataArray, data_input2: xarray.DataArray, time_dim: str = 'time', n_modes=10, standardize=False, use_coslat=False, use_weights=False, weights1=None, weights2=None, n_pca_modes=None, solver='auto', **solver_kwargs)

       
       


.. py:function:: calc_MCA_analysis(model: xeofs.models.mca.MCA, correction=None, alpha=0.05, mini_summary=False)


.. py:function:: get_MCA_projection(model: xeofs.models.mca.MCA, **kwargs)

       
       


.. py:function:: save_MCA_model(model: xeofs.models.mca.MCA, path: save_MCA_model.str)

       
       


.. py:function:: load_MCA_model(path: str)

       
       



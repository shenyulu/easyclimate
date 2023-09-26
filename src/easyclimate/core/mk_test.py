"""
Mann-Kendall Test 
"""
import numpy as np
import xarray as xr
import pymannkendall as mk

from .utility import generate_datatree_dispatcher

@generate_datatree_dispatcher
def original_test(data_input, dim, alpha = 0.05, returns_type = 'dataset_returns'):
        
    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.original_test(data, alpha = alpha)

            trend = mk_test_result.trend
            if trend == 'increasing':
                trend = 1
            elif trend == 'decreasing':
                trend = -1
            elif trend == 'no trend':
                trend = 0
            else:
                raise ValueError('Error `trend` type.')
            
            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims = [["parameter"]],
            output_dtypes=["float64"],
            dask = "parallelized",
            vectorize=True,
            dask_gufunc_kwargs = {"output_sizes": {"parameter": 1}},
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars = {'trend': mk_test_dataarray[...,0],
                        'h':  mk_test_dataarray[...,1],
                        'p':  mk_test_dataarray[...,2],
                        'z':  mk_test_dataarray[...,3],
                        'Tau':  mk_test_dataarray[...,4],
                        's':  mk_test_dataarray[...,5],
                        'var_s':  mk_test_dataarray[...,6],
                        'slope':  mk_test_dataarray[...,7],
                        'intercept':  mk_test_dataarray[...,8],
            }
        )

        mk_test_dataset['trend'].attrs['Description'] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset['h'].attrs['Description'] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset['p'].attrs['Description'] = "p-value of the significance test"
        mk_test_dataset['z'].attrs['Description'] = "normalized test statistics"
        mk_test_dataset['Tau'].attrs['Description'] = "Kendall Tau"
        mk_test_dataset['s'].attrs['Description'] = "Mann-Kendal's score"
        mk_test_dataset['var_s'].attrs['Description'] = "Variance S"
        mk_test_dataset['slope'].attrs['Description'] = "Theil-Sen estimator/slope"
        mk_test_dataset['intercept'].attrs['Description'] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs['Description'] = "original_test"
        return mk_test_dataset
    
    return _top_test(data_input, dim, alpha)

@generate_datatree_dispatcher
def hamed_rao_modification_test(data_input, dim, alpha = 0.05, lag = None, returns_type = 'dataset_returns'):
        
    def _top_test(data_input, dim, alpha, lag = lag):
        def mk_test(data):
            mk_test_result = mk.hamed_rao_modification_test(data, alpha = alpha, lag = lag)

            trend = mk_test_result.trend
            if trend == 'increasing':
                trend = 1
            elif trend == 'decreasing':
                trend = -1
            elif trend == 'no trend':
                trend = 0
            else:
                raise ValueError('Error `trend` type.')
            
            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims = [["parameter"]],
            output_dtypes=["float64"],
            dask = "parallelized",
            vectorize=True,
            dask_gufunc_kwargs = {"output_sizes": {"parameter": 1}},
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars = {'trend': mk_test_dataarray[...,0],
                        'h':  mk_test_dataarray[...,1],
                        'p':  mk_test_dataarray[...,2],
                        'z':  mk_test_dataarray[...,3],
                        'Tau':  mk_test_dataarray[...,4],
                        's':  mk_test_dataarray[...,5],
                        'var_s':  mk_test_dataarray[...,6],
                        'slope':  mk_test_dataarray[...,7],
                        'intercept':  mk_test_dataarray[...,8],
            }
        )

        mk_test_dataset['trend'].attrs['Description'] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset['h'].attrs['Description'] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset['p'].attrs['Description'] = "p-value of the significance test"
        mk_test_dataset['z'].attrs['Description'] = "normalized test statistics"
        mk_test_dataset['Tau'].attrs['Description'] = "Kendall Tau"
        mk_test_dataset['s'].attrs['Description'] = "Mann-Kendal's score"
        mk_test_dataset['var_s'].attrs['Description'] = "Variance S"
        mk_test_dataset['slope'].attrs['Description'] = "Theil-Sen estimator/slope"
        mk_test_dataset['intercept'].attrs['Description'] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs['Description'] = "hamed_rao_modification_test"
        return mk_test_dataset
    
    return _top_test(data_input, dim, alpha, lag = lag)

@generate_datatree_dispatcher
def yue_wang_modification_test(data_input, dim, alpha = 0.05, lag = None, returns_type = 'dataset_returns'):
        
    def _top_test(data_input, dim, alpha, lag = lag):
        def mk_test(data):
            mk_test_result = mk.yue_wang_modification_test(data, alpha = alpha, lag = lag)

            trend = mk_test_result.trend
            if trend == 'increasing':
                trend = 1
            elif trend == 'decreasing':
                trend = -1
            elif trend == 'no trend':
                trend = 0
            else:
                raise ValueError('Error `trend` type.')
            
            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims = [["parameter"]],
            output_dtypes=["float64"],
            dask = "parallelized",
            vectorize=True,
            dask_gufunc_kwargs = {"output_sizes": {"parameter": 1}},
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars = {'trend': mk_test_dataarray[...,0],
                        'h':  mk_test_dataarray[...,1],
                        'p':  mk_test_dataarray[...,2],
                        'z':  mk_test_dataarray[...,3],
                        'Tau':  mk_test_dataarray[...,4],
                        's':  mk_test_dataarray[...,5],
                        'var_s':  mk_test_dataarray[...,6],
                        'slope':  mk_test_dataarray[...,7],
                        'intercept':  mk_test_dataarray[...,8],
            }
        )

        mk_test_dataset['trend'].attrs['Description'] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset['h'].attrs['Description'] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset['p'].attrs['Description'] = "p-value of the significance test"
        mk_test_dataset['z'].attrs['Description'] = "normalized test statistics"
        mk_test_dataset['Tau'].attrs['Description'] = "Kendall Tau"
        mk_test_dataset['s'].attrs['Description'] = "Mann-Kendal's score"
        mk_test_dataset['var_s'].attrs['Description'] = "Variance S"
        mk_test_dataset['slope'].attrs['Description'] = "Theil-Sen estimator/slope"
        mk_test_dataset['intercept'].attrs['Description'] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs['Description'] = "yue_wang_modification_test"
        return mk_test_dataset
    
    return _top_test(data_input, dim, alpha, lag = lag)

@generate_datatree_dispatcher
def pre_whitening_modification_test(data_input, dim, alpha = 0.05, returns_type = 'dataset_returns'):
        
    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.pre_whitening_modification_test(data, alpha = alpha)

            trend = mk_test_result.trend
            if trend == 'increasing':
                trend = 1
            elif trend == 'decreasing':
                trend = -1
            elif trend == 'no trend':
                trend = 0
            else:
                raise ValueError('Error `trend` type.')
            
            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims = [["parameter"]],
            output_dtypes=["float64"],
            dask = "parallelized",
            vectorize=True,
            dask_gufunc_kwargs = {"output_sizes": {"parameter": 1}},
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars = {'trend': mk_test_dataarray[...,0],
                        'h':  mk_test_dataarray[...,1],
                        'p':  mk_test_dataarray[...,2],
                        'z':  mk_test_dataarray[...,3],
                        'Tau':  mk_test_dataarray[...,4],
                        's':  mk_test_dataarray[...,5],
                        'var_s':  mk_test_dataarray[...,6],
                        'slope':  mk_test_dataarray[...,7],
                        'intercept':  mk_test_dataarray[...,8],
            }
        )

        mk_test_dataset['trend'].attrs['Description'] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset['h'].attrs['Description'] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset['p'].attrs['Description'] = "p-value of the significance test"
        mk_test_dataset['z'].attrs['Description'] = "normalized test statistics"
        mk_test_dataset['Tau'].attrs['Description'] = "Kendall Tau"
        mk_test_dataset['s'].attrs['Description'] = "Mann-Kendal's score"
        mk_test_dataset['var_s'].attrs['Description'] = "Variance S"
        mk_test_dataset['slope'].attrs['Description'] = "Theil-Sen estimator/slope"
        mk_test_dataset['intercept'].attrs['Description'] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs['Description'] = "pre_whitening_modification_test"
        return mk_test_dataset
    
    return _top_test(data_input, dim, alpha)

@generate_datatree_dispatcher
def trend_free_pre_whitening_modification_test(data_input, dim, alpha = 0.05, returns_type = 'dataset_returns'):
        
    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.trend_free_pre_whitening_modification_test(data, alpha = alpha)

            trend = mk_test_result.trend
            if trend == 'increasing':
                trend = 1
            elif trend == 'decreasing':
                trend = -1
            elif trend == 'no trend':
                trend = 0
            else:
                raise ValueError('Error `trend` type.')
            
            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims = [["parameter"]],
            output_dtypes=["float64"],
            dask = "parallelized",
            vectorize=True,
            dask_gufunc_kwargs = {"output_sizes": {"parameter": 1}},
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars = {'trend': mk_test_dataarray[...,0],
                        'h':  mk_test_dataarray[...,1],
                        'p':  mk_test_dataarray[...,2],
                        'z':  mk_test_dataarray[...,3],
                        'Tau':  mk_test_dataarray[...,4],
                        's':  mk_test_dataarray[...,5],
                        'var_s':  mk_test_dataarray[...,6],
                        'slope':  mk_test_dataarray[...,7],
                        'intercept':  mk_test_dataarray[...,8],
            }
        )

        mk_test_dataset['trend'].attrs['Description'] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset['h'].attrs['Description'] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset['p'].attrs['Description'] = "p-value of the significance test"
        mk_test_dataset['z'].attrs['Description'] = "normalized test statistics"
        mk_test_dataset['Tau'].attrs['Description'] = "Kendall Tau"
        mk_test_dataset['s'].attrs['Description'] = "Mann-Kendal's score"
        mk_test_dataset['var_s'].attrs['Description'] = "Variance S"
        mk_test_dataset['slope'].attrs['Description'] = "Theil-Sen estimator/slope"
        mk_test_dataset['intercept'].attrs['Description'] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs['Description'] = "trend_free_pre_whitening_modification_test"
        return mk_test_dataset
    
    return _top_test(data_input, dim, alpha)

@generate_datatree_dispatcher
def seasonal_test(data_input, dim, alpha = 0.05, period = 12, returns_type = 'dataset_returns'):
        
    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.seasonal_test(data, alpha = alpha, period = period)

            trend = mk_test_result.trend
            if trend == 'increasing':
                trend = 1
            elif trend == 'decreasing':
                trend = -1
            elif trend == 'no trend':
                trend = 0
            else:
                raise ValueError('Error `trend` type.')
            
            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims = [["parameter"]],
            output_dtypes=["float64"],
            dask = "parallelized",
            vectorize=True,
            dask_gufunc_kwargs = {"output_sizes": {"parameter": 1}},
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars = {'trend': mk_test_dataarray[...,0],
                        'h':  mk_test_dataarray[...,1],
                        'p':  mk_test_dataarray[...,2],
                        'z':  mk_test_dataarray[...,3],
                        'Tau':  mk_test_dataarray[...,4],
                        's':  mk_test_dataarray[...,5],
                        'var_s':  mk_test_dataarray[...,6],
                        'slope':  mk_test_dataarray[...,7],
                        'intercept':  mk_test_dataarray[...,8],
            }
        )

        mk_test_dataset['trend'].attrs['Description'] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset['h'].attrs['Description'] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset['p'].attrs['Description'] = "p-value of the significance test"
        mk_test_dataset['z'].attrs['Description'] = "normalized test statistics"
        mk_test_dataset['Tau'].attrs['Description'] = "Kendall Tau"
        mk_test_dataset['s'].attrs['Description'] = "Mann-Kendal's score"
        mk_test_dataset['var_s'].attrs['Description'] = "Variance S"
        mk_test_dataset['slope'].attrs['Description'] = "Theil-Sen estimator/slope"
        mk_test_dataset['intercept'].attrs['Description'] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs['Description'] = "seasonal_test"
        return mk_test_dataset
    
    return _top_test(data_input, dim, alpha, period)

@generate_datatree_dispatcher
def regional_test(data_input, dim, alpha = 0.05, returns_type = 'dataset_returns'):
        
    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.regional_test(data, alpha = alpha)

            trend = mk_test_result.trend
            if trend == 'increasing':
                trend = 1
            elif trend == 'decreasing':
                trend = -1
            elif trend == 'no trend':
                trend = 0
            else:
                raise ValueError('Error `trend` type.')
            
            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims = [["parameter"]],
            output_dtypes=["float64"],
            dask = "parallelized",
            vectorize=True,
            dask_gufunc_kwargs = {"output_sizes": {"parameter": 1}},
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars = {'trend': mk_test_dataarray[...,0],
                        'h':  mk_test_dataarray[...,1],
                        'p':  mk_test_dataarray[...,2],
                        'z':  mk_test_dataarray[...,3],
                        'Tau':  mk_test_dataarray[...,4],
                        's':  mk_test_dataarray[...,5],
                        'var_s':  mk_test_dataarray[...,6],
                        'slope':  mk_test_dataarray[...,7],
                        'intercept':  mk_test_dataarray[...,8],
            }
        )

        mk_test_dataset['trend'].attrs['Description'] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset['h'].attrs['Description'] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset['p'].attrs['Description'] = "p-value of the significance test"
        mk_test_dataset['z'].attrs['Description'] = "normalized test statistics"
        mk_test_dataset['Tau'].attrs['Description'] = "Kendall Tau"
        mk_test_dataset['s'].attrs['Description'] = "Mann-Kendal's score"
        mk_test_dataset['var_s'].attrs['Description'] = "Variance S"
        mk_test_dataset['slope'].attrs['Description'] = "Theil-Sen estimator/slope"
        mk_test_dataset['intercept'].attrs['Description'] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs['Description'] = "regional_test"
        return mk_test_dataset
    
    return _top_test(data_input, dim, alpha)

@generate_datatree_dispatcher
def correlated_seasonal_test(data_input, dim, alpha = 0.05, period = 12, returns_type = 'dataset_returns'):
        
    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.correlated_seasonal_test(data, alpha = alpha, period = period)

            trend = mk_test_result.trend
            if trend == 'increasing':
                trend = 1
            elif trend == 'decreasing':
                trend = -1
            elif trend == 'no trend':
                trend = 0
            else:
                raise ValueError('Error `trend` type.')
            
            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims = [["parameter"]],
            output_dtypes=["float64"],
            dask = "parallelized",
            vectorize=True,
            dask_gufunc_kwargs = {"output_sizes": {"parameter": 1}},
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars = {'trend': mk_test_dataarray[...,0],
                        'h':  mk_test_dataarray[...,1],
                        'p':  mk_test_dataarray[...,2],
                        'z':  mk_test_dataarray[...,3],
                        'Tau':  mk_test_dataarray[...,4],
                        's':  mk_test_dataarray[...,5],
                        'var_s':  mk_test_dataarray[...,6],
                        'slope':  mk_test_dataarray[...,7],
                        'intercept':  mk_test_dataarray[...,8],
            }
        )

        mk_test_dataset['trend'].attrs['Description'] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset['h'].attrs['Description'] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset['p'].attrs['Description'] = "p-value of the significance test"
        mk_test_dataset['z'].attrs['Description'] = "normalized test statistics"
        mk_test_dataset['Tau'].attrs['Description'] = "Kendall Tau"
        mk_test_dataset['s'].attrs['Description'] = "Mann-Kendal's score"
        mk_test_dataset['var_s'].attrs['Description'] = "Variance S"
        mk_test_dataset['slope'].attrs['Description'] = "Theil-Sen estimator/slope"
        mk_test_dataset['intercept'].attrs['Description'] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs['Description'] = "correlated_seasonal_test"
        return mk_test_dataset
    
    return _top_test(data_input, dim, alpha, period)

@generate_datatree_dispatcher
def partial_test(data_input, dim, alpha = 0.05, returns_type = 'dataset_returns'):
        
    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.partial_test(data, alpha = alpha)

            trend = mk_test_result.trend
            if trend == 'increasing':
                trend = 1
            elif trend == 'decreasing':
                trend = -1
            elif trend == 'no trend':
                trend = 0
            else:
                raise ValueError('Error `trend` type.')
            
            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims = [["parameter"]],
            output_dtypes=["float64"],
            dask = "parallelized",
            vectorize=True,
            dask_gufunc_kwargs = {"output_sizes": {"parameter": 1}},
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars = {'trend': mk_test_dataarray[...,0],
                        'h':  mk_test_dataarray[...,1],
                        'p':  mk_test_dataarray[...,2],
                        'z':  mk_test_dataarray[...,3],
                        'Tau':  mk_test_dataarray[...,4],
                        's':  mk_test_dataarray[...,5],
                        'var_s':  mk_test_dataarray[...,6],
                        'slope':  mk_test_dataarray[...,7],
                        'intercept':  mk_test_dataarray[...,8],
            }
        )

        mk_test_dataset['trend'].attrs['Description'] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset['h'].attrs['Description'] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset['p'].attrs['Description'] = "p-value of the significance test"
        mk_test_dataset['z'].attrs['Description'] = "normalized test statistics"
        mk_test_dataset['Tau'].attrs['Description'] = "Kendall Tau"
        mk_test_dataset['s'].attrs['Description'] = "Mann-Kendal's score"
        mk_test_dataset['var_s'].attrs['Description'] = "Variance S"
        mk_test_dataset['slope'].attrs['Description'] = "Theil-Sen estimator/slope"
        mk_test_dataset['intercept'].attrs['Description'] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs['Description'] = "partial_test"
        return mk_test_dataset
    
    return _top_test(data_input, dim, alpha)

@generate_datatree_dispatcher
def sens_slope(data_input, dim, alpha = 0.05, returns_type = 'dataset_returns'):
        
    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.sens_slope(data, alpha = alpha)

            trend = mk_test_result.trend
            if trend == 'increasing':
                trend = 1
            elif trend == 'decreasing':
                trend = -1
            elif trend == 'no trend':
                trend = 0
            else:
                raise ValueError('Error `trend` type.')
            
            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims = [["parameter"]],
            output_dtypes=["float64"],
            dask = "parallelized",
            vectorize=True,
            dask_gufunc_kwargs = {"output_sizes": {"parameter": 1}},
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars = {'trend': mk_test_dataarray[...,0],
                        'h':  mk_test_dataarray[...,1],
                        'p':  mk_test_dataarray[...,2],
                        'z':  mk_test_dataarray[...,3],
                        'Tau':  mk_test_dataarray[...,4],
                        's':  mk_test_dataarray[...,5],
                        'var_s':  mk_test_dataarray[...,6],
                        'slope':  mk_test_dataarray[...,7],
                        'intercept':  mk_test_dataarray[...,8],
            }
        )

        mk_test_dataset['trend'].attrs['Description'] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset['h'].attrs['Description'] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset['p'].attrs['Description'] = "p-value of the significance test"
        mk_test_dataset['z'].attrs['Description'] = "normalized test statistics"
        mk_test_dataset['Tau'].attrs['Description'] = "Kendall Tau"
        mk_test_dataset['s'].attrs['Description'] = "Mann-Kendal's score"
        mk_test_dataset['var_s'].attrs['Description'] = "Variance S"
        mk_test_dataset['slope'].attrs['Description'] = "Theil-Sen estimator/slope"
        mk_test_dataset['intercept'].attrs['Description'] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs['Description'] = "sens_slope"
        return mk_test_dataset
    
    return _top_test(data_input, dim, alpha)

@generate_datatree_dispatcher
def seasonal_sens_slope(data_input, dim, alpha = 0.05, period = 12, returns_type = 'dataset_returns'):
        
    def _top_test(data_input, dim, alpha):
        def mk_test(data):
            mk_test_result = mk.seasonal_sens_slope(data, alpha = alpha, period = period)

            trend = mk_test_result.trend
            if trend == 'increasing':
                trend = 1
            elif trend == 'decreasing':
                trend = -1
            elif trend == 'no trend':
                trend = 0
            else:
                raise ValueError('Error `trend` type.')
            
            h = mk_test_result.h
            p = mk_test_result.p
            z = mk_test_result.z
            Tau = mk_test_result.Tau
            s = mk_test_result.s
            var_s = mk_test_result.var_s
            slope = mk_test_result.slope
            intercept = mk_test_result.intercept
            return np.array([trend, h, p, z, Tau, s, var_s, slope, intercept])

        # Use xarray apply_ufunc to create DataArray
        mk_test_dataarray = xr.apply_ufunc(
            mk_test,
            data_input,
            input_core_dims=[[dim]],
            output_core_dims = [["parameter"]],
            output_dtypes=["float64"],
            dask = "parallelized",
            vectorize=True,
            dask_gufunc_kwargs = {"output_sizes": {"parameter": 1}},
        )

        # Transform DataArray to Dataset
        mk_test_dataset = xr.Dataset(
            data_vars = {'trend': mk_test_dataarray[...,0],
                        'h':  mk_test_dataarray[...,1],
                        'p':  mk_test_dataarray[...,2],
                        'z':  mk_test_dataarray[...,3],
                        'Tau':  mk_test_dataarray[...,4],
                        's':  mk_test_dataarray[...,5],
                        'var_s':  mk_test_dataarray[...,6],
                        'slope':  mk_test_dataarray[...,7],
                        'intercept':  mk_test_dataarray[...,8],
            }
        )

        mk_test_dataset['trend'].attrs['Description'] = "tells the trend (increasing: 1, decreasing: -1 or no trend: 0)."
        mk_test_dataset['h'].attrs['Description'] = "True (if trend is present) or False (if the trend is absence)."
        mk_test_dataset['p'].attrs['Description'] = "p-value of the significance test"
        mk_test_dataset['z'].attrs['Description'] = "normalized test statistics"
        mk_test_dataset['Tau'].attrs['Description'] = "Kendall Tau"
        mk_test_dataset['s'].attrs['Description'] = "Mann-Kendal's score"
        mk_test_dataset['var_s'].attrs['Description'] = "Variance S"
        mk_test_dataset['slope'].attrs['Description'] = "Theil-Sen estimator/slope"
        mk_test_dataset['intercept'].attrs['Description'] = "intercept of Kendall-Theil Robust Line, for seasonal test, full period cycle consider as unit time step"
        mk_test_dataset.attrs['Description'] = "seasonal_sens_slope"
        return mk_test_dataset
    
    return _top_test(data_input, dim, alpha, period)
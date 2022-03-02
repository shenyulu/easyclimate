import numpy as np
import scipy.stats
import xarray as xr

def linregress(data_input, dim = 'time'):
    x = np.arange(0, data_input[dim].shape[0])

    def linregress_scipy(data):
        LinregressResult = scipy.stats.linregress(x, data)
        slope = LinregressResult.slope
        intercept = LinregressResult.intercept
        rvalue = LinregressResult.rvalue
        pvalue = LinregressResult.pvalue
        stderr = LinregressResult.stderr
        intercept_stderr = LinregressResult.intercept_stderr
        return np.array([slope, intercept, rvalue, pvalue, stderr, intercept_stderr])

    LinregressResult_dataarray = xr.apply_ufunc(
        linregress_scipy,
        data_input,
        input_core_dims=[[dim]],
        output_core_dims = [["parameter"]],
        output_dtypes=["float64"],
        dask = "parallelized",
        vectorize=True,
        dask_gufunc_kwargs = {"output_sizes": {"parameter": 6}},
    )

    return xr.Dataset(
        data_vars = {'slope': LinregressResult_dataarray[:,:,0],
                     'intercept':  LinregressResult_dataarray[:,:,1],
                     'rvalue':  LinregressResult_dataarray[:,:,2],
                     'pvalue':  LinregressResult_dataarray[:,:,3],
                     'stderr':  LinregressResult_dataarray[:,:,4],
                     'intercept_stderr':  LinregressResult_dataarray[:,:,5],
        }
    )
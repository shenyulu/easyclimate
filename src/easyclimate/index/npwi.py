import xarray as xr
import numpy as np
# https://journals.ametsoc.org/view/journals/clim/17/11/1520-0442_2004_017_2241_gumoar_2.0.co_2.xml


# NPVI_op = NPVI_cal(qvi_op)
# monsoon_region_op = monsoonregion(qvi_op)
# monsoon_onset_op = monsoon_onsetdate_cal(NPVI_op)
# monsoon_detreat_op = monsoon_detreatdate_cal(NPVI_op, monsoon_onset_op)

# monsoon_onset_op_draw = monsoon_onset_op.where(monsoon_region_op).compute()
# monsoon_detreat_op_draw = monsoon_detreat_op.where(monsoon_region_op).compute()

def calc_index_NPWI(precipitable_water_data):
    """
    https://journals.ametsoc.org/view/journals/clim/17/11/1520-0442_2004_017_2241_gumoar_2.0.co_2.xml
    """
    # 计算PWmax, PWmin
    PWmax = precipitable_water_data.groupby('time.year').max(dim = 'time').mean(dim = 'year')
    PWmin = precipitable_water_data.groupby('time.year').min(dim = 'time').mean(dim = 'year')

    # 计算NPWI
    NPWI = (precipitable_water_data - PWmin)/(PWmax - PWmin)
    return NPWI

def monsoon_onsetdate_cal(NPWI, thresh = 0.618, rollingday = 3, n = 7):
    '''
    NPWI: 三维数组
    '''
    # 找到NPWI 大于thresh 的部分
    NPWI_overthresh = (NPWI > thresh)

    # 寻找到符合连续3 天满足大于thresh 的部分
    NPWI_overthresh_conti3 = NPWI_overthresh.rolling(time = rollingday, center=True).sum()     # NPWI_overthresh_conti3 = NPWI_overthresh.rolling(time = 3, center=True).sum()
    NPWI_overthresh_over3 = (NPWI_overthresh_conti3 > (rollingday-1)).astype(int)              # NPWI_overthresh_over3 = (NPWI_overthresh_conti3 > 2).astype(int)
    
    # 9 个网格单元中需满足n个网格单元满足上述条件
    monsoonindex = NPWI_overthresh_over3.rolling({'lon': 3, 'lat': 3}, center=True).sum()
    monsoonindex_start = (monsoonindex >= n).astype(int)

    # 获取日期在各年的第几日，并创建数组
    lonshape = NPWI['lon'].shape[0]
    latshape = NPWI['lat'].shape[0]

    NPWI_day = NPWI['time'].to_pandas().dt.day_of_year.values
    NPWI_timelength = NPWI_day.shape[0]
    timearray = np.broadcast_to(NPWI_day, ((latshape, lonshape,NPWI_timelength)))
    dateofyear_dataarray = xr.DataArray(timearray,
                                        coords={'lat': NPWI['lat'], 'lon': NPWI['lon'], 'time': NPWI['time']}
                                        ).transpose('time', 'lat', 'lon')
    
    # 找到具体日期，并对取各年中日期的最小值，并对诸年取中位数
    dateofyear_NPWI = dateofyear_dataarray.where(monsoonindex_start > 0, drop = True).squeeze()
    return dateofyear_NPWI.groupby('time.year').min().median(dim = 'year')

def monsoon_detreatdate_cal(data, monsoon_onset, thresh = 0.618, rollingday = 3, n = 7):
    '''
    NPWI: 三维数组
    '''
    # 忽略1,2,3月
    # months = data.time.dt.month
    # month_idxs = months.isin([7, 8, 9, 10, 11, 12])
    NPWI = data #.isel(time = month_idxs)
    # 找到NPWI 大于thresh 的部分
    NPWI_overthresh = (NPWI < thresh)

    # 寻找到符合连续3 天满足大于thresh 的部分
    NPWI_overthresh_conti3 = NPWI_overthresh.rolling(time = rollingday, center=True).sum()     # NPWI_overthresh_conti3 = NPWI_overthresh.rolling(time = 3, center=True).sum()
    NPWI_overthresh_over3 = (NPWI_overthresh_conti3 > (rollingday-1)).astype(int)              # NPWI_overthresh_over3 = (NPWI_overthresh_conti3 > 2).astype(int)
    
    # 9 个网格单元中需满足n个网格单元满足上述条件
    monsoonindex = NPWI_overthresh_over3.rolling({'lon': 3, 'lat': 3}, center=True).sum()
    monsoonindex_start = (monsoonindex >= n).astype(int)

    # 获取日期在各年的第几日，并创建数组
    lonshape = NPWI['lon'].shape[0]
    latshape = NPWI['lat'].shape[0]

    NPWI_day = NPWI['time'].to_pandas().dt.day_of_year.values
    NPWI_timelength = NPWI_day.shape[0]
    timearray = np.broadcast_to(NPWI_day, ((latshape, lonshape,NPWI_timelength)))
    dateofyear_dataarray = xr.DataArray(timearray,
                                        coords={'lat': NPWI['lat'], 'lon': NPWI['lon'], 'time': NPWI['time']}
                                        ).transpose('time', 'lat', 'lon')
    dateofyear_dataarray_overonset = dateofyear_dataarray.where((dateofyear_dataarray > monsoon_onset))
    
    # 找到具体日期，并对取各年中日期的最小值，并对诸年取中位数
    dateofyear_NPWI = dateofyear_dataarray_overonset.where(monsoonindex_start > 0, drop = True).squeeze()
    return dateofyear_NPWI.groupby('time.year').min().median(dim = 'year')


def monsoonregion(PW):
    PW_averagemonth = PW.groupby('time.month').mean(dim = 'time')
    PWw = PW_averagemonth.sel(month = [6, 7, 8]).max(dim = 'month')
    PWc = PW_averagemonth.sel(month = [1, 2, 12]).max(dim = 'month')
    monsoonregionmask = ((PWw - PWc) > 12)
    return monsoonregionmask
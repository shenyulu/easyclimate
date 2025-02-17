"""
pytest for filter/redfitpy.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
from pathlib import Path
from .const_define import TEST_DATA_PATH

sst_data = xr.open_dataset(str(Path(TEST_DATA_PATH, "test_input_nino34.nc"))).sst

# def test_calc_redfit():
#     result_data = ecl.filter.calc_redfit(sst_data)

#     result_data1 = result_data.gxx[20:30].data
#     # result_data2 = result_data.gxx_corr[20:30].data
#     result_data3 = result_data.gred_th[20:30].data
#     result_data4 = result_data.gred[20:30].data
#     result_data5 = result_data.corrFac[20:30].data
#     result_data6 = result_data.chi2_80[20:30].data
#     result_data7 = result_data.chi2_90[20:30].data
#     result_data8 = result_data.chi2_95[20:30].data
#     result_data9 = result_data.chi2_99[20:30].data

#     refer_data1 = np.array([ 2.6100595 ,  1.3072749 ,  7.6913514 ,  1.6619556 ,  5.278525  ,
#         3.2057655 ,  0.66721755, 17.080835  , 15.754893  ,  4.182542  ])
#     # refer_data2 = np.array([ 2.0198553,  1.0682539,  6.2757936,  1.3678471,  4.170987 ,
#     #     2.510872 ,  0.5208666, 13.633644 , 12.198356 ,  3.278343 ])
#     refer_data3 = np.array([4.374933 , 3.9981651, 3.6673093, 3.3753517, 3.1165347, 2.8861244,
#        2.6801665, 2.4953728, 2.3289793, 2.1786544])
#     refer_data4 = np.array([5.4194856, 5.2379966, 4.6552963, 4.325109 , 3.8132083, 3.6384985,
#        3.162897 , 3.1781259, 3.0590422, 2.7382824])
#     refer_data5 = np.array([1.2387586, 1.3101001, 1.2694038, 1.2813802, 1.2235411, 1.2606866,
#        1.1801122, 1.2736077, 1.313469 , 1.2568686])
#     refer_data6 = np.array([7.0427766, 6.436255 , 5.9036427, 5.4336486, 5.0170045, 4.6460896,
#        4.314538 , 4.0170565, 3.749196 , 3.5072029])
#     refer_data7 = np.array([10.071838 ,  9.204455 ,  8.442769 ,  7.7706327,  7.1747923,
#         6.644348 ,  6.1701984,  5.744772 ,  5.361706 ,  5.0156326])
#     refer_data8 = np.array([13.109244 , 11.980281 , 10.98889  , 10.114055 ,  9.338524 ,
#         8.648112 ,  8.030971 ,  7.4772468,  6.9786577,  6.5282173])
#     refer_data9 = np.array([20.143677 , 18.40891  , 16.885538 , 15.5412655, 14.349585 ,
#        13.288696 , 12.340397 , 11.489544 , 10.723412 , 10.031265 ])

#     assert np.isclose(result_data1, refer_data1, atol=0.01).all()
#     # assert np.isclose(result_data2, refer_data2, atol=0.01).all()
#     assert np.isclose(result_data3, refer_data3, atol=0.01).all()
#     # assert np.isclose(result_data4, refer_data4, atol=0.01).all()
#     # assert np.isclose(result_data5, refer_data5, atol=0.01).all()
#     # assert np.isclose(result_data6, refer_data6, atol=0.01).all()
#     assert np.isclose(result_data7, refer_data7, atol=0.01).all()
#     assert np.isclose(result_data8, refer_data8, atol=0.01).all()
#     assert np.isclose(result_data9, refer_data9, atol=0.01).all()

# def test_calc_redfit_cross():
#     gxx, gyy, gxy, cxy, phxy = ecl.filter.calc_redfit_cross(sst_data.data[:20], sst_data.data[:20])

#     # ------------------------------------
#     result_gxx_data1 = gxx.gxx[5:].data
#     result_gxx_data2 = gxx.gxx_corr[5:].data
#     result_gxx_data3 = gxx.gred_th[5:].data
#     result_gxx_data4 = gxx.gred[5:].data
#     result_gxx_data5 = gxx.corrFac[5:].data
#     result_gxx_data6 = gxx.chi2_90[5:].data
#     result_gxx_data7 = gxx.chi2_95[5:].data
#     result_gxx_data8 = gxx.chi2_99[5:].data
#     result_gxx_data9 = gxx.ci90[5:].data
#     result_gxx_data10 = gxx.ci95[5:].data
#     result_gxx_data11= gxx.ci99[5:].data

#     refer_gxx_data1 = np.array([0.49393213, 0.04329514, 0.10684587, 0.03644858, 0.0762934 ,
#        0.00138481])
#     refer_gxx_data2 = np.array([0.02960541, 0.00239743, 0.00609645, 0.00213968, 0.0042054 ,
#        0.00017695])
#     refer_gxx_data3 = np.array([0.0501421 , 0.03834052, 0.03162566, 0.0277664 , 0.02574903,
#        0.02512014])
#     refer_gxx_data4 = np.array([0.8365632 , 0.692392  , 0.5542689 , 0.47298875, 0.46713322,
#        0.1965879 ])
#     refer_gxx_data5 = np.array([16.683847 , 18.059011 , 17.525925 , 17.034569 , 18.141779 ,
#         7.8259077])
#     refer_gxx_data6 = np.array([0.11543563, 0.0882664 , 0.07280763, 0.06392298, 0.05927864,
#        0.05783083])
#     refer_gxx_data7 = np.array([0.15024804, 0.11488526, 0.09476454, 0.0832005 , 0.07715555,
#        0.07527111])
#     refer_gxx_data8 = np.array([0.23087126, 0.17653279, 0.14561526, 0.12784596, 0.11855727,
#        0.11566166])
#     refer_gxx_data9 = np.array([0.1179627 , 0.09038644, 0.07605662, 0.0689447 , 0.06384886,
#        0.07242694])
#     refer_gxx_data10 = np.array([0.1542198 , 0.12090861, 0.10014802, 0.08840656, 0.08543347,
#        0.09888934])
#     refer_gxx_data11 = np.array([0.23585491, 0.21016547, 0.16784589, 0.15234832, 0.13953124,
#        0.16913408])

#     # ------------------------------------
#     result_gyy_data1 = gyy.gyy[5:].data
#     result_gyy_data2 = gyy.gyy_corr[5:].data
#     result_gyy_data3 = gyy.gred_th[5:].data
#     result_gyy_data4 = gyy.gred[5:].data
#     result_gyy_data5 = gyy.corrFac[5:].data
#     result_gyy_data6 = gyy.chi2_90[5:].data
#     result_gyy_data7 = gyy.chi2_95[5:].data
#     result_gyy_data8 = gyy.chi2_99[5:].data
#     result_gyy_data9 = gyy.ci90[5:].data
#     result_gyy_data10 = gyy.ci95[5:].data
#     result_gyy_data11= gyy.ci99[5:].data

#     refer_gyy_data1 = np.array([0.49393213, 0.04329514, 0.10684587, 0.03644858, 0.0762934 ,
#        0.00138481])
#     refer_gyy_data2 = np.array([0.02958369, 0.00250113, 0.00608998, 0.00213394, 0.00431682,
#        0.00016488])
#     refer_gyy_data3 = np.array([0.0501421 , 0.03834052, 0.03162566, 0.0277664 , 0.02574903,
#        0.02512014])
#     refer_gyy_data4 = np.array([0.8371773 , 0.66368324, 0.55485713, 0.47426134, 0.45507622,
#        0.21098267])
#     refer_gyy_data5 = np.array([16.696095, 17.310228, 17.544523, 17.080402, 17.673529,  8.398944])
#     refer_gyy_data6 = np.array([0.11543563, 0.0882664 , 0.07280763, 0.06392298, 0.05927864,
#        0.05783083])
#     refer_gyy_data7 = np.array([0.15024804, 0.11488526, 0.09476454, 0.0832005 , 0.07715555,
#        0.07527111])
#     refer_gyy_data8 = np.array([0.23087126, 0.17653279, 0.14561526, 0.12784596, 0.11855727,
#        0.11566166])
#     refer_gyy_data9 = np.array([0.12371099, 0.09811879, 0.07792045, 0.06498447, 0.06117768,
#        0.06679989])
#     refer_gyy_data10 = np.array([0.15686867, 0.12784587, 0.10917611, 0.0917956 , 0.08398537,
#        0.09970414])
#     refer_gyy_data11 = np.array([0.2289167 , 0.19085781, 0.16979153, 0.17090872, 0.15318893,
#        0.18701449])

#     # ------------------------------------
#     result_gxy_data1 = gxy.gxy[5:].data

#     refer_gxy_data1 = np.array([0.49393213, 0.04329514, 0.10684587, 0.03644858, 0.0762934 ,
#        0.00138481])

#     # ------------------------------------
#     result_cxy_data1 = cxy.cxy[5:].data
#     result_cxy_data2 = cxy.csig[5:].data
#     result_cxy_data3 = cxy.csig_mc[5:].data
#     result_cxy_data4 = cxy.ci90[5:].data
#     result_cxy_data5 = cxy.ci95[5:].data
#     result_cxy_data6= cxy.ci99[5:].data

#     refer_cxy_data1 = np.array([1., 1., 1., 1., 1., 1.])
#     refer_cxy_data2 = np.array([0., 0., 0., 0., 0., 0.])
#     refer_cxy_data3 = np.array([1., 1., 1., 1., 1., 1.])
#     refer_cxy_data4 = np.array([1., 1., 1., 1., 1., 1.])
#     refer_cxy_data5 = np.array([1., 1., 1., 1., 1., 1.])
#     refer_cxy_data6 = np.array([1., 1., 1., 1., 1., 1.])

#     # ------------------------------------
#     result_phxy_data1 = phxy.phxy[5:].data
#     result_phxy_data2 = phxy.ephi_lower[5:].data
#     result_phxy_data3 = phxy.ephi_upper[5:].data
#     result_phxy_data4 = phxy.ephi_mc_lower[5:].data
#     result_phxy_data5 = phxy.ephi_mc_upper[5:].data

#     refer_phxy_data1 = np.array([0., 0., 0., 0., 0., 0.])
#     refer_phxy_data2 = np.array([0., 0., 0., 0., 0., 0.])
#     refer_phxy_data3 = np.array([0., 0., 0., 0., 0., 0.])
#     refer_phxy_data4 = np.array([0., 0., 0., 0., 0., 0.])
#     refer_phxy_data5 = np.array([0., 0., 0., 0., 0., 0.])

#     assert np.isclose(result_gxx_data1, refer_gxx_data1, atol=0.01).all()
#     assert np.isclose(result_gxx_data2, refer_gxx_data2, atol=0.01).all()
#     assert np.isclose(result_gxx_data3, refer_gxx_data3, atol=0.01).all()
#     assert np.isclose(result_gxx_data4, refer_gxx_data4, atol=0.01).all()
#     assert np.isclose(result_gxx_data5, refer_gxx_data5, atol=0.01).all()
#     assert np.isclose(result_gxx_data6, refer_gxx_data6, atol=0.01).all()
#     assert np.isclose(result_gxx_data7, refer_gxx_data7, atol=0.01).all()
#     assert np.isclose(result_gxx_data8, refer_gxx_data8, atol=0.01).all()
#     assert np.isclose(result_gxx_data9, refer_gxx_data9, atol=0.01).all()
#     assert np.isclose(result_gxx_data10, refer_gxx_data10, atol=0.01).all()
#     assert np.isclose(result_gxx_data11, refer_gxx_data11, atol=0.01).all()

#     assert np.isclose(result_gyy_data1, refer_gyy_data1, atol=0.01).all()
#     assert np.isclose(result_gyy_data2, refer_gyy_data2, atol=0.01).all()
#     assert np.isclose(result_gyy_data3, refer_gyy_data3, atol=0.01).all()
#     assert np.isclose(result_gyy_data4, refer_gyy_data4, atol=0.01).all()
#     assert np.isclose(result_gyy_data5, refer_gyy_data5, atol=0.01).all()
#     assert np.isclose(result_gyy_data6, refer_gyy_data6, atol=0.01).all()
#     assert np.isclose(result_gyy_data7, refer_gyy_data7, atol=0.01).all()
#     assert np.isclose(result_gyy_data8, refer_gyy_data8, atol=0.01).all()
#     assert np.isclose(result_gyy_data9, refer_gyy_data9, atol=0.01).all()
#     assert np.isclose(result_gyy_data10, refer_gyy_data10, atol=0.01).all()
#     assert np.isclose(result_gyy_data11, refer_gyy_data11, atol=0.01).all()

#     assert np.isclose(result_gxy_data1, refer_gxy_data1, atol=0.01).all()

#     assert np.isclose(result_cxy_data1, refer_cxy_data1, atol=0.01).all()
#     assert np.isclose(result_cxy_data2, refer_cxy_data2, atol=0.01).all()
#     assert np.isclose(result_cxy_data3, refer_cxy_data3, atol=0.01).all()
#     assert np.isclose(result_cxy_data4, refer_cxy_data4, atol=0.01).all()
#     assert np.isclose(result_cxy_data5, refer_cxy_data5, atol=0.01).all()
#     assert np.isclose(result_cxy_data6, refer_cxy_data6, atol=0.01).all()

#     assert np.isclose(result_phxy_data1, refer_phxy_data1, atol=0.01).all()
#     assert np.isclose(result_phxy_data2, refer_phxy_data2, atol=0.01).all()
#     assert np.isclose(result_phxy_data3, refer_phxy_data3, atol=0.01).all()
#     assert np.isclose(result_phxy_data4, refer_phxy_data4, atol=0.01).all()
#     assert np.isclose(result_phxy_data5, refer_phxy_data5, atol=0.01).all()

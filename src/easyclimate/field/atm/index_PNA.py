"""
PNA Index

https://psl.noaa.gov/enso/dashboard.html
"""

# from xmca.xarray import xMCA
# from ..core.variability import remove_seasonal_cycle_mean

# def calc_index_PNA(z500_monthly_data, time_dim, lat_dim = 'lat', remove_seasonal_cycle = True, save_analysis_path = None, load_analysis_path = None):
#     """

#     """
#     # Extracting Northern Hemisphere Data
#     z500_monthly_data = z500_monthly_data.sortby(lat_dim).sel(lat = slice(0, 90))

#     if load_analysis_path != None:
#         mca = xMCA()
#         mca.load_analysis(load_analysis_path)
#     else:
#         if remove_seasonal_cycle == True:
#             # Remove seasonal cycle
#             z500_monthly_data_noseason_cycle = remove_seasonal_cycle_mean(z500_monthly_data, time_dim)
#         else:
#             z500_monthly_data_noseason_cycle = z500_monthly_data

#         # calculate empirical orthogonal function
#         mca = xMCA(z500_monthly_data_noseason_cycle)
#         mca.apply_coslat()
#         mca.solve()
#         mca.rotate(n_rot=2, power=1)

#     if save_analysis_path != None:
#         mca.save_analysis(save_analysis_path)

#     pattern_PNA = mca.eofs()['left'].sel(mode = 2)
#     pattern_PNA.name = 'PNA_pattern'

#     pc_PNA = mca.pcs()['left'].sel(mode = 2)
#     pc_PNA.name = 'PNA_index'

#     explained_variance_PNA = mca.explained_variance().sel(mode = 2)

#     return [pc_PNA, pattern_PNA, explained_variance_PNA]

# Changelog♻️

## v2025.12.0

*Released on: 2025/11/29*

* [feat: add ts on normalized funcs](https://github.com/shenyulu/easyclimate/pull/108/commits/df423a3b00bd8b20e3099c8091c6280bce6e6d75)

> ⚠️ From this version, we change the functions name as follows:
>
> - "normalize_zscore" ➡️ "timeseries_normalize_zscore"
> - "normalize_minmax" ➡️ "timeseries_normalize_minmax"
> - "normalize_robust" ➡️ "timeseries_normalize_robust"
> - "normalize_mean" ➡️ "timeseries_normalize_mean"
>
> and add function:
>
> - "calc_precip_anomaly_percentage"

* [fix: use get_axis_num instead](https://github.com/shenyulu/easyclimate/pull/108/commits/60bc622b68299c1bcf17573a13e9dea8b0647188)
* [feat: delete climatological macro](https://github.com/shenyulu/easyclimate/pull/108/commits/42b83e09047635bd2a8a419c8d613b231be95081)

> ⚠️ From this version, we change the functions name as follows:
>
> - "calc_yearly_mean" ➡️ "calc_yearly_climatological_mean"
> - "calc_yearly_sum" ➡️ "calc_yearly_climatological_sum"
> - "calc_yearly_std" ➡️ "calc_yearly_climatological_std"
> - "calc_yearly_var" ➡️ "calc_yearly_climatological_var"
> - "calc_yearly_max" ➡️ "calc_yearly_climatological_max"
> - "calc_yearly_min" ➡️ "calc_yearly_climatological_min"

* [feat: add sample units for q](https://github.com/shenyulu/easyclimate/pull/108/commits/416128a1531d64c8cdb08ebe2338e2e64f4fd7e0)
* [feat: seperate units module](https://github.com/shenyulu/easyclimate/pull/108/commits/cce8c9b66956a5bed14b0277f547b60cae1485cb)
* [fix: radius of the earth](https://github.com/shenyulu/easyclimate/pull/108/commits/8099f2299f7345826374c22c8a2e0be694dfd64e)
* [feat: add monthstat](https://github.com/shenyulu/easyclimate/pull/108/commits/171de4110c68cc3c2c8bd1913ef094d42806dcf8)

> ⚠️ From this version, we add the functions name as follows:
>
> - "calc_monthly_mean"
> - "calc_monthly_sum"
> - "calc_monthly_std"
> - "calc_monthly_var"
> - "calc_monthly_max"
> - "calc_monthly_min"

* [feat & fix: update diff](https://github.com/shenyulu/easyclimate/pull/108/commits/12f8e868bc342c440a47af1870954e182f240997)

> ⚠️ From this version, we change the functions name as follows:
>
> - "calc_lon_gradient" ➡️ "calc_dx_gradient"
> - "calc_lat_gradient" ➡️ "calc_dy_gradient"
>
> add functions:
>
> - "calc_dlon_radian_gradient"
> - "calc_dlon_degree_gradient"
> - "calc_dlat_radian_gradient"
> - "calc_dlat_degree_gradient"
> - "calc_shear_stretch_deform"
>
> and fix functions:
>
> - "calc_delta_pressure"
> - "calc_p_gradient"
> - "calc_p_integral"
> - "calc_top2surface_integral"
> - "calc_divergence"
> - "calc_vorticity"
> - "calc_geostrophic_wind"
> - "calc_geostrophic_wind_vorticity"
> - "calc_water_flux_top2surface_integral"
> - "calc_divergence_watervaporflux"
> - "calc_divergence_watervaporflux_top2surface_integral"
> - "calc_u_advection"
> - "calc_v_advection"
> - "calc_p_advection"

* [feat: add egr example](https://github.com/shenyulu/easyclimate/pull/108/commits/7d8dc200ec2f613d32d9d4afb4e0b446d0f941f2)
* [fix: fix waf calc](https://github.com/shenyulu/easyclimate/pull/108/commits/036f1ba797372790a085fea185c93f98364974ac)

> ⚠️ From this version, we fix the functions name as follows:
>
> - "calc_Plumb_wave_activity_horizontal_flux"
> - "calc_TN_wave_activity_horizontal_flux"

* [feat: add fortran macro](https://github.com/shenyulu/easyclimate/pull/108/commits/7b3a7d933c4e6cc67569b22df33862a775e4ebf2)
* [feat: add remove_sst_trend](https://github.com/shenyulu/easyclimate/pull/108/commits/d5c8b0e584560a3c79c40e13a66376753a656e2b)
* [feat: add add_lon_cyclic_lonarray](https://github.com/shenyulu/easyclimate/pull/108/commits/31b21e1ed46ed6d562be127851101b32646fef58)
* [fix: fix calc_seasonal_mean](https://github.com/shenyulu/easyclimate/pull/108/commits/0d3a5268e533dd478e9e17d68c268ae7196289cd)

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.11.0...v2025.12.0

## v2025.11.0

*Released on: 2025/10/23*

* [fix: EGR error](https://github.com/shenyulu/easyclimate/pull/105/commits/ab0979beac3dbbd36bd16174df23a869252f24da)

> ⚠️ From this version:
>
> - We addressed the issue in ``calc_eady_growth_rate`` where there was an unnecessary multiplication by the gravitational acceleration g. This extra part was confirmed to be incorrect🛠️, so there is now an order of magnitude difference between the current version and the previous version.

* [feat: add time_range for normalize](https://github.com/shenyulu/easyclimate/pull/105/commits/d27d962e8d78b7668caecdfc25fe797f2665d328)
* [feat: add calc_eof_projection_coefficient](https://github.com/shenyulu/easyclimate/pull/105/commits/61a62184d04a360ad1c5a71f1cfb40ece3cf9cb4)
* [feat: add func](https://github.com/shenyulu/easyclimate/pull/105/commits/934a46e0b9b7bd0cfae06550b77d7a3a40984fb1)

  - ``calc_windmask_ttestSpatialPattern_spatial``
  - ``calc_non_centered_corr``
  - ``calc_pattern_corr``

* [fix: fix codacy.instructions](https://github.com/shenyulu/easyclimate/pull/105/commits/d5dc7d06eda628da186bcbd802e90566343f5f03)

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.10.0...v2025.11.0

## v2025.10.0

*Released on: 2025/9/24*

* [fix: dp error](https://github.com/shenyulu/easyclimate/pull/103/commits/4068441335625c862239a8b1261f832e4e2b2afe)
* [feat: add mini*gallery](https://github.com/shenyulu/easyclimate/pull/103/commits/5ca9dab4e611370cdecf5d13b36955cbd6a5c061)
* [feat: add show_versions](https://github.com/shenyulu/easyclimate/pull/103/commits/468a5eeb64bb4879641cc546dad86b5307cca1b2)
* [feat: add windspeed](https://github.com/shenyulu/easyclimate/pull/103/commits/cc510f264e71ffefbac3e1f2167cf604c06d94bd)
* [fix: reg error](https://github.com/shenyulu/easyclimate/pull/103/commits/c37e67e50f79b213abccdb09eb2f05f60aa7ee04)
* [feat: add track and axisymmetric](https://github.com/shenyulu/easyclimate/pull/103/commits/023ea84117170d003e6301228c6315a85b2aaa35)
* [feat: add contributing](https://github.com/shenyulu/easyclimate/pull/103/commits/bf379dbb8aba9eb3fa905e91ebf9ee4c6bc42ca2)
* [fix: system info](https://github.com/shenyulu/easyclimate/pull/103/commits/3b45793bbd18d25a7d09ecbc51fe6b504039e3ef)
* [fix: para error in color](https://github.com/shenyulu/easyclimate/pull/103/commits/6f4f534500a4fca9d2ed969d9ec938b0076db34c)
* [feat: add linux install docs](https://github.com/shenyulu/easyclimate/pull/103/commits/592537ca8ec30246b902231641b9b1f734f7326e)
* [feat: change coastlines_linewidths to coastlines_kwargs](https://github.com/shenyulu/easyclimate/pull/103/commits/d06827f69bd5473b341067182795d8b4d1bd4195)

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.9.0...v2025.10.0

## v2025.9.0

*Released on: 2025/8/21*

* [docs: add uv](https://github.com/shenyulu/easyclimate/pull/102/commits/a7a801fbaa7ff3b47dee9dc219d6309a84f96cb2)
* [docs: add video](https://github.com/shenyulu/easyclimate/pull/102/commits/1ec9d296b0f8c56adb69560edd1b9300a0d354c0)
* [update: black version](https://github.com/shenyulu/easyclimate/pull/102/commits/6cd9019faa831f18b3d24cf19665911d44cdf9e1)
* [update: black exclude](https://github.com/shenyulu/easyclimate/pull/102/commits/ea9d3a235369ae2411ae3987b96f729ef2c7a9fa)
* [fix: build process error](https://github.com/shenyulu/easyclimate/pull/102/commits/264f6e92e88bc700f01dca48c30b8184178ed778)
* [fix: remove_low_frequency_signal](https://github.com/shenyulu/easyclimate/pull/102/commits/a526fe7e397cbf292081fe4c71302c850c1e1e86)
* [fix: fixed dims](https://github.com/shenyulu/easyclimate/pull/102/commits/f60b9a3b5d15deb998c420f26adfb0f38b3db680)
* [fix: explicit normalized](https://github.com/shenyulu/easyclimate/pull/102/commits/39127b776d6fd4e068896442b2a526468b583734)
* [feat: add normalized](https://github.com/shenyulu/easyclimate/pull/102/commits/ea13c89069944dcfe3aacad6fc51cf04f1a7f241)
* [fix: fix bar and line](https://github.com/shenyulu/easyclimate/pull/102/commits/29d6a57e12fde5a3e9a5b3e47e2e49dd2eff4b05)
* [fix: corr fix](https://github.com/shenyulu/easyclimate/pull/102/commits/9b1d6a7d980b4562821b42287997ba9745046f4d)
* [fix: fix pdo](https://github.com/shenyulu/easyclimate/pull/102/commits/3b403ff87e8fa3321cf9d1a753888bd9bb9d7be4)
* [text: fix file](https://github.com/shenyulu/easyclimate/pull/102/commits/a9e18f436ab065b336c1abf75b69a69c08a23c4f)
* [feat: add calc_timeseries_correlations](https://github.com/shenyulu/easyclimate/pull/102/commits/d50a012a88fda692e0d5e521700882e1727ad205)
* [feat: use rich instead](https://github.com/shenyulu/easyclimate/pull/102/commits/c5c6a029748aa8b4e769152a0322d42a41ebc9f6)
* [fix: fix draw_significant_area_contourf](https://github.com/shenyulu/easyclimate/pull/102/commits/523ee935bb77063cdbeaa5f5957814c14bb42587)
* [fix: fix cgt error](https://github.com/shenyulu/easyclimate/pull/102/commits/488c66acd3daa83af27e4218c4ca8b0c7544b829)
* [feat: add calc_leadlag_corr_spatial](https://github.com/shenyulu/easyclimate/pull/102/commits/3ebf91f2aeaec7300f47fe645dd1c4a1eef382f6)
* [update: update version](https://github.com/shenyulu/easyclimate/pull/102/commits/d4088cedb7b535a3d67148e85f741f453a26bdab)
* [fix: wk_spectra docs](https://github.com/shenyulu/easyclimate/pull/102/commits/d5f949de877f398139d058700ec53b9831758bd5)

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.8.0...v2025.9.0

### What's Changed

## v2025.8.0

*Released on: 2025/7/24*

### What's Changed


* [feat: add wet bulk](https://github.com/shenyulu/easyclimate/pull/100/commits/792dbe411933569c6493da9b4f957d3e0dd18e07)
* [feat: add calc_lat_weight_lin_rood](https://github.com/shenyulu/easyclimate/pull/100/commits/af185298fee706f08afdb9b244e4adcef1044d0c)
* [feat: add transfer_relative_humidity_2_dewpoint](https://github.com/shenyulu/easyclimate/pull/100/commits/736616c0a52e9e10182a8ebcc8bdce9a3267a814)
* [feat: add condensation and convection](https://github.com/shenyulu/easyclimate/pull/100/commits/7b003eec652f5fcf02aa2a164702182ea9d1881b)
* [fix: icon error win32](https://github.com/shenyulu/easyclimate/pull/100/commits/a1697b6bb52555d293ba9d57e6c43769687fe040)
* [feat: add vinth2p_dp](https://github.com/shenyulu/easyclimate/pull/100/commits/827604b9060d7783c8e8732240d2880390b33ae0)
* [feat: add enthalpy](https://github.com/shenyulu/easyclimate/pull/100/commits/24b6d9219659c17e3a26bccff19abb65b056a3ea)
* [fix: add mixing_ratio_data_units](https://github.com/shenyulu/easyclimate/pull/100/commits/e20857e0899ed3081cda780aa2371e0da96079b8)
* [fix: q unit error](https://github.com/shenyulu/easyclimate/pull/100/commits/435b7f7918a4bf8362c7852a51b0cc6b8ce255b1)
* [feat: add latent heat water](https://github.com/shenyulu/easyclimate/pull/100/commits/614f1a90ddd7a9971ab494c75a8883836000e1ed)
* [feat: add vinth2p](https://github.com/shenyulu/easyclimate/pull/100/commits/8f799791c8b323955175a91adae3dc993827fe10)
* [feat: add calc_enthalpy](https://github.com/shenyulu/easyclimate/pull/100/commits/004e899ce687e6137e2365f537a97aabed168c0b)
* [feat: add calc_relative_angular_momentum](https://github.com/shenyulu/easyclimate/pull/100/commits/688f58f4974ad0d4531b672ef009f9061b144ac0)
* [feat: add daily std var signal](https://github.com/shenyulu/easyclimate/pull/100/commits/c04565b895177c78a41de28362f97182953e2e3b)
* [feat: add pytest parallel support](https://github.com/shenyulu/easyclimate/pull/100/commits/eb8d99238766d57f15d2ed032b908acbd0388044)
* [feat: add vinth2p](https://github.com/shenyulu/easyclimate/pull/100/commits/ad21fe275feddc4b75194e28ce812b5b701732e5)
* [feat: add spectrum](https://github.com/shenyulu/easyclimate/pull/100/commits/ad61ecb28db27f2ef2fd5b706a6ab7f94868500e)
* [feat: center files in backend](https://github.com/shenyulu/easyclimate/pull/100/commits/d5ef0bec3311e15c8cdeacddd33b515d61cda1b9)
* [feat: add calc_lead_lag_correlation_coefficients](https://github.com/shenyulu/easyclimate/pull/100/commits/0dcc68fe43dd6507a62da4e4396be6ecfbc4ae40)
* [feat: add bsiso](https://github.com/shenyulu/easyclimate/pull/100/commits/06ad3e24fef8f06b47edbd7b517dcd835be2fa04)
* [fix: fix datanode plain text](https://github.com/shenyulu/easyclimate/pull/100/commits/30a8df61f060e76fad1101ed3838a1d3d9012250)
* [feat: migrate datatree to datanode](https://github.com/shenyulu/easyclimate/pull/100/commits/042c56bd86cb847f17b338c0f52dda48c18b6377)
* [feat: add da and bbo](https://github.com/shenyulu/easyclimate/pull/100/commits/0a386caa2f50ad42ddcc79d93707a59bbcdcb55b)
* [feat: add clean_extra_coords](https://github.com/shenyulu/easyclimate/pull/100/commits/29facafa653dd2afe65765a677d4cfd32e2af034)
* [feat: add interp_mesh2point_withtime](https://github.com/shenyulu/easyclimate/pull/100/commits/adb5ae532f7810903afa0351f88e6dd14fcb0e9b)


**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.7.0...v2025.8.0

## v2025.7.0

*Released on: 2025/6/21*

### What's Changed


* [feat: seperate physics part](https://github.com/shenyulu/easyclimate/pull/99/commits/f273ae7c7b52cb16171293b1a2462ec8c3df7879)
* [feat: add mesh2point](https://github.com/shenyulu/easyclimate/pull/99/commits/4b537b56fcf7d933f654e3351d93952263fc5fa0)
* [feat: add emd](https://github.com/shenyulu/easyclimate/pull/99/commits/be566d107bb3bce99ec3c9b2cc80b12a50b81617)
* [feat: add AO support](https://github.com/shenyulu/easyclimate/pull/99/commits/137c8fcd374ec3cb687dc32d198ebd2d3fef874c)
* [feat: add spharm support](https://github.com/shenyulu/easyclimate/pull/99/commits/2ed3dec189f11787a47bc4fabd0515532012b45f)
* [feat: add statellite support](https://github.com/shenyulu/easyclimate/pull/99/commits/f494de19e8cb823cf6e8b2cb06d7c85a13f539a5)
* [feat: add plot of coriolis coefficient](https://github.com/shenyulu/easyclimate/pull/99/commits/c2d3a09fd7f358615f2b75edf9c6c494013bdd23)
* [feat: add gaussian filter](https://github.com/shenyulu/easyclimate/pull/99/commits/f2471c7e5548b174585be373a3923b17df22f60e)

> ⚠️ From this version:
>
> - We removed `calc_potential_intensity_Bister_Emanuel_2002` code from `easyclimate.field.mesoscale.potential_intensity.py`, please use `easyclimate.field.typhoon.potential_intensity.py` instead 🛠️.
> - Seperate `easyclimate.core.diagnosis` to the part `easyclimate.physics`.

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.5.0...v2025.7.0

## v2025.5.0

*Released on: 2025/5/9*

### What's Changed


* [feat: add smooth_daily_annual_cycle gallery](https://github.com/shenyulu/easyclimate/pull/97/commits/e8fb01a110da712b9c5095af3413cf4c10377816)
* [add: add bar and line](https://github.com/shenyulu/easyclimate/pull/97/commits/21b80f8d06e1c99ace985e35e36177d5c27547f1)
  - `bar_plot_with_threshold`
  - `line_plot_with_threshold`
* [add: diagnosis func](https://github.com/shenyulu/easyclimate/pull/97/commits/478ba1ebe5ee77c509a50bad5be8f3c6ae11941c)
  - `calc_potential_temperature`
  - `calc_potential_temperature_vertical`
  - `calc_equivalent_potential_temperature`
  - `calc_equivalent_potential_temperature_davies_jones2009`
  - `calc_wet_bulb_potential_temperature_iteration`
  - `calc_wet_bulb_potential_temperature_davies_jones2008`
  - `calc_wet_bulb_temperature_stull2011`
  - `calc_wet_bulb_temperature_sadeghi2013`
  - `calc_lifting_condensation_level_bolton1980`
  - `calc_lifting_condensation_level_Bohren_Albrecht2023`
  - `calc_moist_adiabatic_lapse_rate`
  - `transfer_dewpoint_2_mixing_ratio`
* [fix: Solve Invalid value for attr solver_kwargs: dict](https://github.com/shenyulu/easyclimate/pull/97/commits/40f6fedd141958e563db711f4b9c9e452f91d472)
* [fix: units fix](https://github.com/shenyulu/easyclimate/pull/97/commits/dae2c447a28971d63bd2285da710671a0866ac23)
* [fix: change calc_potential_temperature_vertical from calc_potential_temperature](https://github.com/shenyulu/easyclimate/pull/97/commits/06eefc3a44d546153bff44b4263588c1497591aa)
* [feat: datanode enhance](https://github.com/shenyulu/easyclimate/pull/97/commits/9f94994352062f2e602fdc5aca305cc16e9f42f3)
* [feat: add datatime steps](https://github.com/shenyulu/easyclimate/pull/97/commits/82dfb5d73f7239877f63d0aaf0649491e92dd4b4)
  - `datetime_to_numeric`
  - `numeric_to_datetime`
  - `calculate_time_steps`
* [fix: fix time-invariant dependence](https://github.com/shenyulu/easyclimate/pull/97/commits/a8841bc9930a3adbc3afb8b1b23fd6230c6135a3)
* [docs: add example for tcpv](https://github.com/shenyulu/easyclimate/pull/97/commits/02471c1ce5be93c479afaafc63f07515fa6e0583)
* [feat: add emd func](https://github.com/shenyulu/easyclimate/pull/97/commits/b28219cd66c4f03dcb5b318c2c194525c222dfbb)
  - `filter_emd`
  - `filter_eemd`
* [fix: list errors](https://github.com/shenyulu/easyclimate/pull/97/commits/87002ab5575234a3f8cb0009c8f6153c5f779805)
* [feat: add open_datanode](https://github.com/shenyulu/easyclimate/pull/97/commits/ac6e274e101afcddc4af0c20878a96a1a3818abf)
  - `open_datanode`
* [feat: add meof support](https://github.com/shenyulu/easyclimate/pull/97/commits/c2400ae8d00878329fc5cf69f10909533fdf8e83)
  - `easyclimate.eof.get_EOF_model`
* [fix: path support linux](https://github.com/shenyulu/easyclimate/pull/97/commits/761fbf0cb208c3ef90c8fce88cfe281b7f0093e6)
* [fix: fix snowballstemmer version](https://github.com/shenyulu/easyclimate/pull/97/commits/288572e9d266cd363532582d320559cc15a8f726)
* [docs: use zarr 2](https://github.com/shenyulu/easyclimate/pull/97/commits/f42e596d2faff5a7eb4fccb889abdcdd60098177)


**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.4.0...v2025.5.0

## v2025.4.0

*Released on: 2025/4/12*

### What's Changed


* [feat: function enhancement](https://github.com/shenyulu/easyclimate/pull/92/commits/c689d06d5c55a2ff175d77525d93469448626e1e "feat: function enhancement")
  - `transfer_data_temperature_units`
  - `check_deprecation_status`
  - `smooth_daily_annual_cycle`
  - `calc_daily_annual_cycle_mean`
  - `remove_smooth_daily_annual_cycle_mean`
  - `draw_mjo_phase_space_basemap`
  - `draw_mjo_phase_space`
  - `remove_dominant_signals`
  - `decompose_symasym`
  - `calc_spectral_coefficients`
  - `draw_wk_anti_analysis`
  - `draw_wk_sym_analysis`
* [feat: add dependency](https://github.com/shenyulu/easyclimate/pull/92/commits/867f6ece202db26860bc7030ba0b69705e15b0ea "feat: add dependency")
* [feat: add check deprecations](https://github.com/shenyulu/easyclimate/pull/92/commits/b74f65e96a6dd49aaef8e9b5788994706ae66c4d "feat: add check deprecations")
* [feat: add calc_corr_spatial](https://github.com/shenyulu/easyclimate/pull/92/commits/6dd4c2df55ee18e057eaed5191d26fd8e58559f0 "feat: add calc_corr_spatial")
* [feat: add validate_dataarrays](https://github.com/shenyulu/easyclimate/pull/92/commits/fc81428f7688120bb17c57937add3a03c1005954 "feat: add validate_dataarrays")
* [fix: add default args](https://github.com/shenyulu/easyclimate/pull/92/commits/839e1e351997c669d9ec1a93a8075d7a0bf2942e "fix: add default args")
* [fix: auto repair lon error](https://github.com/shenyulu/easyclimate/pull/92/commits/96690699e7202062ab8e755f2602935b88eeafc0 "fix: auto repair lon error")
* [feat: add barnes_filter](https://github.com/shenyulu/easyclimate/pull/92/commits/39be22854cd4bf1ef7ba1616f123353609da7203 "feat: add barnes_filter")
* [feat: add gallery](https://github.com/shenyulu/easyclimate/pull/92/commits/edc00b01680160a5e38e1508b85fb9fcc397f4fc "feat: add gallery")
* [feat: add spatial pcf](https://github.com/shenyulu/easyclimate/pull/92/commits/ce74ea028b39a53bfb2de6282e7cacc72e9b3d76 "feat: add spatial pcf")
* [update: docs update](https://github.com/shenyulu/easyclimate/pull/92/commits/bebf4a2ec6f8703821475af9401cf269b77fda20 "update: docs update")

> ⚠️ From this version, we migrated `calc_potential_intensity_Bister_Emanuel_2002` code from `easyclimate.field.mesoscale.potential_intensity.py` (Removal Version `v2025.7.0`) to `easyclimate.field.typhoon.potential_intensity.py` 🛠️.

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.3.0...v2025.4.0

## v2025.3.0

*Released on: 2025/3/16*

### What's Changed

* [fix: fast-barnes-py error](https://github.com/shenyulu/easyclimate/pull/91/commits/cff9df9b6c20ae847cf7de5a55ef9418ca6c853e "fix: fast-barnes-py error")
* [feat: add wrf-python docs](https://github.com/shenyulu/easyclimate/pull/91/commits/fd53fcc6a75d0aa7487515c8a8ea13166ce2dfbe "feat: add wrf-python docs")
* [feat: add test_transfer_data_units](https://github.com/shenyulu/easyclimate/pull/91/commits/5946c39ebc71ea18346aea2317c0e0deb32b99b9 "feat: add test_transfer_data_units")
* [feat: add human index](https://github.com/shenyulu/easyclimate/pull/91/commits/68fc2187411fcf22651fcdf35c347a47a4568b7f "feat: add human index")
* [update: move backend](https://github.com/shenyulu/easyclimate/pull/91/commits/b3a96f328dd0f325f603df1c81a24d7f67d7a38b "update: move backend")
* [fix: redfit check](https://github.com/shenyulu/easyclimate/pull/91/commits/f1a662216e4846910221478f8fbb5eddbd845ed0 "fix: redfit check")
* [feat: add aerobulk](https://github.com/shenyulu/easyclimate/pull/91/commits/341817604291703991946dcff0d20461dc7520a6 "feat: add aerobulk")

> ⚠️ From this version, we uniformly migrated `fortran-related` code to the [easyclimate-backend](https://github.com/shenyulu/easyclimate-backend) repository 🛠️.

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.2.0...v2025.3.0

## v2025.2.0

*Released on: 2025/2/17*

### What's Changed

* [feat: add Lanczos filter](https://github.com/shenyulu/easyclimate/commit/54e407c7d6dfa1ec35aaf0458e87c69d70d98d1e)

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.1.0...v2025.2.0

## v2025.1.0

*Released on: 2025/1/9*

### What's Changed

* [feat: add curly vector support](https://github.com/shenyulu/easyclimate/commit/ca89e9be939fab698dc07a87943b7ad323220f35)
* [fix: fix Collection in matplotlib](https://github.com/shenyulu/easyclimate/commit/1cb5e719303d4bdb78d4856201e941c901b8f991)

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.12.0...v2025.1.0

## v2024.12.0

*Released on: 2024/12/12*

### What's Changed

* [feat: add interp1d_vertical_pressure2altitude](https://github.com/shenyulu/easyclimate/pull/86/commits/0137a7427dd9bd44746f11fdf03b28ee07d713b0)
* [feat: add docs seaborn rely](https://github.com/shenyulu/easyclimate/pull/86/commits/fca5fec9a860b32786ceffac80dfccd4cab97a2e)
* [fix: rtd github download](https://github.com/shenyulu/easyclimate/pull/86/commits/cef99c25a18cf017623926d4bec135c745fb8a1b)
* [feat: support Chinese and Japanese readme](https://github.com/shenyulu/easyclimate/pull/87/commits/0742ed89673d06c452e7e407a1401afb905a8697)
* [update: update index page](https://github.com/shenyulu/easyclimate/pull/87/commits/485563e916aa52430cf4077bcb0e7a251ae6f902)
* [fix: fix fast-barnes-py version temperately](https://github.com/shenyulu/easyclimate/pull/87/commits/277e687c652e53e59f1cc8699d074de194caf231)

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.11.0...v2024.12.0

## v2024.11.0

*Released on: 2024/11/17*

### What's Changed

* [fix: modify `specific_humidity_units` to `specific_humidity_data_units`](https://github.com/shenyulu/easyclimate/pull/85/commits/67ed21523b17105c7e353dae85d9ba8b7e5c36e9)
* [fix: rename file in the `field`](https://github.com/shenyulu/easyclimate/pull/85/commits/e7defd7e398e5836bde1a015ed52c9c2e703ff42)
* [fix: fix error DataTree immigration](https://github.com/shenyulu/easyclimate/pull/85/commits/a4fa73d2bbb333fac5a157857429d6e17cae88d1)
* [feat: add Potential intensity of TC](https://github.com/shenyulu/easyclimate/pull/85/commits/e7defd7e398e5836bde1a015ed52c9c2e703ff42)

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.10.0...v2024.11.0

## v2024.10.0

*Released on: 2024/10/09*

### What's Changed

* [feat: add quick_draw_rectangular_box](https://github.com/shenyulu/easyclimate/pull/77/commits/727ad88dce8902d72bee249f7c485d3860108a75)
* [feat: Modify the units conversion module](https://github.com/shenyulu/easyclimate/pull/77/commits/18740c3ffd78dcf55b1503b1f9284993dad47535)
* [feat: Enhance unit availability](https://github.com/shenyulu/easyclimate/pull/77/commits/470e4451dcb9c67b5d1af8e06de4ed6c1968721f)
* [feat: add water vapor parameter calculation function](https://github.com/shenyulu/easyclimate/pull/77/commits/5692972b44fcc712da6512976ca43da84739308b)
* [fix: fix error introduced by xeofs 3.0.0](https://github.com/shenyulu/easyclimate/pull/77/commits/930444e720dae3a000bd262fb67e59cc91ce8380)
* [fix: Type expression](https://github.com/shenyulu/easyclimate/pull/77/commits/125a161232a3879e6aab13046540685e61dbd29c)
* [update: fix eof error for xeof&gt;3.0.0](https://github.com/shenyulu/easyclimate/pull/77/commits/e6a32dabf54fdf2b7f0df5d27d37c466c34386ea)
* [test: add test for water vapor calculation](https://github.com/shenyulu/easyclimate/pull/77/commits/7314b5f38b3c605d38bb15626e496e521727ac09)
* [fix: accuracy problem on linux](https://github.com/shenyulu/easyclimate/pull/77/commits/fa042c132701c9eb5c707d57514cf29bbc024490)
* [depend: update dependency](https://github.com/shenyulu/easyclimate/pull/77/commits/d1486b0d298a29328bfbfe1b748d3122ed002975)
* [update: update version of pre*commit*config](https://github.com/shenyulu/easyclimate/pull/77/commits/ca25504503ae70a0cfd989e00f7e7cef96d96d17)
* [fix: test dependency for netcdf](https://github.com/shenyulu/easyclimate/pull/77/commits/a4f56b283ed3086f9fd54d13c61651d01bddccd7)

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.08.01...v2024.10.0

## v2024.08.01

*Released on: 2024/08/02*

### What's Changed

* Update requirements text by @shenyulu in https://github.com/shenyulu/easyclimate/pull/51
* Update sphinx from 6.2.1 to 7.3.7 in the requirements.txt by @shenyulu in https://github.com/shenyulu/easyclimate/pull/52
* update sphinx version from 6.2.1 to 7.3.7 by @shenyulu in https://github.com/shenyulu/easyclimate/pull/53
* build(deps): bump sphinx-autoapi from 3.1.1 to 3.1.2 by @dependabot in https://github.com/shenyulu/easyclimate/pull/54
* chore: necessary update (2024.08.01) by @shenyulu in https://github.com/shenyulu/easyclimate/pull/66

  * [fix: update `calc_NPWI_monsoon_retreat` from `calc_NPWI_monsoon_detreat`](https://github.com/shenyulu/easyclimate/pull/66/commits/26d041ce5986104bdd7f0a6c30d195c22b9c8833)
  * [fix: change name `calc_detrend_spatial` from `calc_detrend_data`](https://github.com/shenyulu/easyclimate/pull/66/commits/dcb54394191faf7ae91765f941ab9a157d81ffb2)
  * [test: update calc_NPWI_monsoon_retreat from `calc_NPWI_monsoon_detreat](https://github.com/shenyulu/easyclimate/pull/66/commits/deb0c8c799c758c61ac3cfa741fa3d34c68ca63d)
  * [feat: add quick_draw_spatial_basemap for quick drawing](https://github.com/shenyulu/easyclimate/pull/66/commits/370a26f2660610ee1312e376c25c7f13ceb0e921)
  * [depend: Update dependency requirement](https://github.com/shenyulu/easyclimate/pull/66/commits/39052c8a623ec94b760763c5f6203079ec1e90be)
  * [test: add test for quick_draw_spatial_basemap](https://github.com/shenyulu/easyclimate/pull/66/commits/db3e0738153c5a14704d51bab83c9a38057a7c10)
  * [fix: File reference modification](https://github.com/shenyulu/easyclimate/pull/66/commits/d8067e94b0fe98c2ae001148f268952f14c0d41a)

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.05.01...v2024.08.01

## v2024.05.01

*Released on: 2024/05/19*

### What's Changed

* Add the calculation of virtual temperature and interpolating from the model levels to pressure levels by @shenyulu in https://github.com/shenyulu/easyclimate/pull/27
* [pre-commit.ci] pre-commit autoupdate by @pre-commit-ci in https://github.com/shenyulu/easyclimate/pull/29
* Bump pytest from 8.0.2 to 8.1.1 by @dependabot in https://github.com/shenyulu/easyclimate/pull/28
* [pre-commit.ci] pre-commit autoupdate by @pre-commit-ci in https://github.com/shenyulu/easyclimate/pull/30
* update testdata for the change of pyspharm-syl by @shenyulu in https://github.com/shenyulu/easyclimate/pull/33
* Remove the pin of xarray by @shenyulu in https://github.com/shenyulu/easyclimate/pull/34
* chore: necessary update (2024.04.01) by @shenyulu in https://github.com/shenyulu/easyclimate/pull/36

  * [fix: calculation of ENSO index](https://github.com/shenyulu/easyclimate/pull/36/commits/ed57c3941a926229166b634bec467162096a58bb)
  * [fix: dropping variables using drop is deprecated, use drop_vars](https://github.com/shenyulu/easyclimate/pull/36/commits/f8dc30f9fe46a1bccc0cdb51c525fc2018d2901b)
  * [refactor, test: rename populate_monmean2everymon and `populate_daym…](https://github.com/shenyulu/easyclimate/pull/36/commits/a82c636611690d571683554d93897991a5512332)
  * [docs: add more details for calc_daily_climatological_anomaly](https://github.com/shenyulu/easyclimate/pull/36/commits/0ef7ca145c74222e08fbd591181a26e0a00ab3d8)
  * [feat: supports the removal of seasonal cycles from the climate averag…](https://github.com/shenyulu/easyclimate/pull/36/commits/4930f0d40f5856fb4dc8582e64901145cd45336a)
* chore: necessary update (2024.05.01) by @shenyulu in https://github.com/shenyulu/easyclimate/pull/43

  * [feat: The calculation of the PNA index](https://github.com/shenyulu/easyclimate/pull/43/commits/6570c97e62d32567f1d19f85637dbd7ca79150a0)
  * [feat: The calculation of the NAO index](https://github.com/shenyulu/easyclimate/pull/43/commits/da59768f9ce8c113f718230f46f8aa1360b5e664)
  * [refactor: rename field.atm to field.teleconnection](https://github.com/shenyulu/easyclimate/pull/43/commits/2d4019aa1643710c778985e04e8509473dab2816)
  * [feat: add teleconnection index](https://github.com/shenyulu/easyclimate/pull/43/commits/890af93ecff65c78aa9c48334a48004dec1cef4b)
  * [feat: add calc_seasonal_mean](https://github.com/shenyulu/easyclimate/pull/43/commits/085e966bf517df60378ade3c24c9615f5df42c2e)
  * [feat: add IOD index](https://github.com/shenyulu/easyclimate/pull/43/commits/08d278a5b27ca074a06d5213218f086f60b6b63c)
  * [refactor, test: use SST data to calculate ENSO index](https://github.com/shenyulu/easyclimate/pull/43/commits/777c2492800d3ef596228ff977c1264bc4935cf6)
  * [fix: Change the file name to lowercase](https://github.com/shenyulu/easyclimate/pull/43/commits/03ab712adc190af098d1efbb0813bb2071efb22e)
  * [feat: add IOBM index](https://github.com/shenyulu/easyclimate/pull/43/commits/a938a74e91f5a92076b368bcb2c7f2ed246effeb)
  * [fix: change module name `land_atm_interaction` to `land`](https://github.com/shenyulu/easyclimate/pull/43/commits/30903458cf71af2adf66dbe1593372fc1af519b0)

### New Contributors

* @pre-commit-ci made their first contribution in https://github.com/shenyulu/easyclimate/pull/29

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.03.01...v2024.05.01

## v2024.03.01

*Released on: 2024/03/05*

### What's Changed

* pin xarray version <=2023.12.0 by @shenyulu in https://github.com/shenyulu/easyclimate/pull/13
* modify eof functions and add test and add create .pre-commit-config.yaml by @shenyulu in https://github.com/shenyulu/easyclimate/pull/16
* update requirement by @shenyulu in https://github.com/shenyulu/easyclimate/pull/23
* remove the pin of geocat.viz by @shenyulu in https://github.com/shenyulu/easyclimate/pull/24
* change binder example location by @shenyulu in https://github.com/shenyulu/easyclimate/pull/25

## v2024.01.01

*Released on: 2024/01/04*

### What's Changed

* Publish version 2024.01.01 by @shenyulu in https://github.com/shenyulu/easyclimate/pull/9
* Bump pytest from 7.4.3 to 7.4.4 by @dependabot in https://github.com/shenyulu/easyclimate/pull/8

## v2023.12.01

*Released on: 2023/12/10*

### What's Changed

* Bump sphinx from 6.2.1 to 7.2.6 by @dependabot in https://github.com/shenyulu/easyclimate/pull/2
* Publish version 2023.12.01 by @shenyulu in https://github.com/shenyulu/easyclimate/pull/3
* Add pytest (codecov to 44.37%)
* Improve requirements

### New Contributors

* @shenyulu made their first contribution in https://github.com/shenyulu/easyclimate/pull/3

## v0.0.1

*Released on: 2022/05/02*

### What's Changed

* First release

### New Contributors

* shenyulu

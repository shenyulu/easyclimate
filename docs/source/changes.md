# Changelog

## Version 2024.10.0

*Released on: 2024/10/09*

### What's Changed

* [feat: add quick_draw_rectangular_box](https://github.com/shenyulu/easyclimate/pull/77/commits/727ad88dce8902d72bee249f7c485d3860108a75)
* [feat: Modify the units conversion module](https://github.com/shenyulu/easyclimate/pull/77/commits/18740c3ffd78dcf55b1503b1f9284993dad47535)
* [feat: Enhance unit availability](https://github.com/shenyulu/easyclimate/pull/77/commits/470e4451dcb9c67b5d1af8e06de4ed6c1968721f)
* [feat: add water vapor parameter calculation function](https://github.com/shenyulu/easyclimate/pull/77/commits/5692972b44fcc712da6512976ca43da84739308b)
* [fix: fix error introduced by xeofs 3.0.0](https://github.com/shenyulu/easyclimate/pull/77/commits/930444e720dae3a000bd262fb67e59cc91ce8380)
* [fix: Type expression](https://github.com/shenyulu/easyclimate/pull/77/commits/125a161232a3879e6aab13046540685e61dbd29c)
* [update: fix eof error for xeof>3.0.0](https://github.com/shenyulu/easyclimate/pull/77/commits/e6a32dabf54fdf2b7f0df5d27d37c466c34386ea)
* [test: add test for water vapor calculation](https://github.com/shenyulu/easyclimate/pull/77/commits/7314b5f38b3c605d38bb15626e496e521727ac09)
* [fix: accuracy problem on linux](https://github.com/shenyulu/easyclimate/pull/77/commits/fa042c132701c9eb5c707d57514cf29bbc024490)
* [depend: update dependency](https://github.com/shenyulu/easyclimate/pull/77/commits/d1486b0d298a29328bfbfe1b748d3122ed002975)
* [update: update version of pre*commit*config](https://github.com/shenyulu/easyclimate/pull/77/commits/ca25504503ae70a0cfd989e00f7e7cef96d96d17)
* [fix: test dependency for netcdf](https://github.com/shenyulu/easyclimate/pull/77/commits/a4f56b283ed3086f9fd54d13c61651d01bddccd7)

**Full Changelog**: https://github.com/shenyulu/easyclimate/compare/v2024.08.01...v2024.10.0

## Version 2024.08.01

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

## Version 2024.05.01

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

## Version 2024.03.01

*Released on: 2024/03/05*

### What's Changed
* pin xarray version <=2023.12.0 by @shenyulu in https://github.com/shenyulu/easyclimate/pull/13
* modify eof functions and add test and add create .pre-commit-config.yaml by @shenyulu in https://github.com/shenyulu/easyclimate/pull/16
* update requirement by @shenyulu in https://github.com/shenyulu/easyclimate/pull/23
* remove the pin of geocat.viz by @shenyulu in https://github.com/shenyulu/easyclimate/pull/24
* change binder example location by @shenyulu in https://github.com/shenyulu/easyclimate/pull/25

## Version 2024.01.01

*Released on: 2024/01/04*

### What's Changed
* Publish version 2024.01.01 by @shenyulu in https://github.com/shenyulu/easyclimate/pull/9
* Bump pytest from 7.4.3 to 7.4.4 by @dependabot in https://github.com/shenyulu/easyclimate/pull/8

## Version 2023.12.01

*Released on: 2023/12/10*

### What's Changed
* Bump sphinx from 6.2.1 to 7.2.6 by @dependabot in https://github.com/shenyulu/easyclimate/pull/2
* Publish version 2023.12.01 by @shenyulu in https://github.com/shenyulu/easyclimate/pull/3
* Add pytest (codecov to 44.37%)
* Improve requirements

### New Contributors
* @shenyulu made their first contribution in https://github.com/shenyulu/easyclimate/pull/3

## Version 0.0.1

*Released on: 2022/05/02*

### What's Changed
* First release

### New Contributors
* shenyulu

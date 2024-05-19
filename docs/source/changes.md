# Changelog

## Version 2024.05.01

*Released on: 2024/05/19*

## What's Changed
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


## New Contributors
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

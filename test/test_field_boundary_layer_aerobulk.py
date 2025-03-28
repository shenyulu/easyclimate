"""
pytest for field.detection.aerobulk.py
"""

import pytest
import easyclimate as ecl
import xarray as xr
import numpy as np


def test_calc_turbulent_fluxes_without_skin_correction():
    tmp = ecl.field.boundary_layer.calc_turbulent_fluxes_without_skin_correction(
        sst_data=xr.DataArray(np.array([273.15 + 22.0, 273.15 + 22.0])),
        sst_data_units="degC",
        absolute_temperature_data=xr.DataArray(
            np.array([273.15 + 20.0, 273.15 + 20.0])
        ),
        absolute_temperature_data_units="degC",
        specific_humidity_data=xr.DataArray(np.array([0.012, 0.012])),
        specific_humidity_data_units="g/kg",
        zonal_wind_speed_data=xr.DataArray(np.array([4, 4])),
        zonal_wind_speed_data_units="m/s",
        meridional_wind_speed_data=xr.DataArray(np.array([9, 9])),
        meridional_wind_speed_data_units="m/s",
        mean_sea_level_pressure_data=xr.DataArray(np.array([101000.0, 101000.0])),
        mean_sea_level_pressure_data_units="Pa",
        algorithm="ncar",
    )
    result_data = tmp["ql"].data[0]
    refer_data = np.float64([4082.766776880013])
    assert np.isclose(result_data, refer_data).all()


def test_calc_turbulent_fluxes_skin_correction():
    tmp = ecl.field.boundary_layer.calc_turbulent_fluxes_skin_correction(
        sst_data=xr.DataArray(
            np.array([273.15 + 22.0, 273.15 + 22.0]),
        ),
        sst_data_units="degC",
        absolute_temperature_data=xr.DataArray(
            np.array([273.15 + 20.0, 273.15 + 20.0])
        ),
        absolute_temperature_data_units="degC",
        specific_humidity_data=xr.DataArray(np.array([0.012, 0.012])),
        specific_humidity_data_units="g/kg",
        zonal_wind_speed_data=xr.DataArray(np.array([4, 4])),
        zonal_wind_speed_data_units="m/s",
        meridional_wind_speed_data=xr.DataArray(np.array([9, 9])),
        meridional_wind_speed_data_units="m/s",
        mean_sea_level_pressure_data=xr.DataArray(np.array([101000.0, 101000.0])),
        mean_sea_level_pressure_data_units="Pa",
        downwelling_shortwave_radiation=xr.DataArray(np.array([0, 0])),
        downwelling_shortwave_radiation_units="W/m^2",
        downwelling_longwave_radiation=xr.DataArray(np.array([350.0, 350.0])),
        downwelling_longwave_radiation_units="W/m^2",
        algorithm="coare3p6",
    )
    result_data = tmp["t_s"].data[0]
    refer_data = np.float64([566.6892344670376])
    assert np.isclose(result_data, refer_data).all()

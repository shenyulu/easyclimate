"""
Potential intensity of TC
"""

import xarray as xr
from typing import Literal
from ..typhoon.potential_intensity import (
    calc_potential_intensity_Bister_Emanuel_2002 as newfunc,
)
from ...core.utility import deprecated

__all__ = [
    "calc_potential_intensity_Bister_Emanuel_2002",
]


# 使用示例
@deprecated(
    version="2025.4.0",
    removal_version="2025.7.0",
    replacement="easyclimate.field.typhoon.calc_potential_intensity_Bister_Emanuel_2002",
)
def calc_potential_intensity_Bister_Emanuel_2002(
    sst_data: xr.DataArray,
    sst_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    surface_pressure_data: xr.DataArray,
    surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    temperature_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    specific_humidity_data: xr.DataArray,
    specific_humidity_data_units: str,
    vertical_dim: str,
    vertical_dim_units: str,
    CKCD: float = 0.9,
    ascent_flag: bool = False,
    diss_flag: bool = True,
    V_reduc: float = 0.8,
    ptop: float = 50,
    miss_handle: bool = True,
) -> xr.Dataset:
    """
    Calculate potential intensity of TC (tropical cyclone) according to the Bister and Emanuel (2002) algorithm.

    .. warning::

        Please use `easyclimate.field.typhoon.calc_potential_intensity_Bister_Emanuel_2002` instead.

    This function calculates the maximum wind speed and mimimum central pressure achievable
    in tropical cyclones, given a sounding and a sea surface temperature.

    From Bister and Emanuel (1998) EQN. 21, PI may be computed directly via:

    .. math::

        V_{max}^{2} = \\frac{C_k}{C_D}(\\frac{T_{s} - T_{0}}{T_{0}})(h_0^* - h^*),

    where :math:`C_k` and :math:`C_D` are the enthalpy and momentum surface exchange coefficients,
    respectively; :math:`T_{s}` is the sea surface temperature; :math:`T_{0}` is the mean outflow temperature;
    :math:`h_0^*` is the saturation moist static energy at the sea surface;
    and :math:`h^*` is the moist static energy of the free troposphere.
    The ratio :math:`\\frac{C_k}{C_D}` is an uncertain quantity typically taken
    to be a constant (default is 0.9, see Emanuel 2003 and references therein).

    Building on this definition, one can extract TC efficiency
    and disequilibrium, and decompose the terms to determine their relative contributions to potential intensity.

    The efficiency of TC PI is the Carnot efficiency. Typical values range between 50-70% in the tropics.

    Each term in the PI equation may decomposed by taking the natural logarithm of both sides, arriving at (Wing et al. 2015; EQN. 2):

    .. math::
        2*\\log(V_{max}) = \\log(\\frac{C_k}{C_D}) + \\log(\\frac{T_{s} - T_{0}}{T_{0}}) + \\log(h_0^* - h^*).

    Note that the units of everything input to the functions (and particularly the temperatures) must match.

    Parameters
    ----------
    sst_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The sea surface temperature data.
    sst_data_units: :py:class:`str <str>`.
        The unit corresponding to `sst_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    surface_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Mean surface sea level pressure.
    surface_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `surface_pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The Specific humidity of air.
    specific_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    CKCD: :py:class:`float <float>`, default 0.9.
        Ratio of :math:`C_k` to :math:`C_D` (unitless number), i.e. the ratio of the exchange coefficients of enthalpy and
        momentum flux (e.g. see Bister and Emanuel 1998, EQN. 17-18). More discussion on :math:`\\frac{C_k}{C_D}` is found in Emanuel (2003).
        Default is 0.9 based on e.g. Wing et al. (2015).
    ascent_flag: :py:class:`bool <bool>`, default False.
        Adjustable constant fraction (unitless fraction) for buoyancy of displaced parcels,
        where `True` is Reversible ascent (default) and `False` is Pseudo-adiabatic ascent.
    V_reduc: :py:class:`float <float>`, default 0.8.
        Adjustable constant fraction (unitless fraction) for reduction of gradient winds to 10-m winds
        see Emanuel (2000) and Powell (1980).
    ptop: :py:class:`float <float>`, default 50 **hPa**.
        Pressure below which sounding is ignored (**hPa**).
    miss_handle: :py:class:`bool <bool>`, default True.
        Flag that determines how missing (NaN) values are handled in CAPE calculation.
        - If `False` (BE02 default), NaN values in profile are ignored and PI is still calcuated.
        - If `True`, given NaN values PI will be set to missing (with `IFLAG=3` in CAPE calculation).

        .. note::
            If any missing values are between the lowest valid level and ptop
            then PI will automatically be set to missing (with `IFLAG=3` in CAPE calculation)

    Returns
    -------
    - vmax: The maximum surface wind speed (m/s) reduced to reflect surface drag via :math:`V_{\\text{reduc}}`.
    - pmin: The minimum central pressure (hPa)
    - ifl: A flag value: A value of 1 means OK; a value of 0 indicates no convergence; a value of 2
      means that the CAPE routine failed to converge; a value of 3  means the CAPE routine failed due to
      missing data in the inputs.
    - t0: The outflow temperature (K)
    - otl: The outflow temperature level (hPa), defined as the level of neutral bouyancy
      where the outflow temperature is found, i.e. where buoyancy is actually equal
      to zero under the condition of an air parcel that is saturated at sea level pressure.
    - eff: Tropical cyclone efficiency.
    - diseq: Thermodynamic disequilibrium.
    - lnpi: Natural :math:`\\log(\\text{Potential Intensity})`
    - lneff: Natural :math:`\\log(\\text{Tropical Cyclone Efficiency})`
    - lndiseq: Natural :math:`\\log(\\text{Thermodynamic Disequilibrium})`
    - lnCKCD: Natural :math:`\\log(C_k/C_D)`

    Reference
    --------------
    - https://github.com/dgilford/tcpyPI

    - Bister, M., Emanuel, K.A. Dissipative heating and hurricane intensity. Meteorl. Atmos. Phys. 65, 233–240 (1998). https://doi.org/10.1007/BF01030791
    - Bister, M., and K. A. Emanuel, Low frequency variability of tropical cyclone potential intensity, 1, Interannual to interdecadal variability, J. Geophys. Res., 107(D24), 4801, https://doi.org/10.1029/2001JD000776, 2002.
    - Emanuel, K.: A Statistical Analysis of Tropical Cyclone Intensity, Mon. Weather Rev., 128, 1139–1152, https://doi.org/10.1175/1520-0493(2000)128<1139:ASAOTC>2.0.CO;2, 2000.
    - Emanuel, K.: Tropical Cyclones, Annu. Rev. Earth Pl. Sc., 31, 75–104, https://doi.org/10.1146/annurev.earth.31.100901.141259, 2003.
    - Gilford, D. M.: pyPI (v1.3): Tropical Cyclone Potential Intensity Calculations in Python, Geosci. Model Dev., 14, 2351–2369, https://doi.org/10.5194/gmd-14-2351-2021, 2021.
    - Powell, M. D.: Evaluations of Diagnostic Marine Boundary-Layer Models Applied to Hurricanes, Mon. Weather Rev., 108, 757–766, https://doi.org/10.1175/1520-0493(1980)108<0757:EODMBL>2.0.CO;2, 1980.
    - Wing, A. A., Emanuel, K., and Solomon, S.: On the factors affecting trends and variability in tropical cyclone potential intensity, Geophys. Res. Lett., 42, 8669–8677, https://doi.org/10.1002/2015GL066145, 2015.
    """
    result = newfunc(
        sst_data,
        sst_data_units,
        surface_pressure_data,
        surface_pressure_data_units,
        temperature_data,
        temperature_data_units,
        specific_humidity_data,
        specific_humidity_data_units,
        vertical_dim,
        vertical_dim_units,
        CKCD,
        ascent_flag,
        diss_flag,
        V_reduc,
        ptop,
        miss_handle,
    )
    return result

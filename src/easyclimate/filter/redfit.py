"""
Red-noise spectra estimating
"""

import numpy as np
import xarray as xr
import pandas as pd
import os
from ctypes import c_char
from typing import Literal

from ..backend import _ecl_redfit, _ecl_redfit_x

__all__ = ["calc_redfit", "calc_redfit_cross"]


def cfg_generate(
    nmlfile_name,
    fnin="./red_15.dat",
    fnout="test.red",
    nsim=1000,
    mctest=False,
    rhopre=-99.0,
    ofac=1.0,
    hifac=1.0,
    n50=1,
    iwin="rectangular",
):
    """
    The cfg file generate generator.

    Parameter:

    fnin:   Input filename with time series data
    fnout:  Results are written to this file (ASCII format)
    nsim:   Number of Monte-Carlo simulations (1000–2000 should be o.k. in most cases)
    mctest: Toggle calculation of false-alarm levels based on Monte-Carlo simulation,
            if set to T : perform Monte-Carlo test,
            if set to F : skip Monte-Carlo test (default).
    ofac:   Oversampling factor for Lomb-Scargle Fourier transform (typical values: 2.0– 4.0)
    hifac:  Max. frequency to analyze is set to hifac * <fNyq> (default = 1.0)
    n50:    Number of WOSA segments (with 50 % overlap)
    rhopre: Prescibed value for ρ; unused if < 0 (default = -99.0)
    iwin:   Window-type identifier used to suppress sidelobes in spectral analysis:
            (0: Rectangular, 1: Welch, 2: Hanning, 3: Triangular, 4: Blackman-Harris).

    Parameters ofac, hifac, n50 and window type are identical to the SPECTRUM program
    (see Schulz and Stattegger, 1997 for further details).
    Except mctest, hifac and rhopre all parameters must be specified.
    """

    # Transform boolean type to string
    if mctest == True:
        mctest_str = "T"
    elif mctest == False:
        mctest_str = "F"

    if iwin == "rectangular":
        iwin_integer = 0
    elif iwin == "welch":
        iwin_integer = 1
    elif iwin == "hanning":
        iwin_integer = 2
    elif iwin == "triangular":
        iwin_integer = 3
    elif iwin == "blackmanharris":
        iwin_integer = 4
    else:
        raise ValueError(
            'Parameter iwin should be "rectangular", "welch", "hanning", "triangular", or "blackmanharris"'
        )

    with open(nmlfile_name, mode="w", encoding="utf-8") as f:
        f.writelines(
            [
                "&cfg\n"
                "    fnin = '" + fnin + "',\n"
                "   fnout = '" + fnout + "',\n"
                "    nsim = " + str(nsim) + ",\n"
                "  mctest = " + str(mctest_str) + ",\n"
                "  rhopre = " + str(rhopre) + ",\n"
                "    ofac = " + str(ofac) + ",\n"
                "   hifac = " + str(hifac) + ",\n"
                "     n50 = " + str(n50) + ",\n"
                "    iwin = " + str(iwin_integer) + "\n"
                "/\n"
                "\n"
                "\n"
                "! iwin = 0: Rectangular\n"
                "!        1: Welch\n"
                "!        2: Hanning\n"
                "!        3: Triangular\n"
                "!        4: Blackman-Harris"
            ]
        )


def cfg_generate_cross(
    nmlfile_name,
    fnin_1="./x.dat",
    fnin_2="./y.dat",
    fnout="REDFIT-X-result",
    x_sign=False,
    y_sign=False,
    nsim=1000,
    mctest=True,
    mctest_phi=True,
    rhopre_1=-999.0,
    rhopre_2=-999.0,
    ofac=4.0,
    hifac=1.0,
    n50=8,
    alpha=0.05,
    iwin="rectangular",
):
    """
    The cfg file cross redfit generate generator.

    Parameter:

    fnin:   Input filename with time series data
    fnout:  Results are written to this file (ASCII format)
    nsim:   Number of Monte-Carlo simulations (1000–2000 should be o.k. in most cases)
    mctest: Toggle calculation of false-alarm levels based on Monte-Carlo simulation,
            if set to T : perform Monte-Carlo test,
            if set to F : skip Monte-Carlo test (default).
    ofac:   Oversampling factor for Lomb-Scargle Fourier transform (typical values: 2.0– 4.0)
    hifac:  Max. frequency to analyze is set to hifac * <fNyq> (default = 1.0)
    n50:    Number of WOSA segments (with 50 % overlap)
    rhopre: Prescibed value for ρ; unused if < 0 (default = -99.0)
    iwin:   Window-type identifier used to suppress sidelobes in spectral analysis:
            (0: Rectangular, 1: Welch, 2: Hanning, 3: Triangular, 4: Blackman-Harris).

    Parameters ofac, hifac, n50 and window type are identical to the SPECTRUM program
    (see Schulz and Stattegger, 1997 for further details).
    Except mctest, hifac and rhopre all parameters must be specified.
    """

    # Transform boolean type to string
    if mctest == True:
        mctest_str = "T"
    elif mctest == False:
        mctest_str = "F"

    if x_sign == True:
        x_sign_str = "T"
    elif x_sign == False:
        x_sign_str = "F"
    if y_sign == True:
        y_sign_str = "T"
    elif y_sign == False:
        y_sign_str = "F"

    if mctest_phi == True:
        mctest_phi_str = "T"
    elif mctest_phi == False:
        mctest_phi_str = "F"

    if iwin == "rectangular":
        iwin_integer = 0
    elif iwin == "welch":
        iwin_integer = 1
    elif iwin == "hanning":
        iwin_integer = 2
    elif iwin == "triangular":
        iwin_integer = 3
    elif iwin == "blackmanharris":
        iwin_integer = 4
    else:
        raise ValueError(
            'Parameter iwin should be "rectangular", "welch", "hanning", "triangular", or "blackmanharris"'
        )

    with open(nmlfile_name, mode="w", encoding="utf-8") as f:
        f.writelines(
            [
                "&cfg\n"
                "    fnin(1) = '" + str(fnin_1) + "',\n"
                "    fnin(2) = '" + str(fnin_2) + "',\n"
                "   fnout = '" + str(fnout) + "',\n"
                "    x_sign = " + str(x_sign_str) + ",\n"
                "    y_sign = " + str(y_sign_str) + ",\n"
                "    nsim = " + str(nsim) + ",\n"
                "  mctest = " + str(mctest_str) + ",\n"
                "  mctest_phi = " + str(mctest_phi_str) + ",\n"
                "  rhopre(1) = " + str(rhopre_1) + ",\n"
                "  rhopre(2) = " + str(rhopre_2) + ",\n"
                "    ofac = " + str(ofac) + ",\n"
                "   hifac = " + str(hifac) + ",\n"
                "     n50 = " + str(n50) + ",\n"
                "     alpha = " + str(alpha) + ",\n"
                "    iwin = " + str(iwin_integer) + "\n"
                "/\n"
            ]
        )


def gen_random_string(length):
    """
    A string that generates a random password, containing only letters and numbers.
    You can specify the number of bits in a string.
    From: https://blog.csdn.net/xiong_xin/article/details/122461373

    E.g.,
    gen_random_string(15)

    ---------------------
    '7W492927406483w'
    """
    import random
    import string

    # 随机生成字母和数字的位数
    # The number of randomly generated letters and numbers
    numcount = random.randint(1, length - 1)
    lettercount = length - numcount

    # 随机抽样生成数字序列
    # Random sampling generates a sequence of numbers
    numlist = [random.choice(string.digits) for _ in range(numcount)]

    # 随机抽样生成字母序列
    # Random sampling generates a sequence of letters
    letterlist = [random.choice(string.ascii_letters) for _ in range(lettercount)]

    # 合并字母数字序列
    # Merge alphanumeric sequences
    alllist = numlist + letterlist

    # 乱序
    # out-of-order
    result = random.shuffle(alllist)

    # 生成目标结果字符串
    # Generate the target result string
    result = "".join([i for i in alllist])

    return result


def calc_redfit(
    data: xr.DataArray,
    timearray: np.array = None,
    nsim: int = 1000,
    mctest: bool = False,
    rhopre: float = -99.0,
    ofac: float = 1.0,
    hifac: float = 1.0,
    n50: int = 1,
    iwin: Literal[
        "rectangular", "welch", "hanning", "triangular", "blackmanharris"
    ] = "rectangular",
):
    """
    Estimating red-noise spectra directly from unevenly spaced paleoclimatic time series.

    Parameters
    ----------
    data: :py:class:`xarray.DataArray<xarray.DataArray>`
        Input time series data
    timearray: :py:class:`numpy.array<numpy.array>`
        Time series data array
    nsim: :py:class:`int<int>`
        Number of Monte-Carlo simulations (1000-2000 should be o.k. in most cases)
    mctest: :py:class:`bool<bool>`
        Toggle calculation of false-alarm levels based on Monte-Carlo simulation,
        if set to `True` : perform Monte-Carlo test,
        if set to `False` : skip Monte-Carlo test (default).
    rhopre: :py:class:`float<float>`
        Prescibed value for :math:`\\rho`; unused if < 0 (default = -99.0)
    ofac: :py:class:`float<float>`
        Oversampling factor for Lomb-Scargle Fourier transform (typical values: 2.0-4.0)
    hifac: :py:class:`float<float>`
        Max. frequency to analyze is set to hifac * <fNyq> (default = 1.0)
    n50: :py:class:`int<int>`
        Number of WOSA segments (with 50 % overlap)
    iwin: {"rectangular", "welch", "hanning", "triangular", "blackmanharris"}
        Window-type identifier used to suppress sidelobes in spectral analysis:
        ({"rectangular", "welch", "hanning", "triangular", "blackmanharris"}, optional)

    .. caution::
        Parameters `ofac`, `hifac`, `n50` and window type are identical to the SPECTRUM program
        (see Schulz and Stattegger, 1997 for further details).
        Except mctest, hifac and rhopre all parameters must be specified.

    Returns
    -------
    The red-noise spectra (:py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso::
        - Schulz, M., & Mudelsee, M. (2002). REDFIT: estimating red-noise spectra directly from unevenly spaced paleoclimatic time series [Software]. Computers & Geosciences, 28(3), 421-426. https://doi.org/10.1016/S0098-3004(01)00044-9
        - https://www.marum.de/Prof.-Dr.-michael-schulz/Michael-Schulz-Software.html

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_redfit.py
    """

    # Generate datainput
    # 生成数据
    # -------------------------

    random_str = gen_random_string(15)

    data_size = data.shape[0]
    data_array = data.data
    timestep_array = np.arange(data_size)

    if timearray is not None:
        timestep_array = timearray

    data_pd = pd.DataFrame({"time_step": timestep_array, "value": data_array})
    data_pd.to_csv("tmp_data_" + random_str + ".tmp", sep=" ", header=None, index=False)

    # Generate cfg file
    # 生成配置文件
    # -------------------------

    cfg_generate(
        "tmp_cfg_" + random_str + ".cfg",
        fnin="tmp_data_" + random_str + ".tmp",
        fnout="tmp_summary_" + random_str + ".tmp",
        nsim=nsim,
        mctest=mctest,
        rhopre=rhopre,
        ofac=ofac,
        hifac=hifac,
        n50=n50,
        iwin=iwin,
    )

    # Python字符串处理
    input_str = "tmp_cfg_" + random_str + ".cfg"
    byte_str = input_str.encode("utf-8")  # 编码为字节

    # 调整长度为80，用空字符填充
    padded = byte_str.ljust(80, b"\0")[:80]

    # 创建C字符数组
    cfg_array = (c_char * 80)(*padded)

    # 调用Fortran子例程
    _ecl_redfit.run_redfit(cfg_array)

    # Read data
    # 读取生成数据
    # -------------------------
    # 浮点数数据
    save_real = np.fromfile("tmp_real.output", dtype=np.float32)
    (
        ofac,
        hifac,
        idum,
        varx,
        avgdt,
        rho,
        tau,
        dof,
        sixdB_Bandwidth,
        cfal,
        faccrit,
        rcritlo_90,
        rcrithi_90,
        rcritlo_95,
        rcrithi_95,
        rcritlo_98,
        rcrithi_98,
    ) = save_real
    # 整型数据
    save_int = np.fromfile("tmp_int.output", dtype=np.int32)
    n50, iwin, nsim, rcnt, ntime, nout = save_int
    # 数组数据
    if mctest == True:
        save_array_raw = np.fromfile("tmp_array.output", dtype=np.float32)
        save_array = save_array_raw.reshape((14, nout))
    else:
        save_array_raw = np.fromfile("tmp_array.output", dtype=np.float32)
        save_array = save_array_raw.reshape((10, nout))

    # 生成 Dataset
    gxx = xr.DataArray(
        save_array[1, :],
        name="gxx",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "spectrum of input data"},
    )
    gxx_corr = xr.DataArray(
        save_array[2, :],
        name="gxx_corr",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "bias-corrected spectrum of input data"},
    )
    gred_th = xr.DataArray(
        save_array[3, :],
        name="gred_th",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "theoretical AR(1) spectrum"},
    )
    gred = xr.DataArray(
        save_array[4, :],
        name="gred",
        coords=[("freq", save_array[0, :])],
        attrs={
            "Description": "average spectrum of Nsim AR(1) time series (uncorrected)"
        },
    )
    corrfac = xr.DataArray(
        save_array[5, :],
        name="corrFac",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "Gxx / Gxx_corr"},
    )
    chi2_80 = xr.DataArray(
        save_array[6, :],
        name="chi2_80",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "80-% false-alarm level (Chi^2)"},
    )
    chi2_90 = xr.DataArray(
        save_array[7, :],
        name="chi2_90",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "90-% false-alarm level (Chi^2)"},
    )
    chi2_95 = xr.DataArray(
        save_array[8, :],
        name="chi2_95",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "95-% false-alarm level (Chi^2)"},
    )
    chi2_99 = xr.DataArray(
        save_array[9, :],
        name="chi2_99",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "99-% false-alarm level (Chi^2)"},
    )

    if mctest == True:
        ci80 = xr.DataArray(
            save_array[10, :],
            name="ci80",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "80%-MC = 80-% false-alarm level (MC)"},
        )
        ci90 = xr.DataArray(
            save_array[11, :],
            name="ci90",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "90%-MC = 90-% false-alarm level (MC)"},
        )
        ci95 = xr.DataArray(
            save_array[12, :],
            name="ci95",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "95%-MC = 95-% false-alarm level (MC)"},
        )
        ci99 = xr.DataArray(
            save_array[13, :],
            name="ci99",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "99%-MC = 99-% false-alarm level (MC)"},
        )

    # Merge dataset
    if mctest == True:
        save_dataset = xr.merge(
            [
                gxx,
                gxx_corr,
                gred_th,
                gred,
                corrfac,
                chi2_80,
                chi2_90,
                chi2_95,
                chi2_99,
                ci80,
                ci90,
                ci95,
                ci99,
            ]
        )
    else:
        save_dataset = xr.merge(
            [gxx, gxx_corr, gred_th, gred, corrfac, chi2_80, chi2_90, chi2_95, chi2_99]
        )

    # Add attribution
    save_dataset.attrs["Input"] = (
        "OFAC = "
        + str(ofac)
        + ", HIFAC = "
        + str(hifac)
        + ", n50 = "
        + str(n50)
        + ", Iwin = "
        + str(iwin)
        + ", Nsim = "
        + str(nsim)
    )
    save_dataset.attrs["Initial values"] = (
        "idum = "
        + str(idum)
        + ", Data variance (from data spectrum) = "
        + str(varx)
        + ", Avg. dt = "
        + str(avgdt)
    )

    if rhopre < 0:
        save_dataset.attrs["Results"] = (
            "Avg. autocorr. coeff., rho = "
            + str(rho)
            + ", Avg. tau = "
            + str(tau)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
            + ", Critical false-alarm level (Thomson, 1990) = "
            + str(cfal)
            + ",  ==> corresponding scaling factor for red noise = "
            + str(faccrit)
        )
    else:
        save_dataset.attrs["Results"] = (
            "PRESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre)
            + ", Avg. tau = "
            + str(tau)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
            + ", Critical false-alarm level (Thomson, 1990) = "
            + str(cfal)
            + ",  ==> corresponding scaling factor for red noise = "
            + str(faccrit)
        )

    if (iwin == 0) & (ofac == 1.0) & (n50 == 1):
        save_dataset.attrs["Equality of theoretical and data spectrum"] = (
            "90-% acceptance region: rcritlo = "
            + str(rcritlo_90)
            + ", rcrithi = "
            + str(rcrithi_90)
            + "; 95-% acceptance region: rcritlo = "
            + str(rcritlo_95)
            + ", rcrithi = "
            + str(rcrithi_95)
            + "; 98-% acceptance region: rcritlo = "
            + str(rcritlo_98)
            + ", rcrithi = "
            + str(rcrithi_98)
            + "; r_test = "
            + str(rcnt)
        )
    else:
        if iwin != 0:
            print("[Warning] Test requires iwin = 'rectangular'!")
        elif ofac != 1.0:
            print("[Warning] Test requires OFAC = 1.0!")
        elif n50 != 1:
            print("[Warning] Test requires N50 = 1!")
    save_dataset.attrs["Elapsed time"] = str(ntime) + " [s]"

    save_dataset.attrs["Description"] = (
        "Estimating red-noise spectra directly from unevenly spaced paleoclimatic time series."
    )
    save_dataset.attrs["About"] = (
        "Michael Schulz, Manfred Mudelsee => https://www.marum.de/Prof.-Dr.-michael-schulz/Michael-Schulz-Software.html"
    )
    save_dataset.attrs["Reference"] = (
        "Schulz, M. and Mudelsee, M. (2002) REDFIT: Estimating red-noise spectra directly from unevenly spaced paleoclimatic time series. Computers and Geosciences, 28, 421-426. https://doi.org/10.1016/S0098-3004(01)00044-9"
    )
    save_dataset.attrs["Python platform"] = (
        "Shen yulu => https://github.com/shenyulu/easyclimate-backend"
    )
    # Add coordinate period (reciprocal of frequency)
    save_dataset = save_dataset.assign_coords({"period": 1 / save_dataset.freq})

    # Clean up temporary files
    # 清理临时文件
    # -------------------------
    os.remove("tmp_cfg_" + random_str + ".cfg")
    os.remove("tmp_data_" + random_str + ".tmp")
    os.remove("tmp_real.output")
    os.remove("tmp_int.output")
    os.remove("tmp_array.output")
    os.remove("tmp_summary_" + random_str + ".tmp")

    return save_dataset


def calc_redfit_cross(
    data_x: xr.DataArray,
    data_y: xr.DataArray,
    timearray_x: np.array = None,
    timearray_y: np.array = None,
    x_sign: bool = False,
    y_sign: bool = False,
    nsim: int = 1000,
    mctest: bool = True,
    mctest_phi: bool = True,
    rhopre_1: float = -999.0,
    rhopre_2: float = -999.0,
    ofac: float = 1.0,
    hifac: float = 1.0,
    n50: int = 1,
    alpha: float = 0.05,
    iwin: Literal[
        "rectangular", "welch", "hanning", "triangular", "blackmanharris"
    ] = "rectangular",
):
    """
    Estimating red-noise spectra directly from unevenly spaced paleoclimatic time series.

    Parameters
    ----------
    data_x::py:class:`xarray.DataArray<xarray.DataArray>`
        First input time series data
    data_y: :py:class:`xarray.DataArray<xarray.DataArray>`
        Second input time series data
    timearray_x: :py:class:`numpy.array<numpy.array>`
        First time series data array
    timearray_y: :py:class:`numpy.array<numpy.array>`
        Second time series data array
    x_sign: :py:class:`bool<bool>`
        Change the sign of the first time series:
        if `True`: The sign of the data is changed
        if `False`: The sign of the data is not changed (default)
    y_sign: :py:class:`bool<bool>`
        Change the sign of the second time series:
        if `True`: The sign of the data is changed
        if `False`: The sign of the data is not changed (default)
    nsim: :py:class:`int<int>`
        Number of Monte Carlo simulations (1000-2000 is recommended)
    mctest: :py:class:`bool<bool>`
        Estimate the significance of auto and coherency spectrum with Monte Carlo simulations
        if `True`: perform Monte Carlo simulations
        if `False`: do not perform Monte Carlo simulations
    mctest_phi: :py:class:`bool<bool>`
        Estimate Monte Carlo confidence interval for the phase spectrum
        if `True`: perform Monte Carlo simulations (mctest needs to be true as well)
        if `False`: do not perform Monte Carlo simulations
    rhopre_1: :py:class:`float<float>`
        Prescribed value for :math:`\\rho` for the first time series, not used if :math:`\\rho < 0` (default = -999.0).
    rhopre_2: :py:class:`float<float>`
        Prescribed value for :math:`\\rho` for the second time series, not used if :math:`\\rho< 0` (default = -999.0).
    ofac: :py:class:`float<float>`
        Oversampling factor for Lomb-Scargle Fourier transform (typical values: 2.0-4.0).
    hifac: :py:class:`float<float>`
        Max. frequency to analyze is set to hifac * <fNyq> (default = 1.0).
    n50: :py:class:`int<int>`
        Number of WOSA segments (with 50 % overlap)
    alpha: :py:class:`float<float>`
        Significance level (Note: only 0.01, 0.05 [default], or 0.1 are allowed).
    iwin: {"rectangular", "welch", "hanning", "triangular", "blackmanharris"}
        Window-type identifier used to suppress sidelobes in spectral analysis:
        ({"rectangular", "welch", "hanning", "triangular", "blackmanharris"}, optional).

    .. caution::
        Parameters ofac, hifac, n50 and window type are identical to the SPECTRUM program
        (see Schulz and Stattegger, 1997 for further details).
        Except mctest, hifac, rhopre(1) and rhopre(2) all parameters must be specified.

    .. seealso::
        - Schulz, M., & Mudelsee, M. (2002). REDFIT: estimating red-noise spectra directly from unevenly spaced paleoclimatic time series [Software]. Computers & Geosciences, 28(3), 421-426. https://doi.org/10.1016/S0098-3004(01)00044-9
        - https://www.marum.de/Prof.-Dr.-michael-schulz/Michael-Schulz-Software.html
    """

    # Generate datainput
    # 生成数据
    # -------------------------

    random_str = gen_random_string(15)

    #  -------------------------
    def createdata(data, timearray, str_signal):
        data_size = data.shape[0]
        data_array = data.data
        timestep_array = np.arange(data_size)

        if timearray is not None:
            timestep_array = timearray

        data_pd = pd.DataFrame({"time_step": timestep_array, "value": data_array})
        data_pd.to_csv(
            "tmp_data_" + str_signal + "_" + random_str + ".tmp",
            sep=" ",
            header=None,
            index=False,
        )

    createdata(data_x, timearray_x, "x")
    createdata(data_y, timearray_y, "y")

    if mctest == True and mctest_phi == True:
        print(
            "[Warning] Monte Carlo simulations need more time to calculate, please wait!"
        )
    elif mctest == False and mctest_phi == True:
        error = f"""When mctest is `False`, mctest_phi should be `True`"""
        raise ValueError(error)

    # Generate cfg file
    # 生成配置文件
    # -------------------------

    cfg_generate_cross(
        "tmp_cfg_" + random_str + ".cfg",
        fnin_1="tmp_data_x_" + random_str + ".tmp",
        fnin_2="tmp_data_y_" + random_str + ".tmp",
        fnout="tmp_cross_summary_" + random_str + ".tmp",
        x_sign=x_sign,
        y_sign=y_sign,
        nsim=nsim,
        mctest=mctest,
        mctest_phi=mctest_phi,
        rhopre_1=rhopre_1,
        rhopre_2=rhopre_2,
        ofac=ofac,
        hifac=hifac,
        n50=n50,
        alpha=alpha,
        iwin=iwin,
    )

    # Python字符串处理
    input_str = "tmp_cfg_" + random_str + ".cfg"
    byte_str = input_str.encode("utf-8")  # 编码为字节

    # 调整长度为80，用空字符填充
    padded = byte_str.ljust(80, b"\0")[:80]

    # 创建C字符数组
    cfg_array = (c_char * 80)(*padded)

    # 调用Fortran子例程
    _ecl_redfit_x.run_redfit_x(cfg_array)

    # Read data -gxx
    # 读取生成数据 gxx
    # -------------------------
    # 浮点数数据
    save_real = np.fromfile("tmp_real_gxx.output", dtype=np.float32)
    (
        ofac,
        hifac,
        varx,
        avgdt,
        rhox,
        rhopre,
        taux,
        dof,
        sixdB_Bandwidth,
        cfal,
        faccrit,
    ) = save_real
    # 整型数据
    save_int = np.fromfile("tmp_int_gxx.output", dtype=np.int32)
    n50, iwin, nsim, ntime, nout = save_int
    # 数组数据
    if mctest == True:
        save_array_raw = np.fromfile("tmp_array_gxx.output", dtype=np.float32)
        save_array = save_array_raw.reshape((12, nout))
    else:
        save_array_raw = np.fromfile("tmp_array_gxx.output", dtype=np.float32)
        save_array = save_array_raw.reshape((9, nout))

    # 生成 Dataset
    gxx = xr.DataArray(
        save_array[1, :],
        name="gxx",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "spectrum of input data"},
    )
    gxx_corr = xr.DataArray(
        save_array[2, :],
        name="gxx_corr",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "bias-corrected spectrum of input data"},
    )
    gred_th = xr.DataArray(
        save_array[3, :],
        name="gred_th",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "theoretical AR(1) spectrum"},
    )
    gred = xr.DataArray(
        save_array[4, :],
        name="gred",
        coords=[("freq", save_array[0, :])],
        attrs={
            "Description": "average spectrum of Nsim AR(1) time series (uncorrected)"
        },
    )
    corrfac = xr.DataArray(
        save_array[5, :],
        name="corrFac",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "Gxx / Gxx_corr"},
    )
    chi2_90 = xr.DataArray(
        save_array[6, :],
        name="chi2_90",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "90-% false-alarm level (Chi^2)"},
    )
    chi2_95 = xr.DataArray(
        save_array[7, :],
        name="chi2_95",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "95-% false-alarm level (Chi^2)"},
    )
    chi2_99 = xr.DataArray(
        save_array[8, :],
        name="chi2_99",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "99-% false-alarm level (Chi^2)"},
    )

    if mctest == True:
        ci90 = xr.DataArray(
            save_array[9, :],
            name="ci90",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "90%-MC = 90-% false-alarm level (MC)"},
        )
        ci95 = xr.DataArray(
            save_array[10, :],
            name="ci95",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "95%-MC = 95-% false-alarm level (MC)"},
        )
        ci99 = xr.DataArray(
            save_array[11, :],
            name="ci99",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "99%-MC = 99-% false-alarm level (MC)"},
        )

    # Merge dataset
    if mctest == True:
        save_dataset_gxx = xr.merge(
            [
                gxx,
                gxx_corr,
                gred_th,
                gred,
                corrfac,
                chi2_90,
                chi2_95,
                chi2_99,
                ci90,
                ci95,
                ci99,
            ]
        )
    else:
        save_dataset_gxx = xr.merge(
            [gxx, gxx_corr, gred_th, gred, corrfac, chi2_90, chi2_95, chi2_99]
        )

    # Add attribution
    save_dataset_gxx.attrs["Input"] = (
        "OFAC = "
        + str(ofac)
        + ", HIFAC = "
        + str(hifac)
        + ", n50 = "
        + str(n50)
        + ", Iwin = "
        + str(iwin)
        + ", Nsim = "
        + str(nsim)
    )
    save_dataset_gxx.attrs["Initial values"] = (
        "Data variance (from data spectrum) = "
        + str(varx)
        + ", Avg. dt = "
        + str(avgdt)
    )

    if rhopre < 0:
        save_dataset_gxx.attrs["Results"] = (
            "Avg. autocorr. coeff., rho = "
            + str(rhox)
            + ", Avg. tau = "
            + str(taux)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
            + ", Critical false-alarm level (Thomson, 1990) = "
            + str(cfal)
            + ",  ==> corresponding scaling factor for red noise = "
            + str(faccrit)
        )
    else:
        save_dataset_gxx.attrs["Results"] = (
            "PRESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre)
            + ", Avg. tau = "
            + str(taux)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
            + ", Critical false-alarm level (Thomson, 1990) = "
            + str(cfal)
            + ",  ==> corresponding scaling factor for red noise = "
            + str(faccrit)
        )

    save_dataset_gxx.attrs["Elapsed time"] = str(ntime) + " [s]"

    save_dataset_gxx.attrs["Description"] = (
        "(Autospectrum for the 1st time series) Cross-spectral analysis of unevenly spaced paleoclimate time series."
    )
    save_dataset_gxx.attrs["About"] = (
        "Michael Schulz, Manfred Mudelsee => https://www.marum.de/Prof.-Dr.-michael-schulz/Michael-Schulz-Software.html"
    )
    save_dataset_gxx.attrs["Reference"] = (
        "Olafsdottir, K. B., Schulz, M. and Mudelsee, M. (2016): REDFIT-X: Cross-spectral analysis of unevenly spaced paleoclimate time series. Computers and Geosciences, 91, 11-18. https://doi.org/10.1016/S0098-3004(01)00044-9"
    )
    save_dataset_gxx.attrs["Python platform"] = (
        "Shen yulu => https://github.com/shenyulu/easyclimate-backend"
    )

    # Add coordinate period (reciprocal of frequency)
    save_dataset_gxx = save_dataset_gxx.assign_coords(
        {"period": 1 / save_dataset_gxx.freq}
    )

    # Read data -gyy
    # 读取生成数据 gyy
    # -------------------------
    # 浮点数数据
    save_real = np.fromfile("tmp_real_gyy.output", dtype=np.float32)
    (
        ofac,
        hifac,
        vary,
        avgdt,
        rhoy,
        rhopre,
        tauy,
        dof,
        sixdB_Bandwidth,
        cfal,
        faccrit,
    ) = save_real
    # 整型数据
    save_int = np.fromfile("tmp_int_gyy.output", dtype=np.int32)
    n50, iwin, nsim, ntime, nout = save_int
    # 数组数据
    if mctest == True:
        save_array_raw = np.fromfile("tmp_array_gyy.output", dtype=np.float32)
        save_array = save_array_raw.reshape((12, nout))
    else:
        save_array_raw = np.fromfile("tmp_array_gyy.output", dtype=np.float32)
        save_array = save_array_raw.reshape((9, nout))

    # 生成 Dataset
    gyy = xr.DataArray(
        save_array[1, :],
        name="gyy",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "spectrum of input data"},
    )
    gyy_corr = xr.DataArray(
        save_array[2, :],
        name="gyy_corr",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "bias-corrected spectrum of input data"},
    )
    gred_th = xr.DataArray(
        save_array[3, :],
        name="gred_th",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "theoretical AR(1) spectrum"},
    )
    gred = xr.DataArray(
        save_array[4, :],
        name="gred",
        coords=[("freq", save_array[0, :])],
        attrs={
            "Description": "average spectrum of Nsim AR(1) time series (uncorrected)"
        },
    )
    corrfac = xr.DataArray(
        save_array[5, :],
        name="corrFac",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "Gxx / Gxx_corr"},
    )
    chi2_90 = xr.DataArray(
        save_array[6, :],
        name="chi2_90",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "90-% false-alarm level (Chi^2)"},
    )
    chi2_95 = xr.DataArray(
        save_array[7, :],
        name="chi2_95",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "95-% false-alarm level (Chi^2)"},
    )
    chi2_99 = xr.DataArray(
        save_array[8, :],
        name="chi2_99",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "99-% false-alarm level (Chi^2)"},
    )

    if mctest == True:
        ci90 = xr.DataArray(
            save_array[9, :],
            name="ci90",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "90%-MC = 90-% false-alarm level (MC)"},
        )
        ci95 = xr.DataArray(
            save_array[10, :],
            name="ci95",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "95%-MC = 95-% false-alarm level (MC)"},
        )
        ci99 = xr.DataArray(
            save_array[11, :],
            name="ci99",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "99%-MC = 99-% false-alarm level (MC)"},
        )

    # Merge dataset
    if mctest == True:
        save_dataset_gyy = xr.merge(
            [
                gyy,
                gyy_corr,
                gred_th,
                gred,
                corrfac,
                chi2_90,
                chi2_95,
                chi2_99,
                ci90,
                ci95,
                ci99,
            ]
        )
    else:
        save_dataset_gyy = xr.merge(
            [gyy, gyy_corr, gred_th, gred, corrfac, chi2_90, chi2_95, chi2_99]
        )

    # Add attribution
    save_dataset_gyy.attrs["Input"] = (
        "OFAC = "
        + str(ofac)
        + ", HIFAC = "
        + str(hifac)
        + ", n50 = "
        + str(n50)
        + ", Iwin = "
        + str(iwin)
        + ", Nsim = "
        + str(nsim)
    )
    save_dataset_gyy.attrs["Initial values"] = (
        "Data variance (from data spectrum) = "
        + str(vary)
        + ", Avg. dt = "
        + str(avgdt)
    )

    if rhopre < 0:
        save_dataset_gyy.attrs["Results"] = (
            "Avg. autocorr. coeff., rho = "
            + str(rhoy)
            + ", Avg. tau = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
            + ", Critical false-alarm level (Thomson, 1990) = "
            + str(cfal)
            + ",  ==> corresponding scaling factor for red noise = "
            + str(faccrit)
        )
    else:
        save_dataset_gyy.attrs["Results"] = (
            "PRESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre)
            + ", Avg. tau = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
            + ", Critical false-alarm level (Thomson, 1990) = "
            + str(cfal)
            + ",  ==> corresponding scaling factor for red noise = "
            + str(faccrit)
        )

    save_dataset_gyy.attrs["Elapsed time"] = str(ntime) + " [s]"

    save_dataset_gyy.attrs["Description"] = (
        "(Autospectrum for the 2nd time series) Cross-spectral analysis of unevenly spaced paleoclimate time series."
    )
    save_dataset_gyy.attrs["About"] = (
        "Michael Schulz, Manfred Mudelsee => https://www.marum.de/Prof.-Dr.-michael-schulz/Michael-Schulz-Software.html"
    )
    save_dataset_gyy.attrs["Reference"] = (
        "Olafsdottir, K. B., Schulz, M. and Mudelsee, M. (2016): REDFIT-X: Cross-spectral analysis of unevenly spaced paleoclimate time series. Computers and Geosciences, 91, 11-18."
    )
    save_dataset_gyy.attrs["Python platform"] = (
        "Shen yulu => https://github.com/shenyulu/easyclimate-backend"
    )
    # Add coordinate period (reciprocal of frequency)
    save_dataset_gyy = save_dataset_gyy.assign_coords(
        {"period": 1 / save_dataset_gyy.freq}
    )

    # Read data -gxy
    # 读取生成数据 gxy
    # -------------------------
    # 浮点数数据
    save_real = np.fromfile("tmp_real_gxy.output", dtype=np.float32)
    (
        ofac,
        hifac,
        alpha,
        avgdty,
        rhox,
        rhopre_1,
        rhoy,
        rhopre_2,
        taux,
        tauy,
        dof,
        sixdB_Bandwidth,
    ) = save_real
    # 整型数据
    save_int = np.fromfile("tmp_int_gxy.output", dtype=np.int32)
    n50, iwin, nsim, ntime, nout = save_int
    # 数组数据
    save_array_raw = np.fromfile("tmp_array_gxy.output", dtype=np.float32)
    save_array = save_array_raw.reshape((2, nout))

    # 生成 Dataset
    gxy = xr.DataArray(
        save_array[1, :],
        name="gxy",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "cross-spectrum of input data"},
    )

    save_dataset_gxy = gxy.to_dataset()

    # Add attribution
    save_dataset_gxy.attrs["Description"] = (
        "(Cross-spectrum) Cross-spectral analysis of unevenly spaced paleoclimate time series."
    )

    save_dataset_gxy.attrs["Input"] = (
        "OFAC = "
        + str(ofac)
        + ", HIFAC = "
        + str(hifac)
        + ", n50 = "
        + str(n50)
        + ", Iwin = "
        + str(iwin)
        + ", Nsim = "
        + str(nsim)
        + ", Level of signif. = "
        + str(alpha)
    )
    save_dataset_gxy.attrs["Initial values"] = "Avg. dtxy = " + str(avgdty)

    if (rhopre_1 < 0) and (rhopre_2 < 0):
        save_dataset_gxy.attrs["Results"] = (
            "Avg. autocorr. coeff., rhox = "
            + str(rhox)
            + ", Avg. autocorr. coeff., rhoy = "
            + str(rhoy)
            + ", Avg. taux = "
            + str(taux)
            + ", Avg. tauy = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
        )
    elif (rhopre_1 >= 0) and (rhopre_2 < 0):
        save_dataset_gxy.attrs["Results"] = (
            "RESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre_1)
            + ", Avg. autocorr. coeff., rhoy = "
            + str(rhoy)
            + ", Avg. taux = "
            + str(taux)
            + ", Avg. tauy = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
        )
    elif (rhopre_1 < 0) and (rhopre_2 >= 0):
        save_dataset_gxy.attrs["Results"] = (
            "Avg. autocorr. coeff., rhox = "
            + str(rhox)
            + ", RESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre_2)
            + ", Avg. taux = "
            + str(taux)
            + ", Avg. tauy = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
        )
    elif (rhopre_1 >= 0) and (rhopre_2 >= 0):
        save_dataset_gxy.attrs["Results"] = (
            "RESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre_1)
            + ", RESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre_2)
            + ", Avg. taux = "
            + str(taux)
            + ", Avg. tauy = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
        )
    else:
        raise ValueError("There are some internal mistake in data process.")

    save_dataset_gxy.attrs["Elapsed time"] = str(ntime) + " [s]"
    save_dataset_gxy.attrs["About"] = (
        "Michael Schulz, Manfred Mudelsee => https://www.marum.de/Prof.-Dr.-michael-schulz/Michael-Schulz-Software.html"
    )
    save_dataset_gxy.attrs["Reference"] = (
        "Olafsdottir, K. B., Schulz, M. and Mudelsee, M. (2016): REDFIT-X: Cross-spectral analysis of unevenly spaced paleoclimate time series. Computers and Geosciences, 91, 11-18."
    )
    save_dataset_gxy.attrs["Python platform"] = (
        "Shen yulu => https://github.com/shenyulu/easyclimate-backend"
    )
    # Add coordinate period (reciprocal of frequency)
    save_dataset_gxy = save_dataset_gxy.assign_coords(
        {"period": 1 / save_dataset_gxy.freq}
    )

    # Read data -cxy
    # 读取生成数据 cxy
    # -------------------------
    # 浮点数数据
    save_real = np.fromfile("tmp_real_cxy.output", dtype=np.float32)
    (
        ofac,
        hifac,
        alpha,
        avgdty,
        rhox,
        rhopre_1,
        rhoy,
        rhopre_2,
        taux,
        tauy,
        dof,
        sixdB_Bandwidth,
        csig_value,
        csig_mc_value,
    ) = save_real
    # 整型数据
    save_int = np.fromfile("tmp_int_cxy.output", dtype=np.int32)
    n50, iwin, nsim, ntime, nout = save_int
    # 数组数据
    if mctest == True:
        save_array_raw = np.fromfile("tmp_array_cxy.output", dtype=np.float32)
        save_array = save_array_raw.reshape((5, nout))
    else:
        save_array_raw = np.fromfile("tmp_array_cxy.output", dtype=np.float32)
        save_array = save_array_raw.reshape((2, nout))

    # 生成 Dataset
    cxy = xr.DataArray(
        save_array[1, :],
        name="cxy",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "Coherency-spectrum of input data"},
    )
    csig = xr.DataArray(
        np.full((nout), csig_value),
        name="csig",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "Theoretical False alarm level"},
    )
    csig_mc = xr.DataArray(
        np.full((nout), csig_mc_value),
        name="csig_mc",
        coords=[("freq", save_array[0, :])],
        attrs={
            "Description": "MC-CSig = Mean Monte Carlo false-alarm level (for alpha = "
            + str(alpha)
        },
    )

    if mctest == True:
        ci90 = xr.DataArray(
            save_array[2, :],
            name="ci90",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "90%-MC = 90-% false-alarm level (MC)"},
        )
        ci95 = xr.DataArray(
            save_array[3, :],
            name="ci95",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "95%-MC = 95-% false-alarm level (MC)"},
        )
        ci99 = xr.DataArray(
            save_array[4, :],
            name="ci99",
            coords=[("freq", save_array[0, :])],
            attrs={"Description": "99%-MC = 99-% false-alarm level (MC)"},
        )

    # Merge dataset
    if mctest == True:
        save_dataset_cxy = xr.merge([cxy, csig, csig_mc, ci90, ci95, ci99])
    else:
        save_dataset_cxy = xr.merge([cxy, csig, csig_mc])

    # Add attribution
    save_dataset_cxy.attrs["Description"] = (
        "(Coherency spectrum) Cross-spectral analysis of unevenly spaced paleoclimate time series."
    )

    save_dataset_cxy.attrs["Input"] = (
        "OFAC = "
        + str(ofac)
        + ", HIFAC = "
        + str(hifac)
        + ", n50 = "
        + str(n50)
        + ", Iwin = "
        + str(iwin)
        + ", Nsim = "
        + str(nsim)
        + ", Level of signif. = "
        + str(alpha)
    )
    save_dataset_cxy.attrs["Initial values"] = "Avg. dtxy = " + str(avgdty)

    if (rhopre_1 < 0) and (rhopre_2 < 0):
        save_dataset_cxy.attrs["Results"] = (
            "Avg. autocorr. coeff., rhox = "
            + str(rhox)
            + ", Avg. autocorr. coeff., rhoy = "
            + str(rhoy)
            + ", Avg. taux = "
            + str(taux)
            + ", Avg. tauy = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
        )
    elif (rhopre_1 >= 0) and (rhopre_2 < 0):
        save_dataset_cxy.attrs["Results"] = (
            "RESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre_1)
            + ", Avg. autocorr. coeff., rhoy = "
            + str(rhoy)
            + ", Avg. taux = "
            + str(taux)
            + ", Avg. tauy = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
        )
    elif (rhopre_1 < 0) and (rhopre_2 >= 0):
        save_dataset_cxy.attrs["Results"] = (
            "Avg. autocorr. coeff., rhox = "
            + str(rhox)
            + ", RESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre_2)
            + ", Avg. taux = "
            + str(taux)
            + ", Avg. tauy = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
        )
    elif (rhopre_1 >= 0) and (rhopre_2 >= 0):
        save_dataset_cxy.attrs["Results"] = (
            "RESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre_1)
            + ", RESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre_2)
            + ", Avg. taux = "
            + str(taux)
            + ", Avg. tauy = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
        )
    else:
        raise ValueError("There are some internal mistake in data process.")

    save_dataset_cxy.attrs["Elapsed time"] = str(ntime) + " [s]"
    save_dataset_cxy.attrs["About"] = (
        "Michael Schulz, Manfred Mudelsee => https://www.marum.de/Prof.-Dr.-michael-schulz/Michael-Schulz-Software.html"
    )
    save_dataset_cxy.attrs["Reference"] = (
        "Olafsdottir, K. B., Schulz, M. and Mudelsee, M. (2016): REDFIT-X: Cross-spectral analysis of unevenly spaced paleoclimate time series. Computers and Geosciences, 91, 11-18."
    )
    save_dataset_cxy.attrs["Python platform"] = (
        "Shen yulu => https://github.com/shenyulu/easyclimate-backend"
    )
    # Add coordinate period (reciprocal of frequency)
    save_dataset_cxy = save_dataset_cxy.assign_coords(
        {"period": 1 / save_dataset_cxy.freq}
    )

    # Read data -phxy
    # 读取生成数据 phxy
    # -------------------------
    # 浮点数数据
    save_real = np.fromfile("tmp_real_phxy.output", dtype=np.float32)
    (
        ofac,
        hifac,
        alpha,
        avgdty,
        rhox,
        rhopre_1,
        rhoy,
        rhopre_2,
        taux,
        tauy,
        dof,
        sixdB_Bandwidth,
    ) = save_real
    # 整型数据
    save_int = np.fromfile("tmp_int_phxy.output", dtype=np.int32)
    n50, iwin, nsim, ntime, nout = save_int
    # 数组数据
    if mctest == True and mctest_phi == True:
        save_array_raw = np.fromfile("tmp_array_phxy.output", dtype=np.float32)
        save_array = save_array_raw.reshape((6, nout))
    else:
        save_array_raw = np.fromfile("tmp_array_phxy.output", dtype=np.float32)
        save_array = save_array_raw.reshape((4, nout))

    # 生成 Dataset
    phxy = xr.DataArray(
        save_array[1, :],
        name="phxy",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "Coherency-spectrum of input data"},
    )
    ephi_lower = xr.DataArray(
        save_array[1, :],
        name="ephi_lower",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "CI-low = Theoretical Confidence Interval - lower"},
    )
    ephi_upper = xr.DataArray(
        save_array[1, :],
        name="ephi_upper",
        coords=[("freq", save_array[0, :])],
        attrs={"Description": "CI-up =  Theoretical Confidence Interval - upper"},
    )

    if mctest == True and mctest_phi == True:
        ephi_mc_lower = xr.DataArray(
            save_array[1, :],
            name="ephi_mc_lower",
            coords=[("freq", save_array[0, :])],
            attrs={
                "Description": "CI-mc-low = Monte Carlo Confidence Interval - lower   (percentiles)"
            },
        )
        ephi_mc_upper = xr.DataArray(
            save_array[1, :],
            name="ephi_mc_upper",
            coords=[("freq", save_array[0, :])],
            attrs={
                "Description": "CI-mc-up =  Monte Carlo Confidence Interval - upper   (percentiles)"
            },
        )

    # Merge dataset
    if mctest == True and mctest_phi == True:
        save_dataset_phxy = xr.merge(
            [phxy, ephi_lower, ephi_upper, ephi_mc_lower, ephi_mc_upper]
        )
    else:
        save_dataset_phxy = xr.merge([phxy, ephi_lower, ephi_upper])

    # Add attribution
    save_dataset_phxy.attrs["Description"] = (
        "(Phase spectrum) Cross-spectral analysis of unevenly spaced paleoclimate time series."
    )

    save_dataset_phxy.attrs["Input"] = (
        "OFAC = "
        + str(ofac)
        + ", HIFAC = "
        + str(hifac)
        + ", n50 = "
        + str(n50)
        + ", Iwin = "
        + str(iwin)
        + ", Nsim = "
        + str(nsim)
        + ", Level of signif. = "
        + str(alpha)
    )
    save_dataset_phxy.attrs["Initial values"] = "Avg. dtxy = " + str(avgdty)

    if (rhopre_1 < 0) and (rhopre_2 < 0):
        save_dataset_phxy.attrs["Results"] = (
            "Avg. autocorr. coeff., rhox = "
            + str(rhox)
            + ", Avg. autocorr. coeff., rhoy = "
            + str(rhoy)
            + ", Avg. taux = "
            + str(taux)
            + ", Avg. tauy = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
        )
    elif (rhopre_1 >= 0) and (rhopre_2 < 0):
        save_dataset_phxy.attrs["Results"] = (
            "RESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre_1)
            + ", Avg. autocorr. coeff., rhoy = "
            + str(rhoy)
            + ", Avg. taux = "
            + str(taux)
            + ", Avg. tauy = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
        )
    elif (rhopre_1 < 0) and (rhopre_2 >= 0):
        save_dataset_phxy.attrs["Results"] = (
            "Avg. autocorr. coeff., rhox = "
            + str(rhox)
            + ", RESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre_2)
            + ", Avg. taux = "
            + str(taux)
            + ", Avg. tauy = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
        )
    elif (rhopre_1 >= 0) and (rhopre_2 >= 0):
        save_dataset_phxy.attrs["Results"] = (
            "RESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre_1)
            + ", RESCRIBED avg. autocorr. coeff., rho = "
            + str(rhopre_2)
            + ", Avg. taux = "
            + str(taux)
            + ", Avg. tauy = "
            + str(tauy)
            + ", Degrees of freedom = "
            + str(dof)
            + ", 6-dB Bandwidth = "
            + str(sixdB_Bandwidth)
        )
    else:
        raise ValueError("There are some internal mistake in data process.")

    save_dataset_phxy.attrs["Elapsed time"] = str(ntime) + " [s]"
    save_dataset_phxy.attrs["About"] = (
        "Michael Schulz, Manfred Mudelsee => https://www.marum.de/Prof.-Dr.-michael-schulz/Michael-Schulz-Software.html"
    )
    save_dataset_phxy.attrs["Reference"] = (
        "Olafsdottir, K. B., Schulz, M. and Mudelsee, M. (2016): REDFIT-X: Cross-spectral analysis of unevenly spaced paleoclimate time series. Computers and Geosciences, 91, 11-18."
    )
    save_dataset_phxy.attrs["Python platform"] = (
        "Shen yulu => https://github.com/shenyulu/easyclimate-backend"
    )
    # Add coordinate period (reciprocal of frequency)
    save_dataset_phxy = save_dataset_phxy.assign_coords(
        {"period": 1 / save_dataset_phxy.freq}
    )

    # Clean up temporary files
    # 清理临时文件
    # -------------------------
    os.remove("tmp_real_gxx.output")
    os.remove("tmp_int_gxx.output")
    os.remove("tmp_array_gxx.output")

    os.remove("tmp_real_gyy.output")
    os.remove("tmp_int_gyy.output")
    os.remove("tmp_array_gyy.output")

    os.remove("tmp_real_gxy.output")
    os.remove("tmp_int_gxy.output")
    os.remove("tmp_array_gxy.output")

    os.remove("tmp_real_cxy.output")
    os.remove("tmp_int_cxy.output")
    os.remove("tmp_array_cxy.output")

    os.remove("tmp_real_phxy.output")
    os.remove("tmp_int_phxy.output")
    os.remove("tmp_array_phxy.output")

    os.remove("tmp_data_x_" + random_str + ".tmp")
    os.remove("tmp_data_y_" + random_str + ".tmp")
    os.remove("tmp_cfg_" + random_str + ".cfg")

    os.remove("tmp_cross_summary_" + random_str + ".tmp" + ".cxy")
    os.remove("tmp_cross_summary_" + random_str + ".tmp" + ".gxx")
    os.remove("tmp_cross_summary_" + random_str + ".tmp" + ".gxy")
    os.remove("tmp_cross_summary_" + random_str + ".tmp" + ".gyy")
    os.remove("tmp_cross_summary_" + random_str + ".tmp" + ".phxy")

    return (
        save_dataset_gxx,
        save_dataset_gyy,
        save_dataset_gxy,
        save_dataset_cxy,
        save_dataset_phxy,
    )

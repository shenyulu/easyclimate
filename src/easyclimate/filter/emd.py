"""
Empirical Mode Decomposition and Ensemble Empirical Mode Decomposition (EEMD)
"""

from __future__ import annotations
import xarray as xr
import numpy as np
from typing import Literal
from PyEMD import EMD, EEMD
from ..core.utility import datetime_to_numeric

__all__ = ["filter_emd", "filter_eemd"]


def filter_emd(
    input_data: xr.DataArray,
    time_step: Literal["ns", "us", "ms", "s", "m", "h", "D", "W", "M", "Y"],
    time_array=None,
    time_dim: str = "time",
    spline_kind: Literal[
        "cubic", "akima", "pchip", "cubic_hermite", "slinear", "quadratic", "linear"
    ] = "cubic",
    nbsym: int = 2,
    max_iteration: int = 1000,
    energy_ratio_thr: float = 0.2,
    std_thr: float = 0.2,
    svar_thr: float = 0.001,
    total_power_thr: float = 0.005,
    range_thr: float = 0.001,
    extrema_detection: Literal["simple", "parabol"] = "simple",
    max_imf: int = -1,
    dtype=np.float64,
):
    """
    Empirical Mode Decomposition

    Method of decomposing signal into Intrinsic Mode Functions (IMFs) based on algorithm presented in Huang et al (1998).

    Algorithm was validated with Rilling et al (2003). Matlab’s version from 3.2007.

    Threshold which control the goodness of the decomposition:

    - ``std_thr``: Test for the proto-IMF how variance changes between siftings.
    - ``svar_thr``: Test for the proto-IMF how energy changes between siftings.
    - ``total_power_thr``: Test for the whole decomp how much of energy is solved.
    - ``range_thr``: Test for the whole decomp whether the difference is tiny.

    Parameters
    ----------
    input_data: :py:class:`xarray.DataArray <xarray.DataArray>`
        Input signal data to be decomposed.
    time_step: :py:class:`str <str>`
        Time step unit for datetime conversion (e.g., 's', 'ms', 'us', 'ns').
    time_array: :py:class:`array-like <array-like>`, optional
        Custom time array for the input signal. If None, uses input_data's time dimension.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    spline_kind: ``Literal["cubic", "akima", "pchip", "cubic_hermite", "slinear", "quadratic", "linear"]``, default: "cubic"
        Type of spline used for envelope interpolation.
    nbsym: :py:class:`int <int>`, default: 2
        Number of points to add at signal boundaries for mirroring.
    max_iteration: :py:class:`int <int>`, default: 1000
        Maximum number of iterations per single sifting in EMD.
    energy_ratio_thr: :py:class:`float <float>`, default: 0.2
        Threshold value on energy ratio per IMF check.
    std_thr: :py:class:`float <float>`, default: 0.2
        Threshold value on standard deviation per IMF check.
    svar_thr: :py:class:`float <float>`, default: 0.001
        Threshold value on scaled variance per IMF check.
    total_power_thr: :py:class:`float <float>`, default: 0.005
        Threshold value on total power per EMD decomposition.
    range_thr: :py:class:`float <float>`, default: 0.001
        Threshold for amplitude range (after scaling) per EMD decomposition.
    extrema_detection: ``Literal["simple", "parabol"]``, default: "simple"
        Method used to finding extrema.
    max_imf: :py:class:`int <int>`, default: -1
        IMF number to which decomposition should be performed. Negative value means all.
    dtype: :py:class:`numpy.dtype <numpy.dtype>`, default: np.float64
        Data type used for calculations.


    Reference
    --------------
    - Huang Norden E., Shen Zheng, Long Steven R., Wu Manli C., Shih Hsing H., Zheng Quanan, Yen Nai-Chyuan, Tung Chi Chao and Liu Henry H. 1998 The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysisProc. R. Soc. Lond. A.454903–995 http://doi.org/10.1098/rspa.1998.0193
    - Gabriel Rilling, Patrick Flandrin, Paulo Gonçalves. On empirical mode decomposition and its algorithms. IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing NSIP-03, Jun 2003, Grado, Italy. https://inria.hal.science/inria-00570628v1
    - Colominas, M. A., Schlotthauer, G., and Torres, M. E. (2014). Improved complete ensemble EMD: A suitable tool for biomedical signal processing. Biomedical Signal Processing and Control, 14, 19-29. https://doi.org/10.1016/j.bspc.2014.06.009
    - Hsuan, R. (2014). Ensemble Empirical Mode Decomposition Parameters Optimization for Spectral Distance Measurement in Hyperspectral Remote Sensing Data. Remote Sens. 6(3), 2069-2083. http://doi.org/10.3390/rs6032069. http://www.mdpi.com/2072-4292/6/3/2069
    - Kim, D., and HS. Uh (2009). EMD: A Package for Empirical Mode Decomposition and Hilbert Spectrum. https://journal.r-project.org/archive/2009-1/RJournal_2009-1_Kim+Oh.pdf
    - Lambert et al. Empirical Mode Decomposition. https://www.clear.rice.edu/elec301/Projects02/empiricalMode/
    - Meta Trader. Introduction to the Empirical Mode Decomposition Method. https://www.mql5.com/en/articles/439
    - Salisbury, J.I. and Wimbush, M. (2002). Using modern time series analysis techniques to predict ENSO events from the SOI time series. Nonlinear Processes in Geophysics, 9, 341-345. http://www.nonlin-processes-geophys.net/9/341/2002/npg-9-341-2002.pdf
    - Torres, M. E., Colominas, M. A., Schlotthauer, G., & Flandrin, P. (2011). A complete ensemble empirical mode decomposition with adaptive noise. ICASSP, 4144-4147. http://doi.org/10.1109/ICASSP.2011.5947265.
    - Wang, T., M. Zhang, Q. Yu, and H. Zhang (2012). Comparing the applications of EMD and EEMD on time-frequency analysis of seismic signal. J. Appl. Geophys., 83, 29-34. http://doi.org/10.1016/j.jappgeo.2012.05.002.
    - Wu, Z., & Huang, N. E. (2009). Ensemble empirical mode decomposition: a noise-assisted data analysis method. Advances in Adaptive Data Analysis, 01(01), 1-41. https://doi.org/10.1142/S1793536909000047
    - Wu, Z, et al (2015). Fast multidimensional ensemble empirical mode decomposition for the analysis of big spatio-temporal datasets. Philos Trans A Math Phys Eng Sci, 374(2065), 20150197. http://doi.org/10.1098/rsta.2015.0197. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4792406/
    - Wu, Y. and Shen, BW (2016). An Evaluation of the Parallel Ensemble Empirical Mode Decomposition Method in Revealing the Role of Downscaling Processes Associated with African Easterly Waves in Tropical Cyclone Genesis. http://doi.org/10.1175/JTECH-D-15-0257.1. http://journals.ametsoc.org/doi/abs/10.1175/JTECH-D-15-0257.1

    .. seealso::

        - https://pyemd.readthedocs.io/
        - https://www.ncl.ucar.edu/Applications/eemd.shtml

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_emd.py
    """
    if time_array is None:
        t = datetime_to_numeric(input_data[time_dim].data, unit=time_step)
        t1 = input_data[time_dim].data
    else:
        t = time_array
        t1 = t

    s = input_data.data

    IMF = EMD(
        spline_kind=spline_kind,
        nbsym=nbsym,
        MAX_ITERATION=max_iteration,
        energy_ratio_thr=energy_ratio_thr,
        std_thr=std_thr,
        svar_thr=svar_thr,
        total_power_thr=total_power_thr,
        range_thr=range_thr,
        extrema_detection=extrema_detection,
        DTYPE=dtype,
    ).emd(s, t, max_imf=max_imf)
    emd_size = IMF.shape[0]

    result = xr.Dataset()
    result["input"] = input_data
    for iter in np.arange(0, emd_size):
        result["imf" + str(iter)] = xr.DataArray(
            IMF[iter], dims=time_dim, coords={time_dim: t1}
        )

    return result


def filter_eemd(
    input_data: xr.DataArray,
    time_step,
    time_array=None,
    time_dim: str = "time",
    noise_seed: None | int = None,
    trials: int = 100,
    noise_width: float = 0.05,
    parallel: bool = False,
    processes: None | int = None,
    separate_trends: bool = False,
    spline_kind: Literal[
        "cubic", "akima", "pchip", "cubic_hermite", "slinear", "quadratic", "linear"
    ] = "cubic",
    nbsym: int = 2,
    max_iteration: int = 1000,
    energy_ratio_thr: float = 0.2,
    std_thr: float = 0.2,
    svar_thr: float = 0.001,
    total_power_thr: float = 0.005,
    range_thr: float = 0.001,
    extrema_detection: Literal["simple", "parabol"] = "parabol",
    dtype=np.float64,
):
    """
    Ensemble Empirical Mode Decomposition (EEMD)

    Ensemble empirical mode decomposition (EEMD) is noise-assisted technique (Wu & Huang, 2009),
    which is meant to be more robust than simple Empirical Mode Decomposition (EMD).
    The robustness is checked by performing many decompositions on signals slightly perturbed from their initial position.
    In the grand average over all IMF results the noise will cancel each other out and the result is pure decomposition.

    Parameters
    ----------
    input_data: :py:class:`xarray.DataArray <xarray.DataArray>`
        Input signal data to be decomposed.
    time_step: :py:class:`str <str>`
        Time step unit for datetime conversion (e.g., 's', 'ms', 'us', 'ns').
    time_array: :py:class:`array-like <array-like>`, optional
        Custom time array for the input signal. If None, uses input_data's time dimension.
    time_dim: :py:class:`str <str>`, default: "time"
        The time coordinate dimension name.
    noise_seed: :py:class:`int <int>` or None, default: None
        Set seed for noise generation.

        .. warning::

            Given the nature of EEMD, each time you decompose a signal you will obtain a different set of components. That’s the expected consequence of adding noise which is going to be random. To make the decomposition reproducible, one needs to set a seed for the random number generator used in EEMD.

    trials: :py:class:`int <int>`, default: 100
        Number of trials or EMD performance with added noise.
    noise_width: :py:class:`float <float>`, default: 0.05
        Standard deviation of Gaussian noise (:math:`\\hat\\sigma`).
        It's relative to absolute amplitude of the signal, i.e.
        :math:`\\hat\\sigma = \\sigma\\cdot|\\max(S)-\\min(S)|`, where
        :math:`\\sigma` is noise_width.
    parallel: :py:class:`bool <bool>`, default: False
        Flag whether to use multiprocessing in EEMD execution. Since each EMD(s+noise) is independent this should improve execution speed considerably. Note that it’s disabled by default because it’s the most common problem when EEMD takes too long time to finish. If you set the flag to True, make also sure to set processes to some reasonable value.
    processes: :py:class:`int <int>` or None, default: None
        Number of processes harness when executing in parallel mode. The value should be between 1 and max that depends on your hardware. If None, uses all available cores.
    separate_trends: :py:class:`bool <bool>`, default: False
        Flag whether to isolate trends from each EMD decomposition into a separate component. If ``True``, the resulting EEMD will contain ensemble only from IMFs and the mean residue will be stacked as the last element.
    spline_kind: ``Literal["cubic", "akima", "pchip", "cubic_hermite", "slinear", "quadratic", "linear"]``, default: "cubic"
        Type of spline used for envelope interpolation.
    nbsym: :py:class:`int <int>`, default: 2
        Number of points to add at signal boundaries for mirroring.
    max_iteration: :py:class:`int <int>`, default: 1000
        Maximum number of iterations per single sifting in EMD.
    energy_ratio_thr: :py:class:`float <float>`, default: 0.2
        Threshold value on energy ratio per IMF check.
    std_thr: :py:class:`float <float>`, default: 0.2
        Threshold value on standard deviation per IMF check.
    svar_thr: :py:class:`float <float>`, default: 0.001
        Threshold value on scaled variance per IMF check.
    total_power_thr: :py:class:`float <float>`, default: 0.005
        Threshold value on total power per EMD decomposition.
    range_thr: :py:class:`float <float>`, default: 0.001
        Threshold for amplitude range (after scaling) per EMD decomposition.
    extrema_detection: ``Literal["simple", "parabol"]``, default: "parabol"
        Method used to finding extrema.
    dtype: :py:class:`numpy.dtype <numpy.dtype>`, default: np.float64
        Data type used for calculations.

    Returns
    -------
    :py:class:`xarray.Dataset <xarray.Dataset>`
        Dataset containing the input data and the decomposed eIMFs (Ensemble IMFs).

    Reference
    --------------
    Wu, Z., & Huang, N. E. (2009). Ensemble empirical mode decomposition: a noise-assisted data analysis method. Advances in Adaptive Data Analysis, 01(01), 1-41. https://doi.org/10.1142/S1793536909000047

    .. seealso::
        - :func:`filter_emd` : Standard Empirical Mode Decomposition
        - https://pyemd.readthedocs.io/
        - https://www.ncl.ucar.edu/Applications/eemd.shtml

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_emd.py
    """
    if time_array is None:
        t = datetime_to_numeric(input_data[time_dim].data, unit=time_step)
        t1 = input_data[time_dim].data
    else:
        t = time_array
        t1 = t

    s = input_data.data

    # Assign EEMD to `eemd` variable
    eemd = EEMD(
        trials=trials,
        noise_width=noise_width,
        parallel=parallel,
        processes=processes,
        separate_trends=separate_trends,
    )

    if noise_seed == None:
        pass
    elif isinstance(noise_seed, int):
        eemd.noise_seed(seed=noise_seed)
    else:
        raise ValueError(
            f"The value of noise_seed ({noise_seed}) should be `None` or `int`."
        )

    # Say we want detect extrema using parabolic method
    emd = eemd.EMD
    emd.spline_kind = spline_kind
    emd.nbsym = nbsym
    emd.MAX_ITERATION = max_iteration
    emd.energy_ratio_thr = energy_ratio_thr
    emd.std_thr = std_thr
    emd.svar_thr = svar_thr
    emd.total_power_thr = total_power_thr
    emd.range_thr = range_thr
    emd.extrema_detection = extrema_detection
    emd.DTYPE = dtype

    # Execute EEMD on s
    eIMFs = eemd.eemd(s, t)
    nIMFs = eIMFs.shape[0]

    result = xr.Dataset()
    result["input"] = input_data
    for iter in np.arange(0, nIMFs):
        result["eimf" + str(iter)] = xr.DataArray(
            eIMFs[iter], dims=time_dim, coords={time_dim: t1}
        )

    return result

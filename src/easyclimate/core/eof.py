"""
The analysis of the EOF and MCA

.. note::
    `xeofs`: https://xeofs.readthedocs.io/en/latest/
"""

import xarray as xr
import numpy as np
import xeofs
import warnings
from .variability import remove_seasonal_cycle_mean as remove_seasonal_cycle_mean_func
from .datanode import DataNode
from typing import Literal

import warnings

# -------------------------------------------------------------------
# EOF analysis

__all__ = [
    "get_EOF_model",
    "save_EOF_model",
    "load_EOF_model",
    "calc_EOF_analysis",
    "get_EOF_projection",
    "get_REOF_model",
    "save_REOF_model",
    "load_REOF_model",
    "calc_REOF_analysis",
    "get_REOF_projection",
    "get_MCA_model",
    "save_MCA_model",
    "load_MCA_model",
    "calc_MCA_analysis",
    "get_MCA_projection",
    "calc_eof_projection_coefficient",
]


def reformat_brace(element: xr.DataArray):
    """
    Solve `Invalid value for attr 'solver_kwargs': {}`
    """
    tmp = element.attrs["solver_kwargs"]
    element.attrs["solver_kwargs"] = f"{tmp}"
    return element


def get_EOF_model(
    data_input: xr.DataArray | list,
    lat_dim: str,
    lon_dim: str,
    time_dim: str = "time",
    n_modes: int = 10,
    remove_seasonal_cycle_mean=False,
    center: bool = False,
    standardize: bool = False,
    use_coslat: bool = True,
    random_state: int | None = None,
    solver: Literal["auto", "full", "randomized"] = "auto",
    solver_kwargs: dict = {},
) -> xeofs.single.eof.EOF:
    """
    Build the model of the Empirical Orthogonal Functions (EOF) analysis, more commonly known as Principal Component Analysis (PCA).

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`list <list>`
         The spatio-temporal data to be calculated.
    lat_dim: :py:class:`str <str>`.
        Latitude coordinate dimension name.
    lon_dim: :py:class:`str <str>`.
        Longitude coordinate dimension name.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    n_modes: :py:class:`int<int>`, default `10`.
        Number of modes to calculate.
    remove_seasonal_cycle_mean: :py:class:`bool<bool>`, default `False`.
        Whether to remove seasonal cycle mean of the input data.
        If it is `True`, the function will use :py:func:`easyclimate.remove_seasonal_cycle_mean<easyclimate.remove_seasonal_cycle_mean>` to remove seasonal cycle mean of the input data.
    center: :py:class:`bool<bool>`, default `False`.
        Whether to center the input data.
    standardize: :py:class:`bool<bool>`, default `False`.
        Whether to standardize the input data.
    use_coslat: :py:class:`bool<bool>`, default `True`.
        Whether to use cosine of latitude for scaling.
    random_state: :py:class:`int<int>`, default `None`.
        Seed for the random number generator.
    solver: {"auto", "full", "randomized"}, default: "auto".
        Solver to use for the SVD computation.
    solver_kwargs: :py:class:`dict<dict>`, default `{}`.
        Additional keyword arguments to be passed to the SVD solver.

    Returns
    -------
    :py:class:`xeofs.single.EOF<xeofs.single.EOF>`

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_multieof.py
        ./dynamic_docs/plot_basic_eof.py
    """
    from xeofs.single import EOF
    from xeofs.utils.constants import VALID_LONGITUDE_NAMES
    from xeofs.utils.constants import VALID_LATITUDE_NAMES

    if remove_seasonal_cycle_mean == True:
        if isinstance(data_input, xr.DataArray):
            data_input = remove_seasonal_cycle_mean_func(data_input, dim=time_dim)
        elif isinstance(data_input, list):
            tmp = list()
            for item in data_input:
                tmp_data = remove_seasonal_cycle_mean_func(item, dim=time_dim)
                tmp.append(tmp_data)
            data_input = tmp

    if (lat_dim in VALID_LATITUDE_NAMES) == False:
        warnings.warn(
            f"The lat_dim '{lat_dim}' in the data_input is not valid, we rename it to the 'lat'."
        )
        data_input = data_input.rename({lat_dim, "lat"})
    if (lon_dim in VALID_LONGITUDE_NAMES) == False:
        warnings.warn(
            f"The lon_dim '{lon_dim}' in the data_input is not valid, we rename it to the 'lon'."
        )
        data_input = data_input.rename({lon_dim, "lon"})

    model = EOF(
        n_modes=n_modes,
        center=center,
        standardize=standardize,
        use_coslat=use_coslat,
        random_state=random_state,
        solver=solver,
        solver_kwargs=solver_kwargs,
    )
    model.fit(data_input, dim=time_dim)

    return model


def save_EOF_model(
    model: xeofs.single.eof.EOF,
    path: str,
    overwrite: bool = False,
    save_data: bool = False,
    engine: Literal["zarr", "netcdf4", "h5netcdf"] = "zarr",
    **kwargs,
):
    """
    Save the model.

    Parameters
    ----------
    model: :py:class:`xeofs.single.EOF<xeofs.single.EOF>`
        The model of :py:class:`xeofs.single.EOF<xeofs.single.EOF>` is the results from :py:func:`easyclimate.eof.get_EOF_model <easyclimate.core.eof.get_EOF_model>` or :py:func:`xeofs.single.eof.EOF.fit <xeofs.single.eof.EOF.fit>`.
    path: :py:class:`str <str>`
        Path to save the model.
    overwrite: :py:class:`bool <bool>`, default `False`
        Whether or not to overwrite the existing path if it already exists. Ignored unless `engine = "zarr"`.
    save_data: :py:class:`bool <bool>`, default `False`
        Whether or not to save the full input data along with the fitted components.
    engine: {"zarr", "netcdf4", "h5netcdf"}, default `"zarr"`
        Xarray backend engine to use for writing the saved model.
    **kwargs: :py:class:`dict <dict>`.
        Additional keyword arguments to pass to `xarray.DataTree.to_netcdf()` or `xarray.DataTree.to_zarr()`.
    """
    model.save(
        path=path, overwrite=overwrite, save_data=save_data, engine=engine, **kwargs
    )


def load_EOF_model(
    path: str, engine: Literal["zarr", "netcdf4", "h5netcdf"] = "zarr", **kwargs
) -> xeofs.single.eof.EOF:
    """
    Load a saved EOF model.

    Parameters
    ----------
    path: :py:class:`str <str>`
        Path to the saved model.
    engine: {"zarr", "netcdf4", "h5netcdf"}, default `"zarr"`
        Xarray backend engine to use for reading the saved model.
    **kwargs: :py:class:`dict <dict>`.
        Additional keyword arguments to pass to `open_datatree()`.

    Returns
    -------
    The model of :py:class:`xeofs.single.EOF<xeofs.single.EOF>` is the results from :py:func:`easyclimate.eof.get_EOF_model <easyclimate.core.eof.get_EOF_model>` or :py:func:`xeofs.single.eof.EOF.fit <xeofs.single.eof.EOF.fit>`.
    """
    return xeofs.single.eof.EOF.load(path=path, engine=engine, **kwargs)


def calc_EOF_analysis(
    model: xeofs.single.eof.EOF, PC_normalized: bool = True
) -> xr.Dataset:
    """
    Calculate the results of the EOF model.

    Parameters
    ----------
    model: :py:class:`xeofs.single.EOF<xeofs.single.EOF>`
        The model of :py:class:`xeofs.single.EOF<xeofs.single.EOF>` is the results from :py:func:`easyclimate.eof.get_EOF_model <easyclimate.core.eof.get_EOF_model>` or :py:func:`xeofs.single.eof.EOF.fit <xeofs.single.eof.EOF.fit>`.
    PC_normalized: :py:class:`bool`, default `True`.
        Whether to normalize the scores by the L2 norm (singular values).
    Returns
    -------
    The results of the EOF model :py:class:`xarray.Dataset<xarray.Dataset>`.

    - **EOF: The (EOF) components**: The components in EOF anaylsis are the eigenvectors of the covariance/correlation matrix. Other names include the principal components or EOFs.
    - **PC: The (PC) scores**: The scores in EOF anaylsis are the projection of the data matrix onto the eigenvectors of the covariance matrix (or correlation) matrix. Other names include the principal component (PC) scores or just PCs.
    - **explained_variance**: The explained variance. The explained variance :math:`\\lambda_i` is the variance explained by each mode. It is defined as

    .. math::

        \\lambda_i = \\frac{\\sigma_i^2}{N-1}


    where :math:`\\sigma_i` is the singular value of the :math:`i`-th mode and :math:`N` is the number of samples. Equivalently, :math:`\\lambda_i` is the :math:`i`-th eigenvalue of the covariance matrix.

    - **explained_variance_ratio**: The explained variance ratio. The explained variance ratio :math:`\\gamma_i` is the variance explained by each mode normalized by the total variance. It is defined as

    .. math::

        \\gamma_i = \\frac{\\lambda_i}{\\sum_{j=1}^M \\lambda_j}


    where :math:`\\lambda_i` is the explained variance of the :math:`i`-th mode and :math:`M` is the total number of modes.

    - **singular_values**: The singular values of the Singular Value Decomposition (SVD).

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_multieof.py
        ./dynamic_docs/plot_basic_eof.py
    """
    if isinstance(model.components(), xr.DataArray or xr.Dataset):
        model_output = xr.Dataset()
        # Solve `Invalid value for attr 'solver_kwargs': {}`
        model_output["EOF"] = reformat_brace(model.components())
        model_output["PC"] = reformat_brace(model.scores(normalized=PC_normalized))
        model_output["explained_variance"] = reformat_brace(model.explained_variance())
        model_output["explained_variance_ratio"] = reformat_brace(
            model.explained_variance_ratio()
        )
        model_output["singular_values"] = reformat_brace(model.singular_values())
    elif isinstance(model.components(), list):
        model_output = DataNode(name="root")

        components_length = len(model.components())

        for num in np.arange(components_length):
            model_output[f"EOF/var{num}"] = model.components()[num]

        model_output["PC"] = model.scores(normalized=PC_normalized)
        model_output["explained_variance"] = model.explained_variance()
        model_output["explained_variance_ratio"] = model.explained_variance_ratio()
        model_output["singular_values"] = model.singular_values()
    return model_output


def get_EOF_projection(
    model: xeofs.single.eof.EOF,
    data: xr.DataArray,
    normalized: bool = True,
):
    """
    Project data onto the components.

    Parameters
    ----------
    model: :py:class:`xeofs.single.EOF<xeofs.single.EOF>`
        The model of :py:class:`xeofs.single.EOF<xeofs.single.EOF>` is the results from :py:func:`easyclimate.eof.get_EOF_model <easyclimate.core.eof.get_EOF_model>` or :py:func:`xeofs.single.eof.EOF.fit <xeofs.single.eof.EOF.fit>`.
    data: :py:class:`xarray.DataArray<xarray.DataArray>`
        Data to be transformed.
    normalized: :py:class:`bool<bool>`, default `True`.
        Whether to normalize the scores by the L2 norm.

    Returns
    -------
    projections: :py:class:`xarray.DataArray<xarray.DataArray>`
        Projections of the data onto the components.
    """
    return model.transform(data, normalized=normalized)


# -------------------------------------------------------------------
# Rotate EOF analysis


def get_REOF_model(
    data_input: xr.DataArray,
    lat_dim: str,
    lon_dim: str,
    time_dim: str = "time",
    n_modes: int = 2,
    power: int = 1,
    max_iter: int = None,
    rtol: float = 1e-8,
    remove_seasonal_cycle_mean=False,
    standardize: bool = False,
    use_coslat: bool = True,
    random_state=None,
    solver: Literal["auto", "full", "randomized"] = "auto",
    solver_kwargs={},
) -> xeofs.single.EOFRotator:
    """
    Build the model of the Rotate Empirical Orthogonal Functions (REOF) analysis.

    Parameters
    ----------
    data_input: :py:class:`xarray.DataArray<xarray.DataArray>`
         The spatio-temporal data to be calculated.
    lat_dim: :py:class:`str <str>`.
        Latitude coordinate dimension name.
    lon_dim: :py:class:`str <str>`.
        Longitude coordinate dimension name.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    n_modes: :py:class:`int<int>`, default `10`.
        Number of modes to calculate.
    remove_seasonal_cycle_mean: :py:class:`bool<bool>`, default `False`.
        Whether to remove seasonal cycle mean of the input data.
        If it is `True`, the function will use :py:func:`easyclimate.remove_seasonal_cycle_mean<easyclimate.remove_seasonal_cycle_mean>` to remove seasonal cycle mean of the input data.
    standardize: :py:class:`bool<bool>`, default `False`.
        Whether to standardize the input data.
    use_coslat: :py:class:`bool<bool>`, default `True`.
        Whether to use cosine of latitude for scaling.
    random_state: :py:class:`int<int>`, default `None`.
        Seed for the random number generator.
    solver: {"auto", "full", "randomized"}, default: "auto".
        Solver to use for the SVD computation.
    solver_kwargs: :py:class:`dict<dict>`, default `{}`.
        Additional keyword arguments to be passed to the SVD solver.

    Returns
    -------
    :py:class:`xeofs.single.EOFRotator<xeofs.single.EOFRotator>`

    Reference
    --------------
    Richman, M.B. (1986), Rotation of principal components. J. Climatol., 6: 293-335. https://doi.org/10.1002/joc.3370060305
    """
    from xeofs.single import EOF
    from xeofs.single import EOFRotator
    from xeofs.utils.constants import VALID_LONGITUDE_NAMES
    from xeofs.utils.constants import VALID_LATITUDE_NAMES

    if remove_seasonal_cycle_mean == True:
        data_input = remove_seasonal_cycle_mean_func(data_input, dim=time_dim)

    if (lat_dim in VALID_LATITUDE_NAMES) == False:
        warnings.warn(
            f"The lat_dim '{lat_dim}' in the data_input is not valid, we rename it to the 'lat'."
        )
        data_input = data_input.rename({lat_dim, "lat"})
    if (lon_dim in VALID_LONGITUDE_NAMES) == False:
        warnings.warn(
            f"The lon_dim '{lon_dim}' in the data_input is not valid, we rename it to the 'lon'."
        )
        data_input = data_input.rename({lon_dim, "lon"})

    model = EOF(
        n_modes=n_modes,
        standardize=standardize,
        use_coslat=use_coslat,
        random_state=random_state,
        solver=solver,
        solver_kwargs=solver_kwargs,
    )
    model.fit(data_input, dim=time_dim)
    rotator = EOFRotator(
        n_modes=n_modes,
        power=power,
        max_iter=max_iter,
        rtol=rtol,
    )
    rotator.fit(model)
    return rotator


def save_REOF_model(
    model: xeofs.single.EOFRotator,
    path: str,
    overwrite: bool = False,
    save_data: bool = False,
    engine: Literal["zarr", "netcdf4", "h5netcdf"] = "zarr",
    **kwargs,
):
    """
    Save the model.

    Parameters
    ----------
    model: :py:class:`xeofs.single.EOFRotator <xeofs.single.EOFRotator>`
        The model of :py:class:`xeofs.single.EOFRotator <xeofs.single.EOFRotator>` is the results from :py:func:`easyclimate.eof.get_REOF_model <easyclimate.core.eof.get_REOF_model>` or :py:func:`xeofs.single.EOFRotator.fit <xeofs.single.EOFRotator.fit>`.
    path: :py:class:`str <str>`
        Path to save the model.
    overwrite: :py:class:`bool <bool>`, default `False`
        Whether or not to overwrite the existing path if it already exists. Ignored unless `engine = "zarr"`.
    save_data: :py:class:`bool <bool>`, default `False`
        Whether or not to save the full input data along with the fitted components.
    engine: {"zarr", "netcdf4", "h5netcdf"}, default `"zarr"`
        Xarray backend engine to use for writing the saved model.
    **kwargs: :py:class:`dict <dict>`.
        Additional keyword arguments to pass to `xarray.DataTree.to_netcdf()` or `xarray.DataTree.to_zarr()`.
    """
    model.save(
        path=path, overwrite=overwrite, save_data=save_data, engine=engine, **kwargs
    )


def load_REOF_model(
    path: str, engine: Literal["zarr", "netcdf4", "h5netcdf"] = "zarr", **kwargs
) -> xeofs.single.EOFRotator:
    """
    Load a saved REOF model.

    Parameters
    ----------
    path: :py:class:`str <str>`
        Path to the saved model.
    engine: {"zarr", "netcdf4", "h5netcdf"}, default `"zarr"`
        Xarray backend engine to use for reading the saved model.
    **kwargs: :py:class:`dict <dict>`.
        Additional keyword arguments to pass to `open_datatree()`.

    Returns
    -------
    The model of :py:class:`xeofs.single.EOFRotator <xeofs.single.EOFRotator>` is the results from :py:func:`easyclimate.eof.get_REOF_model <easyclimate.core.eof.get_REOF_model>` or :py:func:`xeofs.single.EOFRotator.fit <xeofs.single.EOFRotator.fit>`.
    """
    return xeofs.single.EOFRotator.load(path=path, engine=engine, **kwargs)


def calc_REOF_analysis(
    model: xeofs.single.EOFRotator, PC_normalized: bool = True
) -> xr.Dataset:
    """
    Calculate the results of the REOF model.

    Parameters
    ----------
    model: :py:class:`xeofs.single.EOFRotator <xeofs.single.EOFRotator>`
        The model of :py:class:`xeofs.single.EOFRotator <xeofs.single.EOFRotator>` is the results from :py:func:`easyclimate.eof.get_REOF_model <easyclimate.core.eof.get_REOF_model>` or :py:func:`xeofs.single.EOFRotator.fit <xeofs.single.EOFRotator.fit>`.
    PC_normalized: :py:class:`bool`, default `True`.
        Whether to normalize the scores by the L2 norm (singular values).

    Returns
    -------
    The results of the EOF model :py:class:`xarray.Dataset<xarray.Dataset>`.

    - **EOF: The (EOF) components**: The components in EOF anaylsis are the eigenvectors of the covariance/correlation matrix. Other names include the principal components or EOFs.
    - **PC: The (PC) scores**: The scores in EOF anaylsis are the projection of the data matrix onto the eigenvectors of the covariance matrix (or correlation) matrix. Other names include the principal component (PC) scores or just PCs.
    - **explained_variance**: The explained variance. The explained variance :math:`\\lambda_i` is the variance explained by each mode. It is defined as

    .. math::

        \\lambda_i = \\frac{\\sigma_i^2}{N-1}

    where :math:`\\sigma_i` is the singular value of the :math:`i`-th mode and :math:`N` is the number of samples. Equivalently, :math:`\\lambda_i` is the :math:`i`-th eigenvalue of the covariance matrix.

    - **explained_variance_ratio**: The explained variance ratio. The explained variance ratio :math:`\\gamma_i` is the variance explained by each mode normalized by the total variance. It is defined as

    .. math::

        \\gamma_i = \\frac{\\lambda_i}{\\sum_{j=1}^M \\lambda_j}

    where :math:`\\lambda_i` is the explained variance of the :math:`i`-th mode and :math:`M` is the total number of modes.


    - **singular_values**: The singular values of the Singular Value Decomposition (SVD).
    """
    model_output = xr.Dataset()
    model_output["EOF"] = model.components()
    model_output["PC"] = model.scores(normalized=PC_normalized)
    model_output["explained_variance"] = model.explained_variance()
    model_output["explained_variance_ratio"] = model.explained_variance_ratio()
    model_output["singular_values"] = model.singular_values()

    return model_output


def get_REOF_projection(
    model: xeofs.single.EOFRotator,
    data: xr.DataArray,
    normalized: bool = True,
):
    """
    Project data onto the components.

    Parameters
    ----------
    model: :py:class:`xeofs.single.EOFRotator <xeofs.single.EOFRotator>`
        The model of :py:class:`xeofs.single.EOFRotator <xeofs.single.EOFRotator>` is the results from :py:func:`easyclimate.eof.get_REOF_model <easyclimate.core.eof.get_REOF_model>` or :py:func:`xeofs.single.EOFRotator.fit <xeofs.single.EOFRotator.fit>`.
    data: :py:class:`xarray.DataArray<xarray.DataArray>`
        Data to be transformed.
    normalized: :py:class:`bool<bool>`, default `True`.
        Whether to normalize the scores by the L2 norm.

    Returns
    -------
    projections: :py:class:`xarray.DataArray<xarray.DataArray>`
        Projections of the data onto the components.
    """
    return model.transform(data, normalized=normalized)


# -------------------------------------------------------------------
# MCA analysis


def get_MCA_model(
    data_left: xr.DataArray,
    data_right: xr.DataArray,
    lat_dim: str,
    lon_dim: str,
    time_dim: str = "time",
    n_modes=10,
    standardize: bool = False,
    use_coslat: bool = False,
    n_pca_modes: int = "auto",
    weights_left: xr.DataArray = None,
    weights_right: xr.DataArray = None,
    random_state: int = None,
    solver: Literal["auto", "full", "randomized"] = "auto",
    solver_kwargs: dict = {},
) -> xeofs.cross.MCA:
    """
    Build the model of the Maximum Covariance Analyis (MCA). MCA is a statistical method that finds patterns of maximum covariance between two datasets.

    .. note::
        MCA is similar to Principal Component Analysis (PCA) and Canonical Correlation Analysis (CCA), but while PCA finds modes of maximum variance and CCA finds modes of maximum correlation, MCA finds modes of maximum covariance.

    Parameters
    ----------
    data_left: :py:class:`xarray.DataArray <xarray.DataArray>`
        Left input data.
    data_right: :py:class:`xarray.DataArray <xarray.DataArray>`
        Right input data.
    lat_dim: :py:class:`str <str>`.
        Latitude coordinate dimension name.
    lon_dim: :py:class:`str <str>`.
        Longitude coordinate dimension name.
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.
    n_modes: :py:class:`int <int>`, default `10`.
        Number of modes to calculate.
    standardize: :py:class:`bool <bool>`, default `False`.
        Whether to standardize the input data.
    use_coslat: :py:class:`bool <bool>`, default `True`.
        Whether to use cosine of latitude for scaling.
    n_pca_modes: :py:class:`int <int>`, default same as `n_modes`, i.e, 'auto'.
        The number of principal components to retain during the PCA preprocessing step applied to both data sets prior to executing MCA.
        If set to None, PCA preprocessing will be bypassed, and the MCA will be performed on the original datasets.
        Specifying an integer value greater than 0 for `n_pca_modes` will trigger the PCA preprocessing,
        retaining only the specified number of principal components. This reduction in dimensionality
        can be especially beneficial when dealing with high-dimensional data, where computing the
        cross-covariance matrix can become computationally intensive or in scenarios where multicollinearity is a concern.
    weights_left: :py:class:`xarray.DataArray <xarray.DataArray>`
        Weights to be applied to the left input data.
    weights_right: :py:class:`xarray.DataArray <xarray.DataArray>`
        Weights to be applied to the right input data.
    random_state: :py:class:`int<int>`, default `None`.
        Seed for the random number generator.
    solver: {"auto", "full", "randomized"}, default: "auto".
        Solver to use for the SVD computation.
    solver_kwargs: :py:class:`dict<dict>`, default `{}`.
        Additional keyword arguments to be passed to the SVD solver.

    Returns
    -------
    :py:class:`xeofs.cross.MCA <xeofs.cross.MCA>`

    Reference
    --------------
    - Bretherton, C. S., Smith, C., & Wallace, J. M. (1992). An Intercomparison of Methods for Finding Coupled Patterns in Climate Data. Journal of Climate, 5(6), 541-560. https://doi.org/10.1175/1520-0442(1992)005<0541:AIOMFF>2.0.CO;2
    - Cherry, S. (1996). Singular Value Decomposition Analysis and Canonical Correlation Analysis. Journal of Climate, 9(9), 2003-2009. https://doi.org/10.1175/1520-0442(1996)009<2003:SVDAAC>2.0.CO;2

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_basic_eof.py
    """
    from xeofs.cross import MCA
    from xeofs.utils.constants import VALID_LONGITUDE_NAMES
    from xeofs.utils.constants import VALID_LATITUDE_NAMES

    time_length_data_left = data_left[time_dim].shape[0]
    time_length_data_right = data_right[time_dim].shape[0]
    if time_length_data_left != time_length_data_right:
        raise ValueError(
            f"The time length of the data_left is {time_length_data_left}, but the time length of the data_right is {time_length_data_right}. It is not equal, please check the time length of the input data."
        )

    if (lat_dim in VALID_LATITUDE_NAMES) == False:
        warnings.warn(
            f"The lat_dim '{lat_dim}' in the data_input is not valid, we rename it to the 'lat'."
        )
        data_left = data_left.rename({lat_dim, "lat"})
        data_right = data_right.rename({lat_dim, "lat"})
    if (lon_dim in VALID_LONGITUDE_NAMES) == False:
        warnings.warn(
            f"The lon_dim '{lon_dim}' in the data_input is not valid, we rename it to the 'lon'."
        )
        data_left = data_left.rename({lon_dim, "lon"})
        data_right = data_right.rename({lon_dim, "lon"})

    if n_pca_modes == "auto":
        n_pca_modes = n_modes

    model = MCA(
        n_modes=n_modes,
        standardize=standardize,
        use_coslat=use_coslat,
        solver=solver,
        n_pca_modes=n_pca_modes,
        random_state=random_state,
        solver_kwargs=solver_kwargs,
    )
    model.fit(
        data_left,
        data_right,
        dim=time_dim,
        weights_X=weights_left,
        weights_Y=weights_right,
    )
    return model


def save_MCA_model(
    model: xeofs.cross.MCA,
    path: str,
    overwrite: bool = False,
    save_data: bool = False,
    engine: Literal["zarr", "netcdf4", "h5netcdf"] = "zarr",
    **kwargs,
):
    """
    Save the model.

    Parameters
    ----------
    model: :py:class:`xeofs.cross.MCA <xeofs.cross.MCA>`
        The model of :py:class:`xeofs.cross.MCA <xeofs.cross.MCA>` is the results from :py:func:`easyclimate.eof.get_MCA_model <easyclimate.core.eof.get_MCA_model>` or :py:func:`xeofs.cross.mca.MCA.fit <xeofs.cross.mca.MCA.fit>`.
    path: :py:class:`str <str>`
        Path to save the model.
    overwrite: :py:class:`bool <bool>`, default `False`
        Whether or not to overwrite the existing path if it already exists. Ignored unless `engine = "zarr"`.
    save_data: :py:class:`bool <bool>`, default `False`
        Whether or not to save the full input data along with the fitted components.
    engine: {"zarr", "netcdf4", "h5netcdf"}, default `"zarr"`
        Xarray backend engine to use for writing the saved model.
    **kwargs: :py:class:`dict <dict>`.
        Additional keyword arguments to pass to `xarray.DataTree.to_netcdf()` or `xarray.DataTree.to_zarr()`.
    """
    model.save(
        path=path, overwrite=overwrite, save_data=save_data, engine=engine, **kwargs
    )


def load_MCA_model(
    path: str, engine: Literal["zarr", "netcdf4", "h5netcdf"] = "zarr", **kwargs
) -> xeofs.cross.MCA:
    """
    Load a saved MCA model.

    Parameters
    ----------
    path: :py:class:`str <str>`
        Path to the saved model.
    engine: {"zarr", "netcdf4", "h5netcdf"}, default `"zarr"`
        Xarray backend engine to use for reading the saved model.
    **kwargs: :py:class:`dict <dict>`.
        Additional keyword arguments to pass to `open_datatree()`.

    Returns
    -------
    The model of :py:class:`xeofs.cross.MCA <xeofs.cross.EOFRotator>` is the results from :py:func:`easyclimate.eof.get_MCA_model <easyclimate.core.eof.get_MCA_model>` or :py:func:`xeofs.cross.mca.MCA.fit <xeofs.cross.mca.MCA.fit>`.
    """
    return xeofs.cross.MCA.load(path=path, engine=engine, **kwargs)


def calc_MCA_analysis(
    model: xeofs.cross.MCA, correction=None, alpha=0.05, PC_normalized: bool = True
) -> DataNode:
    """
    Calculate the results of the EOF model.

    Parameters
    ----------
    model: :py:class:`xeofs.cross.MCA <xeofs.cross.MCA>`
        The model of :py:class:`xeofs.cross.MCA <xeofs.cross.MCA>` is the results from :py:func:`easyclimate.eof.get_MCA_model <easyclimate.core.eof.get_MCA_model>` or :py:func:`xeofs.cross.mca.MCA.fit <xeofs.cross.mca.MCA.fit>`.
    correction: :py:class:`str <str>`, default `None`
        Method to apply a multiple testing correction. If None, no correction is applied. Available methods are:

        - bonferroni : one-step correction
        - sidak : one-step correction
        - holm-sidak : step down method using Sidak adjustments
        - holm : step-down method using Bonferroni adjustments
        - simes-hochberg : step-up method (independent)
        - hommel : closed method based on Simes tests (non-negative)
        - fdr_bh : Benjamini/Hochberg (non-negative) (default)
        - fdr_by : Benjamini/Yekutieli (negative)
        - fdr_tsbh : two stage fdr correction (non-negative)
        - fdr_tsbky : two stage fdr correction (non-negative)

    alpha: :py:class:`float <float>`, default `0.05`
        The desired family-wise error rate. Not used if correction is None.
    PC_normalized: :py:class:`bool`, default `True`.
        Whether to normalize the scores by the L2 norm (singular values).

    Returns
    -------
    The results of the MCA model (:py:class:`easyclimate.DataNode <easyclimate.core.datanode.DataNode>`).

    - **EOF**: The singular vectors of the left and right field.
    - **PC**: The scores of the left and right field. The scores in MCA are the projection of the left and right field onto the left and right singular vector of the cross-covariance matrix.
    - **correlation_coefficients_X**: Get the correlation coefficients for the scores of :math:`X`.

    The correlation coefficients of the scores of :math:`X` are given by:

    .. math::

        c_{x, ij} = \\text{corr} \\left(\\mathbf{r}_{x, i}, \\mathbf{r}_{x, j} \\right)

    where :math:`\\mathbf{r}_{x, i}` and :math:`\\mathbf{r}_{x, j}` are the :math:`i` th and :math:`j` th scores of :math:`X`.

    - **correlation_coefficients_Y**: Get the correlation coefficients for the scores of :math:`Y`.

    The correlation coefficients of the scores of :math:`Y` are given by:

    .. math::

        c_{y, ij} = \\text{corr} \\left(\\mathbf{r}_{y, i}, \\mathbf{r}_{y, j} \\right)

    where :math:`\\mathbf{r}_{y, i}` and :math:`\\mathbf{r}_{y, j}` are the :math:`i` th and :math:`j` th scores of :math:`Y`.
    - **covariance_fraction_CD95**: Get the covariance fraction (CF).

    Cheng and Dunkerton (1995) define the CF as follows:

    .. math::

        CF_i = \\frac{\\sigma_i}{\\sum_{i=1}^{m} \\sigma_i}

    where :math:`m` is the total number of modes and :math:`\\sigma_i` is the :math:`i`-th singular value of the covariance matrix.

    This implementation estimates the sum of singular values from the first n modes,
    therefore one should aim to retain as many modes as possible to get a good estimate of the covariance fraction.

    .. note::

        In MCA, the focus is on maximizing the squared covariance (SC).
        As a result, this quantity is preserved during decomposition - meaning the SC of both datasets
        remains unchanged before and after decomposition. Each mode explains a fraction of the total SC,
        and together, all modes can reconstruct the total SC of the cross-covariance matrix.
        However, the (non-squared) covariance is not invariant in MCA;
        it is not preserved by the individual modes and cannot be reconstructed from them.
        Consequently, the squared covariance fraction (SCF) is invariant in MCA and is typically
        used to assess the relative importance of each mode. In contrast, the convariance fraction (CF) is not invariant.
        Cheng and Dunkerton (1995) introduced the CF to compare the relative importance of modes
        before and after Varimax rotation in MCA. Notably, when the data fields in MCA are identical,
        the CF corresponds to the explained variance ratio in Principal Component Analysis (PCA).

    - **cross_correlation_coefficients**: Get the cross-correlation coefficients.

    The cross-correlation coefficients between the scores of :math:`X` and :math:`Y` are computed as:

    .. math::

        c_{xy, i} = \\text{corr} \\left(\\mathbf{r}_{x, i}, \\mathbf{r}_{y, i} \\right)

    where :math:`\\mathbf{r}_{x, i}` and :math:`\\mathbf{r}_{y, i}` are the :math:`i` th scores of :math:`X` and :math:`Y`.

    .. note::

        When :math:`\\alpha=0`, the cross-correlation coefficients are equivalent to the canonical correlation coefficients.

    - **fraction_variance_X_explained_by_X**: Get the fraction of variance explained (FVE X).

    The FVE X is the fraction of variance in :math:`X` explained by the scores of :math:`X`.

    It is computed as a weighted mean-square error (see equation (15) in Swenson (2015)) :

    .. math::

        FVE_{X|X,i} = 1 - \\frac{\\|\\mathbf{d}_{X,i}\\|_F^2}{\\|X\\|_F^2}

    where :math:`\\mathbf{d}_{X,i}` are the residuals of the input data :math:`X` after reconstruction by the :math:`i` th scores of :math:`X`.

    - **fraction_variance_Y_explained_by_X**: Get the fraction of variance explained (FVE YX).

    The FVE YX is the fraction of variance in :math:`Y` explained by the scores of :math:`X`.
    It is computed as a weighted mean-square error (see equation (15) in Swenson (2015)) :

    .. math::

        FVE_{Y|X,i} = 1 - \\frac{\\|(X^TX)^{-1/2} \\mathbf{d}_{X,i}^T \\mathbf{d}_{Y,i}\\|_F^2}{\\|(X^TX)^{-1/2} X^TY\\|_F^2}

    where :math:`\\mathbf{d}_{X,i}` and :math:`\\mathbf{d}_{Y,i}` are the residuals of the input data :math:`X`
    and :math:`Y` after reconstruction by the :math:`i` th scores of :math:`X` and :math:`Y`, respectively.

    - **fraction_variance_Y_explained_by_Y**: Get the fraction of variance explained (FVE Y).

    The FVE Y is the fraction of variance in :math:`Y` explained by the scores of :math:`Y`.
    It is computed as a weighted mean-square error (see equation (15) in Swenson (2015)) :

    .. math::

        FVE_{Y|Y,i} = 1 - \\frac{\\|\\mathbf{d}_{Y,i}\\|_F^2}{\\|Y\\|_F^2}

    where :math:`\\mathbf{d}_{Y,i}` are the residuals of the input data :math:`Y`
    after reconstruction by the :math:`i` th scores of :math:`Y`.

    - **squared_covariance_fraction**: Get the squared covariance fraction (SCF).

    The SCF is computed as a weighted mean-square error (see equation (15) in Swenson (2015)) :

    .. math::

        SCF_{i} = 1 - \\frac{\\|\\mathbf{d}_{X,i}^T \\mathbf{d}_{Y,i}\\|_F^2}{\\|X^TY\\|_F^2}

    where :math:`\\mathbf{d}_{X,i}` and :math:`\\mathbf{d}_{Y,i}` are the residuals of the input data :math:`X`
    and :math:`Y` after reconstruction by the :math:`i` th scores of :math:`X` and :math:`Y`, respectively.

    - **heterogeneous_patterns**: The heterogeneous patterns of the left and right field.

    The heterogeneous patterns are the correlation coefficients between the input data and the scores of the other field.

    More precisely, the heterogeneous patterns :math:`r_{\\mathrm{het}}` are defined as

    .. math::

        r_{\\mathrm{het}, x} = corr \\left(X, A_y \\right), \\ r_{\\mathrm{het}, y} = corr \\left(Y, A_x \\right)

    where :math:`X` and :math:`Y` are the input data, :math:`A_x` and :math:`A_y` are the scores of the left and right field, respectively.

    - **homogeneous_patterns**: The homogeneous patterns of the left and right field.

    The homogeneous patterns are the correlation coefficients between the input data and the scores.

    More precisely, the homogeneous patterns :math:`r_{\\mathrm{hom}}` are defined as

    .. math::

        r_{\\mathrm{hom}, x} = corr \\left(X, A_x \\right), \\ r_{\\mathrm{hom}, y} = corr \\left(Y, A_y \\right)

    where :math:`X` and :math:`Y` are the input data, :math:`A_x` and :math:`A_y` are the scores of the left and right field, respectively.

    Reference
    --------------
    - Cheng, X., & Dunkerton, T. J. (1995). Orthogonal Rotation of Spatial Patterns Derived from Singular Value Decomposition Analysis. Journal of Climate, 8(11), 2631-2643. https://doi.org/10.1175/1520-0442(1995)008<2631:OROSPD>2.0.CO;2
    - Swenson, E. (2015). Continuum Power CCA: A Unified Approach for Isolating Coupled Modes. Journal of Climate, 28(3), 1016-1030. https://doi.org/10.1175/JCLI-D-14-00451.1

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_basic_eof.py
    """
    model_output = DataNode(name="root")

    # components
    left_components = model.components()[0]
    right_components = model.components()[1]
    left_components.name = "left_EOF"
    right_components.name = "right_EOF"

    model_output["EOF/left_EOF"] = left_components
    model_output["EOF/right_EOF"] = right_components

    # scores
    left_scores = model.scores(normalized=PC_normalized)[0]
    right_scores = model.scores(normalized=PC_normalized)[1]
    left_scores.name = "left_PC"
    right_scores.name = "right_PC"

    model_output["PC/left_PC"] = left_scores
    model_output["PC/right_PC"] = right_scores

    # correlation_coefficients_X
    correlation_coefficients_X = model.correlation_coefficients_X()
    model_output["correlation_coefficients_X"] = correlation_coefficients_X

    # correlation_coefficients_Y
    correlation_coefficients_Y = model.correlation_coefficients_Y()
    model_output["correlation_coefficients_Y"] = correlation_coefficients_Y

    # covariance_fraction_CD95
    covariance_fraction = model.covariance_fraction_CD95()
    model_output["covariance_fraction"] = covariance_fraction

    # cross_correlation_coefficients
    cross_correlation_coefficients = model.cross_correlation_coefficients()
    model_output["cross_correlation_coefficients"] = cross_correlation_coefficients

    # fraction_variance_X_explained_by_X
    fraction_variance_X_explained_by_X = model.fraction_variance_X_explained_by_X()
    model_output["fraction_variance_X_explained_by_X"] = (
        fraction_variance_X_explained_by_X
    )

    # fraction_variance_Y_explained_by_X
    fraction_variance_Y_explained_by_X = model.fraction_variance_Y_explained_by_X()
    model_output["fraction_variance_Y_explained_by_X"] = (
        fraction_variance_Y_explained_by_X
    )

    # fraction_variance_Y_explained_by_Y
    fraction_variance_Y_explained_by_Y = model.fraction_variance_Y_explained_by_Y()
    model_output["fraction_variance_Y_explained_by_Y"] = (
        fraction_variance_Y_explained_by_Y
    )

    # squared_covariance_fraction
    squared_covariance_fraction = model.squared_covariance_fraction()
    model_output["squared_covariance_fraction"] = squared_covariance_fraction

    # heterogeneous_patterns
    tmp = model.heterogeneous_patterns(correction=correction, alpha=alpha)
    left_heterogeneous_patterns = tmp[0][0]
    right_heterogeneous_patterns = tmp[0][1]
    pvalues_of_left_heterogeneous_patterns = tmp[1][0]
    pvalues_of_right_heterogeneous_patterns = tmp[1][1]

    model_output["heterogeneous_patterns/left_heterogeneous_patterns"] = (
        left_heterogeneous_patterns
    )
    model_output["heterogeneous_patterns/right_heterogeneous_patterns"] = (
        right_heterogeneous_patterns
    )
    model_output["heterogeneous_patterns/pvalues_of_left_heterogeneous_patterns"] = (
        pvalues_of_left_heterogeneous_patterns
    )
    model_output["heterogeneous_patterns/pvalues_of_right_heterogeneous_patterns"] = (
        pvalues_of_right_heterogeneous_patterns
    )

    # homogeneous_patterns
    tmp = model.homogeneous_patterns(correction=correction, alpha=alpha)
    left_homogeneous_patterns = tmp[0][0]
    right_homogeneous_patterns = tmp[0][1]
    pvalues_of_left_homogeneous_patterns = tmp[1][0]
    pvalues_of_right_homogeneous_patterns = tmp[1][1]

    model_output["homogeneous_patterns/left_homogeneous_patterns"] = (
        left_homogeneous_patterns
    )
    model_output["homogeneous_patterns/right_homogeneous_patterns"] = (
        right_homogeneous_patterns
    )
    model_output["homogeneous_patterns/pvalues_of_left_homogeneous_patterns"] = (
        pvalues_of_left_homogeneous_patterns
    )
    model_output["homogeneous_patterns/pvalues_of_right_homogeneous_patterns"] = (
        pvalues_of_right_homogeneous_patterns
    )

    return model_output


def get_MCA_projection(
    model: xeofs.cross.mca.MCA,
    data_left: xr.DataArray | xr.Dataset,
    data_right: xr.DataArray | xr.Dataset,
    normalized: bool = True,
) -> DataNode:
    """
    Get the expansion coefficients of "unseen" data. The expansion coefficients are obtained by projecting data onto the singular vectors.

    Parameters
    ----------
    model: :py:class:`xeofs.cross.MCA <xeofs.cross.MCA>`
        The model of :py:class:`xeofs.cross.MCA <xeofs.cross.MCA>` is the results from :py:func:`easyclimate.eof.get_MCA_model <easyclimate.core.eof.get_MCA_model>` or :py:meth:`xeofs.cross.mca.MCA.fit <xeofs.cross.mca.MCA.fit>`.
    data_left: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        Left input data. Must be provided if ``data_right`` is not provided.
    data_right: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        Right input data. Must be provided if ``data_left`` is not provided.
    normalized: :py:class:`bool`, default `False`.
        Whether to return L2 normalized scores.

    Returns
    -------
    scores: :py:class:`easyclimate.DataNode <easyclimate.core.datanode.DataNode>`

    - **scores1**: Left scores.
    - **scores2**: Right scores.
    """
    scores1, scores2 = model.transform(
        data1=data_left, data2=data_right, normalized=normalized
    )
    scores1.name = "scores1"
    scores2.name = "scores2"

    projection_output = DataNode(name="root")
    projection_output["scores1"] = scores1
    projection_output["scores2"] = scores2
    return projection_output


def calc_eof_projection_coefficient(
    data_field: xr.DataArray,
    eof_mode: xr.DataArray,
    time_dim: str = "time",
):
    """
    Compute the EOF projection coefficients for projecting a data field onto an EOF mode over common
    spatial dimensions. This is useful in EOF (Empirical Orthogonal Function) analysis for climate or
    geophysical data, where the field is projected onto spatial modes to obtain time-varying coefficients.

    The mathematical foundation is based on the decomposition (:math:`\\mathbf{X} = \\mathbf{V} \\mathbf{T}`),
    solving for the coefficients (:math:`\\mathbf{T} = \\mathbf{V}^{-1} \\mathbf{X}`). For a single normalized
    EOF mode (:math:`\\mathbf{V}`), this simplifies to the projection:

    .. math::

        t = \\frac{\\sum (x \\cdot v)}{\\sum v^2}

    where the summation is over the stacked spatial (pattern) dimensions, and (:math:`\\mathbf{X}`) is the
    data field (potentially with a 'time' dimension), (:math:`\\mathbf{V}`) is the EOF mode (spatial pattern).

    The spatial pattern dimensions are automatically detected as the intersection of the input dimensions,
    excluding 'time' (if present). Both inputs are stacked along these pattern dimensions into a temporary
    'pattern' dimension, and the projection is computed along it. NaN values are filled with 0 before computation.

    - If `data_field` lacks 'time', the result is a scalar.
    - If `data_field` has 'time' and `eof_mode` does not, the result preserves the 'time' dimension.
    - Broadcasting occurs automatically for compatible shapes.

    Parameters
    ----------
    data_field : :py:class:`xarray.DataArray`
        The input data field to project (e.g., time series of spatial fields :math:`\\mathbf{X}`).
    eof_mode : :py:class:`xarray.DataArray`
        The EOF spatial mode :math:`\\mathbf{V}` (must have compatible spatial dimensions).
    time_dim: :py:class:`str <str>`, default: `time`.
        The time coordinate dimension name.

    Returns
    -------
    coefficients : :py:class:`xarray.DataArray` or scalar
        The EOF projection coefficients (:math:`\\mathbf{T}`). Dimensions match the non-spatial dimensions of `data_field`
        (e.g., 'time' if present).

    .. note::
        - Assumes inputs have compatible shapes and the only differing dimension is 'time' in `data_field`.
        - NaNs are filled with 0 to avoid propagation; adjust if needed.
        - For zero-norm cases in the denominator, the result is set to 0.

    Examples
    --------
    Scalar projection for a single spatial field:

    >>> import xarray as xr
    >>> import numpy as np
    >>> import easyclimate as ecl
    >>> # Create a random number generator with a fixed seed.
    >>> rng = np.random.default_rng(42)
    >>> field = xr.DataArray(rng.random((2, 3)), dims=['lat', 'lon'])
    >>> eof_v = xr.DataArray(rng.random((2, 3)), dims=['lat', 'lon'])
    >>> coeff = ecl.eof.calc_eof_projection_coefficient(field, eof_v)
    >>> print(coeff)
    <xarray.DataArray 'eof_projection' ()> Size: 8B
    array(0.95208032)
    Attributes:
        long_name:  EOF Projection Coefficient
        units:

    Time series projection:

    >>> # Create a random number generator with a fixed seed.
    >>> rng = np.random.default_rng(42)
    >>> time = xr.DataArray(np.arange(4), dims=['time'])
    >>> timed_field = xr.DataArray(rng.random((4, 2, 3)), dims=['time', 'lat', 'lon'])
    >>> coeff_time = calc_eof_projection_coefficient(timed_field, eof_v)
    >>> print(coeff_time)
    <xarray.DataArray 'eof_projection' (time: 4)> Size: 32B
    array([0.95208032, 1.        , 0.64684219, 1.06549741])
    Dimensions without coordinates: time
    Attributes:
        long_name:  EOF Projection Coefficient
        units:
    """
    # Ensure inputs are DataArray
    if not isinstance(data_field, xr.DataArray) or not isinstance(
        eof_mode, xr.DataArray
    ):
        raise TypeError("Inputs must be xarray.DataArray objects")

    # Detect common spatial dimensions (exclude 'time')
    dims_field = set(data_field.dims)
    dims_mode = set(eof_mode.dims)
    common_dims = dims_field.intersection(dims_mode)
    pattern_dims_set = common_dims - {time_dim}
    if not pattern_dims_set:
        raise ValueError("No common spatial (non-time) dimensions found between inputs")
    pattern_dims = sorted(pattern_dims_set)  # Sorted for consistent stacking

    # Fill NaNs with 0
    data_field_filled = data_field.fillna(0)
    eof_mode_filled = eof_mode.fillna(0)

    # Stack both inputs along pattern dimensions
    field_stacked = data_field_filled.stack(pattern=tuple(pattern_dims))
    mode_stacked = eof_mode_filled.stack(pattern=tuple(pattern_dims))

    # Compute projection: dot product along pattern / norm^2 of mode
    numerator = (field_stacked * mode_stacked).sum(dim="pattern")
    denominator = (mode_stacked**2).sum()

    # Avoid division by zero
    coefficients = xr.where(denominator != 0, numerator / denominator, 0)

    # Add necessary name and attributes
    coefficients.name = "eof_projection_coefficient"
    coefficients.attrs = {
        "long_name": "EOF Projection Coefficient",
        "units": data_field.attrs.get(
            "units", ""
        ),  # Inherit units from data_field if available
    }

    return coefficients

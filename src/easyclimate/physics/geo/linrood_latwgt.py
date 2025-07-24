"""
Latitudes and weights in the Lin-Rood model

This module provides functions to calculate the latitudes and weights used in the Lin-Rood model,
which is commonly employed in atmospheric modeling and semi-Lagrangian transport schemes.
"""

import numpy as np
from numpy import ndarray
from typing import Tuple

__all__ = ["calc_lat_weight_lin_rood"]


def calc_lat_weight_lin_rood(nlat: int) -> Tuple[ndarray, ndarray]:
    """
    Calculate the latitudes and weights used by the Lin-Rood model.

    The Lin-Rood model requires a specific distribution of latitudes and corresponding weights
    for numerical integration on a spherical grid. This function generates these values based
    on the number of desired latitudes.

    Parameters
    ----------
    nlat : :py:class:`int <int>`
        Number of latitudes. Must be at least 2 to define a valid grid (from pole to pole).

    Returns
    -------
    Tuple[ndarray, ndarray]
        A tuple containing two numpy arrays:
        - lat : ndarray
            Array of latitudes in degrees, ranging from -90 (South Pole) to 90 (North Pole).
        - weight : ndarray
            Array of weights corresponding to each latitude, used for numerical integration.

    .. tip::

        The weights are computed such that they are suitable for use in the Lin-Rood semi-Lagrangian
        transport scheme. The latitudes are uniformly spaced between the poles.

    References
    ----------
    - Lin, S., & Rood, R. B. (1996). Multidimensional Flux-Form Semi-Lagrangian Transport Schemes. Monthly Weather Review, 124(9), 2046-2070. https://journals.ametsoc.org/view/journals/mwre/124/9/1520-0493_1996_124_2046_mffslt_2_0_co_2.xml
    - Lin, S.-J. and Rood, R.B. (1997), An explicit flux-form semi-lagrangian shallow-water model on the sphere. Q.J.R. Meteorol. Soc., 123: 2477-2498. https://doi.org/10.1002/qj.49712354416

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Built-in/linrood_latwgt.shtml
    """
    if nlat < 2:
        raise ValueError(
            "nlat must be at least 2 to define a valid grid from pole to pole."
        )

    def linroodwt(nlat: int) -> ndarray:
        """
        Calculate Lin-Rood weights for numerical integration.

        Parameters
        ----------
        nlat : int
            Number of latitudes.

        Returns
        -------
        ndarray
            Array of weights, with values derived from the sine of latitude intervals.
        """
        weight = np.zeros(nlat)
        sine, cosp, sinp, cose = setrig(nlat)

        # Compute weights for intermediate latitudes
        for nl in range(1, nlat - 1):
            weight[nl] = sine[nl + 1] - sine[nl]

        # Handle poles
        weight[0] = 1.0 + sine[1]  # South Pole weight
        weight[-1] = 1.0 - sine[-1]  # North Pole weight

        return weight

    def setrig(jm: int) -> Tuple[ndarray, ndarray, ndarray, ndarray]:
        """
        Setup trigonometric quantities needed for weight calculations.

        Parameters
        ----------
        jm : int
            Number of latitudes.

        Returns
        -------
        Tuple[ndarray, ndarray, ndarray, ndarray]
            A tuple containing four numpy arrays:
            - sine : Sine values at latitude midpoints.
            - cosp : Cosine values at latitude points.
            - sinp : Sine values at latitude points.
            - cose : Cosine values at latitude edges.
        """
        sine = np.zeros(jm)
        cosp = np.zeros(jm)
        sinp = np.zeros(jm)
        cose = np.zeros(jm)

        jm1 = jm - 1
        pi = np.pi

        # Calculate sine of latitude midpoints
        for j in range(1, jm):
            ph5 = -0.5 * pi + (j - 0.5) * (pi / jm1)
            sine[j] = np.sin(ph5)

        # Handle poles for cosp
        cosp[0] = 0.0  # South Pole
        cosp[-1] = 0.0  # North Pole

        # Calculate cosp for intermediate latitudes
        dp = pi / jm1
        for j in range(1, jm - 1):
            cosp[j] = (sine[j + 1] - sine[j]) / dp

        # Calculate cose (cosine at edges)
        for j in range(1, jm):
            cose[j] = 0.5 * (cosp[j - 1] + cosp[j])
        cose[0] = cose[1]  # South Pole edge

        # Handle poles for sinp
        sinp[0] = -1.0  # South Pole
        sinp[-1] = 1.0  # North Pole

        # Calculate sinp for intermediate latitudes
        for j in range(1, jm - 1):
            sinp[j] = 0.5 * (sine[j] + sine[j + 1])

        return sine, cosp, sinp, cose

    # Calculate latitudes (uniformly spaced from -90 to 90 degrees)
    lat = np.zeros(nlat)
    dlat = 180.0 / (nlat - 1)
    lat[0] = -90.0
    for nl in range(1, nlat - 1):
        lat[nl] = lat[nl - 1] + dlat
    lat[-1] = 90.0

    # Calculate weights
    weight = linroodwt(nlat)

    return lat, weight

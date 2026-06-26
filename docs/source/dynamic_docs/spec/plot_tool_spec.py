# -*- coding: utf-8 -*-
"""
Tools for Spherical Harmonics
==============================================


This page introduces low-level helpers for spherical-harmonic grids and
triangular spectral indexing. Each function is available in both the standard
backend and a Rust-backed ``*_rs`` variant where applicable.

"""
import easyclimate as ecl

# %%
# Gaussian latitudes and quadrature weights are commonly used by spherical-harmonic transforms on Gaussian grids. The input ``nlat`` is the number of latitude points.
ecl.spec.calc_gaussian_latitudes(72)
ecl.spec.calc_gaussian_latitudes_rs(72)

# %%
# Geodesic points provide nearly even point locations on the sphere for an icosahedral geodesic. The argument ``m`` controls the number of points along one geodesic-triangle edge.
ecl.spec.calc_geodesic_points(10)
ecl.spec.calc_geodesic_points_rs(10)

# %%
# Spherical-harmonic coefficient indices describe the packed triangular order used for complex spectral coefficients. ``ntrunc`` is the triangular truncation limit, such as T42.
ecl.spec.calc_spherical_harmonic_coefficients(42)
ecl.spec.calc_spherical_harmonic_coefficients_rs(42)

# %%
# Associated Legendre functions are evaluated at one latitude for the selected triangular truncation. They are the basis functions used internally by the transforms.
ecl.spec.calc_legendre_functions(45.0, 42)
ecl.spec.calc_legendre_functions_rs(45.0, 42)

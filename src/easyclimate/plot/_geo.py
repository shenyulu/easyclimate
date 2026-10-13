"""Geographic extent handling shared by native mesh plots."""

import cartopy.crs as ccrs
import numpy as np


def _central_longitude(crs, default=0.0):
    """Read the central meridian in both old and new Cartopy projections."""
    params = getattr(crs, "proj4_params", {})
    # Cartopy 0.26 switched PlateCarree from eqc/lon_0 to latlong/pm.
    return float(params.get("lon_0", params.get("pm", default)))


def _set_cartopy_extent(ax, lon_min, lon_max, lat_min, lat_max):
    """Keep continuous longitude windows on cylindrical map projections."""
    lon_min, lon_max = float(lon_min), float(lon_max)
    while lon_max <= lon_min:
        lon_max += 360.0
    lat_min = max(-90.0, float(lat_min))
    lat_max = min(90.0, float(lat_max))
    projection = ax.projection

    if isinstance(projection, (ccrs.PlateCarree, ccrs.Mercator)):
        span = lon_max - lon_min
        if span >= 360.0 - 1e-6:
            x_min, x_max = projection.x_limits
        else:
            center = _central_longitude(projection)
            start = (lon_min - center + 180.0) % 360.0 - 180.0
            end = start + span
            # A window crossing the projection seam needs both map edges.
            if end > 180.0 + 1e-6:
                x_min, x_max = projection.x_limits
            else:
                x_min = projection.transform_point(
                    center + start, 0, ccrs.PlateCarree()
                )[0]
                x_max = projection.transform_point(center + end, 0, ccrs.PlateCarree())[
                    0
                ]
                if np.isclose(start, -180.0, atol=1e-8, rtol=0):
                    x_min = projection.x_limits[0]
                if np.isclose(end, 180.0, atol=1e-8, rtol=0):
                    x_max = projection.x_limits[1]
        center = _central_longitude(projection)
        y_min = projection.transform_point(center, lat_min, ccrs.PlateCarree())[1]
        y_max = projection.transform_point(center, lat_max, ccrs.PlateCarree())[1]
        if lat_min == -90.0:
            y_min = projection.y_limits[0]
        if lat_max == 90.0:
            y_max = projection.y_limits[1]
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
    elif lon_max - lon_min >= 360.0 - 1e-6:
        if lat_min == -90 and lat_max == 90:
            ax.set_global()
        elif isinstance(projection, (ccrs.NorthPolarStereo, ccrs.SouthPolarStereo)):
            outer_lat = (
                lat_min if isinstance(projection, ccrs.NorthPolarStereo) else lat_max
            )
            x, y = projection.transform_point(0, outer_lat, ccrs.PlateCarree())
            radius = np.hypot(x, y)
            ax.set_xlim(
                max(-radius, projection.x_limits[0]),
                min(radius, projection.x_limits[1]),
            )
            ax.set_ylim(
                max(-radius, projection.y_limits[0]),
                min(radius, projection.y_limits[1]),
            )
        else:
            lons, lats = np.meshgrid(
                np.linspace(-180, 180, 721), np.linspace(lat_min, lat_max, 181)
            )
            points = projection.transform_points(ccrs.PlateCarree(), lons, lats)
            valid = np.all(np.isfinite(points[..., :2]), axis=-1)
            if not np.any(valid):
                raise ValueError(
                    "The geographic window is outside the projection domain."
                )
            xy = points[valid, :2]
            ax.set_xlim(np.min(xy[:, 0]), np.max(xy[:, 0]))
            ax.set_ylim(np.min(xy[:, 1]), np.max(xy[:, 1]))
    else:
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

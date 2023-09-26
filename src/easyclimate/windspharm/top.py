"""
top interface

https://ajdawson.github.io/windspharm/latest/api/windspharm.xarray.html
"""

from .xarray import *
import xarray as xr

def calc_wind_speed(u, v, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)
    return w.magnitude()

def calc_relative_vorticity_and_horizontal_divergence(u, v, truncation=None, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)
    vrt, div = w.vrtdiv(truncation=truncation)

    data = xr.Dataset()
    data['vrt'] = vrt
    data['div'] = div
    return data

def calc_relative_vorticity(u, v, truncation=None, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)
    vrt = w.vorticity(truncation=truncation)
    return vrt

def calc_divergence(u, v, truncation=None, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)
    div = w.divergence(truncation=truncation)
    return div

def calc_planetary_vorticity(u, v, omega=7.2921150, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)
    pvrt = w.planetaryvorticity(omega=omega)
    return pvrt

def calc_absolute_vorticity(u, v, truncation=None, omega=7.2921150, R = 6371200.0, legfunc='stored'):
    """

    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)
    avrt = w.absolutevorticity(omega=omega, truncation=truncation)
    return avrt

def calc_streamfunction_and_velocity_potential(u, v, truncation=None, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)
    sf, vp = w.sfvp(truncation=truncation)

    data = xr.Dataset()
    data['stream'] = sf
    data['pv'] = vp
    return data

def calc_streamfunction(u, v, truncation=None, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)
    sf = w.streamfunction(truncation=truncation)
    return sf

def calc_velocity_potential(u, v, truncation=None, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)
    vp = w.velocitypotential(truncation=truncation)
    return vp

def calc_helmholtz(u, v, truncation=None, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)
    uchi, vchi, upsi, vpsi = w.helmholtz(truncation=truncation)

    data = xr.Dataset()
    data['uchi'] = uchi
    data['vchi'] = vchi
    data['upsi'] = upsi
    data['vpsi'] = vpsi
    return data

def calc_irrotational_component(u, v, truncation=None, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)
    uchi, vchi = w.irrotationalcomponent(truncation=truncation)

    data = xr.Dataset()
    data['uchi'] = uchi
    data['vchi'] = vchi
    return data

def calc_nondivergent_component(u, v, truncation=None, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)
    upsi, vpsi = w.nondivergentcomponent(truncation=truncation)

    data = xr.Dataset()
    data['upsi'] = upsi
    data['vpsi'] = vpsi
    return data

def calc_rossby_wave_source(u, v, truncation=None, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(u, v, rsphere=R, legfunc=legfunc)

    eta = w.absolutevorticity()
    div = w.divergence(truncation=truncation)
    uchi, vchi = w.irrotationalcomponent(truncation=truncation)
    etax, etay = w.gradient(eta, truncation=truncation)
    etax.attrs['units'] = 'm**-1 s**-1'
    etay.attrs['units'] = 'm**-1 s**-1'

    S = eta * -1. * div - (uchi * etax + vchi * etay)
    return S

def calc_gradient(data, truncation=None, R = 6371200.0, legfunc='stored'):
    """
    
    """
    w = VectorWind(data, data, rsphere=R, legfunc=legfunc)
    data_zonal, data_meridional = w.gradient(data, truncation=truncation)

    data = xr.Dataset()
    data['zonal_gradient'] = data_zonal
    data['meridional_gradient'] = data_meridional
    return data
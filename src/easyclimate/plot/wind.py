def barbs(ds, x, y, ax, u, v, cbar_kwargs = {}, **kwargs):
    """Barbs plot of Dataset variables.

    Wraps :py:func:`matplotlib:matplotlib.pyplot.barbs`.
    """
    import matplotlib.pyplot as plt
    from xarray.core.alignment import broadcast

    if x is None or y is None or u is None or v is None:
        raise ValueError("Must specify x, y, u, v for quiver plots.")

    x, y, u, v = broadcast(ds[x], ds[y], ds[u], ds[v])

    args = [x.values, y.values, u.values, v.values]
    hueexist = "hue" in kwargs

    if "hue" in kwargs or "cmap_params" in kwargs or "hue_style" in kwargs:
        hue = kwargs.pop("hue")

        if hue:
            args.append(ds[hue].values)
            
    kwargs.setdefault("pivot", "middle")
    hdl = ax.barbs(*args, **kwargs)

    if hueexist:
        hdl1 = plt.colorbar(hdl, ax = ax, **cbar_kwargs)
    else:
        hdl1 = hdl
   
    return hdl1
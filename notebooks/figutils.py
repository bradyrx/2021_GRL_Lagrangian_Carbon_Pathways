import numpy as np
import gsw
import PyCO2SYS as pyco2
import xarray as xr
from xhistogram.xarray import histogram as xhist

ANT_PEN_LAT = -63.3833
REGION_NAMES = ["drake", "campbell", "crozet", "kerguelan"]
BOUNDS = {
    "drake": [-71, -25, ANT_PEN_LAT, -50],
    "campbell": [145, 180, -63, -47],
    "crozet": [22, 55, -60, -45],
    "kerguelan": [68, 100, -60, -45],
}


def autocorr_by_hand(x, lag):
    """Computes the autocorrelation coefficient.

    See:
    https://stackoverflow.com/questions/36038927/
    whats-the-difference-between-pandas-acf-and-statsmodel-acf
    
    x : time series to autocorrelate
    lag : int of how many lags to evaluate at
    """
    # Slice the relevant subseries based on the lag
    y1 = x[: (len(x) - lag)]
    y2 = x[lag:]
    # Subtract the subseries means
    sum_product = np.sum((y1 - np.mean(y1)) * (y2 - np.mean(y2)))
    # Normalize with the subseries stds
    return sum_product / ((len(x) - lag) * np.std(y1) * np.std(y2))


def compute_idx_of_first_200m_crossing(z):
    """Find first time particle upwells across 200 m.

    Generally used after subsetting the last 1000 m
    crossing.
    
    z : zLevelParticle
    """
    currentDepth = z
    previousDepth = np.roll(z, 1)
    previousDepth[0] = 999
    cond = (currentDepth >= -200) & (previousDepth < -200)
    idx = cond.argmax()
    return idx


def compute_idx_of_last_1000m_crossing(z):
    """Find index of final time particle upwells across 1000 m.
    
    z : zLevelParticle
    """
    currentDepth = z
    previousDepth = np.roll(z, 1)
    previousDepth[0] = 999  # So we're not dealing with a nan here.
    cond = (currentDepth >= -1000) & (previousDepth < -1000)
    idx = (
        len(cond) - np.flip(cond).argmax() - 1
    )  # Finds last location that condition is true.
    return idx


def deep_tracer_origin(x, y, z, tracer):
    """Finds x, y, z origin and content of tracer before its last 1000 m crossing.

    Memory time is decided by the e-folding time.

        x (ndarray): lonParticle
        y (ndarray): latParticle
        z (ndarray): zLevelDepth
        tracer (ndarray): e.g., particleDIC
    """
    # For when I input single-particle DataArrays.
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    tracer = np.asarray(tracer)

    # find last crossing to 1000m and trim to this length
    idx = compute_idx_of_last_1000m_crossing(z)
    x_subset = x[0:idx]
    y_subset = y[0:idx]
    z_subset = z[0:idx]
    tracer_subset = tracer[0:idx]

    # compute ACF, must be at least a year long.
    auto = np.asarray(
        [autocorr_by_hand(tracer_subset, i) for i in range(len(tracer_subset))]
    )
    e_folding = (auto <= (1 / np.e)).argmax()
    # Go backward this many steps.
    subset_idx = len(z_subset) - e_folding - 1
    # If there is no zero crossing, argmax will return exactly zero. Which can't be the
    # case ever since at zero, ACF == 1 by definition.
    if e_folding == 0:
        x_origin = np.nan
        y_origin = np.nan
        z_origin = np.nan
        tracer_origin = np.nan
    else:
        x_origin = x_subset[subset_idx]
        y_origin = y_subset[subset_idx]
        z_origin = z_subset[subset_idx]
        tracer_origin = tracer_subset[subset_idx]
    return np.array([x_origin, y_origin, z_origin, tracer_origin])


def derive_variables(ds):
    """Derives a bunch of useful variables from the particle output.

    * dbar : depth * -1
    * in situ temperature : in situ temperature
    * rho : in situ density
    * sigma0 : potential density
    * pCO2 : in situ pCO2
    * pco2_theta : potential pCO2
    """
    ds["dbar"] = ("time", ds["zLevelParticle"] * -1)

    # in situ temperature
    insitu = gsw.pt_from_t(ds.particleSalinity, ds.particleTemperature, 0, ds.dbar)
    ds["insitu_temp"] = ("time", insitu)

    # in situ density
    ds["rho"] = (
        "time",
        gsw.density.rho(ds.particleSalinity, ds.particleTemperature, ds.dbar),
    )

    # potential density relative to surface.
    ds["sigma0"] = (
        "time",
        gsw.density.sigma0(ds.particleSalinity, ds.particleTemperature) + 1000,
    )

    # in situ pCO2
    # based on in situ density. converts from mmol/m3 to umol/kg.
    conversion = 1000 * (1 / ds.rho)

    result = pyco2.CO2SYS_nd(
        ds.particleALK * conversion,  # in situ Alk, umol/kg
        ds.particleDIC * conversion,  # in situ DIC, umol/kg
        1,
        2,
        salinity=ds.particleSalinity,  # in situ salinity
        temperature=ds.insitu_temp,  # in situ temperature
        pressure=ds.dbar,  # pressure (depth) in db/m (where alk/dic were taken)
        total_silicate=ds.particleSiO3 * conversion,  # in situ SiO3 in umol/kg
        total_phosphate=ds.particlePO4 * conversion,  # in situ PO4 in umol/kg
    )["pCO2"]
    ds["pCO2"] = ("time", result)

    # potential pCO2 -- defining as the pCO2 it'd be at the ambient temperature of the
    # particle's mixed layer crossing.
    ambient = find_200m_ambient_temperature(
        ds.zLevelParticle.values, ds.insitu_temp.values
    )
    ds["pco2_theta"] = ("time", ds.pCO2 * (1 + 0.0423 * (ambient - ds.insitu_temp)))
    return ds


def find_200m_ambient_temperature(z, t):
    """
    Finds the ambient temperature at the first 200m crossing point
    following the particle's last 1000 m crossing. I.e., for its
    mixed layer crossing.
    """
    idx_a = compute_idx_of_last_1000m_crossing(z)
    z_subset = z[idx_a::]
    idx_b = compute_idx_of_first_200m_crossing(z_subset)
    return t[idx_a + idx_b]


def mixed_layer_tracer_origin(x, y, z, tracer):
    """Finds x, y, z origin and content of tracer before its mixed layer
    crossing (following its last 1000m crossing).

    Memory time is decided by the e-folding time.
    """
    # For when I input single-particle DataArrays.
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    tracer = np.asarray(tracer)

    idx_a = compute_idx_of_last_1000m_crossing(z)
    z_subset = z[idx_a::]
    idx_b = compute_idx_of_first_200m_crossing(z_subset)
    x_subset = x[idx_a : idx_a + idx_b + 1]
    y_subset = y[idx_a : idx_a + idx_b + 1]
    z_subset = z[idx_a : idx_a + idx_b + 1]
    tracer_subset = tracer[idx_a : idx_a + idx_b + 1]

    auto = np.asarray(
        [autocorr_by_hand(tracer_subset, i) for i in range(len(tracer_subset))]
    )
    e_folding = (auto <= (1 / np.e)).argmax()
    # Go backward this many steps
    subset_idx = len(z_subset) - e_folding - 1
    # If there is no zero crossing, argmax will return exactly zero. Which can't be the
    # case ever since at zero, ACF == 1 by definition.
    if e_folding == 0:
        x_origin = np.nan
        y_origin = np.nan
        z_origin = np.nan
        tracer_origin = np.nan
    else:
        x_origin = x_subset[subset_idx]
        y_origin = y_subset[subset_idx]
        z_origin = z_subset[subset_idx]
        tracer_origin = tracer_subset[subset_idx]
    return np.array([x_origin, y_origin, z_origin, tracer_origin])


def plot_box(
    ax, x0, x1, y0, y1, color="black", linewidth=1.5, linestyle="--", **kwargs
):
    def _box(ax, x0, x1, y0, y1, **kwargs):
        ax.plot([x0, x1], [y0, y1], **kwargs)

    _box(
        ax,
        x0,
        x1,
        y0,
        y0,
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
        **kwargs
    )
    _box(
        ax,
        x1,
        x1,
        y0,
        y1,
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
        **kwargs
    )
    _box(
        ax,
        x0,
        x1,
        y1,
        y1,
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
        **kwargs
    )
    _box(
        ax,
        x0,
        x0,
        y0,
        y1,
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
        **kwargs
    )


def return_xhist_crossings(lon_crossing, lat_crossing, res):
    """Uses `xhistogram` to create map of N crossings into 1000m

    Args:
        lon_crossing (xr.DataArray): Lon. location of final crossing into
            1000m (degrees)
        lat_crossing (xr.DataArray): Lat. location of final crossing into
            1000m (degrees)
        res (int, float): Resolution of grid for bin counting.

    Returns:
        xr.DataArray of N crossings at each grid cell location.
    """
    lon_crossing = lon_crossing.rename("x")
    lat_crossing = lat_crossing.rename("y")
    # I add one to the lat/lon because `xhist` reduces it by one.
    grid_x = np.linspace(0, 360, int((360 / res)) + 1)
    grid_y = np.linspace(-90, 90, int((180 / res)) + 1)
    histogram = xhist(lat_crossing, lon_crossing, bins=[grid_y, grid_x])
    histogram = histogram.rename({"y_bin": "lat", "x_bin": "lon"})
    return histogram

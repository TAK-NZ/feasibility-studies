"""OpenDrift run helpers for the hindcast: one run carries many drifters (and many windage values).

All runs: OpenDrift 1.14.12, readers = opendrift.readers.reader_netCDF_CF_generic, time step 900 s, output every 6 h,
Runge-Kutta advection, no vertical motion, 48 h long, fixed random seed (the `seed` argument of the model).
Positions are returned at the leads 6, 12, 24 and 48 h. Deactivated elements keep their last position.
"""
import logging
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import xarray as xr
from opendrift.models.leeway import Leeway
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic as rr

import config as C
from drift_lib import geod_km, lon180

logging.getLogger("opendrift").setLevel(logging.CRITICAL)
_READERS = {}
LEADS = list(C.LEAD_IDX)
IDX = [C.LEAD_IDX[h] for h in LEADS]


def get_reader(path):
    path = str(path)
    if path not in _READERS:
        if len(_READERS) > 12:
            _READERS.clear()
        _READERS[path] = rr.Reader(path)
    return _READERS[path]


def tag_of(t0):
    return pd.Timestamp(t0).strftime("%Y%m%dT%H")


def forcing_paths(block, current, wind, t0, stokes):
    """Reader file list: current file, Stokes file (if used or required), wind file."""
    p = [C.FORCING_DIR / str(block) / f"{current}.nc"]
    if stokes:
        p.append(C.FORCING_DIR / str(block) / "stokes.nc")
    p.append(C.FORCING_DIR / "winds" / f"{tag_of(t0)}_{wind}.nc")
    return p


def wind_available(t0, wind):
    return (C.FORCING_DIR / "winds" / f"{tag_of(t0)}_{wind}.nc").exists()


def _extract(o, n):
    """lon, lat (n, 4), status code (n,), count of elements not active at 48 h."""
    r = o.result
    lon = r.lon.values.astype(float)
    lat = r.lat.values.astype(float)
    st = r.status.values.astype(float)
    order = np.argsort(r.trajectory.values)
    lon, lat, st = lon[order], lat[order], st[order]
    assert lon.shape == (n, 9), lon.shape
    # forward fill positions of deactivated elements
    for a in (lon, lat, st):
        df = pd.DataFrame(a).ffill(axis=1)
        a[:] = df.values
    cats = list(o.status_categories)
    code = st[:, -1]
    name = np.array([cats[int(c)] if np.isfinite(c) else "none" for c in code])
    return lon180(lon[:, IDX]), lat[:, IDX], name


def run_oceandrift(t0, block, lon, lat, current, wind, stokes, windage, cdf=None,
                   cur_unc=0.0, wind_unc=0.0, seed=C.SIM_SEED):
    """Seed n elements at (lon, lat) at time t0. windage and cdf are per element (array or scalar)."""
    n = len(lon)
    o = OceanDrift(loglevel=50, seed=seed)
    for p in forcing_paths(block, current, wind, t0, stokes):
        o.add_reader(get_reader(p))
    o.set_config("general:use_auto_landmask", True)
    o.set_config("drift:vertical_mixing", False)
    o.set_config("drift:vertical_advection", False)
    o.set_config("drift:advection_scheme", "runge-kutta")
    o.set_config("drift:stokes_drift", bool(stokes))
    o.set_config("seed:ocean_only", False)          # coastal seeds stay where the drifter was; stranding handles land
    o.set_config("general:time_step_minutes", C.TIME_STEP_S / 60)
    o.set_config("drift:current_uncertainty", float(cur_unc))
    o.set_config("drift:wind_uncertainty", float(wind_unc))
    kw = dict(wind_drift_factor=np.broadcast_to(np.asarray(windage, float), (n,)).copy())
    if cdf is not None:
        kw["current_drift_factor"] = np.broadcast_to(np.asarray(cdf, float), (n,)).copy()
    o.seed_elements(lon=np.asarray(lon, float), lat=np.asarray(lat, float), time=pd.Timestamp(t0).to_pydatetime(),
                    number=n, origin_marker=np.zeros(n, dtype=int), **kw)
    o.run(duration=timedelta(hours=C.WINDOW_HOURS), time_step=C.TIME_STEP_S,
          time_step_output=timedelta(hours=C.OUTPUT_STEP_H))
    return _extract(o, n)


def run_leeway(t0, block, lon, lat, current, wind, object_type=1, cur_unc=0.0, wind_unc=0.0, seed=C.SIM_SEED):
    n = len(lon)
    o = Leeway(loglevel=50, seed=seed)
    for p in forcing_paths(block, current, wind, t0, True):
        o.add_reader(get_reader(p))
    o.set_config("general:use_auto_landmask", True)
    o.set_config("drift:advection_scheme", "runge-kutta")
    o.set_config("seed:ocean_only", False)
    o.set_config("general:time_step_minutes", C.TIME_STEP_S / 60)
    o.set_config("drift:current_uncertainty", float(cur_unc))
    o.set_config("drift:wind_uncertainty", float(wind_unc))
    o.seed_elements(lon=np.asarray(lon, float), lat=np.asarray(lat, float), time=pd.Timestamp(t0).to_pydatetime(),
                    number=n, object_type=object_type)
    o.run(duration=timedelta(hours=C.WINDOW_HOURS), time_step=C.TIME_STEP_S,
          time_step_output=timedelta(hours=C.OUTPUT_STEP_H))
    return _extract(o, n)


# ------------------------------------------------------------------ values at the datum (baseline inputs)
def sample_field(path, uname, vname, lon, lat, t0):
    """Bilinear (space) value of a 2-D vector field at time t0 (exact time step in the files). Land (NaN) cells
    fall back to the mean of the valid cells within 2 grid cells; returns (u, v, filled flag)."""
    ds = xr.open_dataset(path)
    lon360 = np.mod(np.asarray(lon, float), 360.0)
    pts = dict(latitude=xr.DataArray(np.asarray(lat, float), dims="p"), longitude=xr.DataArray(lon360, dims="p"))
    sl = ds.sel(time=np.datetime64(pd.Timestamp(t0)), method="nearest")
    ip = sl.interp(**pts)
    u, v = ip[uname].values.astype(float), ip[vname].values.astype(float)
    filled = ~np.isfinite(u)
    for i in np.where(filled)[0]:
        ii = int(np.abs(sl.latitude.values - lat[i]).argmin())
        jj = int(np.abs(sl.longitude.values - lon360[i]).argmin())
        sub_u = sl[uname].values[max(ii - 2, 0):ii + 3, max(jj - 2, 0):jj + 3]
        sub_v = sl[vname].values[max(ii - 2, 0):ii + 3, max(jj - 2, 0):jj + 3]
        u[i] = np.nanmean(sub_u) if np.isfinite(sub_u).any() else 0.0
        v[i] = np.nanmean(sub_v) if np.isfinite(sub_v).any() else 0.0
    ds.close()
    return u, v, filled


def advect_constant(lon, lat, u, v, hours):
    """Position after `hours` of constant velocity (u east, v north, m/s) on the WGS84 ellipsoid."""
    from drift_lib import GEOD
    spd = np.hypot(u, v)
    az = np.degrees(np.arctan2(u, v))
    lo, la, _ = GEOD.fwd(np.asarray(lon, float), np.asarray(lat, float), az, spd * hours * 3600.0)
    return lon180(lo), la


def obs_at_leads(w):
    """Observed lon, lat (n, 4) at the leads from a windows table."""
    lon = np.column_stack([w[f"lon_{C.LEAD_IDX[h]}"].values for h in LEADS])
    lat = np.column_stack([w[f"lat_{C.LEAD_IDX[h]}"].values for h in LEADS])
    return lon, lat


def err_km(w, lon, lat):
    olon, olat = obs_at_leads(w)
    return geod_km(olon, olat, lon, lat)

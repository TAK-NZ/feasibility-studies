"""Build OpenDrift-ready forcing files, the final windows table and the Eulerian-skill diagnostic.

Reads  data/cmems/*.nc, data/winds/{gfs,ecmwf}/*.nc, results/windows_selected.csv, data/gdp/gdp_6h_box.csv
Writes data/forcing/<block>/{cur_tide,cur_notide,cur_6h_surf,cur_6h_15m,cur_6h_15m_tide,stokes}.nc
       data/forcing/winds/<t0>_{gfs,ecmwf,both}.nc
       results/windows.csv         one row per selected window: split, stratum, season, wind class, drifter id ...
       results/forcing_check.csv   model current vs drifter ve/vn per class and current variant

Current variants (D1 in work/plan.md):
  cur_tide        hourly uo+utide, vo+vtide at 0.494 m   (utotal is NOT used, it contains Stokes drift)
  cur_notide      hourly uo, vo at 0.494 m
  cur_6h_surf     6-hourly instantaneous uo, vo at 0.494 m
  cur_6h_15m      6-hourly instantaneous uo, vo at 15.81 m (drogue depth of an SVP drifter)
  cur_6h_15m_tide cur_6h_15m interpolated to hourly + hourly tide at 0.494 m (assumes the tide is depth-uniform)
Stokes drift (vsdx, vsdy) comes from the hourly merged-uv file for every variant.
Wind "both" = arithmetic mean of GFS and ECMWF u and v after ECMWF is bilinearly interpolated to the GFS grid (D5).
Usage: python scripts/s05_build_forcing.py [--workers 4]
"""
import argparse
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import xarray as xr

import config as C

FORC = C.DATA / "forcing"
ENC = dict(zlib=True, complevel=3, dtype="float32")
UN = "m s-1"
SN = dict(
    u="eastward_sea_water_velocity", v="northward_sea_water_velocity",
    sx="sea_surface_wave_stokes_drift_x_velocity", sy="sea_surface_wave_stokes_drift_y_velocity",
)
CUR_VARIANTS = ["cur_tide", "cur_notide", "cur_6h_surf", "cur_6h_15m", "cur_6h_15m_tide"]
WIND_VARIANTS = ["gfs", "ecmwf", "both"]


def _cur_ds(u, v):
    return xr.Dataset({
        "uo": u.assign_attrs(standard_name=SN["u"], units=UN),
        "vo": v.assign_attrs(standard_name=SN["v"], units=UN)})


def _write(ds, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.nc")
    ds.to_netcdf(tmp, encoding={v: ENC for v in ds.data_vars})
    tmp.rename(path)


def build_block(block):
    out = FORC / block
    names = ["cur_tide", "cur_notide", "cur_6h_surf", "cur_6h_15m", "cur_6h_15m_tide", "stokes"]
    if all((out / f"{n}.nc").exists() for n in names):
        return block, "skipped"
    h = xr.open_dataset(C.CMEMS_DIR / f"uv1h_{block}.nc").isel(depth=0, drop=True).load()
    s = xr.open_dataset(C.CMEMS_DIR / f"cur6h_s_{block}.nc").isel(depth=0, drop=True).load()
    d = xr.open_dataset(C.CMEMS_DIR / f"cur6h_15_{block}.nc").isel(depth=0, drop=True).load()
    _write(_cur_ds(h.uo + h.utide, h.vo + h.vtide), out / "cur_tide.nc")
    _write(_cur_ds(h.uo, h.vo), out / "cur_notide.nc")
    _write(_cur_ds(s.uo, s.vo), out / "cur_6h_surf.nc")
    _write(_cur_ds(d.uo, d.vo), out / "cur_6h_15m.nc")
    for x in (s, d):   # float32 coordinates differ in the last bits: align on the hourly file's grid
        assert np.allclose(x.longitude, h.longitude, atol=1e-3) and np.allclose(x.latitude, h.latitude, atol=1e-3)
    d = d.assign_coords(latitude=h.latitude, longitude=h.longitude)
    s = s.assign_coords(latitude=h.latitude, longitude=h.longitude)
    di = d.interp(time=h.time)   # linear in time, as OpenDrift would do
    _write(_cur_ds(di.uo + h.utide, di.vo + h.vtide), out / "cur_6h_15m_tide.nc")
    st = xr.Dataset({
        "vsdx": h.vsdx.assign_attrs(standard_name=SN["sx"], units=UN),
        "vsdy": h.vsdy.assign_attrs(standard_name=SN["sy"], units=UN)})
    _write(st, out / "stokes.nc")
    return block, "built"


def build_winds(t0):
    """t0 = 'YYYYMMDDTHH'. Writes gfs/ecmwf/both files with x_wind/y_wind standard names."""
    outdir = FORC / "winds"
    outdir.mkdir(parents=True, exist_ok=True)
    g = xr.open_dataset(C.WIND_DIR / "gfs" / f"{t0}.nc").load()
    status = []

    def std(ds):
        return xr.Dataset({
            "u10": ds.u10.assign_attrs(standard_name="eastward_wind", units=UN),
            "v10": ds.v10.assign_attrs(standard_name="northward_wind", units=UN)})

    if not (outdir / f"{t0}_gfs.nc").exists():
        _write(std(g), outdir / f"{t0}_gfs.nc")
    ef = C.WIND_DIR / "ecmwf" / f"{t0}.nc"
    if ef.exists():
        if not (outdir / f"{t0}_ecmwf.nc").exists() or not (outdir / f"{t0}_both.nc").exists():
            e = xr.open_dataset(ef).load()
            _write(std(e), outdir / f"{t0}_ecmwf.nc")
            ei = e.interp(latitude=g.latitude, longitude=g.longitude, kwargs={"fill_value": "extrapolate"})
            both = (g + ei) / 2.0     # mean of u and v separately
            _write(std(both), outdir / f"{t0}_both.nc")
            status.append("ecmwf+both")
    else:
        status.append("no ecmwf")
    return t0, ",".join(status)


def wind_class(ms):
    return "<5" if ms < 5 else ("5-10" if ms < 10 else ">=10")


def write_windows():
    w = pd.read_csv(C.RESULTS / "windows_selected.csv", dtype={"drifter_id": str, "block_id": str})
    w["t0"] = pd.to_datetime(w.t0)
    w["t0_tag"] = w.t0.dt.strftime("%Y%m%dT%H")
    w["stratum"] = w.zone + "_" + w.drogue_class
    w["has_ecmwf"] = [(C.WIND_DIR / "ecmwf" / f"{t}.nc").exists() for t in w.t0_tag]
    # mean GFS forecast wind speed at the datum over the 17 valid times (t0 .. t0+48 h)
    ws = np.full(len(w), np.nan)
    for tag, idx in w.groupby("t0_tag").groups.items():
        g = xr.open_dataset(C.WIND_DIR / "gfs" / f"{tag}.nc")
        sel = w.loc[idx]
        lon = np.mod(sel.lon_0.values, 360.0)
        ii = np.abs(g.latitude.values[:, None] - sel.lat_0.values[None]).argmin(0)
        jj = np.abs(g.longitude.values[:, None] - lon[None]).argmin(0)
        u = g.u10.values[:, ii, jj]
        v = g.v10.values[:, ii, jj]
        ws[w.index.get_indexer(idx)] = np.hypot(u, v).mean(axis=0)
    w["gfs_wind_ms"] = ws
    w["wind_class"] = [wind_class(x) for x in ws]
    # a drifter must be in one split only
    nsplit = w.groupby("drifter_id").split.nunique()
    assert (nsplit == 1).all(), "a drifter appears in both tuning and test"
    cols = ["window_id", "drifter_id", "split", "stratum", "zone", "drogue_class", "season", "wind_class",
            "gfs_wind_ms", "t0", "block_id", "has_ecmwf", "min_land_km"] + \
           [c for c in w.columns if c.startswith(("lon_", "lat_"))] + ["drogue_lost_date", "ve_pre", "vn_pre"]
    w = w[cols]
    w.to_csv(C.RESULTS / "windows.csv", index=False)
    print(f"results/windows.csv: {len(w)} windows, {w.drifter_id.nunique()} drifters")
    print(w.groupby(["split", "stratum"]).agg(windows=("window_id", "size"), drifters=("drifter_id", "nunique")))
    print(w.groupby(["split", "wind_class"]).size().unstack())
    return w


def forcing_check(w):
    """Model current at the observed positions versus the drifter velocity ve, vn (6-hourly rows)."""
    g = pd.read_csv(C.GDP_CSV, dtype={"ID": str}, usecols=["ID", "time", "lon360", "latitude", "ve", "vn"])
    g["time"] = pd.to_datetime(g.time, utc=True).dt.tz_localize(None)
    for v in ("ve", "vn"):
        g.loc[g[v] <= -999990, v] = np.nan
    recs = []
    for _, r in w.iterrows():
        for k in range(1, C.N_OBS):   # skip k=0 (ve/vn missing at the start of tracks)
            recs.append((r.window_id, r.drifter_id, r.block_id, r.drogue_class, r.split,
                         r.t0 + pd.Timedelta(hours=6 * k), r[f"lon_{k}"], r[f"lat_{k}"]))
    o = pd.DataFrame(recs, columns=["window_id", "drifter_id", "block_id", "drogue_class", "split", "time", "lon", "lat"])
    o = o.merge(g[["ID", "time", "ve", "vn"]].rename(columns={"ID": "drifter_id"}), on=["drifter_id", "time"], how="left")
    o["lon"] = np.mod(o.lon, 360.0)
    o = o.drop_duplicates(["drifter_id", "time"]).dropna(subset=["ve", "vn"])
    for var in CUR_VARIANTS:
        o[f"u_{var}"] = np.nan
        o[f"v_{var}"] = np.nan
    for blk, idx in o.groupby("block_id").groups.items():
        sel = o.loc[idx]
        pts = dict(time=xr.DataArray(sel.time.values, dims="p"), latitude=xr.DataArray(sel.lat.values, dims="p"),
                   longitude=xr.DataArray(sel.lon.values, dims="p"))
        for var in CUR_VARIANTS:
            ds = xr.open_dataset(FORC / blk / f"{var}.nc")
            ip = ds.interp(**pts)
            o.loc[idx, f"u_{var}"] = ip.uo.values
            o.loc[idx, f"v_{var}"] = ip.vo.values
    rows = []
    for (cls, split), grp in o.groupby(["drogue_class", "split"]):
        for var in CUR_VARIANTS:
            m = grp.dropna(subset=[f"u_{var}", f"v_{var}"])
            du, dv = m[f"u_{var}"] - m.ve, m[f"v_{var}"] - m.vn
            # vector correlation (Kundu): sum of dot products over the root of the sum of squares
            num = (m.ve * m[f"u_{var}"] + m.vn * m[f"v_{var}"]).sum()
            den = np.sqrt((m.ve ** 2 + m.vn ** 2).sum() * (m[f"u_{var}"] ** 2 + m[f"v_{var}"] ** 2).sum())
            rows.append(dict(drogue_class=cls, split=split, current=var, n_obs=len(m), n_drifters=m.drifter_id.nunique(),
                             bias_u=du.mean(), bias_v=dv.mean(), rmse_vec=float(np.sqrt((du ** 2 + dv ** 2).mean())),
                             speed_obs=float(np.hypot(m.ve, m.vn).mean()), speed_model=float(np.hypot(m[f"u_{var}"], m[f"v_{var}"]).mean()),
                             vector_corr=float(num / den)))
    fc = pd.DataFrame(rows)
    fc.to_csv(C.RESULTS / "forcing_check.csv", index=False)
    print(fc.round(3).to_string())
    return fc


def smoke_test(block):
    from opendrift.readers import reader_netCDF_CF_generic as rr
    need = {"cur": ["x_sea_water_velocity", "y_sea_water_velocity"],
            "stokes": ["sea_surface_wave_stokes_drift_x_velocity", "sea_surface_wave_stokes_drift_y_velocity"]}
    for var in CUR_VARIANTS:
        r = rr.Reader(str(FORC / block / f"{var}.nc"))
        assert all(x in r.variables for x in need["cur"]), (var, r.variables)
        assert r.covers_positions(np.array([-179.5]), np.array([-40.0]))[0].size == 1, var
    r = rr.Reader(str(FORC / block / "stokes.nc"))
    assert all(x in r.variables for x in need["stokes"]), r.variables
    for wv in WIND_VARIANTS:
        t0 = pd.Timestamp(block).strftime("%Y%m%dT00")
        f = FORC / "winds" / f"{t0}_{wv}.nc"
        if f.exists():
            r = rr.Reader(str(f))
            assert "x_wind" in r.variables and "y_wind" in r.variables, (wv, r.variables)
    print("smoke test OK on block", block)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    blocks = pd.read_csv(C.RESULTS / "blocks.csv", dtype={"block_id": str}).block_id.tolist()
    with ProcessPoolExecutor(a.workers) as ex:
        for b, st in ex.map(build_block, blocks):
            print("block", b, st)
    tags = sorted(p.stem for p in (C.WIND_DIR / "gfs").glob("*.nc"))
    with ProcessPoolExecutor(a.workers) as ex:
        for t, st in ex.map(build_winds, tags):
            if st:
                print("wind", t, st)
    smoke_test(blocks[0])
    w = write_windows()
    forcing_check(w)

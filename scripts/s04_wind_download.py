"""Download forecast-mode 10 m winds (GFS 0.25 deg, ECMWF IFS open data) for every seed time.

For seed time t0 (00Z) the cycle t0-12 h (12Z of the previous day) is used, forecast hours 12..60 step 3
(17 valid times covering t0 .. t0+48 h). Fields are fetched by HTTP byte range from the free AWS buckets
(no login), decoded with eccodes, cut to the forcing box, stored with longitude 0..360.
  GFS:   s3 noaa-gfs-bdp-pds        (.idx text index)
  ECMWF: s3 ecmwf-forecasts (CC-BY-4.0)  (.index JSON-lines index; path layout changes over time)
Output: data/winds/<model>/<t0>.nc and data/winds/manifest.csv
Usage: python scripts/s04_wind_download.py
"""
import json
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import eccodes
import numpy as np
import pandas as pd
import requests
import xarray as xr

import config as C

H = pd.Timedelta(hours=1)
_local = threading.local()
_prefix_cache = {}
_lock = threading.Lock()


def session():
    if not hasattr(_local, "s"):
        _local.s = requests.Session()
    return _local.s


def get(url, rng=None, tries=8):
    """GET with retries; back off on 503 SlowDown and transient errors. Returns bytes or None for 404."""
    hdr = {"Range": f"bytes={rng[0]}-{rng[1]}"} if rng else {}
    for i in range(tries):
        try:
            r = session().get(url, headers=hdr, timeout=60)
            if r.status_code == 404:
                return None
            if r.status_code in (200, 206):
                return r.content
            if r.status_code not in (503, 500, 502, 504, 429):
                r.raise_for_status()
        except (requests.ConnectionError, requests.Timeout):
            pass
        time.sleep(min(60, 2 ** i))
    raise RuntimeError(f"giving up on {url}")


_decode_lock = threading.Lock()   # eccodes/HDF5 are not thread-safe: use processes (ProcessPoolExecutor) for parallelism


def decode(buf):
    """All GRIB messages in buf -> dict shortName -> (lats, lons, values2d)."""
    with _decode_lock:
        return _decode(buf)


def _decode(buf):
    out = {}
    pos = 0
    while pos < len(buf):
        assert buf[pos:pos + 4] == b"GRIB", "unexpected bytes between GRIB messages"
        length = int.from_bytes(buf[pos + 8:pos + 16], "big")   # GRIB2 total message length
        h = eccodes.codes_new_from_message(buf[pos:pos + length])
        pos += length
        try:
            name = eccodes.codes_get(h, "shortName")
            ni, nj = eccodes.codes_get(h, "Ni"), eccodes.codes_get(h, "Nj")
            lat = eccodes.codes_get_array(h, "latitudes").reshape(nj, ni)[:, 0]
            lon = eccodes.codes_get_array(h, "longitudes").reshape(nj, ni)[0, :]
            val = eccodes.codes_get_values(h).reshape(nj, ni)
            out[name] = (lat, lon, val)
        finally:
            eccodes.codes_release(h)
    return out


def cut(lat, lon, val):
    """Cut to the forcing box, longitude 0..360 ascending, latitude ascending."""
    lon = np.where(lon < 0, lon + 360.0, lon)
    jo, io_ = np.argsort(lat), np.argsort(lon)
    lat, lon, val = lat[jo], lon[io_], val[jo][:, io_]
    m = 1.0
    jm = (lat >= C.FORCING_LAT[0] - m) & (lat <= C.FORCING_LAT[1] + m)
    im = (lon >= C.FORCING_LON[0] - m) & (lon <= C.FORCING_LON[1] + m)
    return lat[jm], lon[im], val[np.ix_(jm, im)].astype("float32")


# ---------------- GFS
def gfs_field(cycle, fh):
    d, hh = cycle.strftime("%Y%m%d"), cycle.strftime("%H")
    base = f"{C.GFS_BASE}/gfs.{d}/{hh}/atmos/gfs.t{hh}z.pgrb2.0p25.f{fh:03d}"
    idx = get(base + ".idx")
    if idx is None:
        return None, "missing"
    lines = idx.decode().strip().split("\n")
    offs = [int(l.split(":")[1]) for l in lines]
    iu = next(i for i, l in enumerate(lines) if ":UGRD:10 m above ground:" in l)
    iv = next(i for i, l in enumerate(lines) if ":VGRD:10 m above ground:" in l)
    assert iv == iu + 1, "UGRD/VGRD lines not adjacent"
    end = offs[iv + 1] - 1 if iv + 1 < len(offs) else ""
    buf = get(base, (offs[iu], end if end != "" else ""))
    return decode(buf), "gfs"


# ---------------- ECMWF
def ecmwf_prefix(cycle):
    d, hh = cycle.strftime("%Y%m%d"), cycle.strftime("%H")
    key = (d, hh)
    with _lock:
        if key in _prefix_cache:
            return _prefix_cache[key]
    found = None
    for p in C.ECMWF_PREFIXES:
        u = f"{C.ECMWF_BASE}/{d}/{hh}z/{p}/oper/{d}{hh}0000-12h-oper-fc.index"
        if get(u) is not None:
            found = p
            break
    with _lock:
        _prefix_cache[key] = found
    return found


def ecmwf_field(cycle, fh):
    d, hh = cycle.strftime("%Y%m%d"), cycle.strftime("%H")
    p = ecmwf_prefix(cycle)
    if p is None:
        return None, "missing"
    base = f"{C.ECMWF_BASE}/{d}/{hh}z/{p}/oper/{d}{hh}0000-{fh}h-oper-fc"
    idx = get(base + ".index")
    if idx is None:
        return None, "missing"
    recs = [json.loads(l) for l in idx.decode().strip().split("\n")]
    out = {}
    for par in ("10u", "10v"):
        r = next(x for x in recs if x["param"] == par and x["levtype"] == "sfc")
        buf = get(base + ".grib2", (r["_offset"], r["_offset"] + r["_length"] - 1))
        out.update(decode(buf))
    return out, p


def fetch_one(model, t0):
    out = C.WIND_DIR / model / f"{t0:%Y%m%dT%H}.nc"
    if out.exists():
        try:
            d = xr.open_dataset(out)
            d.load()
            ok = ("u10" in d and "v10" in d and d.sizes["time"] == len(C.WIND_FHOURS)
                  and bool(np.isfinite(d.latitude.values).all() and np.isfinite(d.longitude.values).all()))
            d.close()
            if not ok:
                raise OSError("incomplete file")
        except OSError:   # truncated by an earlier crash: fetch again
            out.unlink()
    if out.exists():
        with xr.open_dataset(out) as d:
            return dict(model=model, t0=str(t0), status="ok", n_steps=d.sizes["time"], variant=d.attrs.get("path_variant"),
                        res_deg=d.attrs.get("res_deg"), size_mb=round(out.stat().st_size / 1e6, 3), cached=True)
    cycle = t0 - C.WIND_CYCLE_OFFSET_H * H
    fn = gfs_field if model == "gfs" else ecmwf_field
    us, vs, times, variant, grid = [], [], [], None, None
    for fh in C.WIND_FHOURS:
        res, var = fn(cycle, fh)
        if res is None:
            return dict(model=model, t0=str(t0), status=f"missing f{fh:03d} ({var})", n_steps=len(times))
        (la, lo, u), (_, _, v) = res["10u"], res["10v"]
        la, lo, u = cut(la, lo, u)
        _, _, v = cut(res["10v"][0], res["10v"][1], res["10v"][2])
        us.append(u); vs.append(v); times.append(cycle + fh * H)
        variant, grid = var, (la, lo)
    la, lo = grid
    ds = xr.Dataset(
        {"u10": (("time", "latitude", "longitude"), np.stack(us), {"standard_name": "eastward_wind", "units": "m s-1"}),
         "v10": (("time", "latitude", "longitude"), np.stack(vs), {"standard_name": "northward_wind", "units": "m s-1"})},
        coords={"time": np.array(times, dtype="datetime64[ns]"),
                "latitude": ("latitude", la, {"units": "degrees_north", "standard_name": "latitude"}),
                "longitude": ("longitude", lo, {"units": "degrees_east", "standard_name": "longitude"})},
        attrs={"source": model, "cycle": str(cycle), "path_variant": str(variant),
               "res_deg": float(round(abs(lo[1] - lo[0]), 3)), "fhours": " ".join(map(str, C.WIND_FHOURS))})
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".nc.tmp")
    ds.to_netcdf(tmp, encoding={v: {"zlib": True, "complevel": 4} for v in ds.data_vars})
    tmp.replace(out)
    return dict(model=model, t0=str(t0), status="ok", n_steps=len(times), variant=variant,
                res_deg=ds.attrs["res_deg"], size_mb=round(out.stat().st_size / 1e6, 3), cached=False)


def sanity(t0):
    g = xr.open_dataset(C.WIND_DIR / "gfs" / f"{t0:%Y%m%dT%H}.nc")
    e = xr.open_dataset(C.WIND_DIR / "ecmwf" / f"{t0:%Y%m%dT%H}.nc")
    ei = e.interp(latitude=g.latitude, longitude=g.longitude)
    sg, se = np.hypot(g.u10, g.v10), np.hypot(ei.u10, ei.v10)
    dsp = float(abs(sg - se).mean())
    pt = dict(latitude=-40.0, longitude=170.0, method="nearest")
    gu, gv = float(g.u10.isel(time=4).sel(**pt)), float(g.v10.isel(time=4).sel(**pt))
    eu, ev = float(e.u10.isel(time=4).sel(**pt)), float(e.v10.isel(time=4).sel(**pt))
    ang = lambda u, v: (270 - np.degrees(np.arctan2(v, u))) % 360   # meteorological "from" direction
    print(f"sanity t0={t0}: domain-mean |speed diff|={dsp:.2f} m/s; at (-40,170) valid {str(g.time.values[4])[:16]}: "
          f"GFS {np.hypot(gu, gv):.1f} m/s from {ang(gu, gv):.0f} deg, ECMWF {np.hypot(eu, ev):.1f} m/s from {ang(eu, ev):.0f} deg")
    assert dsp < 3.0


def main():
    b = pd.read_csv(C.RESULTS / "blocks.csv", parse_dates=["block_start"])
    t0s = sorted({r.block_start + o * H for r in b.itertuples() for o in C.SEED_OFFSETS_H})
    tasks = [(m, t) for t in t0s for m in ("gfs", "ecmwf")]
    print(f"{len(t0s)} seed times, {len(tasks)} (model, t0) files, {len(C.WIND_FHOURS)} steps each")
    C.WIND_DIR.mkdir(parents=True, exist_ok=True)
    rows, t_start = [], time.time()
    with ProcessPoolExecutor(C.WIND_WORKERS) as ex:
        futs = {ex.submit(fetch_one, m, t): (m, t) for m, t in tasks}
        for i, f in enumerate(as_completed(futs), 1):
            r = f.result()
            rows.append(r)
            if i % 10 == 0 or r["status"] != "ok":
                print(f"[{i}/{len(tasks)}] {r['model']} {r['t0']} {r['status']} (elapsed {time.time() - t_start:.0f}s)", flush=True)
    m = pd.DataFrame(rows).sort_values(["model", "t0"])
    m.to_csv(C.WIND_DIR / "manifest.csv", index=False)
    print(m.groupby(["model", "status"]).size())
    print(m[m.status == "ok"].groupby(["model", "variant", "res_deg"]).size())
    print("total size GB:", round(m.size_mb.sum() / 1000, 3))
    both = sorted(set(m[(m.model == "gfs") & (m.status == "ok")].t0) & set(m[(m.model == "ecmwf") & (m.status == "ok")].t0))
    for t in (both[0], both[len(both) // 2], both[-1]):
        sanity(pd.Timestamp(t))


if __name__ == "__main__":
    main()

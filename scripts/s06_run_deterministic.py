"""Deterministic OpenDrift runs for every selected window (tuning and test together) plus the hand-calculation baselines.

One OpenDrift run per (seed time, class, current, wind, Stokes). All drifters of the class that have a window at that seed
time, times every windage value of the grid, are carried as separate elements (zero uncertainty, so the elements do
not interact).
  lost class     OceanDrift, wind_drift_factor = windage grid, surface currents (cur_tide, cur_notide, cur_6h_surf),
                 wind (gfs, ecmwf, both), Stokes drift on/off          -> 3 x 3 x 2 x 6 = 108 configurations
  drogued class  OceanDrift, wind_drift_factor = 0, Stokes off, currents cur_tide, cur_notide, cur_6h_surf,
                 cur_6h_15m, cur_6h_15m_tide                           -> 5 configurations
Baselines (also written down in work/baseline.md). Everything uses only what an operator has at t0:
  B1 (the hand calculation)   position(t) = datum + V * t on the WGS84 ellipsoid, V constant, evaluated once at the datum at t0:
        V = U_cur + k * U_wind
        U_cur  lost windows:    hourly Eulerian current uo, vo at 0.494 m (no tide, no Stokes), bilinear at the datum, time t0
               drogued windows: 6-hourly current uo, vo at 15.81 m (the drogue depth), bilinear at the datum, time t0
        U_wind forecast-mode 10 m wind valid at t0 (cycle t0 - 12 h, f012), bilinear at the datum; source = gfs, ecmwf or both
        k      lost windows: windage fraction from config.WINDAGE_GRID, chosen on the tuning set (lowest 24 h median B1 error) in s07;
               drogued windows: 0
  expanding circle: radius R_p(t) = a_p * t, a_p (km/h) = p-quantile of the tuning-set 24 h B1 separation divided by 24 h, p = 0.5, 0.95,
        per class (computed in s07, applied in s08). Hit = observed position inside the circle around the B1 position.
  B0  stay at the datum.   B3  datum + the drifter's own mean velocity of the 12 h before t0 (ve_pre, vn_pre) x t. B3 is not
  available to an operator (it needs the drifter's own track); it shows how much of the skill is plain persistence.
Rows for all windage values and all wind sources are written, so s07 can choose k and the wind source.
Windows without an ECMWF file (seed time 2023-05-03) get no ecmwf/both runs.

Writes data/runs/det/<seed>.parquet (resumable), data/runs/det_errors.parquet (all rows) and results/baseline_inputs.csv.
Usage: python scripts/s06_run_deterministic.py [--benchmark] [--workers 7]
"""
import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import config as C
import runner as R

DET = C.RUNS_DIR / "det"


def load_windows():
    w = pd.read_csv(C.RESULTS / "windows.csv", dtype={"drifter_id": str, "block_id": str}, parse_dates=["t0"])
    w["tag"] = w.t0.map(R.tag_of)
    return w


def rows_for(w, method, current, wind, stokes, windage, lon, lat, status):
    """Long-format rows for elements ordered [windage-major, window-minor]; lon, lat shape (n_elem, 4)."""
    n = len(w)
    nwd = len(windage)
    olon, olat = R.obs_at_leads(w)
    out = []
    for j, wd in enumerate(windage):
        sl = slice(j * n, (j + 1) * n) if lon.shape[0] == n * nwd else slice(0, n)
        e = R.geod_km(olon, olat, lon[sl], lat[sl])
        for k, h in enumerate(R.LEADS):
            out.append(pd.DataFrame(dict(
                window_id=w.window_id.values, method=method, current=current, wind=wind, stokes=stokes, windage=wd,
                lead_h=h, lon=lon[sl][:, k], lat=lat[sl][:, k], err_km=e[:, k],
                status=status[sl] if status is not None else "baseline")))
    return pd.concat(out, ignore_index=True)


def baseline_rows(w, t0, inputs):
    """B0, B1, B3 rows and the drift-vector inputs for the windows w of one seed time."""
    blocks = []
    lon0, lat0 = w.lon_0.values, w.lat_0.values
    zero = np.zeros(len(w))
    # B0
    L = np.repeat(lon0[:, None], 4, 1)
    A = np.repeat(lat0[:, None], 4, 1)
    blocks.append(rows_for(w, "B0", "none", "none", False, [0.0], L, A, None))
    # B3 (persistence of the drifter's own velocity)
    u3, v3 = w.ve_pre.fillna(0.0).values, w.vn_pre.fillna(0.0).values
    lo = np.column_stack([R.advect_constant(lon0, lat0, u3, v3, h)[0] for h in R.LEADS])
    la = np.column_stack([R.advect_constant(lon0, lat0, u3, v3, h)[1] for h in R.LEADS])
    b3 = rows_for(w, "B3", "none", "none", False, [0.0], lo, la, None)
    b3.loc[np.tile(w.ve_pre.isna().values, 4), "err_km"] = np.nan
    blocks.append(b3)
    # B1: current at the datum at t0 (lost: hourly Eulerian 0.494 m; drogued: 6-hourly 15.81 m) + windage x wind at t0
    blk = str(w.block_id.iloc[0])
    for cls in ("lost", "drogued"):
        m = (w.drogue_class == cls).values
        if not m.any():
            continue
        ws = w[m]
        if cls == "lost":
            cu, cv, cf = R.sample_field(C.FORCING_DIR / blk / "cur_notide.nc", "uo", "vo", ws.lon_0.values, ws.lat_0.values, t0)
            cur_name = "cur_notide"
        else:
            cu, cv, cf = R.sample_field(C.FORCING_DIR / blk / "cur_6h_15m.nc", "uo", "vo", ws.lon_0.values, ws.lat_0.values, t0)
            cur_name = "cur_6h_15m"
        inp = pd.DataFrame(dict(window_id=ws.window_id.values, drogue_class=cls, cur_variant=cur_name, cur_u=cu, cur_v=cv,
                                cur_filled=cf))
        winds = C.WINDS if cls == "lost" else ["gfs"]
        for wv in winds:
            if not R.wind_available(t0, wv):
                continue
            wu, wvv, _ = R.sample_field(C.FORCING_DIR / "winds" / f"{R.tag_of(t0)}_{wv}.nc", "u10", "v10",
                                        ws.lon_0.values, ws.lat_0.values, t0)
            inp[f"wind_u_{wv}"], inp[f"wind_v_{wv}"] = wu, wvv
            grid = C.WINDAGE_GRID if cls == "lost" else [0.0]
            for wd in grid:
                u, v = cu + wd * wu, cv + wd * wvv
                lo = np.column_stack([R.advect_constant(ws.lon_0.values, ws.lat_0.values, u, v, h)[0] for h in R.LEADS])
                la = np.column_stack([R.advect_constant(ws.lon_0.values, ws.lat_0.values, u, v, h)[1] for h in R.LEADS])
                blocks.append(rows_for(ws, "B1", cur_name, wv if cls == "lost" else "none", False, [wd], lo, la, None))
        inputs.append(inp)
    return pd.concat(blocks, ignore_index=True)


def run_seed(tag):
    out = DET / f"{tag}.parquet"
    inp_out = DET / f"{tag}_inputs.parquet"
    if out.exists() and inp_out.exists():
        return tag, "cached", 0.0, 0
    t_start = time.time()
    w_all = load_windows()
    w = w_all[w_all.tag == tag].reset_index(drop=True)
    t0 = w.t0.iloc[0]
    block = str(w.block_id.iloc[0])
    res = []
    nrun = 0
    inputs = []
    res.append(baseline_rows(w, t0, inputs))
    wd_grid = np.array(C.WINDAGE_GRID)
    for cls in ("lost", "drogued"):
        ws = w[w.drogue_class == cls].reset_index(drop=True)
        if ws.empty:
            continue
        n = len(ws)
        if cls == "lost":
            combos = [(c, wv, s) for c in C.LOST_CURRENTS for wv in C.WINDS for s in C.STOKES]
            grid = wd_grid
        else:
            combos = [(c, "gfs", False) for c in C.DROGUED_CURRENTS]
            grid = np.array([0.0])
        for c, wv, s in combos:
            if not R.wind_available(t0, wv):
                continue
            lon, lat, st = R.run_oceandrift(t0, block, np.tile(ws.lon_0.values, len(grid)), np.tile(ws.lat_0.values, len(grid)),
                                            c, wv, s, np.repeat(grid, n))
            nrun += 1
            res.append(rows_for(ws, "OD", c, wv if cls == "lost" else "none", s, list(grid), lon, lat, st))
    df = pd.concat(res, ignore_index=True)
    for col in ("method", "current", "wind", "status"):
        df[col] = df[col].astype("category")
    df.to_parquet(out)
    pd.concat(inputs, ignore_index=True).to_parquet(inp_out)
    return tag, "ran", time.time() - t_start, nrun


def benchmark(tag):
    """Time one seed time, check a drogued displacement against the forcing velocity and the reader reuse."""
    w = load_windows()
    ws = w[(w.tag == tag) & (w.drogue_class == "drogued")].reset_index(drop=True)
    t0 = ws.t0.iloc[0]
    block = str(ws.block_id.iloc[0])
    lines = []
    t = time.time()
    lon, lat, st = R.run_oceandrift(t0, block, ws.lon_0.values, ws.lat_0.values, "cur_notide", "gfs", False, 0.0)
    dt1 = time.time() - t
    t = time.time()
    lon2, lat2, _ = R.run_oceandrift(t0, block, ws.lon_0.values, ws.lat_0.values, "cur_notide", "gfs", False, 0.0)
    dt2 = time.time() - t
    R._READERS.clear()
    lon3, lat3, _ = R.run_oceandrift(t0, block, ws.lon_0.values, ws.lat_0.values, "cur_notide", "gfs", False, 0.0)
    lines.append(f"seed {tag}: {len(ws)} drogued drifters, 1 run: first {dt1:.2f} s, second {dt2:.2f} s")
    lines.append(f"reused vs fresh readers: max position difference {np.abs(lon - lon3).max():.2e} deg (lon), {np.abs(lat - lat3).max():.2e} (lat)")
    # displacement check: model 6 h displacement vs forcing velocity at the datum x 6 h
    cu, cv, _ = R.sample_field(C.FORCING_DIR / block / "cur_notide.nc", "uo", "vo", ws.lon_0.values, ws.lat_0.values, t0)
    pu, pv = R.advect_constant(ws.lon_0.values, ws.lat_0.values, cu, cv, 6)
    d_pred = R.geod_km(ws.lon_0.values, ws.lat_0.values, pu, pv)
    d_model = R.geod_km(ws.lon_0.values, ws.lat_0.values, lon[:, 0], lat[:, 0])
    ok = d_pred > 1.0
    lines.append(f"6 h displacement, model / (velocity at datum x 6 h): median {np.median(d_model[ok] / d_pred[ok]):.2f}, "
                 f"correlation {np.corrcoef(d_model, d_pred)[0, 1]:.3f} (n={len(ws)})")
    # full seed time with 1 and 6 windage values
    t = time.time()
    n = len(ws)
    R.run_oceandrift(t0, block, np.tile(ws.lon_0.values, 6), np.tile(ws.lat_0.values, 6), "cur_notide", "gfs", True, np.repeat(C.WINDAGE_GRID, n))
    lines.append(f"one run with {n}x6 elements: {time.time() - t:.2f} s")
    txt = "\n".join(lines)
    print(txt)
    (C.RESULTS / "benchmark.txt").write_text(txt + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", action="store_true")
    ap.add_argument("--workers", type=int, default=7)
    a = ap.parse_args()
    DET.mkdir(parents=True, exist_ok=True)
    w = load_windows()
    tags = sorted(w.tag.unique())
    if a.benchmark:
        benchmark(tags[len(tags) // 2])
        sys.exit(0)
    t_start = time.time()
    done = 0
    with ProcessPoolExecutor(a.workers) as ex:
        for tag, st, secs, nrun in ex.map(run_seed, tags):
            done += 1
            if done % 10 == 0 or done == len(tags):
                print(f"[{done}/{len(tags)}] {tag} {st} {secs:.0f}s {nrun} runs (elapsed {time.time() - t_start:.0f}s)", flush=True)
    det = pd.concat([pd.read_parquet(DET / f"{t}.parquet") for t in tags], ignore_index=True)
    det.to_parquet(C.RUNS_DIR / "det_errors.parquet")
    inp = pd.concat([pd.read_parquet(DET / f"{t}_inputs.parquet") for t in tags], ignore_index=True)
    inp.to_csv(C.RESULTS / "baseline_inputs.csv", index=False)
    print("rows", len(det), "windows", det.window_id.nunique())
    print(det.groupby(["method"], observed=True).size())
    od = det[(det.method == "OD") & (det.lead_h == 48)]
    print("OD elements not active at 48 h:", (od.status != "active").mean().round(4), od.status.value_counts().to_dict())
    print("NaN errors:", int(det.err_km.isna().sum()), "(B3 without a pre-t0 velocity:", int(((det.method == "B3") & det.err_km.isna()).sum()), ")")

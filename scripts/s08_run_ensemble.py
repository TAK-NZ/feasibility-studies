"""Ensemble runs, uncertainty tuning (tuning set only) and the test-set ensemble evaluation.

Part 1  (tune-runs)  Level-1 grid of native OpenDrift uncertainty on the OFFSHORE TUNING windows, 300 particles per drifter,
                     for the model chosen by s07 (per class). Level-1 keys: drift:current_uncertainty, drift:wind_uncertainty
                     (lost class), per-element windage spread (lost class), per-element current_drift_factor spread (drogued).
Part 2  (tune-eval)  Level 2 = post hoc persistent velocity offset b ~ N(0, sigma_b^2) per component per particle (outside OpenDrift),
                     displacement b x t added at every lead. Contour hit rates of the 50 and 95 percent highest-density contours
                     (drift_lib.hdr_levels) for every Level-1 x Level-2 combination; selection rule D8: the smallest mean 95 percent
                     contour area at 24 h among settings whose 24 h 95 percent hit rate is at least 0.90 (config.HIT_TARGET_TUNE);
                     if none reaches it, the highest hit rate. Area is computed on a fixed subsample of tuning windows.
                     Result written to results/tuning_choices.json and appended to work/tuning.md.
Part 3  (test-runs)  The chosen setting, 500 particles per drifter, on all TEST windows (all zones), plus the Leeway reference
                     (lost class, object class 1 = PIW unknown state, chosen current and wind, native uncertainty only).
Part 4               results/ens_metrics.csv, one row per window x method x lead (methods: od_final, od_native, leeway_native,
                     leeway_offset, B1_circle) with centroid error, 50/95 hit flags and areas. Raw particle displacements are in
                     data/runs/ens/test_<seed>.npz.
Seeds are fixed (config.SIM_SEED); the offsets use a per-window generator so they do not depend on the order of the runs.
Usage: python scripts/s08_run_ensemble.py [--workers 7] [--part tune-runs|tune-eval|test|all]
"""
import argparse
import json
import time
import zlib
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import config as C
import runner as R
from drift_lib import hdr_areas, hdr_levels, kde_density, local_km

ENS = C.RUNS_DIR / "ens"
LEADS = list(C.LEADS_H)
KM_PER_MS_H = 3.6          # km travelled in 1 h at 1 m/s


def crc(s):
    return zlib.crc32(str(s).encode()) & 0x7FFFFFFF


def load_windows():
    w = pd.read_csv(C.RESULTS / "windows.csv", dtype={"drifter_id": str, "block_id": str}, parse_dates=["t0"])
    w["tag"] = w.t0.map(R.tag_of)
    return w


def choices():
    return json.loads((C.RESULTS / "tuning_choices.json").read_text())


def l1_settings(cls):
    """Level-1 grid as a list of dicts with an id string."""
    out = []
    if cls == "lost":
        for cu in C.CUR_UNC_GRID:
            for wu in C.WIND_UNC_GRID:
                for rw in C.WINDAGE_REL_SIGMA_GRID:
                    out.append(dict(id=f"cu{cu}_wu{wu}_rw{rw}", cur_unc=cu, wind_unc=wu, rel_windage=rw, rel_cdf=0.0))
    else:
        for cu in C.CUR_UNC_GRID:
            for rc in C.CDF_REL_SIGMA_GRID:
                out.append(dict(id=f"cu{cu}_rc{rc}", cur_unc=cu, wind_unc=0.0, rel_windage=0.0, rel_cdf=rc))
    return out


def model_of(cls, ch):
    m = ch["choices"][cls]
    if cls == "lost":
        return dict(current=m["current"], wind=m["wind"], stokes=m["stokes"], windage=m["windage"])
    return dict(current=m["current"], wind="gfs", stokes=False, windage=0.0)


def run_ensemble(ws, cls, model, s, n_p, seed):
    """Run all windows ws (one seed time) with n_p particles each. Returns x, y (n, n_p, 4) in km relative to the datum
    and the status code (0 active, 1 stranded, 2 other) at 48 h (n, n_p)."""
    n = len(ws)
    rng = np.random.default_rng(seed)
    lon = np.repeat(ws.lon_0.values, n_p)
    lat = np.repeat(ws.lat_0.values, n_p)
    wd = model["windage"] * np.clip(1 + s["rel_windage"] * rng.standard_normal(n * n_p), 0, None) if cls == "lost" \
        else np.zeros(n * n_p)
    cdf = np.clip(1 + s["rel_cdf"] * rng.standard_normal(n * n_p), 0, None) if s["rel_cdf"] > 0 else None
    t0 = ws.t0.iloc[0]
    blk = str(ws.block_id.iloc[0])
    lo, la, st = R.run_oceandrift(t0, blk, lon, lat, model["current"], model["wind"], model["stokes"], wd, cdf,
                                  s["cur_unc"], s["wind_unc"], seed=seed)
    return to_local(ws, lo, la, st, n_p)


def to_local(ws, lo, la, st, n_p):
    n = len(ws)
    lo = lo.reshape(n, n_p, 4)
    la = la.reshape(n, n_p, 4)
    x, y = local_km(lo, la, ws.lon_0.values[:, None, None], ws.lat_0.values[:, None, None])
    code = np.where(st == "active", 0, np.where(st == "stranded", 1, 2)).astype(np.int8).reshape(n, n_p)
    return x.astype(np.float32), y.astype(np.float32), code


def obs_local(ws):
    olon, olat = R.obs_at_leads(ws)
    x, y = local_km(olon, olat, ws.lon_0.values[:, None], ws.lat_0.values[:, None])
    return x, y


def offsets(window_id, n_p):
    """Standard normal draws (n_p, 2) for the persistent velocity offset of a window, independent of the run order."""
    return np.random.default_rng([C.SIM_SEED, crc(window_id)]).standard_normal((n_p, 2))


def displaced(x, y, z, sigma_b, lead_h):
    """x, y (n_p,) displacements at the lead plus the offset z*sigma_b (m/s) times the lead."""
    k = sigma_b * KM_PER_MS_H * lead_h
    return np.vstack([x + z[:, 0] * k, y + z[:, 1] * k])


# ------------------------------------------------------------------ part 1
def tune_runs_seed(tag):
    path = ENS / f"tune_{tag}.npz"
    if path.exists():
        return tag, "cached", 0
    ch = choices()
    w = load_windows()
    w = w[(w.tag == tag) & (w.split == "tuning") & (w.zone == "offshore")]
    out, nrun = {}, 0
    for cls in ("lost", "drogued"):
        ws = w[w.drogue_class == cls].reset_index(drop=True)
        model = model_of(cls, ch)
        if ws.empty or (cls == "lost" and not R.wind_available(ws.t0.iloc[0], model["wind"])):
            continue
        out[f"{cls}|wid"] = np.array(ws.window_id.tolist(), dtype="U40")
        ox, oy = obs_local(ws)
        out[f"{cls}|obs_x"], out[f"{cls}|obs_y"] = ox.astype(np.float32), oy.astype(np.float32)
        for s in l1_settings(cls):
            x, y, code = run_ensemble(ws, cls, model, s, C.N_PER_DRIFTER_TUNE, C.SIM_SEED + crc(tag + s["id"]) % 100000)
            out[f"{cls}|{s['id']}|x"], out[f"{cls}|{s['id']}|y"], out[f"{cls}|{s['id']}|code"] = x, y, code
            nrun += 1
    np.savez_compressed(path, **out)
    return tag, "ran", nrun


# ------------------------------------------------------------------ part 2
def area_subset(n, k=C.N_AREA_SUBSAMPLE):
    step = max(1, int(np.ceil(n / k)))
    return step


def tune_eval_seed(tag):
    d = np.load(ENS / f"tune_{tag}.npz")
    rows = []
    for cls in ("lost", "drogued"):
        if f"{cls}|wid" not in d:
            continue
        wids = d[f"{cls}|wid"]
        ox, oy = d[f"{cls}|obs_x"], d[f"{cls}|obs_y"]
        for s in l1_settings(cls):
            X, Y = d[f"{cls}|{s['id']}|x"], d[f"{cls}|{s['id']}|y"]
            for i, wid in enumerate(wids):
                z = offsets(wid, X.shape[1])
                for sb in C.SIGMA_B_GRID:
                    for k, lead in enumerate(LEADS):
                        smp = displaced(X[i, :, k], Y[i, :, k], z, sb, lead)
                        H, lev = hdr_levels(smp)
                        dens = float(kde_density(smp, np.array([[ox[i, k]], [oy[i, k]]]), H)[0])
                        row = dict(cls=cls, l1=s["id"], sigma_b=sb, window_id=wid, lead_h=lead,
                                   hit50=dens >= lev[0.5], hit95=dens >= lev[0.95], area95=np.nan, area50=np.nan)
                        rows.append(row)
    return pd.DataFrame(rows)


def tune_area_seed(args):
    """Area of the 50/95 percent contours at 24 h for a fixed subsample of windows (list of window ids)."""
    tag, sub = args
    d = np.load(ENS / f"tune_{tag}.npz")
    rows = []
    k = LEADS.index(24)
    for cls in ("lost", "drogued"):
        if f"{cls}|wid" not in d:
            continue
        wids = list(d[f"{cls}|wid"])
        for s in l1_settings(cls):
            X, Y = d[f"{cls}|{s['id']}|x"], d[f"{cls}|{s['id']}|y"]
            for i, wid in enumerate(wids):
                if wid not in sub:
                    continue
                z = offsets(wid, X.shape[1])
                for sb in C.SIGMA_B_GRID:
                    a = hdr_areas(displaced(X[i, :, k], Y[i, :, k], z, sb, 24))
                    rows.append(dict(cls=cls, l1=s["id"], sigma_b=sb, window_id=wid, lead_h=24, area50=a[0.5], area95=a[0.95]))
    return pd.DataFrame(rows)


def tune_eval(workers):
    tags = sorted(p.stem[5:] for p in ENS.glob("tune_*.npz"))
    t = time.time()
    with ProcessPoolExecutor(workers) as ex:
        hits = pd.concat(ex.map(tune_eval_seed, tags), ignore_index=True)
    print(f"hit rates: {len(hits)} rows in {time.time() - t:.0f}s")
    w = load_windows()
    sub = set()
    for cls in ("lost", "drogued"):
        ids = sorted(hits[hits.cls == cls].window_id.unique())
        sub |= set(ids[::area_subset(len(ids))])
    t = time.time()
    with ProcessPoolExecutor(workers) as ex:
        areas = pd.concat(ex.map(tune_area_seed, [(tg, sub) for tg in tags]), ignore_index=True)
    print(f"areas: {len(areas)} rows in {time.time() - t:.0f}s")
    hits = hits.drop(columns=["area95", "area50"]).merge(areas, on=["cls", "l1", "sigma_b", "window_id", "lead_h"], how="left")
    hits.to_parquet(ENS / "tune_hits.parquet")
    g = hits.groupby(["cls", "l1", "sigma_b", "lead_h"]).agg(
        n_windows=("hit95", "size"), hit50=("hit50", "mean"), hit95=("hit95", "mean"),
        mean_area50_km2=("area50", "mean"), mean_area95_km2=("area95", "mean"), n_area=("area95", "count")).reset_index()
    g.to_csv(C.RESULTS / "tables_tuning_uncertainty.csv", index=False)
    ch = choices()
    ch["uncertainty"] = {}
    lines = ["", "## Part B: uncertainty tuning (tuning set only, offshore, 24 h)", "",
             f"Particles per drifter: {C.N_PER_DRIFTER_TUNE}. Level 1 = native OpenDrift keys; Level 2 = post hoc persistent velocity offset "
             f"(sigma_b in {C.SIGMA_B_GRID} m/s per component). Selection rule: smallest mean 95% contour area at 24 h among settings with a 24 h "
             f"95% hit rate of at least {C.HIT_TARGET_TUNE}; if none reaches it, the highest hit rate. Areas use a fixed subsample of tuning windows "
             f"(every k-th window of the sorted ids, about {C.N_AREA_SUBSAMPLE} per class; the count is `n_area` in results/tables_tuning_uncertainty.csv).", ""]
    for cls in ("lost", "drogued"):
        s24 = g[(g.cls == cls) & (g.lead_h == 24)].copy()
        ok = s24[s24.hit95 >= C.HIT_TARGET_TUNE]
        if len(ok):
            pick = ok.sort_values(["mean_area95_km2", "sigma_b"]).iloc[0]
            rule = "smallest area with hit95 >= target"
        else:
            pick = s24.sort_values("hit95", ascending=False).iloc[0]
            rule = "no setting reached the target: highest hit95"
        nat = s24[s24.sigma_b == 0.0].sort_values("hit95", ascending=False).iloc[0]
        s_l1 = [x for x in l1_settings(cls) if x["id"] == pick.l1][0]
        ch["uncertainty"][cls] = dict(rule=rule, l1=s_l1, sigma_b=float(pick.sigma_b), tuning_hit50_24h=float(pick.hit50),
                                      tuning_hit95_24h=float(pick.hit95), tuning_area95_24h_km2=float(pick.mean_area95_km2),
                                      tuning_area50_24h_km2=float(pick.mean_area50_km2), n_windows=int(pick.n_windows),
                                      native_only_best=dict(l1=nat.l1, hit50=float(nat.hit50), hit95=float(nat.hit95),
                                                            area95=float(nat.mean_area95_km2)))
        u = ch["uncertainty"][cls]
        lines += [f"### {cls}: {len(s24[['l1']].drop_duplicates())} Level-1 settings x {len(C.SIGMA_B_GRID)} offsets, {int(pick.n_windows)} windows",
                  f"Selected ({rule}): Level 1 `{pick.l1}`, sigma_b = {pick.sigma_b} m/s; tuning-set 24 h hit rate 50%: {pick.hit50:.3f}, 95%: {pick.hit95:.3f}, "
                  f"mean 95% area {pick.mean_area95_km2:.0f} km2.",
                  f"Native only (sigma_b = 0), best Level-1 setting `{nat.l1}`: hit 50%: {nat.hit50:.3f}, 95%: {nat.hit95:.3f}.", "",
                  "Hit rates at 24 h by sigma_b for the selected Level-1 setting:", "", "| sigma_b | hit 50% | hit 95% | mean 95% area km2 |", "|---|---|---|---|"]
        for r in s24[s24.l1 == pick.l1].sort_values("sigma_b").itertuples():
            lines.append(f"| {r.sigma_b} | {r.hit50:.3f} | {r.hit95:.3f} | {r.mean_area95_km2:.0f} |")
        lines.append("")
    ts = ch["tuning_set"]
    lines += ["## Part C: what was tuned, and on what data", "",
              f"Data for every tuned quantity below: the tuning split only (drifters assigned by `s02`, 40 percent of drifters), offshore windows "
              f"(more than 50 km from land): lost class {ts['lost']['windows']} windows from {ts['lost']['drifters']} drifters "
              f"({ts['lost_common_wind_windows']} of them have both wind sources), drogued class {ts['drogued']['windows']} windows from "
              f"{ts['drogued']['drifters']} drifters. Test-set windows and drifters were not used for any of it.", "",
              "| quantity | how chosen | metric on the tuning set |", "|---|---|---|",
              "| lost class: windage fraction of the OceanDrift model | grid `config.WINDAGE_GRID`, per variant | 24 h median separation error |",
              "| lost class: current, wind source, Stokes on/off | best of 18 variants; ties within 2 percent resolved by the 48 h median | 24 h median separation error |",
              "| drogued class: current variant | best of 5; same tie rule | 24 h median separation error |",
              "| baseline: windage k, wind source (same as the model), circle slopes a_0.5 and a_0.95 | grid / quantile of the 24 h baseline error | 24 h median, 24 h quantiles |",
              f"| uncertainty: Level-1 settings and sigma_b (post hoc offset), per class | rule above | 24 h 95 percent hit rate and 95 percent area |", "",
              "Fixed by hand, not tuned: time step 900 s, output step 6 h, Runge-Kutta advection, 300 (tuning) and 500 (test) particles per drifter, "
              "the 3 km kernel floor of the contour estimator, the grids themselves. The sigma_b grid was refined once (0.125 m/s added) after the first "
              "tuning pass showed the 0.1 to 0.15 step was coarse (tuning data only, before any test-set run). Leeway was not tuned (object class 1).", ""]
    (C.RESULTS / "tuning_choices.json").write_text(json.dumps(ch, indent=2, default=str))
    txt = (C.WORK / "tuning.md").read_text()
    txt = txt.split("\n## Part B")[0].rstrip("\n") + "\n" + "\n".join(lines) + "\n"
    (C.WORK / "tuning.md").write_text(txt)
    print("\n".join(lines))


# ------------------------------------------------------------------ part 3 and 4
def contour_rows(ws, i, method, X, Y, code, ox, oy, sigma_b, with_area=True):
    rows = []
    wid = ws.window_id.iloc[i]
    z = offsets(wid, X.shape[1])
    for k, lead in enumerate(LEADS):
        smp = displaced(X[i, :, k], Y[i, :, k], z, sigma_b, lead)
        H, lev = hdr_levels(smp)
        dens = float(kde_density(smp, np.array([[ox[i, k]], [oy[i, k]]]), H)[0])
        a = hdr_areas(smp) if with_area else {0.5: np.nan, 0.95: np.nan}
        cx, cy = smp.mean(1)
        rows.append(dict(window_id=wid, method=method, lead_h=lead, obs_x=float(ox[i, k]), obs_y=float(oy[i, k]), cx=float(cx), cy=float(cy),
                         centroid_err_km=float(np.hypot(cx - ox[i, k], cy - oy[i, k])), hit50=bool(dens >= lev[0.5]),
                         hit95=bool(dens >= lev[0.95]), area50_km2=a[0.5], area95_km2=a[0.95],
                         n_particles=X.shape[1], frac_active=float((code[i] == 0).mean()), frac_stranded=float((code[i] == 1).mean()),
                         frac_other_deactivated=float((code[i] == 2).mean()), sigma_b=sigma_b))
    return rows


def test_seed(tag):
    path = ENS / f"test_{tag}.parquet"
    if path.exists():
        return tag, "cached", 0
    ch = choices()
    w = load_windows()
    w = w[(w.tag == tag) & (w.split == "test")]
    det = None
    rows, particles, nrun = [], {}, 0
    for cls in ("lost", "drogued"):
        ws = w[w.drogue_class == cls].reset_index(drop=True)
        model = model_of(cls, ch)
        if ws.empty or (cls == "lost" and not R.wind_available(ws.t0.iloc[0], model["wind"])):
            continue
        u = ch["uncertainty"][cls]
        ox, oy = obs_local(ws)
        x, y, code = run_ensemble(ws, cls, model, u["l1"], C.N_PER_DRIFTER_TEST, C.SIM_SEED + crc(tag + "test" + cls) % 100000)
        nrun += 1
        particles[f"{cls}|wid"] = np.array(ws.window_id.tolist(), dtype="U40")
        particles[f"{cls}|od|x"], particles[f"{cls}|od|y"], particles[f"{cls}|od|code"] = x, y, code
        for i in range(len(ws)):
            rows += contour_rows(ws, i, "od_final", x, y, code, ox, oy, u["sigma_b"])
            rows += contour_rows(ws, i, "od_native", x, y, code, ox, oy, 0.0, with_area=True)
        if cls == "lost":
            lo, la, st = R.run_leeway(ws.t0.iloc[0], str(ws.block_id.iloc[0]), np.repeat(ws.lon_0.values, C.N_PER_DRIFTER_TEST),
                                      np.repeat(ws.lat_0.values, C.N_PER_DRIFTER_TEST), model["current"], model["wind"], 1,
                                      u["l1"]["cur_unc"], u["l1"]["wind_unc"], seed=C.SIM_SEED + crc(tag + "leeway") % 100000)
            nrun += 1
            lx, ly, lcode = to_local(ws, lo, la, st, C.N_PER_DRIFTER_TEST)
            particles["lost|leeway|x"], particles["lost|leeway|y"], particles["lost|leeway|code"] = lx, ly, lcode
            for i in range(len(ws)):
                rows += contour_rows(ws, i, "leeway_native", lx, ly, lcode, ox, oy, 0.0)
                rows += contour_rows(ws, i, "leeway_offset", lx, ly, lcode, ox, oy, u["sigma_b"])
    np.savez_compressed(ENS / f"test_{tag}.npz", **particles)
    pd.DataFrame(rows).to_parquet(path)
    return tag, "ran", nrun


def baseline_circle_rows(ch):
    """Baseline expanding-circle hits and areas on the test windows (radius = a_p km/h x lead)."""
    w = load_windows()
    w = w[w.split == "test"]
    det = pd.read_parquet(C.RUNS_DIR / "det_errors.parquet", filters=[("method", "==", "B1")])
    for c in ("wind", "method"):
        det[c] = det[c].astype(str)
    out = []
    for cls in ("lost", "drogued"):
        b = ch["baseline"][cls]
        d = det[(det.wind == b["wind"]) & np.isclose(det.windage, b["windage"]) & det.window_id.isin(w[w.drogue_class == cls].window_id)]
        d = d.merge(w[["window_id", "lon_0", "lat_0"]], on="window_id")
        for lead in LEADS:
            dl = d[d.lead_h == lead].copy()
            ww = w.set_index("window_id").loc[dl.window_id]
            ox, oy = local_km(ww[f"lon_{C.LEAD_IDX[lead]}"].values, ww[f"lat_{C.LEAD_IDX[lead]}"].values, ww.lon_0.values, ww.lat_0.values)
            bx, by = local_km(dl.lon.values, dl.lat.values, dl.lon_0.values, dl.lat_0.values)
            dist = np.hypot(ox - bx, oy - by)
            out.append(pd.DataFrame(dict(
                window_id=dl.window_id.values, method="B1_circle", lead_h=lead, obs_x=ox, obs_y=oy, cx=bx, cy=by, centroid_err_km=dist,
                hit50=dist <= b["a50_km_per_h"] * lead, hit95=dist <= b["a95_km_per_h"] * lead,
                area50_km2=np.pi * (b["a50_km_per_h"] * lead) ** 2, area95_km2=np.pi * (b["a95_km_per_h"] * lead) ** 2,
                n_particles=0, frac_active=np.nan, frac_stranded=np.nan, frac_other_deactivated=np.nan, sigma_b=np.nan)))
    return pd.concat(out, ignore_index=True)


def test_all(workers):
    w = load_windows()
    tags = sorted(w[w.split == "test"].tag.unique())
    t = time.time()
    done = 0
    with ProcessPoolExecutor(workers) as ex:
        for tag, st, n in ex.map(test_seed, tags):
            done += 1
            if done % 10 == 0 or done == len(tags):
                print(f"[{done}/{len(tags)}] {tag} {st} ({time.time() - t:.0f}s)", flush=True)
    ens = pd.concat([pd.read_parquet(ENS / f"test_{tg}.parquet") for tg in tags if (ENS / f"test_{tg}.parquet").exists()], ignore_index=True)
    ch = choices()
    ens = pd.concat([ens, baseline_circle_rows(ch)], ignore_index=True)
    meta = w[["window_id", "drifter_id", "split", "stratum", "zone", "drogue_class", "season", "wind_class", "t0", "has_ecmwf"]]
    ens = ens.merge(meta, on="window_id", how="left")
    ens.to_csv(C.RESULTS / "ens_metrics.csv", index=False, float_format="%.5g")
    print("ens_metrics rows", len(ens))
    print(ens.groupby(["method", "drogue_class"]).window_id.nunique().unstack())
    off = ens[(ens.zone == "offshore") & (ens.lead_h == 24)]
    print("offshore test, 24 h hit rates and areas:")
    print(off.groupby(["drogue_class", "method"]).agg(n=("window_id", "nunique"), hit50=("hit50", "mean"), hit95=("hit95", "mean"),
                                                       area95=("area95_km2", "median"), centroid_err_median=("centroid_err_km", "median")).round(3))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--part", default="all", choices=["tune-runs", "tune-eval", "test", "all"])
    a = ap.parse_args()
    ENS.mkdir(parents=True, exist_ok=True)
    w = load_windows()
    if a.part in ("tune-runs", "all"):
        tags = sorted(w[(w.split == "tuning") & (w.zone == "offshore")].tag.unique())
        t = time.time()
        done = 0
        with ProcessPoolExecutor(a.workers) as ex:
            for tag, st, n in ex.map(tune_runs_seed, tags):
                done += 1
                if done % 10 == 0 or done == len(tags):
                    print(f"tune-runs [{done}/{len(tags)}] {tag} {st} ({time.time() - t:.0f}s)", flush=True)
    if a.part in ("tune-eval", "all"):
        tune_eval(a.workers)
    if a.part in ("test", "all"):
        test_all(a.workers)

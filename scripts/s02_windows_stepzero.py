"""Step zero: count usable drifters, tracks and 48 h windows at every filter step, then
choose the 96 h time blocks and the tuning/test split used by all later steps.

Reads data/gdp/gdp_6h_box.csv. Writes
  results/step0_counts.csv, work/step0_counts.md,
  results/blocks.csv, results/windows_selected.csv, results/drifters_split.csv
Usage: python scripts/s02_windows_stepzero.py
"""
import sys
import warnings

import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view as swv

import config as C
from drift_lib import LandDistance, geod_km

H = pd.Timedelta(hours=1)
N = C.N_OBS
rows_out = []   # funnel rows


def add(step, stratum, n_drifters=None, n_tracks=None, n_windows=None, n_days=None, note=""):
    rows_out.append(dict(step=step, stratum=stratum, n_drifters=n_drifters, n_tracks=n_tracks,
                         n_windows=n_windows, n_drifter_days=n_days, note=note))
    print(f"{step:<34} {stratum:<34} drifters={n_drifters} tracks={n_tracks} windows={n_windows} days={n_days}")


def load():
    df = pd.read_csv(C.GDP_CSV, dtype={"ID": str}, low_memory=False)
    for c in ("time", "drogue_lost_date", "start_date", "end_date"):
        df[c] = pd.to_datetime(df[c], utc=True).dt.tz_localize(None)
    for v in ("ve", "vn"):
        df.loc[df[v] <= -999990, v] = np.nan
    return df


def tracks_of(g):
    """Contiguous 6-hourly runs. Returns list of (start_idx, n_points) on the 6 h index."""
    idx = ((g.time - g.time.iloc[0]) / (6 * H)).round().astype(int).values
    brk = np.where(np.diff(idx) != 1)[0] + 1
    runs = np.split(np.arange(len(idx)), brk)
    return [len(r) for r in runs]


def count_tracks(df, min_pts=N):
    n = 0
    for _, g in df.groupby("ID"):
        n += sum(1 for L in tracks_of(g) if L >= min_pts)
    return n


def build_windows(df, ld):
    """All complete 9-point windows (6 h stride) for every drifter, with positions and QC fields."""
    recs = []
    pos_lon, pos_lat = [], []
    for did, g in df.groupby("ID"):
        g = g.sort_values("time")
        t_ref = g.time.iloc[0]
        idx = ((g.time - t_ref) / (6 * H)).round().astype(int).values
        n = idx[-1] + 1
        lon = np.full(n, np.nan); lat = np.full(n, np.nan); dl = np.full(n, np.nan)
        ve = np.full(n, np.nan); vn = np.full(n, np.nan)
        lon[idx] = g.longitude.values; lat[idx] = g.latitude.values
        dl[idx] = g.dland.values; ve[idx] = g.ve.values; vn[idx] = g.vn.values
        valid = ~np.isnan(lon)
        if n < N:
            continue
        ok = swv(valid, N).all(axis=1)
        starts = np.where(ok)[0]
        if len(starts) == 0:
            continue
        step_km = np.full(n, np.nan)
        step_km[:-1] = geod_km(lon[:-1], lat[:-1], lon[1:], lat[1:])
        ii = starts[:, None] + np.arange(N)[None, :]
        smax = swv(step_km, N - 1).max(axis=1)[starts]
        dmin = swv(dl, N).min(axis=1)[starts]
        # mean drifter velocity over the 12 h before t0 (points t0-12h, t0-6h, t0), NaN-aware
        pre = np.clip(starts[:, None] + np.arange(-2, 1)[None, :], 0, n - 1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN slices give NaN, as intended
            pve = np.nanmean(np.where(starts[:, None] - 2 >= 0, ve[pre], np.nan), axis=1)
            pvn = np.nanmean(np.where(starts[:, None] - 2 >= 0, vn[pre], np.nan), axis=1)
        for k, s in enumerate(starts):
            recs.append((did, t_ref + s * 6 * H, smax[k], dmin[k], pve[k], pvn[k]))
        pos_lon.append(lon[ii]); pos_lat.append(lat[ii])
    w = pd.DataFrame(recs, columns=["drifter_id", "t0", "max_step_km", "min_land_km", "ve_pre", "vn_pre"])
    LON = np.vstack(pos_lon); LAT = np.vstack(pos_lat)
    for k in range(N):
        w[f"lon_{k}"] = LON[:, k]; w[f"lat_{k}"] = LAT[:, k]
    return w


def season(month):
    return {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA",
            9: "SON", 10: "SON", 11: "SON"}[month]


def summarize(w, label, by=("drogue_class", "zone")):
    """Add class x zone rows to the funnel."""
    for (c, z), g in w.groupby(list(by)):
        add(label, f"{c} / {z}", g.drifter_id.nunique(), None, len(g))


def main():
    df = load()
    sd, ed = pd.Timestamp(C.SEED_START), pd.Timestamp(C.SEED_END)
    add("1 wide box + time (raw)", "all", df.ID.nunique(), count_tracks(df), len(df), round(len(df) / 4))
    print("buoy types:", df.drop_duplicates("ID").typebuoy.value_counts().to_dict())

    # QC: interpolation error, end_date
    q = df[(df.err_lat < C.ERR_DEG_MAX) & (df.err_lon < C.ERR_DEG_MAX)]
    add("2 position error < 0.05 deg", "all", q.ID.nunique(), count_tracks(q), len(q), round(len(q) / 4))
    q = q[q.end_date.isna() | (q.time <= q.end_date)]
    add("3 time <= end_date (last good)", "all", q.ID.nunique(), count_tracks(q), len(q), round(len(q) / 4))
    q = q[(q.lon360 >= C.OBS_LON[0]) & (q.lon360 <= C.OBS_LON[1]) &
          (q.latitude >= C.OBS_LAT[0]) & (q.latitude <= C.OBS_LAT[1])]
    add("4 observation box 160-185E 55-25S", "all", q.ID.nunique(), count_tracks(q), len(q), round(len(q) / 4))
    t_hi = ed + pd.Timedelta(hours=C.WINDOW_HOURS)
    q_box = q
    q = q[(q.time >= sd) & (q.time <= t_hi)].copy()
    add("5 seed period 2023-01-21..2025-06-16", "all", q.ID.nunique(), count_tracks(q), len(q), round(len(q) / 4))

    ld = LandDistance()
    print(f"land raster: land fraction {ld.land_fraction:.3f}, boundary cells {ld.n_boundary}")
    # rows from 12 h before the first seed are kept so the pre-t0 drifter velocity exists for early windows
    qw = q_box[(q_box.time >= sd - 12 * H) & (q_box.time <= t_hi)].copy()
    qw["dland"] = ld.km(qw.lon360.values, qw.latitude.values)
    q = qw

    w = build_windows(q, ld)
    w = w[(w.t0 >= sd) & (w.t0 <= ed)].reset_index(drop=True)
    add("6 candidate 48 h windows (6 h stride)", "all", w.drifter_id.nunique(), None, len(w))
    w = w[w.max_step_km <= C.MAX_SPEED_MS * 6 * 3.6].reset_index(drop=True)
    add("7 after max 2.5 m/s speed check", "all", w.drifter_id.nunique(), None, len(w))

    # drogue class
    dl = q.drop_duplicates("ID").set_index("ID").drogue_lost_date
    assert q.groupby("ID").drogue_lost_date.nunique(dropna=False).max() == 1
    tl = w.drifter_id.map(dl)
    t0 = w.t0; t1 = t0 + C.WINDOW_HOURS * H; buf = C.DROGUE_BUFFER_H * H
    cls = np.where(tl.isna(), "drogued",
          np.where(tl == pd.Timestamp("1970-01-01"), "excl_uncertain",
          np.where(tl <= t0 - buf, "lost",
          np.where(tl >= t1 + buf, "drogued", "excl_straddle"))))
    w["drogue_class"] = cls
    w["drogue_lost_date"] = tl.values
    for c, g in w.groupby("drogue_class"):
        add("8 drogue state", c, g.drifter_id.nunique(), None, len(g))
    w = w[w.drogue_class.isin(["drogued", "lost"])].reset_index(drop=True)
    w["zone"] = np.where(w.min_land_km > C.OFFSHORE_KM, "offshore",
                np.where(w.min_land_km < C.COASTAL_KM, "coastal", "mid"))
    summarize(w, "9 class x zone (6 h stride)")
    w["season"] = w.t0.dt.month.map(season)

    # non-overlapping (48 h stride) and 00Z-only counts for comparison
    w = w.sort_values(["drifter_id", "t0"]).reset_index(drop=True)
    keep = np.zeros(len(w), bool)
    for did, g in w.groupby("drifter_id"):
        last = None
        for i, t in zip(g.index, g.t0):
            if last is None or t >= last + C.WINDOW_HOURS * H:
                keep[i] = True; last = t
    w["nonoverlap"] = keep
    summarize(w[w.nonoverlap], "10 non-overlapping (48 h stride)")
    w["is00z"] = w.t0.dt.hour == 0
    summarize(w[w.is00z], "11 seed at 00Z only")

    # ----- blocks (D4)
    w["id_t0"] = w.t0.dt.strftime("%Y%m%dT%H")
    usable = w[w.is00z & w.zone.isin(["offshore", "coastal"])]
    bs_all = pd.date_range(sd, ed, freq="D")
    cnt = usable.groupby("t0").size()
    cnt_c = usable[usable.zone == "coastal"].groupby("t0").size()
    cand = pd.DataFrame({"bs": bs_all})
    cand["n_usable"] = [int(cnt.get(b, 0)) + int(cnt.get(b + 48 * H, 0)) for b in bs_all]
    cand["n_coastal"] = [int(cnt_c.get(b, 0)) + int(cnt_c.get(b + 48 * H, 0)) for b in bs_all]
    # the block must fit in the data: second window ends at bs + 96 h
    cand = cand[cand.bs + C.BLOCK_HOURS * H <= pd.Timestamp(C.GDP_TIME_MAX.rstrip("Z"))]
    rng = np.random.default_rng(C.RNG_SEED)
    chosen = []   # (bs, kind)

    def free(b):
        return all(abs((b - c[0]) / H) >= C.BLOCK_HOURS for c in chosen)

    short_months = []
    for ym, g in cand.groupby(cand.bs.dt.to_period("M")):
        elig = g[g.n_usable >= C.MIN_WINDOWS_PER_BLOCK].bs.tolist()
        rng.shuffle(elig)
        got = 0
        for b in elig:
            if got >= C.BLOCKS_PER_MONTH:
                break
            if free(b):
                chosen.append((b, "random")); got += 1
        if got < C.BLOCKS_PER_MONTH:
            short_months.append((str(ym), got))
    n_random = len(chosen)
    for _, r in cand.sort_values(["n_coastal", "bs"], ascending=[False, True]).iterrows():
        if sum(1 for c in chosen if c[1] == "coastal") >= C.COASTAL_EXTRA_BLOCKS or r.n_coastal < 1:
            break
        if free(r.bs):
            chosen.append((r.bs, "coastal"))
    blocks = pd.DataFrame(chosen, columns=["block_start", "kind"]).sort_values("block_start").reset_index(drop=True)
    blocks["block_id"] = blocks.block_start.dt.strftime("%Y%m%d")
    blocks["block_end"] = blocks.block_start + C.BLOCK_HOURS * H
    bmap = {}
    for _, b in blocks.iterrows():
        for off in C.SEED_OFFSETS_H:
            bmap[b.block_start + off * H] = b.block_id
    sel = w[w.is00z & w.t0.isin(bmap.keys())].copy()
    sel["block_id"] = sel.t0.map(bmap)
    blocks["n_windows"] = blocks.block_id.map(sel.groupby("block_id").size()).fillna(0).astype(int)
    blocks["n_offshore"] = blocks.block_id.map(sel[sel.zone == "offshore"].groupby("block_id").size()).fillna(0).astype(int)
    blocks["n_coastal"] = blocks.block_id.map(sel[sel.zone == "coastal"].groupby("block_id").size()).fillna(0).astype(int)
    add("12 blocks chosen", f"{len(blocks)} blocks = {n_random} random + {len(blocks) - n_random} coastal-targeted",
        None, None, None,
        note=f"months with fewer than {C.BLOCKS_PER_MONTH} eligible non-overlapping blocks: {short_months}")
    summarize(sel, "13 selected windows (blocks)")

    # ----- tuning / test split by drifter (D7)
    dr = sel.groupby("drifter_id").agg(
        has_coastal=("zone", lambda s: (s == "coastal").any()),
        has_drogued=("drogue_class", lambda s: (s == "drogued").any()),
        n_windows=("t0", "size")).reset_index()
    dr["stratum"] = np.where(dr.has_coastal, "coastal", "offshore") + "_" + np.where(dr.has_drogued, "drogued", "lost")
    rs = np.random.default_rng(C.RNG_SEED)
    dr["split"] = "test"
    for s, g in dr.groupby("stratum"):
        ids = np.array(g.drifter_id.tolist())
        rs.shuffle(ids)
        n_tune = int(round(C.TUNING_FRACTION * len(ids))) if len(ids) >= 2 else 0
        dr.loc[dr.drifter_id.isin(ids[:n_tune]), "split"] = "tuning"
    sel["split"] = sel.drifter_id.map(dr.set_index("drifter_id").split)
    for sp in ("tuning", "test"):
        for (c, z), g in sel[sel.split == sp].groupby(["drogue_class", "zone"]):
            add(f"14 {sp} set", f"{c} / {z}", g.drifter_id.nunique(), None, len(g))

    off = sel[sel.zone == "offshore"]
    for (sp, s), g in off.groupby(["split", "season"]):
        add("15 offshore by season", f"{sp} / {s}", g.drifter_id.nunique(), None, len(g))

    # ----- outputs
    sel["window_id"] = sel.drifter_id + "_" + sel.t0.dt.strftime("%Y%m%dT%H")
    assert sel.window_id.is_unique
    pos = [f"{a}_{k}" for k in range(N) for a in ("lon", "lat")]
    for k in range(N):
        assert sel[f"lon_{k}"].notna().all() and sel[f"lat_{k}"].notna().all()   # 9 observations
    cols = ["window_id", "drifter_id", "t0", "block_id", "lon_0", "lat_0"] + \
           [c for c in pos if c not in ("lon_0", "lat_0")] + \
           ["zone", "min_land_km", "drogue_class", "split", "season", "drogue_lost_date", "ve_pre", "vn_pre"]
    C.RESULTS.mkdir(exist_ok=True)
    sel.sort_values(["t0", "drifter_id"])[cols].to_csv(C.RESULTS / "windows_selected.csv", index=False)
    blocks.to_csv(C.RESULTS / "blocks.csv", index=False)
    dr.to_csv(C.RESULTS / "drifters_split.csv", index=False)
    cdf = pd.DataFrame(rows_out)
    cdf.to_csv(C.RESULTS / "step0_counts.csv", index=False)

    # ----- gate + markdown
    notes = []
    off_test = sel[(sel.split == "test") & (sel.zone == "offshore")]
    tot_off = int((sel.zone == "offshore").sum())
    for c in ("drogued", "lost"):
        g = off_test[off_test.drogue_class == c]
        flag = g.drifter_id.nunique() < C.ANECDOTAL_MIN_DRIFTERS or len(g) < C.ANECDOTAL_MIN_WINDOWS
        notes.append(f"- offshore {c}, test set: {len(g)} windows, {g.drifter_id.nunique()} drifters -> "
                     f"{'ANECDOTAL' if flag else 'enough for interpretation (but windows are correlated)'}")
    for label, cz in (("all selected", sel[sel.zone == "coastal"]),
                      ("test set only (used for evaluation)", sel[(sel.zone == "coastal") & (sel.split == "test")])):
        flag = cz.drifter_id.nunique() < C.ANECDOTAL_MIN_DRIFTERS or len(cz) < C.ANECDOTAL_MIN_WINDOWS
        notes.append(f"- coastal (<{C.COASTAL_KM:.0f} km), {label}: {len(cz)} windows, {cz.drifter_id.nunique()} drifters -> "
                     f"{'ANECDOTAL: any coastal result is anecdotal' if flag else 'above the anecdotal thresholds'}")
    gate = tot_off >= C.OFFSHORE_MIN_TOTAL_WINDOWS
    md = ["# Step zero: usable sample", "",
          f"Generated by `scripts/s02_windows_stepzero.py` from `{C.GDP_CSV.relative_to(C.ROOT)}`.",
          "Windows are 48 h long (9 six-hourly observations). `6 h stride` means every possible start time; "
          "`blocks` are the 96 h time blocks that get forcing data (2 seed times each, 48 h apart).", "",
          "| step | stratum | drifters | tracks | windows | drifter-days |", "|---|---|---|---|---|---|"]
    for r in rows_out:
        f = lambda x: "" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:,}"
        md.append(f"| {r['step']} | {r['stratum']} | {f(r['n_drifters'])} | {f(r['n_tracks'])} | "
                  f"{f(r['n_windows'])} | {f(r['n_drifter_days'])} |")
    md += ["", "For steps 1-5 the `windows` column holds the number of six-hourly observation rows, not windows.",
           "", "A track is a run of at least 9 consecutive six-hourly observations after the filters of that step.",
           "", "## Sample adequacy", *notes, "",
           f"Selected offshore windows in total: {tot_off} (gate: at least {C.OFFSHORE_MIN_TOTAL_WINDOWS}) -> "
           f"{'gate passed' if gate else 'GATE FAILED'}.",
           f"Anecdotal rule: fewer than {C.ANECDOTAL_MIN_WINDOWS} windows or fewer than {C.ANECDOTAL_MIN_DRIFTERS} drifters.",
           "Windows from one drifter, and windows on the same date, are correlated, so the effective sample size is "
           "closer to the number of drifters than the number of windows.", "",
           "Not counted here: wind-speed classes need the downloaded GFS winds (done in the metrics step). "
           "Seasons are counted in step 15.", ""]
    (C.WORK / "step0_counts.md").write_text("\n".join(md))
    print("\n".join(notes))
    print("blocks:", len(blocks), "selected windows:", len(sel), "offshore:", tot_off)
    if not gate:
        print("GATE FAILED: fewer than 150 offshore windows", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()

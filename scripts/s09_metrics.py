"""Metric tables and pass criteria from the raw simulation results (no new OpenDrift runs).

Reads  results/windows.csv, results/ens_metrics.csv, results/tuning_choices.json, results/tables_tuning_*.csv,
       results/forcing_check.csv, data/runs/det_errors.parquet, data/runs/ens/test_*.npz
Writes results/metrics_*.csv, results/reliability_offshore_test.csv, results/pass_criteria.csv, results/pass_criteria.md

Definitions (work/plan.md D9, D10, work/baseline.md):
  error      separation distance (km) between the observed position and the model position at the lead.
             model_ens    = centroid of the 500-particle ensemble with the tuned persistent offset (the product that carries the contours; HEADLINE)
             model_native = centroid of the same ensemble without the offset (native OpenDrift noise only)
             model_det    = the deterministic single-particle run of the chosen configuration (no noise)
             leeway_*     = OpenDrift Leeway, object class 1, lost class only; B0 stay, B1 hand calculation, B3 persistence (not operator-available)
  hit50/95   observed position inside the 50/95 percent highest-density contour of the particle cloud (circle for B1_circle)
  intervals  95 percent percentile intervals of a cluster bootstrap that resamples whole drifters (B = config.BOOT_B, fixed seed);
             headline numbers are repeated with 96 h blocks (a date cluster) as the cluster, and the wider interval is quoted.
  n_eff      effective number of independent windows = n / (1 + (m - 1) * rho), m = mean windows per drifter, rho = intra-drifter
             correlation of the 24 h model error (one-way ANOVA estimate, floored at 0). It is an approximation, not a bootstrap output.
Test numbers and tuning-set (in-sample) numbers are written to separate files.
Usage: python scripts/s09_metrics.py
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
import config as C
from drift_lib import hdr_levels, kde_density

LEADS = list(C.LEADS_H)
OUT = C.RESULTS
ERR_METHODS = ["model_ens", "model_native", "model_det", "leeway_native", "leeway_offset", "B0", "B1", "B3"]
CONT_METHODS = ["model_ens", "model_native", "leeway_native", "leeway_offset", "B1_circle"]
ENS_NAME = {"od_final": "model_ens", "od_native": "model_native", "leeway_native": "leeway_native",
            "leeway_offset": "leeway_offset", "B1_circle": "B1_circle"}
IMPR_METHODS = ["model_ens", "model_native", "model_det", "leeway_native", "leeway_offset", "B3"]


# ------------------------------------------------------------------ data assembly
def load_all():
    w = pd.read_csv(C.RESULTS / "windows.csv", dtype={"drifter_id": str, "block_id": str})
    ch = json.loads((C.RESULTS / "tuning_choices.json").read_text())
    det = pd.read_parquet(C.RUNS_DIR / "det_errors.parquet")
    for c in ("method", "current", "wind"):
        det[c] = det[c].astype(str)
    return w, ch, det


def det_pick(det, method, current, wind, stokes, windage):
    d = det[(det.method == method) & (det.current == current) & (det.wind == wind) & (det.stokes == stokes) & np.isclose(det.windage, windage)]
    return d.pivot(index="window_id", columns="lead_h", values="err_km")


def chosen_det(det, ch, cls):
    """(OD, B0, B1, B3) wide error tables (index window_id, columns lead) of the chosen model/baseline of a class."""
    c, b = ch["choices"][cls], ch["baseline"][cls]
    if cls == "lost":
        od = det_pick(det, "OD", c["current"], c["wind"], c["stokes"], c["windage"])
        b1 = det_pick(det, "B1", "cur_notide", b["wind"], False, b["windage"])
    else:
        od = det_pick(det, "OD", c["current"], "none", False, 0.0)
        b1 = det_pick(det, "B1", "cur_6h_15m", "none", False, 0.0)
    b0 = det_pick(det, "B0", "none", "none", False, 0.0)
    b3 = det_pick(det, "B3", "none", "none", False, 0.0)
    return od, b0, b1, b3


def build_test_table(w, ch, det):
    """One row per test window that has a model result (2,504), columns kind|method|lead."""
    ens = pd.read_csv(C.RESULTS / "ens_metrics.csv")
    ens["m"] = ens.method.map(ENS_NAME)
    meta_cols = ["window_id", "drifter_id", "block_id", "split", "stratum", "zone", "drogue_class", "season", "wind_class", "t0", "has_ecmwf"]
    base = ens[ens.method == "od_final"].drop_duplicates("window_id")[["window_id"]].merge(w[meta_cols], on="window_id")
    base = base.set_index("window_id")
    cols = {}
    for m in ENS_NAME.values():
        e = ens[ens.m == m]
        for lead in LEADS:
            el = e[e.lead_h == lead].set_index("window_id")
            for kind, src in (("hit50", "hit50"), ("hit95", "hit95"), ("area50", "area50_km2"), ("area95", "area95_km2")):
                cols[f"{kind}|{m}|{lead}"] = el[src].astype(float)
            if m != "B1_circle":
                cols[f"err|{m}|{lead}"] = el["centroid_err_km"].astype(float)
            else:
                cols[f"chk|B1c|{lead}"] = el["centroid_err_km"].astype(float)
    for cls in ("lost", "drogued"):
        od, b0, b1, b3 = chosen_det(det, ch, cls)
        for name, t in (("model_det", od), ("B0", b0), ("B1", b1), ("B3", b3)):
            for lead in LEADS:
                s = t[lead].astype(float)
                key = f"err|{name}|{lead}"
                cols[key] = s if key not in cols else pd.concat([cols[key], s]).groupby(level=0).first()
    T = base.join(pd.DataFrame(cols), how="left")
    # B3 (persistence of the drifter's own 12 h pre-window velocity) is undefined where that velocity is missing; use B0 there
    T["b3_fallback_to_B0"] = T["err|B3|6"].isna()
    for lead in LEADS:
        T[f"err|B3|{lead}"] = T[f"err|B3|{lead}"].fillna(T[f"err|B0|{lead}"])
    # det tables hold both classes' windows only for their own class; the join leaves other-class windows NaN, so fill by class
    return T


# ------------------------------------------------------------------ bootstrap engine
class Spec:
    """Statistic vector over a window table: per (kind, method, lead) a median, p90 or mean, plus improvements."""

    def __init__(self, T):
        self.names, self.ops, self.cols = [], [], []
        for lead in LEADS:
            for m in ERR_METHODS:
                self.add(T, f"err|{m}|{lead}", "err_med|%s|%d" % (m, lead), "median")
                self.add(T, f"err|{m}|{lead}", "err_p90|%s|%d" % (m, lead), "p90")
            for m in CONT_METHODS:
                self.add(T, f"hit50|{m}|{lead}", "hit50|%s|%d" % (m, lead), "mean")
                self.add(T, f"hit95|{m}|{lead}", "hit95|%s|%d" % (m, lead), "mean")
                self.add(T, f"area50|{m}|{lead}", "area50_med|%s|%d" % (m, lead), "median")
                self.add(T, f"area95|{m}|{lead}", "area95_med|%s|%d" % (m, lead), "median")
                self.add(T, f"area95|{m}|{lead}", "area95_mean|%s|%d" % (m, lead), "mean")
        self.ops = np.array(self.ops)
        # improvements: 1 - median(model)/median(baseline)
        self.imp = []
        for lead in LEADS:
            for m in IMPR_METHODS:
                for b in ("B1", "B0"):
                    a, c = f"err_med|{m}|{lead}", f"err_med|{b}|{lead}"
                    if a in self.names and c in self.names:
                        self.imp.append((f"impr_{b}|{m}|{lead}", self.names.index(a), self.names.index(c)))
        self.V = None

    def add(self, T, col, name, op):
        if col in T.columns and T[col].notna().all():
            self.names.append(name)
            self.ops.append(op)
            self.cols.append(col)

    def matrix(self, T):
        return T[self.cols].to_numpy(float)

    def stats(self, V, idx):
        sub = V[idx]
        out = np.empty(len(self.names))
        for op, fn in (("median", lambda a: np.median(a, 0)), ("p90", lambda a: np.quantile(a, 0.9, 0)), ("mean", lambda a: a.mean(0))):
            sel = self.ops == op
            if sel.any():
                out[sel] = fn(sub[:, sel])
        imp = np.array([1.0 - out[a] / out[c] for _, a, c in self.imp])
        return np.concatenate([out, imp])

    @property
    def all_names(self):
        return self.names + [n for n, _, _ in self.imp]


def cluster_groups(clusters):
    uniq, inv = np.unique(np.asarray(clusters), return_inverse=True)
    return [np.where(inv == i)[0] for i in range(len(uniq))]


def boot(spec, T, cluster_col, B=C.BOOT_B, seed=C.BOOT_SEED):
    V = spec.matrix(T)
    groups = cluster_groups(T[cluster_col].to_numpy())
    rng = np.random.default_rng(seed)
    pt = spec.stats(V, np.arange(len(T)))
    draws = np.empty((B, len(pt)))
    for b in range(B):
        pick = rng.integers(0, len(groups), len(groups))
        draws[b] = spec.stats(V, np.concatenate([groups[i] for i in pick]))
    return pt, np.quantile(draws, 0.025, axis=0), np.quantile(draws, 0.975, axis=0), draws


def icc_neff(T, col="err|model_ens|24", cluster="drifter_id"):
    """ANOVA intra-cluster correlation and the effective number of independent windows."""
    y = T[col].to_numpy(float) if col in T.columns else T["err|model_det|24"].to_numpy(float)
    g = T[cluster].to_numpy()
    uniq, inv = np.unique(g, return_inverse=True)
    k, n = len(uniq), len(y)
    if k < 2 or n <= k:
        return np.nan, float(n)
    ni = np.bincount(inv).astype(float)
    mi = np.bincount(inv, weights=y) / ni
    grand = y.mean()
    msb = (ni * (mi - grand) ** 2).sum() / (k - 1)
    msw = ((y - mi[inv]) ** 2).sum() / (n - k)
    n0 = (n - (ni ** 2).sum() / n) / (k - 1)
    rho = max(0.0, (msb - msw) / (msb + (n0 - 1) * msw))
    m = n / k
    return float(rho), float(n / (1 + (m - 1) * rho))


def sample_info(T):
    rho, neff = icc_neff(T)
    return dict(n_windows=len(T), n_drifters=T.drifter_id.nunique(), n_blocks=T.block_id.nunique(), n_seed_dates=T.t0.nunique(),
                mean_windows_per_drifter=len(T) / max(1, T.drifter_id.nunique()), icc_24h_error=rho, n_eff=neff,
                anecdotal=bool(len(T) < C.ANECDOTAL_MIN_WINDOWS or T.drifter_id.nunique() < C.ANECDOTAL_MIN_DRIFTERS))


def summarise(T, group, by_block=False, B=C.BOOT_B):
    """Wide rows (one per lead) for a group table T. Returns (DataFrame, boot draws dict)."""
    spec = Spec(T)
    pt, lo, hi, _ = boot(spec, T, "drifter_id", B=B)
    if by_block:
        _, lo_b, hi_b, _ = boot(spec, T, "block_id", B=B)
    names = spec.all_names
    info = sample_info(T)
    rows = []
    for lead in LEADS:
        r = dict(group, lead_h=lead, **info)
        for i, n in enumerate(names):
            parts = n.split("|")
            if int(parts[-1]) != lead:
                continue
            key = "|".join(parts[:-1]).replace("|", "_")
            r[key] = pt[i]
            r[key + "_lo"] = lo[i]
            r[key + "_hi"] = hi[i]
            if by_block:
                r[key + "_blk_lo"] = lo_b[i]
                r[key + "_blk_hi"] = hi_b[i]
        rows.append(r)
    return pd.DataFrame(rows)


def groups_of(T, kind):
    """Group tables for a stratum type (always offshore unless the kind is 'zone')."""
    off = T[T.zone == "offshore"]
    out = []
    if kind == "headline":
        out = [("offshore", "all", off), ("offshore", "lost", off[off.drogue_class == "lost"]), ("offshore", "drogued", off[off.drogue_class == "drogued"])]
    elif kind == "zone":
        for z in ("offshore", "mid", "coastal"):
            for cls in ("all", "lost", "drogued"):
                s = T[T.zone == z]
                s = s if cls == "all" else s[s.drogue_class == cls]
                if len(s):
                    out.append((z, cls, s))
    else:
        col = {"season": "season", "wind_class": "wind_class", "ecmwf_grid": "ecmwf_grid"}[kind]
        for v in sorted(off[col].unique()):
            for cls in ("all", "lost", "drogued"):
                s = off[off[col] == v]
                s = s if cls == "all" else s[s.drogue_class == cls]
                if len(s):
                    out.append((v, cls, s))
    return out


def _job(args):
    kind, stratum, cls, T, by_block = args
    return summarise(T, dict(stratum_type=kind, stratum=stratum, drogue_class=cls), by_block=by_block)


# ------------------------------------------------------------------ tuning set (in sample)
def tuning_table(w, ch, det):
    t = w[(w.split == "tuning") & (w.zone == "offshore") & (w.has_ecmwf | (w.drogue_class == "drogued"))].copy().set_index("window_id")
    frames = []
    for cls in ("lost", "drogued"):
        od, b0, b1, b3 = chosen_det(det, ch, cls)
        idx = t.index[t.drogue_class == cls]
        f = t.loc[idx].copy()
        for name, tab in (("model_det", od), ("B0", b0), ("B1", b1), ("B3", b3)):
            for lead in LEADS:
                f[f"err|{name}|{lead}"] = tab[lead].reindex(idx).to_numpy()
        frames.append(f)
    T = pd.concat(frames)
    T["window_id"] = T.index
    for lead in LEADS:   # B3 undefined without a pre-window velocity: use B0 there (same rule as the test table)
        T[f"err|B3|{lead}"] = T[f"err|B3|{lead}"].fillna(T[f"err|B0|{lead}"])
    return T


# ------------------------------------------------------------------ variants
def variant_tables(w, ch, det):
    tv = pd.read_csv(C.RESULTS / "tables_tuning_variants.csv")
    test = w[(w.split == "test") & (w.zone == "offshore") & (w.has_ecmwf | (w.drogue_class == "drogued"))]
    rows, allrows = [], []
    for cls in ("lost", "drogued"):
        ids = test[test.drogue_class == cls]
        T = ids.set_index("window_id")
        od, b0, b1, b3 = chosen_det(det, ch, cls)
        sel = det[(det.method == "OD") & det.window_id.isin(ids.window_id) & det.lead_h.isin([24, 48])]
        if cls == "drogued":
            sel = sel[(sel.wind == "none") & (~sel.stokes) & np.isclose(sel.windage, 0.0)]
        else:
            sel = sel[sel.wind != "none"]
        # all configurations (point estimates, includes configurations never chosen: for context)
        for key, g in sel.groupby(["current", "wind", "stokes", "windage"]):
            p = g.pivot(index="window_id", columns="lead_h", values="err_km").reindex(ids.window_id)
            allrows.append(dict(drogue_class=cls, current=key[0], wind=key[1], stokes=key[2], windage=key[3], n_windows=int(p[24].notna().sum()),
                                med24_test=p[24].median(), med48_test=p[48].median()))
        # tuned variants (windage per variant chosen on the tuning set) with bootstrap
        tvc = tv[tv.cls == cls]
        cols, meta = {}, []
        for r in tvc.itertuples():
            cur, wnd, stk, wd = r.current, (r.wind if cls == "lost" else "none"), (bool(r.stokes) if cls == "lost" else False), (float(r.windage) if cls == "lost" else 0.0)
            g = sel[(sel.current == cur) & (sel.wind == wnd) & (sel.stokes == stk) & np.isclose(sel.windage, wd)]
            p = g.pivot(index="window_id", columns="lead_h", values="err_km").reindex(ids.window_id)
            cols[(cur, wnd, stk, wd)] = p
            meta.append(dict(drogue_class=cls, role="variant with the windage tuned on the tuning set", current=cur, wind=wnd, stokes=stk, windage=wd,
                             chosen_on_tuning=bool(r.selected), tuning_med24=r.med24,
                             tuning_tied_with_best=bool(r.tie)))
        ref24, ref48 = b1[24].reindex(ids.window_id).to_numpy(), b1[48].reindex(ids.window_id).to_numpy()
        chosen_key = [k for k, m in zip(cols, meta) if m["chosen_on_tuning"]][0]
        # best configuration on the TEST set, over all windage values (optimistic: selection used test data)
        ar = pd.DataFrame([a for a in allrows if a["drogue_class"] == cls])
        best = ar.sort_values("med24_test").iloc[0]
        bkey = (best.current, best.wind, bool(best.stokes), float(best.windage))
        g = sel[(sel.current == bkey[0]) & (sel.wind == bkey[1]) & (sel.stokes == bkey[2]) & np.isclose(sel.windage, bkey[3])]
        cols[("TESTBEST",) + bkey] = g.pivot(index="window_id", columns="lead_h", values="err_km").reindex(ids.window_id)
        meta.append(dict(drogue_class=cls, role="best configuration on the TEST set (post hoc, optimistic)", current=bkey[0], wind=bkey[1], stokes=bkey[2], windage=bkey[3], chosen_on_tuning=False,
                         tuning_med24=np.nan, tuning_tied_with_best=False, note="best configuration on the test set (post hoc, optimistic)"))
        keys = list(cols)
        E24 = np.column_stack([cols[k][24].to_numpy() for k in keys])
        E48 = np.column_stack([cols[k][48].to_numpy() for k in keys])
        clusters = ids.drifter_id.to_numpy()
        groups = cluster_groups(clusters)
        ci = keys.index(chosen_key)
        rng = np.random.default_rng(C.BOOT_SEED)

        def stat(idx):
            m24, m48 = np.median(E24[idx], 0), np.median(E48[idx], 0)
            b24, b48 = np.median(ref24[idx]), np.median(ref48[idx])
            return np.concatenate([m24, m48, 1 - m24 / b24, 1 - m48 / b48, m24 - m24[ci]])

        pt = stat(np.arange(len(ids)))
        dr = np.empty((C.BOOT_B, len(pt)))
        for b in range(C.BOOT_B):
            pick = rng.integers(0, len(groups), len(groups))
            dr[b] = stat(np.concatenate([groups[i] for i in pick]))
        lo, hi = np.quantile(dr, 0.025, 0), np.quantile(dr, 0.975, 0)
        K = len(keys)
        for i, m in enumerate(meta):
            m.update(n_windows=len(ids), n_drifters=ids.drifter_id.nunique(), med24_test=pt[i], med24_lo=lo[i], med24_hi=hi[i],
                     med48_test=pt[K + i], impr24_vs_B1=pt[2 * K + i], impr24_lo=lo[2 * K + i], impr24_hi=hi[2 * K + i],
                     impr48_vs_B1=pt[3 * K + i], impr48_lo=lo[3 * K + i], impr48_hi=hi[3 * K + i],
                     med24_minus_chosen=pt[4 * K + i], med24_minus_chosen_lo=lo[4 * K + i], med24_minus_chosen_hi=hi[4 * K + i],
                     B1_med24=float(np.median(ref24)))
            rows.append(m)
    return pd.DataFrame(rows), pd.DataFrame(allrows)


# ------------------------------------------------------------------ reliability
def _rel_job(path):
    sys.path.insert(0, str(C.ROOT / "scripts"))
    import s08_run_ensemble as S8
    ch = json.loads((C.RESULTS / "tuning_choices.json").read_text())
    z = np.load(path, allow_pickle=True)
    w = pd.read_csv(C.RESULTS / "windows.csv", dtype={"drifter_id": str}).set_index("window_id")
    ens = pd.read_csv(C.RESULTS / "ens_metrics.csv", usecols=["window_id", "method", "lead_h", "obs_x", "obs_y", "hit50", "hit95"])
    ens = ens[ens.method == "od_final"].set_index(["window_id", "lead_h"])
    rows = []
    for cls in ("lost", "drogued"):
        if f"{cls}|wid" not in z.files:
            continue
        wids = z[f"{cls}|wid"]
        X, Y = z[f"{cls}|od|x"], z[f"{cls}|od|y"]
        sb = ch["uncertainty"][cls]["sigma_b"]
        for i, wid in enumerate(wids):
            if w.loc[wid, "zone"] != "offshore":
                continue
            zz = S8.offsets(wid, X.shape[1])
            for k, lead in enumerate(LEADS):
                smp = S8.displaced(X[i, :, k], Y[i, :, k], zz, sb, lead)
                H, lev = hdr_levels(smp, ps=tuple(C.RELIABILITY_PS))
                e = ens.loc[(wid, lead)]
                d = float(kde_density(smp, np.array([[e.obs_x], [e.obs_y]]), H)[0])
                r = dict(window_id=wid, drogue_class=cls, lead_h=lead, drifter_id=w.loc[wid, "drifter_id"])
                for p in C.RELIABILITY_PS:
                    r[f"p{p}"] = bool(d >= lev[p])
                r["chk50"] = bool(r["p0.5"]) == bool(e.hit50)
                r["chk95"] = bool(r["p0.95"]) == bool(e.hit95)
                rows.append(r)
    return rows


def reliability(workers=7):
    paths = sorted((C.RUNS_DIR / "ens").glob("test_*.npz"))
    with ProcessPoolExecutor(workers) as ex:
        rows = [r for part in ex.map(_rel_job, paths) for r in part]
    R = pd.DataFrame(rows)
    assert R.chk50.all() and R.chk95.all(), "recomputed hit flags differ from results/ens_metrics.csv"
    out = []
    for cls in ("all", "lost", "drogued"):
        s = R if cls == "all" else R[R.drogue_class == cls]
        groups = cluster_groups(s.drifter_id.to_numpy())
        for lead in LEADS:
            sl = s[s.lead_h == lead].reset_index(drop=True)
            M = sl[[f"p{p}" for p in C.RELIABILITY_PS]].to_numpy(float)
            gl = cluster_groups(sl.drifter_id.to_numpy())
            rng = np.random.default_rng(C.BOOT_SEED)
            dr = np.array([M[np.concatenate([gl[i] for i in rng.integers(0, len(gl), len(gl))])].mean(0) for _ in range(C.BOOT_B)])
            for j, p in enumerate(C.RELIABILITY_PS):
                out.append(dict(drogue_class=cls, lead_h=lead, nominal=p, hit_rate=M[:, j].mean(), lo=np.quantile(dr[:, j], 0.025),
                                hi=np.quantile(dr[:, j], 0.975), n_windows=len(sl), n_drifters=sl.drifter_id.nunique()))
    return pd.DataFrame(out)


# ------------------------------------------------------------------ pass criteria
def verdict_improvement(x):
    return "PASS" if x >= C.PASS_IMPROVEMENT else ("BORDERLINE" if x >= C.BORDER_IMPROVEMENT else "FAIL")


def verdict_hit(x):
    return "PASS" if x >= C.PASS_HIT95 else ("BORDERLINE" if x >= C.BORDER_HIT95 else "FAIL")


def wide(a, b):
    return min(a), max(b)


def main():
    OUT.mkdir(exist_ok=True)
    w, ch, det = load_all()
    T = build_test_table(w, ch, det)
    T["ecmwf_grid"] = np.where(pd.to_datetime(T.t0) < pd.Timestamp("2024-02-01"), "ECMWF 0.4 deg (before 2024-02-01)", "ECMWF 0.25 deg (from 2024-02-01)")
    bad = [c for c in T.columns if c.startswith("err|") and "leeway" not in c and T[c].isna().any()]
    assert not bad, f"missing deterministic or ensemble errors: {bad[:4]}"
    off = T[T.zone == "offshore"]
    print(f"test windows scored: {len(T)} (offshore {len(off)}); drifters {T.drifter_id.nunique()} (offshore {off.drifter_id.nunique()})")
    # consistency with windows.csv (plan step 10 assert)
    wt = w[(w.split == "test") & (w.zone == "offshore") & (w.has_ecmwf | (w.drogue_class == "drogued"))]
    assert len(off) == len(wt), (len(off), len(wt))
    # B1 from the deterministic (geodesic) table agrees with the circle-centre distance in the ensemble file (local planar km)
    dev = max(((T[f"err|B1|{l}"] - T[f"chk|B1c|{l}"]).abs() / T[f"err|B1|{l}"]).max() for l in LEADS)
    assert dev < 0.03, dev
    print(f"B1 error, geodesic vs planar local-km: largest relative difference {dev:.3f}")
    # ---- test-set tables
    jobs = []
    for kind in ("headline", "zone", "season", "wind_class", "ecmwf_grid"):
        for stratum, cls, g in groups_of(T, kind):
            jobs.append((kind, stratum, cls, g, kind == "headline"))
    with ProcessPoolExecutor(7) as ex:
        res = list(ex.map(_job, jobs))
    S = pd.concat(res, ignore_index=True)
    S.insert(0, "set", "test")
    S.to_csv(OUT / "metrics_test_by_stratum.csv", index=False, float_format="%.5g")
    H = S[S.stratum_type == "headline"].copy()
    keep = ["set", "stratum", "drogue_class", "lead_h", "n_windows", "n_drifters", "n_blocks", "n_eff", "icc_24h_error"]
    errc = [c for c in H.columns if c.startswith(("err_med_", "err_p90_", "impr_"))]
    H[keep + errc].dropna(axis=1, how="all").to_csv(OUT / "metrics_test_by_lead.csv", index=False, float_format="%.5g")
    contc = [c for c in H.columns if c.startswith(("hit50_", "hit95_", "area50_", "area95_"))]
    H[keep + contc].dropna(axis=1, how="all").to_csv(OUT / "metrics_test_contours_by_lead.csv", index=False, float_format="%.5g")
    ss = S.drop_duplicates(["stratum_type", "stratum", "drogue_class"])[["set", "stratum_type", "stratum", "drogue_class", "n_windows", "n_drifters", "n_blocks",
                                                                          "n_seed_dates", "mean_windows_per_drifter", "icc_24h_error", "n_eff", "anecdotal"]]
    # ---- tuning set (in-sample)
    TT = tuning_table(w, ch, det)
    tres = []
    for cls in ("all", "lost", "drogued"):
        g = TT if cls == "all" else TT[TT.drogue_class == cls]
        tres.append(summarise(g, dict(stratum_type="headline", stratum="offshore", drogue_class=cls)))
    TS = pd.concat(tres, ignore_index=True)
    TS.insert(0, "set", "tuning_in_sample")
    tcols = [c for c in TS.columns if c.startswith(("err_med_", "err_p90_", "impr_"))]
    TS[keep + tcols].dropna(axis=1, how="all").to_csv(OUT / "metrics_tuning_insample_by_lead.csv", index=False, float_format="%.5g")
    tu = pd.read_csv(C.RESULTS / "tables_tuning_uncertainty.csv")
    sel = []
    for cls in ("lost", "drogued"):
        u = ch["uncertainty"][cls]
        a = tu[(tu.cls == cls) & (tu.l1 == u["l1"]["id"]) & np.isclose(tu.sigma_b, u["sigma_b"])].assign(setting="chosen")
        nat = tu[(tu.cls == cls) & (tu.l1 == u["native_only_best"]["l1"]) & (tu.sigma_b == 0)].assign(setting="native_only_best")
        sel += [a, nat]
    pd.concat(sel).to_csv(OUT / "metrics_tuning_calibration_by_lead.csv", index=False, float_format="%.5g")
    ss = pd.concat([ss, TS.drop_duplicates(["stratum", "drogue_class"])[["set", "stratum_type", "stratum", "drogue_class", "n_windows", "n_drifters", "n_blocks",
                                                                          "n_seed_dates", "mean_windows_per_drifter", "icc_24h_error", "n_eff", "anecdotal"]]])
    ss.to_csv(OUT / "metrics_sample_sizes.csv", index=False, float_format="%.5g")
    # ---- variants
    V, A = variant_tables(w, ch, det)
    V.to_csv(OUT / "metrics_test_variants.csv", index=False, float_format="%.5g")
    A.to_csv(OUT / "metrics_test_all_configs_point.csv", index=False, float_format="%.5g")
    # ---- reliability (recomputes the contours of the saved particle clouds, checks the flags against ens_metrics.csv)
    Rl = reliability()
    Rl.to_csv(OUT / "reliability_offshore_test.csv", index=False, float_format="%.5g")
    # ---- pass criteria
    write_pass_criteria(T, S, TS, V, A, Rl, ch)
    print("done")


def get(S, stratum, cls, lead, key):
    r = S[(S.stratum_type == "headline") & (S.stratum == stratum) & (S.drogue_class == cls) & (S.lead_h == lead)].iloc[0]
    return r[key], r[key + "_lo"], r[key + "_hi"], r.get(key + "_blk_lo", np.nan), r.get(key + "_blk_hi", np.nan), r


def write_pass_criteria(T, S, TS, V, A, Rl, ch):
    rows = []

    def add(**kw):
        rows.append(kw)

    # criterion 1
    for cls, label in (("all", "pooled (headline)"), ("lost", "lost-drogue (the SAR-relevant class)"), ("drogued", "drogued")):
        for mkey, mlabel in (("model_ens", "ensemble centroid (headline model)"), ("model_det", "deterministic particle (sensitivity)"),
                             ("model_native", "ensemble centroid without offset (sensitivity)")):
            v, lo, hi, blo, bhi, r = get(S, "offshore", cls, 24, f"impr_B1_{mkey}")
            med_m = r[f"err_med_{mkey}"]
            med_b = r["err_med_B1"]
            lo_w, hi_w = wide([lo, blo], [hi, bhi])
            headline = (cls == "all" and mkey == "model_ens")
            add(criterion="1 offshore skill vs hand-calculation baseline", item=f"{cls}: {mlabel}", headline=headline,
                lead_h=24, n_windows=int(r.n_windows), n_drifters=int(r.n_drifters), n_eff=r.n_eff,
                model_median_km=med_m, baseline_median_km=med_b, value=v, ci_lo_drifter=lo, ci_hi_drifter=hi, ci_lo_block=blo, ci_hi_block=bhi,
                ci_lo_wider=lo_w, ci_hi_wider=hi_w, threshold=C.PASS_IMPROVEMENT, verdict=verdict_improvement(v),
                ci_includes_threshold=bool(lo_w <= C.PASS_IMPROVEMENT <= hi_w), note="value = 1 - median(model error) / median(B1 error), same windows")
    # 48 h for context
    for cls in ("all", "lost", "drogued"):
        v, lo, hi, blo, bhi, r = get(S, "offshore", cls, 48, "impr_B1_model_ens")
        lo_w, hi_w = wide([lo, blo], [hi, bhi])
        add(criterion="1 offshore skill vs hand-calculation baseline", item=f"{cls}: ensemble centroid, 48 h (context, not the criterion)", headline=False,
            lead_h=48, n_windows=int(r.n_windows), n_drifters=int(r.n_drifters), n_eff=r.n_eff, model_median_km=r["err_med_model_ens"],
            baseline_median_km=r["err_med_B1"], value=v, ci_lo_drifter=lo, ci_hi_drifter=hi, ci_lo_block=blo, ci_hi_block=bhi,
            ci_lo_wider=lo_w, ci_hi_wider=hi_w, threshold=C.PASS_IMPROVEMENT, verdict=verdict_improvement(v),
            ci_includes_threshold=bool(lo_w <= C.PASS_IMPROVEMENT <= hi_w), note="48 h is not the agreed criterion lead")
    # best-on-test oracle (post hoc)
    vv = V[V.role.str.contains("TEST set")]
    for cls in ("lost", "drogued"):
        r = vv[vv.drogue_class == cls].iloc[0]
        add(criterion="1 offshore skill vs hand-calculation baseline", item=f"{cls}: best deterministic configuration chosen ON THE TEST SET (post hoc, optimistic)", headline=False,
            lead_h=24, n_windows=int(r.n_windows), n_drifters=int(r.n_drifters), n_eff=np.nan, model_median_km=r.med24_test, baseline_median_km=r.B1_med24,
            value=r.impr24_vs_B1, ci_lo_drifter=r.impr24_lo, ci_hi_drifter=r.impr24_hi, ci_lo_block=np.nan, ci_hi_block=np.nan, ci_lo_wider=r.impr24_lo,
            ci_hi_wider=r.impr24_hi, threshold=C.PASS_IMPROVEMENT, verdict=verdict_improvement(r.impr24_vs_B1) + " (post hoc)",
            ci_includes_threshold=bool(r.impr24_lo <= C.PASS_IMPROVEMENT <= r.impr24_hi),
            note=f"{r.current}, wind {r.wind}, Stokes {r.stokes}, windage {r.windage}; the interval ignores the selection on test data")
    # criterion 2
    for cls in ("all", "lost", "drogued"):
        v95, lo95, hi95, blo, bhi, r = get(S, "offshore", cls, 24, "hit95_model_ens")
        v50, lo50, hi50, _, _, _ = get(S, "offshore", cls, 24, "hit50_model_ens")
        lo_w, hi_w = wide([lo95, blo], [hi95, bhi])
        add(criterion="2 calibration (95% contour hit rate, after tuning)", item=f"{cls}: 95% contour", headline=(cls == "all"), lead_h=24,
            n_windows=int(r.n_windows), n_drifters=int(r.n_drifters), n_eff=r.n_eff, model_median_km=np.nan, baseline_median_km=np.nan, value=v95,
            ci_lo_drifter=lo95, ci_hi_drifter=hi95, ci_lo_block=blo, ci_hi_block=bhi, ci_lo_wider=lo_w, ci_hi_wider=hi_w, threshold=C.PASS_HIT95,
            verdict=verdict_hit(v95), ci_includes_threshold=bool(lo_w <= C.PASS_HIT95 <= hi_w),
            note=f"tuned on the tuning set to {ch['uncertainty']['lost']['tuning_hit95_24h']:.2f} (lost) / {ch['uncertainty']['drogued']['tuning_hit95_24h']:.2f} (drogued)")
        add(criterion="2 calibration (50% contour hit rate, reported only)", item=f"{cls}: 50% contour", headline=False, lead_h=24,
            n_windows=int(r.n_windows), n_drifters=int(r.n_drifters), n_eff=r.n_eff, model_median_km=np.nan, baseline_median_km=np.nan, value=v50,
            ci_lo_drifter=lo50, ci_hi_drifter=hi50, ci_lo_block=np.nan, ci_hi_block=np.nan, ci_lo_wider=lo50, ci_hi_wider=hi50, threshold=np.nan,
            verdict="REPORTED (nominal 0.50)", ci_includes_threshold=bool(lo50 <= 0.5 <= hi50), note="no pass line was agreed for the 50% contour")
    # native-only calibration (what OpenDrift's own uncertainty keys give without the post hoc offset)
    for cls in ("all", "lost", "drogued"):
        v, lo, hi, _, _, r = get(S, "offshore", cls, 24, "hit95_model_native")
        add(criterion="2 calibration (95% contour hit rate, after tuning)", item=f"{cls}: chosen native OpenDrift keys, without the offset layer", headline=False, lead_h=24,
            n_windows=int(r.n_windows), n_drifters=int(r.n_drifters), n_eff=r.n_eff, model_median_km=np.nan, baseline_median_km=np.nan, value=v,
            ci_lo_drifter=lo, ci_hi_drifter=hi, ci_lo_block=np.nan, ci_hi_block=np.nan, ci_lo_wider=lo, ci_hi_wider=hi, threshold=C.PASS_HIT95,
            verdict=verdict_hit(v) + " (without the offset layer)", ci_includes_threshold=bool(lo <= C.PASS_HIT95 <= hi),
            note="same particles as the headline model with the offset left out; the best native-only setting of the tuning grid reached 0.23 (lost) and 0.38 (drogued) on the tuning set and was not run on the test set")
    # criterion 3: coastal
    for cls in ("all", "lost", "drogued"):
        sub = S[(S.stratum_type == "zone") & (S.stratum == "coastal") & (S.drogue_class == cls) & (S.lead_h == 24)]
        if sub.empty:
            continue
        r = sub.iloc[0]
        add(criterion="3 coastal skill (reported, NOT pass/fail)", item=f"{cls}: coastal (min distance to land under {C.COASTAL_KM:.0f} km)", headline=False, lead_h=24,
            n_windows=int(r.n_windows), n_drifters=int(r.n_drifters), n_eff=r.n_eff, model_median_km=r.err_med_model_ens, baseline_median_km=r.err_med_B1,
            value=r.impr_B1_model_ens, ci_lo_drifter=r.impr_B1_model_ens_lo, ci_hi_drifter=r.impr_B1_model_ens_hi, ci_lo_block=np.nan, ci_hi_block=np.nan,
            ci_lo_wider=r.impr_B1_model_ens_lo, ci_hi_wider=r.impr_B1_model_ens_hi, threshold=np.nan, verdict="ANECDOTAL, NOT PASS/FAIL" if r.anecdotal else "NOT PASS/FAIL",
            ci_includes_threshold=False, note=f"95% contour hit rate {r.hit95_model_ens:.2f}, 50% {r.hit50_model_ens:.2f} (n={int(r.n_windows)}); value = improvement over B1")
    P = pd.DataFrame(rows)
    P.to_csv(OUT / "pass_criteria.csv", index=False, float_format="%.5g")
    write_md(P, T, S, TS, V, A, Rl, ch)


def pct(x, d=1):
    return f"{100 * x:.{d}f}%"


def ci(lo, hi, d=1):
    return f"{100 * lo:.{d}f}% to {100 * hi:.{d}f}%"


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def write_md(P, T, S, TS, V, A, Rl, ch):
    w = pd.read_csv(C.RESULTS / "windows.csv", dtype={"drifter_id": str})
    fc = pd.read_csv(C.RESULTS / "forcing_check.csv")
    g = lambda st, cls, lead, key: S[(S.stratum_type == 'headline') & (S.stratum == st) & (S.drogue_class == cls) & (S.lead_h == lead)].iloc[0][key]
    head = P[P.headline & P.criterion.str.startswith("1")].iloc[0]
    c1 = P[P.criterion.str.startswith("1 ")]
    lost = c1[c1.item.str.startswith("lost: ensemble centroid (headline")].iloc[0]
    drog = c1[c1.item.str.startswith("drogued: ensemble centroid (headline")].iloc[0]
    h95 = P[P.headline & P.criterion.str.startswith("2")].iloc[0]
    r50 = P[P.item.eq("all: 50% contour")].iloc[0]
    L = []
    A_ = L.append
    A_("# Pass criteria, evaluated on the test set (offshore)\n")
    A_("Generated by `scripts/s09_metrics.py` from `results/ens_metrics.csv`, `data/runs/det_errors.parquet` and the saved particle clouds. "
       "Machine-readable version: `results/pass_criteria.csv`. All numbers are test-set numbers (drifters never used for any tuning) unless a section says otherwise. "
       "Intervals are 95% percentile intervals of a cluster bootstrap that resamples whole drifters (2,000 replicates, seed 42); the interval is also "
       "computed with 96 h blocks as clusters, and the wider of the two is quoted as \"CI\".\n")
    A_("## Verdicts\n")
    v1 = head.verdict
    A_(md_table(["Criterion", "Result", "Verdict"], [
        ["1. Offshore: median 24 h separation error at least about 30% below the hand-calculation baseline (pooled, headline)",
         f"{pct(head.value)} lower (CI {ci(head.ci_lo_wider, head.ci_hi_wider)}); {head.model_median_km:.2f} km vs {head.baseline_median_km:.2f} km", f"**{v1}**"],
        ["1b. Same, lost-drogue class only (the SAR-relevant class)", f"{pct(lost.value)} lower (CI {ci(lost.ci_lo_wider, lost.ci_hi_wider)}); {lost.model_median_km:.2f} km vs {lost.baseline_median_km:.2f} km", lost.verdict],
        ["1c. Same, drogued class only", f"{pct(drog.value)} lower (CI {ci(drog.ci_lo_wider, drog.ci_hi_wider)}); {drog.model_median_km:.2f} km vs {drog.baseline_median_km:.2f} km", drog.verdict],
        ["2. Calibration: 95% contour holds at least about 80% of observed endpoints at 24 h (pooled, after tuning)",
         f"{pct(h95.value)} (CI {ci(h95.ci_lo_wider, h95.ci_hi_wider)}); the 50% contour holds {pct(r50.value)} (CI {ci(r50.ci_lo_wider, r50.ci_hi_wider)})", f"**{h95.verdict}**"],
        ["3. Coastal skill", f"{int(P[P.item.str.startswith('all: coastal')].n_windows.iloc[0])} test windows from {int(P[P.item.str.startswith('all: coastal')].n_drifters.iloc[0])} drifters", "reported, not pass/fail"],
    ]))
    A_(f"\nRules, fixed before the test-set numbers were looked at (`scripts/config.py`): criterion 1 is PASS at {pct(C.PASS_IMPROVEMENT, 0)} or more, BORDERLINE from "
       f"{pct(C.BORDER_IMPROVEMENT, 0)} to under {pct(C.PASS_IMPROVEMENT, 0)} (\"about 30%\"), FAIL below {pct(C.BORDER_IMPROVEMENT, 0)}. Criterion 2 is PASS at {pct(C.PASS_HIT95, 0)} or more, BORDERLINE from "
       f"{pct(C.BORDER_HIT95, 0)}, FAIL below. The verdict uses the point estimate; the CI column says whether the threshold is inside the interval.\n")
    A_("## Criterion 1: offshore median 24 h error against the hand-calculation baseline\n")
    A_("The baseline (B1) is the operator's hand calculation: datum, one drift vector from the same current (and, for lost drifters, 0.5% of the forecast wind) evaluated once at the datum, "
       "and an expanding circle; its windage and circle growth rates were tuned on the tuning set (`work/baseline.md`). The model is OpenDrift with the configuration chosen on the tuning set "
       "(`work/tuning.md`). Model error is the distance from the observed position to the centroid of the 500-particle ensemble (the product that carries the contours; this was the headline "
       "definition, chosen before the deterministic-particle and no-offset variants were computed; the three agree within 0.3 percentage points, and all are shown). \"Reduction\" is 1 - median(model error) / median(B1 error) on the same windows.\n")
    rows = []
    for r in c1.itertuples():
        rows.append([r.item, r.lead_h, f"{r.n_windows:,}", r.n_drifters, "" if pd.isna(r.n_eff) else f"{r.n_eff:.0f}", f"{r.model_median_km:.2f}", f"{r.baseline_median_km:.2f}", pct(r.value),
                     ci(r.ci_lo_drifter, r.ci_hi_drifter), "" if pd.isna(r.ci_lo_block) else ci(r.ci_lo_block, r.ci_hi_block), r.verdict])
    A_(md_table(["Row", "Lead h", "Windows", "Drifters", "n_eff*", "Model median km", "B1 median km", "Reduction", "CI by drifter", "CI by 96 h block", "Verdict"], rows))
    A_("\n\\* n_eff is the approximate number of independent windows: windows / (1 + (m - 1) x rho), where m is the mean number of windows per drifter and rho the intra-drifter correlation of the 24 h error "
       "(ANOVA estimate). It does not account for drifters that share weather on the same date; the by-block interval does.\n")
    A_("The row labelled \"best configuration chosen ON THE TEST SET\" picks, among all 108 lost-class (or 5 drogued-class) deterministic configurations, the one with the lowest test-set median. "
       "It is selection on test data, so it is an optimistic ceiling for what choosing among the variants could give, and its interval ignores the selection.\n")
    headline_fail = head.value < C.PASS_IMPROVEMENT
    if headline_fail:
        A_(f"**The model does not clearly beat the baseline offshore at 24 h.** The pooled reduction is {pct(head.value)} (CI {ci(head.ci_lo_wider, head.ci_hi_wider)}), "
           f"{'the whole interval is below ' + pct(C.PASS_IMPROVEMENT, 0) if head.ci_hi_wider < C.PASS_IMPROVEMENT else 'the interval includes ' + pct(C.PASS_IMPROVEMENT, 0)}. "
           f"For lost drifters it is {pct(lost.value)} (CI {ci(lost.ci_lo_wider, lost.ci_hi_wider)}) and for drogued drifters {pct(drog.value)} (CI {ci(drog.ci_lo_wider, drog.ci_hi_wider)}). "
           f"Criterion 1 is therefore not met in the agreed form. At 48 h the pooled reduction is {pct(c1[c1.item.str.startswith('all: ensemble centroid, 48')].value.iloc[0])} "
           f"(CI {ci(*c1[c1.item.str.startswith('all: ensemble centroid, 48')][['ci_lo_wider', 'ci_hi_wider']].iloc[0])}), but 24 h is the agreed lead.\n")
    # --- evidence for what would need to change
    need = 0.70 * head.baseline_median_km
    lv = V[V.role.str.startswith("variant") & (V.drogue_class == "lost")]
    dv = V[V.role.str.startswith("variant") & (V.drogue_class == "drogued")]
    lost_best = V[V.role.str.contains("TEST") & (V.drogue_class == "lost")].iloc[0]
    drog_best = V[V.role.str.contains("TEST") & (V.drogue_class == "drogued")].iloc[0]
    f_test = fc[fc.split == "test"]
    chosen_cur = {"lost": ch["choices"]["lost"]["current"], "drogued": ch["choices"]["drogued"]["current"]}
    fl = f_test[(f_test.drogue_class == "lost") & (f_test.current == chosen_cur["lost"])].iloc[0]
    fd = f_test[(f_test.drogue_class == "drogued") & (f_test.current == chosen_cur["drogued"])].iloc[0]
    b3 = g("offshore", "all", 24, "err_med_B3")
    pd_ = pd.concat([lv, dv])
    n_sig = int(((pd_.med24_minus_chosen_lo > 0) | (pd_.med24_minus_chosen_hi < 0)).sum())
    A_("### If the model does not clearly beat the baseline: what would need to change\n")
    A_("The numbers below are what the data show. The remedies are inferences from them, not things this spike tested.\n")
    A_(f"- **Size of the gap.** A 30% reduction would need a pooled 24 h median of {need:.2f} km. The model gives {head.model_median_km:.2f} km, so it is {head.model_median_km - need:.2f} km short "
       f"(the median error would have to fall by a further {pct(1 - need / head.model_median_km, 0)}).")
    A_(f"- **Forcing choices do not close it.** Across the 18 lost-class variants (currents with and without tide, 6-hourly and hourly, GFS, ECMWF or the mean of both, Stokes drift on or off) with windage tuned on the tuning set, the "
       f"test-set 24 h median ranges from {lv.med24_test.min():.2f} to {lv.med24_test.max():.2f} km; the five drogued-class current variants range from {dv.med24_test.min():.2f} to {dv.med24_test.max():.2f} km. "
       f"Even the best configuration chosen on the test set reaches only {pct(lost_best.impr24_vs_B1)} (lost) and {pct(drog_best.impr24_vs_B1)} (drogued). "
       f"Adding tides, Stokes drift or a second wind source is therefore not where the missing skill is: of the {len(lv) + len(dv) - 2} paired differences against the chosen variant, {n_sig} have an interval that excludes zero (`results/metrics_test_variants.csv`). "
       f"One note on the drogued class: the 15 m current variants (the drogue depth) were best at 24 h on both sets ({dv.set_index('current').loc['cur_6h_15m_tide'].impr24_vs_B1:.1%} and {dv.set_index('current').loc['cur_6h_15m'].impr24_vs_B1:.1%} on the test set), but within the 2% tie band on the tuning set, "
       f"and the pre-set tie rule (lower 48 h median) picked the no-tide surface current, which gives {drog.value:.1%}. The rule was not changed after the test set was seen.")
    A_(f"- **Part of the error is a drift that persists and the analysis misses.** B3, which carries the drifter's own velocity from the 12 h before the window forward, has a pooled 24 h median of {b3:.2f} km "
       f"({pct(1 - b3 / head.baseline_median_km, 0)} below B1 and {pct(1 - b3 / head.model_median_km, 0)} below the model). A real operator does not have that information, but it is consistent with a good part (not all) of the 24 h error being a drift "
       f"that persists for a day or more and is absent from the 1/12 degree analysis. The forcing check points the same way: for the chosen currents the vector RMSE against the drifters' own velocity is "
       f"{fl.rmse_vec:.2f} m/s (lost) and {fd.rmse_vec:.2f} m/s (drogued), against observed mean speeds of {fl.speed_obs:.2f} and {fd.speed_obs:.2f} m/s, with vector correlation {fl.vector_corr:.2f} and {fd.vector_corr:.2f} (`results/forcing_check.csv`, test set).")
    A_("- **What the evidence points to.** (i) Better current forcing is the lever: higher-resolution or regional ocean forcing (the brief's example), or a different or additional ocean model for eddy and boundary-current fields. "
       "(ii) Anything that anchors the cloud to an observed drift, such as a datum marker buoy or a drifter report, would address the persistence-type error that B3 exposes. "
       "(iii) A hindcast with archived forecasts instead of analysis currents would show how much of this is lost live; it can only get worse than the numbers here. None of these was tested in this spike, and (i) in particular is untested: "
       "the data show the 1/12 degree analysis is not enough for the 30% target at 24 h, not that a finer product would be.")
    A_(f"- **The baseline is not equally strong for both classes.** For drogued drifters B1 is no better than leaving the datum where it is (B0 median {g('offshore', 'drogued', 24, 'err_med_B0'):.2f} km vs B1 {g('offshore', 'drogued', 24, 'err_med_B1'):.2f} km at 24 h), "
       f"so the 30% line is against a weak reference there and the model still falls short ({pct(drog.value)}). Against B0 (no calculation at all) the pooled reduction is {pct(g('offshore', 'all', 24, 'impr_B0_model_ens'))} "
       f"(CI {ci(g('offshore', 'all', 24, 'impr_B0_model_ens_lo'), g('offshore', 'all', 24, 'impr_B0_model_ens_hi'))}).")
    A_("\n## Criterion 2: calibration of the contours\n")
    A_("The 95% and 50% contours are highest-density regions of the 500-particle cloud (Gaussian kernel density estimate, 3 km minimum bandwidth). The cloud is the OpenDrift ensemble with the Level-1 uncertainty setting chosen on the tuning set "
       f"plus a persistent velocity offset applied to each particle outside OpenDrift ({ch['uncertainty']['lost']['sigma_b']} m/s for lost, {ch['uncertainty']['drogued']['sigma_b']} m/s for drogued, per component), selected on the tuning set so that the 24 h 95% hit rate was at least 0.90 with the smallest area. "
       "That offset is an extra tuned parameter, not an OpenDrift feature. Without it the chosen native keys hold almost nothing (last three rows of the table).\n")
    rows = []
    for r in P[P.criterion.str.startswith("2")].itertuples():
        rows.append([r.item, f"{r.n_windows:,}", r.n_drifters, f"{r.n_eff:.0f}", pct(r.value), ci(r.ci_lo_drifter, r.ci_hi_drifter), "" if pd.isna(r.ci_lo_block) else ci(r.ci_lo_block, r.ci_hi_block), r.verdict])
    A_(md_table(["Row (24 h, offshore)", "Windows", "Drifters", "n_eff", "Hit rate", "CI by drifter", "CI by 96 h block", "Verdict"], rows))
    A_("\nHit rate by lead (offshore, test; the offset was tuned at 24 h only). `circle` is the hand-calculation baseline's expanding circle, calibrated on the tuning set at 24 h to the same 50% and 95% targets.\n")
    rows = []
    for cls in ("all", "lost", "drogued"):
        for lead in LEADS:
            r = S[(S.stratum_type == "headline") & (S.drogue_class == cls) & (S.lead_h == lead)].iloc[0]
            rows.append([cls, lead, f"{pct(r.hit50_model_ens, 0)} ({pct(r.hit50_model_ens_lo, 0)} to {pct(r.hit50_model_ens_hi, 0)})", f"{pct(r.hit95_model_ens, 0)} ({pct(r.hit95_model_ens_lo, 0)} to {pct(r.hit95_model_ens_hi, 0)})",
                         pct(r.hit50_B1_circle, 0), pct(r.hit95_B1_circle, 0), f"{r.area95_med_model_ens:,.0f}", f"{r.area95_med_B1_circle:,.0f}"])
    A_(md_table(["Class", "Lead h", "Model 50% hit (CI)", "Model 95% hit (CI)", "Circle 50% hit", "Circle 95% hit", "Model median 95% area km2", "Circle 95% area km2"], rows))
    r24 = S[(S.stratum_type == "headline") & (S.drogue_class == "all") & (S.lead_h == 24)].iloc[0]
    rl24 = S[(S.stratum_type == "headline") & (S.drogue_class == "lost") & (S.lead_h == 24)].iloc[0]
    A_(f"\nAt 24 h the model's 95% contour has a median area of {r24.area95_med_model_ens:,.0f} km2 against {r24.area95_med_B1_circle:,.0f} km2 for the baseline circle "
       f"({pct(1 - r24.area95_med_model_ens / r24.area95_med_B1_circle, 0)} smaller) at hit rates of {pct(r24.hit95_model_ens, 0)} and {pct(r24.hit95_B1_circle, 0)}. For lost drifters: {rl24.area95_med_model_ens:,.0f} vs "
       f"{rl24.area95_med_B1_circle:,.0f} km2 at {pct(rl24.hit95_model_ens, 0)} vs {pct(rl24.hit95_B1_circle, 0)}. The two are not matched on hit rate exactly, and the circle's growth rate was fitted to a single 24 h target, "
       "so this is an indication of a smaller search area at similar coverage, not a precise equal-coverage comparison.\n")
    rl = Rl[(Rl.drogue_class == "all") & (Rl.lead_h == 24)].set_index("nominal")
    A_("Reliability at 24 h (all offshore test windows; nominal contour level against the fraction of observed endpoints inside; `results/reliability_offshore_test.csv`, `figures/fig2_calibration.png`): "
       + ", ".join(f"{int(p * 100)}% holds {pct(rl.loc[p].hit_rate, 1)}" for p in (0.1, 0.5, 0.8, 0.9, 0.95, 0.99))
       + ". Up to about the 80% level the contours hold slightly more than nominal; at 90% and above they hold less. The cloud is somewhat too concentrated in the middle and its tails are too light, "
       "so the 95% contour is reliable only because the offset was tuned for it, and the 99% contour would not be.\n")
    h6, h12 = g("offshore", "all", 6, "hit95_model_ens"), g("offshore", "all", 12, "hit95_model_ens")
    A_(f"Short leads are under-dispersed: the pooled 95% contour holds {pct(h6, 1)} at 6 h and {pct(h12, 1)} at 12 h, and for drogued drifters {pct(g('offshore', 'drogued', 6, 'hit95_model_ens'), 1)} at 6 h. "
       "One persistent offset grows linearly from zero, while real short-lead errors do not. The 80% pass line is reached at every lead for the pooled sample (6 h only just), for the lost class at every lead, "
       f"and for the drogued class from 12 h. The nominal 95% is not reached at any lead for the pooled sample (best: {pct(g('offshore', 'all', 48, 'hit95_model_ens'), 1)} at 48 h); only the lost class reaches it, at 48 h ({pct(g('offshore', 'lost', 48, 'hit95_model_ens'), 1)}).\n")
    A_("## Criterion 3: coastal skill (reported, not pass/fail)\n")
    wc = w[w.zone == "coastal"]
    A_(f"Coastal means the minimum distance to land over the 48 h window is under {C.COASTAL_KM:.0f} km. Across the whole study there are {len(wc)} coastal windows from {wc.drifter_id.nunique()} drifters "
       f"({(wc.split == 'tuning').sum()} windows in the tuning set, {(wc.split == 'test').sum()} in the test set). The test-set numbers are below. These are anecdotal (fewer than {C.ANECDOTAL_MIN_WINDOWS} windows) "
       "and are not interpreted. The model's 1/12 degree grid (about 8 km) cannot resolve harbours, estuaries or tidal straits, and a window that is coastal by distance may still be in open water.\n")
    rows = []
    for r in P[P.criterion.str.startswith("3")].itertuples():
        sub = S[(S.stratum_type == "zone") & (S.stratum == "coastal") & (S.drogue_class == r.item.split(":")[0]) & (S.lead_h == 24)].iloc[0]
        rows.append([r.item, r.n_windows, r.n_drifters, f"{r.model_median_km:.2f}", f"{r.baseline_median_km:.2f}", f"{pct(r.value, 0)} ({pct(r.ci_lo_wider, 0)} to {pct(r.ci_hi_wider, 0)})",
                     pct(sub.hit95_model_ens, 0), pct(sub.hit50_model_ens, 0), "anecdotal"])
    A_(md_table(["Coastal group (24 h)", "Windows", "Drifters", "Model median km", "B1 median km", "Reduction (CI)", "95% hit", "50% hit", "Status"], rows))
    A_("\nMid-zone windows (10 to 50 km from land) are in `results/metrics_test_by_stratum.csv` (stratum_type `zone`).\n")
    A_("## Sample sizes behind every table\n")
    ss = pd.read_csv(C.RESULTS / "metrics_sample_sizes.csv")
    ss = ss[(ss.stratum_type.isin(["headline", "zone"])) & ~((ss.stratum_type == "zone") & (ss.stratum == "offshore"))]
    A_(md_table(["Set", "Stratum", "Class", "Windows", "Drifters", "96 h blocks", "Mean windows per drifter", "ICC (24 h error)", "n_eff (approx.)", "Anecdotal"],
                [[r.set, r.stratum, r.drogue_class, f"{r.n_windows:,}", r.n_drifters, r.n_blocks, f"{r.mean_windows_per_drifter:.1f}", f"{r.icc_24h_error:.2f}", f"{r.n_eff:.0f}", "yes" if r.anecdotal else "no"] for r in ss.itertuples()]))
    A_("\nWindows from one drifter overlap in the sense of sharing the same drifter, regime and often the same weather, so the interval widths above come from resampling drifters, not windows. "
       f"The pooled offshore test sample is {int(head.n_windows):,} windows but only about {head.n_eff:.0f} effectively independent ones (75 drifters).\n")
    A_("## Tuning-set numbers (in-sample, not evidence of skill)\n")
    A_("Median 24 h error of the deterministic chosen configuration on the offshore tuning windows, which were used to choose it: "
       + "; ".join(f"{cls}: model {TS[(TS.drogue_class == cls) & (TS.lead_h == 24)].err_med_model_det.iloc[0]:.2f} km vs B1 {TS[(TS.drogue_class == cls) & (TS.lead_h == 24)].err_med_B1.iloc[0]:.2f} km "
                  f"({pct(TS[(TS.drogue_class == cls) & (TS.lead_h == 24)].impr_B1_model_det.iloc[0])}, n={int(TS[(TS.drogue_class == cls) & (TS.lead_h == 24)].n_windows.iloc[0]):,} windows, "
                  f"{int(TS[(TS.drogue_class == cls) & (TS.lead_h == 24)].n_drifters.iloc[0])} drifters)" for cls in ("lost", "drogued"))
       + f". Tuning-set 24 h hit rates of the chosen uncertainty setting: lost 50% {ch['uncertainty']['lost']['tuning_hit50_24h']:.3f}, 95% {ch['uncertainty']['lost']['tuning_hit95_24h']:.3f}; "
       f"drogued 50% {ch['uncertainty']['drogued']['tuning_hit50_24h']:.3f}, 95% {ch['uncertainty']['drogued']['tuning_hit95_24h']:.3f}. Files: `results/metrics_tuning_insample_by_lead.csv`, `results/metrics_tuning_calibration_by_lead.csv`.\n")
    A_("## Reading these results\n")
    A_("- The test uses best-estimate analysis currents, not archived forecasts, so live skill would be lower by an amount this test cannot measure.\n"
       "- Drifters are 35 cm floats (drogued ones follow the 15 m current), not people or boats; the tuned windage of 0.5% does not transfer to a SAR object.\n"
       "- Every window starts from an exactly known position and time. A real datum has neither.\n"
       "- The tuned offset layer means the calibration shown is for settings fitted to this region and period; it was fitted on 40% of drifters and tested on the other 60%.\n"
       "- Strata results (season, wind class, ECMWF grid period) are shown without multiple-comparison correction; differences between strata are mostly within their intervals.\n"
       "- Windows without a pre-window drifter velocity (13 windows of 4,086) use B0 in the B3 reference line only.\n")
    (C.RESULTS / "pass_criteria.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()

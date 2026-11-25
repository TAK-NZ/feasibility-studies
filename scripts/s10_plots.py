"""Figures from the metric tables and the saved particle clouds (no new OpenDrift runs). PNGs go to figures/, each side under 2500 px.
Reads results/metrics_*.csv, results/reliability_offshore_test.csv, results/windows.csv, results/tuning_choices.json,
data/runs/det_errors.parquet, data/runs/ens/test_<tag>.npz. Run s09_metrics.py first.
  fig1_error_vs_lead.png     median and 90th percentile separation error by lead: model, hand calculation (B1), stay (B0), persistence (B3), Leeway
  fig2_calibration.png       hit rate of the 50% / 95% contours by lead, reliability at 24 h, contour area by lead (model vs baseline circle)
  fig3_example_windows.png   observed track, ensemble cloud, 50% / 95% contours and baseline for four example windows, with locator maps
  fig4_strata.png            24 h reduction against B1 and 95% hit rate by stratum (class, zone, season, wind speed, ECMWF grid)
  fig5_variants.png          24 h median error of every forcing variant on the test set against the baseline and the 30% line
The example windows are chosen by rule (the window whose 24 h error is closest to the class median or 90th percentile), not by eye.
Usage: python scripts/s10_plots.py
"""
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402
import config as C  # noqa: E402
from drift_lib import hdr_levels, kde_density, local_km  # noqa: E402

FIG = C.FIGURES
LEADS = list(C.LEADS_H)
COL = {"model_ens": "#1b6ca8", "model_det": "#6baed6", "B1": "#d95f02", "B0": "#7f7f7f", "B3": "#1b9e77", "leeway_offset": "#9467bd", "B1_circle": "#d95f02"}
LAB = {"model_ens": "OpenDrift ensemble centroid", "B1": "B1 hand calculation", "B0": "B0 stay at datum", "B3": "B3 own prior velocity (not operator data)",
       "leeway_offset": "Leeway class 1 (untuned)", "B1_circle": "B1 expanding circle"}


def save(fig, name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / name, dpi=C.FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print("wrote", FIG / name)


def sel(S, st, cls, stype="headline"):
    return S[(S.stratum_type == stype) & (S.stratum == st) & (S.drogue_class == cls)].sort_values("lead_h")


def fig1(S):
    fig, ax = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    for j, (cls, title) in enumerate((("all", "offshore, both classes"), ("lost", "offshore, lost-drogue"), ("drogued", "offshore, drogued"))):
        d = sel(S, "offshore", cls)
        for i, (kind, ylabel) in enumerate((("med", "median separation error (km)"), ("p90", "90th percentile separation error (km)"))):
            a = ax[i, j]
            for m in ("B0", "B1", "B3", "leeway_offset", "model_ens"):
                k = f"err_{kind}_{m}"
                if k not in d.columns or d[k].isna().all():
                    continue
                a.plot(d.lead_h, d[k], marker="o", color=COL[m], label=LAB[m], lw=2.2 if m == "model_ens" else 1.4, ms=4, ls="--" if m in ("B3", "leeway_offset") else "-")
                a.fill_between(d.lead_h, d[k + "_lo"], d[k + "_hi"], color=COL[m], alpha=0.15 if m in ("model_ens", "B1") else 0.07, lw=0)
            a.set_xticks(LEADS)
            a.grid(alpha=0.3)
            if i == 1:
                a.set_xlabel("lead time (h)")
            if j == 0:
                a.set_ylabel(ylabel)
            if i == 0:
                n = int(d.n_windows.iloc[0])
                a.set_title(f"{title}\n{n:,} test windows, {int(d.n_drifters.iloc[0])} drifters, n_eff about {d.n_eff.iloc[0]:.0f}", fontsize=10)
    h, l = ax[0, 1].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.03))
    fig.suptitle("Separation error by lead time, test set (bands: 95% cluster bootstrap by drifter)", fontsize=12)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    save(fig, "fig1_error_vs_lead.png")


def fig2(S, R):
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
    d = sel(S, "offshore", "all")
    a = ax[0]
    for nominal, ls in ((0.5, ":"), (0.95, ":")):
        a.axhline(nominal, color="k", ls=ls, lw=1)
    a.axhline(C.PASS_HIT95, color="r", ls="--", lw=1)
    a.text(48, C.PASS_HIT95 - 0.03, "80% pass line", color="r", ha="right", fontsize=8)
    for lvl, mk in ((95, "o"), (50, "s")):
        for m, c in (("model_ens", COL["model_ens"]), ("B1_circle", COL["B1_circle"])):
            k = f"hit{lvl}_{m}"
            lo, hi = d.get(k + "_lo"), d.get(k + "_hi")
            a.errorbar(d.lead_h + (0.4 if m == "B1_circle" else -0.4), d[k], yerr=[d[k] - lo, hi - d[k]], fmt="-" + mk, color=c, capsize=2, ms=5, lw=1.5,
                       label=f"{LAB[m]} {lvl}%")
    a.set_xticks(LEADS)
    a.set_ylim(0, 1.02)
    a.set_xlabel("lead time (h)")
    a.set_ylabel("fraction of observed positions inside")
    a.set_title("Contour hit rate by lead (dotted: nominal 50% and 95%)", fontsize=10)
    a.legend(fontsize=7, loc="lower right")
    a.grid(alpha=0.3)
    a = ax[1]
    a.plot([0, 1], [0, 1], "k:", lw=1)
    for cls, c in (("all", "k"), ("lost", "#1b6ca8"), ("drogued", "#d95f02")):
        r = R[(R.drogue_class == cls) & (R.lead_h == 24)].sort_values("nominal")
        a.errorbar(r.nominal, r.hit_rate, yerr=[r.hit_rate - r.lo, r.hi - r.hit_rate], fmt="-o", color=c, capsize=2, ms=4, lw=1.4, label=cls)
    a.axhline(C.PASS_HIT95, color="r", ls="--", lw=0.8)
    a.set_xlabel("nominal contour level")
    a.set_ylabel("observed fraction inside")
    a.set_title("Reliability at 24 h (offshore, test)", fontsize=10)
    a.legend(fontsize=8)
    a.set_aspect("equal")
    a.grid(alpha=0.3)
    a = ax[2]
    for cls, c, ls in (("lost", "#1b6ca8", "-"), ("drogued", "#6baed6", "-")):
        dd = sel(S, "offshore", cls)
        a.plot(dd.lead_h, dd.area95_med_model_ens, ls, marker="o", color=c, label=f"model {cls}")
        a.plot(dd.lead_h, dd.area95_med_B1_circle, "--", marker="s", color="#d95f02" if cls == "lost" else "#fdae6b", label=f"circle {cls}")
    a.set_yscale("log")
    a.set_xticks(LEADS)
    a.set_xlabel("lead time (h)")
    a.set_ylabel("median area of the 95% contour (km$^2$)")
    a.set_title("Search area at 95% (model vs baseline circle)", fontsize=10)
    a.legend(fontsize=8)
    a.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save(fig, "fig2_calibration.png")


# ---------------------------------------------------------------- example windows
def land_grid():
    from roaring_landmask import RoaringLandmask
    lon = np.arange(C.FORCING_LON[0], C.FORCING_LON[1] + 1e-6, 0.2)
    lat = np.arange(C.FORCING_LAT[0], C.FORCING_LAT[1] + 1e-6, 0.2)
    g = np.meshgrid(lon, lat)
    land = RoaringLandmask.new().contains_many_par(((g[0].ravel() + 180) % 360) - 180, g[1].ravel()).reshape(g[0].shape)
    return lon, lat, land


def fig3(W, ch, det, E):
    import s08_run_ensemble as S8
    test = E[(E.method == "od_final") & (E.lead_h == 24)]
    picks = []
    for cls, zone, q, label in (("lost", "offshore", 0.5, "lost, offshore, median error"), ("lost", "offshore", 0.9, "lost, offshore, 90th percentile error"),
                                ("drogued", "offshore", 0.5, "drogued, offshore, median error"), ("lost", "coastal", 0.5, "lost, coastal, median error (anecdotal)")):
        t = test[(test.drogue_class == cls) & (test.zone == zone)]
        target = t.centroid_err_km.quantile(q, interpolation="nearest")
        wid = t.iloc[(t.centroid_err_km - target).abs().argmin()].window_id
        picks.append((wid, cls, label))
    lon, lat, land = land_grid()
    fig, ax = plt.subplots(2, 4, figsize=(17, 8.8), gridspec_kw=dict(height_ratios=[1.35, 1]))
    for j, (wid, cls, label) in enumerate(picks):
        w = W.set_index("window_id").loc[wid]
        tag = pd.Timestamp(w.t0).strftime("%Y%m%dT%H")
        z = np.load(C.RUNS_DIR / "ens" / f"test_{tag}.npz", allow_pickle=True)
        i = int(np.where(z[f"{cls}|wid"] == wid)[0][0])
        X, Y = z[f"{cls}|od|x"][i], z[f"{cls}|od|y"][i]
        sb = ch["uncertainty"][cls]["sigma_b"]
        zz = S8.offsets(wid, X.shape[0])
        a = ax[0, j]
        # observed track
        ol = np.array([[w[f"lon_{k}"], w[f"lat_{k}"]] for k in range(9)])
        ox, oy = local_km(ol[:, 0], ol[:, 1], w.lon_0, w.lat_0)
        # baseline positions
        b = ch["baseline"][cls]
        dd = det[(det.window_id == wid) & (det.method == "B1")]
        dd = dd[(dd.wind == b["wind"]) & np.isclose(dd.windage, b["windage"])] if cls == "lost" else dd[dd.current == "cur_6h_15m"]
        bx, by = local_km(dd.sort_values("lead_h").lon.values, dd.sort_values("lead_h").lat.values, w.lon_0, w.lat_0)
        allx, ally = [ox], [oy]
        for k, (lead, col) in enumerate(zip(LEADS, ("#c6dbef", "#9ecae1", "#3182bd", "#08306b"))):
            smp = S8.displaced(X[:, k], Y[:, k], zz, sb, lead)
            if lead in (24, 48):
                a.scatter(smp[0][::3], smp[1][::3], s=3, color=col, alpha=0.5, lw=0, label=f"particles {lead} h")
                H, lev = hdr_levels(smp)
                pad = 3 * np.sqrt(np.diag(H)) + 5
                gx = np.linspace(smp[0].min() - pad[0], smp[0].max() + pad[0], 160)
                gy = np.linspace(smp[1].min() - pad[1], smp[1].max() + pad[1], 160)
                GX, GY = np.meshgrid(gx, gy)
                dens = kde_density(smp, np.vstack([GX.ravel(), GY.ravel()]), H).reshape(GX.shape)
                a.contour(GX, GY, dens, levels=[lev[0.95]], colors=[col], linewidths=1.6)
                a.contour(GX, GY, dens, levels=[lev[0.5]], colors=[col], linewidths=1.0, linestyles="--")
            allx.append(smp[0])
            ally.append(smp[1])
        a.plot(ox, oy, "-o", color="k", ms=4, lw=1.4, label="observed track (6-hourly)")
        a.plot(ox[[0, 4, 8]], oy[[0, 4, 8]], "o", color="k", ms=7, mfc="none")
        a.text(ox[4], oy[4], " 24 h", fontsize=7)
        a.text(ox[8], oy[8], " 48 h", fontsize=7)
        a.plot(bx, by, "-s", color=COL["B1"], ms=4, lw=1.2, label="B1 hand calculation")
        a.add_patch(Circle((bx[2], by[2]), b["a95_km_per_h"] * 24, fill=False, ec=COL["B1"], ls=":", lw=1.2))
        a.plot(0, 0, "*", color="gold", mec="k", ms=12, label="datum")
        e = test.set_index("window_id").loc[wid]
        e48 = E[(E.method == "od_final") & (E.lead_h == 48) & (E.window_id == wid)].iloc[0]
        a.set_title(f"{label}\ndrifter {wid[:15]}, {w.t0}\n24 h error {e.centroid_err_km:.1f} km; 95% contour hit: 24 h {'yes' if e.hit95 else 'NO'}, 48 h {'yes' if e48.hit95 else 'NO'}", fontsize=8)
        a.set_aspect("equal")
        a.set_xlabel("east of datum (km)")
        if j == 0:
            a.set_ylabel("north of datum (km)")
            a.legend(fontsize=6, loc="upper left")
        a.grid(alpha=0.3)
        m = ax[1, j]
        m.pcolormesh(lon, lat, np.where(land, 1.0, np.nan), cmap="Greys", vmin=0, vmax=2, shading="auto")
        m.plot(ol[:, 0] % 360, ol[:, 1], "r-", lw=2)
        m.plot(w.lon_0 % 360, w.lat_0, "r*", ms=8)
        m.set_xlim(C.FORCING_LON)
        m.set_ylim(C.FORCING_LAT)
        m.set_aspect(1 / np.cos(np.deg2rad(40)))
        m.set_title(f"location ({w.lat_0:.1f}, {w.lon_0 % 360:.1f}E), {w.min_land_km:.0f} km from land", fontsize=8)
        m.grid(alpha=0.3)
    fig.suptitle("Example test windows: solid contour = 95%, dashed = 50% (24 h and 48 h clouds); dotted circle = baseline 95% circle at 24 h. "
                 "Windows chosen by rule, not by eye.", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, "fig3_example_windows.png")


def fig4(S):
    S24 = S[(S.lead_h == 24) & S.stratum_type.isin(["zone", "season", "wind_class", "ecmwf_grid"]) & S.drogue_class.isin(["lost", "drogued"])].copy()
    order = ["zone", "season", "wind_class", "ecmwf_grid"]
    S24["o"] = S24.stratum_type.map({k: i for i, k in enumerate(order)})
    S24 = S24.sort_values(["o", "stratum", "drogue_class"], ascending=[True, True, True]).reset_index(drop=True)
    labels = [f"{r.stratum_type.replace('ecmwf_grid', 'ECMWF').replace('wind_class', 'wind m/s')}: {r.stratum.replace(' (before 2024-02-01)', ' <2024-02').replace(' (from 2024-02-01)', ' >=2024-02')} (n={int(r.n_windows):,}, {int(r.n_drifters)} dr.)"
              for r in S24.itertuples()]
    fig, ax = plt.subplots(1, 2, figsize=(15, 10), sharey=True)
    y = np.arange(len(S24))[::-1]
    for a, key, thr, xl in ((ax[0], "impr_B1_model_ens", C.PASS_IMPROVEMENT, "24 h median error reduction vs B1"), (ax[1], "hit95_model_ens", C.PASS_HIT95, "95% contour hit rate at 24 h")):
        for r, yy in zip(S24.itertuples(), y):
            c = "#1b6ca8" if r.drogue_class == "lost" else "#d95f02"
            v, lo, hi = getattr(r, key), getattr(r, key + "_lo"), getattr(r, key + "_hi")
            a.errorbar(v, yy, xerr=[[v - lo], [hi - v]], fmt="o", color=c, mfc="white" if r.anecdotal else c, capsize=2, ms=6)
        a.axvline(thr, color="r", ls="--", lw=1)
        a.set_xlabel(xl)
        a.grid(alpha=0.3)
    ax[0].axvline(0, color="k", lw=0.6)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels(labels, fontsize=7)
    ax[0].plot([], [], "o", color="#1b6ca8", label="lost-drogue")
    ax[0].plot([], [], "o", color="#d95f02", label="drogued")
    ax[0].plot([], [], "o", mfc="white", color="k", label="anecdotal (<30 windows)")
    ax[0].legend(fontsize=8, loc="lower right")
    fig.suptitle("Test set by stratum, 24 h (95% cluster bootstrap by drifter; red dashed = pass line; zone rows include mid and coastal)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, "fig4_strata.png")


def fig5(V, ch):
    fig, ax = plt.subplots(1, 2, figsize=(15, 8), gridspec_kw=dict(width_ratios=[18, 6]))
    for a, cls in zip(ax, ("lost", "drogued")):
        d = V[(V.drogue_class == cls)].copy().reset_index(drop=True)
        d["lab"] = [f"{r.current} / {r.wind if cls == 'lost' else '-'} / Stokes {'on' if r.stokes else 'off'}{' / windage ' + str(r.windage) if cls == 'lost' else ''}" + (" [TEST-BEST, post hoc]" if "TEST" in r.role else "")
                    + (" [chosen on tuning]" if r.chosen_on_tuning else "") for r in d.itertuples()]
        y = np.arange(len(d))[::-1]
        for r, yy in zip(d.itertuples(), y):
            c = "#1b6ca8" if r.chosen_on_tuning else ("#9467bd" if "TEST" in r.role else "#7f7f7f")
            a.errorbar(r.med24_test, yy, xerr=[[r.med24_test - r.med24_lo], [r.med24_hi - r.med24_test]], fmt="o", color=c, capsize=2, ms=5)
        b1 = d.B1_med24.iloc[0]
        a.axvline(b1, color="#d95f02", lw=1.5, label=f"B1 median {b1:.2f} km")
        a.axvline(b1 * (1 - C.PASS_IMPROVEMENT), color="r", ls="--", lw=1, label=f"30% below B1 ({b1 * (1 - C.PASS_IMPROVEMENT):.2f} km)")
        a.set_yticks(y)
        a.set_yticklabels(d.lab, fontsize=7)
        a.set_xlabel("test-set 24 h median error (km), deterministic particle")
        a.set_title(f"{cls}: {int(d.n_windows.iloc[0]):,} windows, {int(d.n_drifters.iloc[0])} drifters", fontsize=10)
        a.legend(fontsize=8, loc="lower right")
        a.grid(alpha=0.3)
    fig.suptitle("Forcing variants, test set, offshore (blue: chosen on the tuning set; purple: best on the test set, post hoc)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, "fig5_variants.png")


def main():
    S = pd.read_csv(C.RESULTS / "metrics_test_by_stratum.csv")
    R = pd.read_csv(C.RESULTS / "reliability_offshore_test.csv")
    V = pd.read_csv(C.RESULTS / "metrics_test_variants.csv")
    ch = json.loads((C.RESULTS / "tuning_choices.json").read_text())
    W = pd.read_csv(C.RESULTS / "windows.csv", dtype={"drifter_id": str})
    E = pd.read_csv(C.RESULTS / "ens_metrics.csv", dtype={"drifter_id": str})
    det = pd.read_parquet(C.RUNS_DIR / "det_errors.parquet", filters=[("method", "==", "B1")])
    for c in ("method", "current", "wind"):
        det[c] = det[c].astype(str)
    fig1(S)
    fig2(S, R)
    fig3(W, ch, det, E)
    fig4(S)
    fig5(V, ch)


if __name__ == "__main__":
    main()

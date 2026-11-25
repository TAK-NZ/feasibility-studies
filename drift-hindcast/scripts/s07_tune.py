"""Tune windage, choose the model variant and calibrate the baseline. TUNING SET ONLY (split == 'tuning', offshore).

Reads data/runs/det_errors.parquet and results/windows.csv. Writes results/tuning_choices.json,
results/tables_tuning_windage.csv, results/tables_tuning_variants.csv and work/tuning.md (part A).

Rules (decided before looking at the test set):
 a) per (current, wind, Stokes) of the lost class: the windage of the grid with the lowest median 24 h error (tie: lower 48 h median).
 b) lost class: the variant with the lowest tuning 24 h median error after (a); a variant within 2 percent of the best is a tie and
    the lower 48 h median wins. Drogued class: the same rule over the 5 current variants (no windage, no Stokes).
    Only windows that have every wind source (ECMWF file present) are used to compare lost-class variants.
 c) baseline B1 (work/baseline.md): the windage with the lowest B1 24 h median error, for the wind source of the chosen model;
    expanding-circle slopes a_0.5 and a_0.95 (km per hour) = the p-quantile of the 24 h separation between the B1 position and the
    observation, divided by 24 h, per class.
If the optimal windage sits at the upper edge of the grid (0 is the physical lower bound, and Stokes-on variants
legitimately choose 0) the script stops and asks for a wider grid (config.WINDAGE_GRID).
Usage: python scripts/s07_tune.py
"""
import json
import sys

import numpy as np
import pandas as pd

import config as C


def med(x):
    return float(np.median(x))


def main():
    w = pd.read_csv(C.RESULTS / "windows.csv", dtype={"drifter_id": str}, parse_dates=["t0"])
    det = pd.read_parquet(C.RUNS_DIR / "det_errors.parquet")
    for c in ("method", "current", "wind", "status"):
        det[c] = det[c].astype(str)
    tw = w[(w.split == "tuning") & (w.zone == "offshore")]
    out = {"tuning_set": {}, "choices": {}, "notes": []}
    lines = ["# Tuning record (part A: windage, variant choice, baseline)", "",
             "Tuned on the **tuning set only**: windows with `split == tuning` and `zone == offshore` (more than 50 km from land). "
             "Tuning and test drifters are disjoint (checked in `s05_build_forcing.py`).", ""]
    tables_w, tables_v = [], []
    for cls in ("lost", "drogued"):
        tc = tw[tw.drogue_class == cls]
        out["tuning_set"][cls] = dict(windows=len(tc), drifters=int(tc.drifter_id.nunique()))
        lines.append(f"- {cls}: {len(tc)} offshore tuning windows from {tc.drifter_id.nunique()} drifters")
    lines.append("")

    # ---------------------------------------------------------------- lost class
    tl = tw[(tw.drogue_class == "lost") & tw.has_ecmwf]
    out["tuning_set"]["lost_common_wind_windows"] = len(tl)
    d = det[(det.method == "OD") & det.window_id.isin(tl.window_id) & det.lead_h.isin([24, 48]) & det.wind.isin(C.WINDS)]
    p = d.groupby(["current", "wind", "stokes", "windage", "lead_h"]).err_km.median().unstack("lead_h").reset_index()
    p.columns = ["current", "wind", "stokes", "windage", "med24", "med48"]
    p["n_windows"] = len(tl)
    tables_w.append(p.assign(cls="lost"))
    best = []
    for (c, wv, s), g in p.groupby(["current", "wind", "stokes"]):
        g = g.sort_values(["med24", "med48"])
        b = g.iloc[0]
        best.append(b)
        if b.windage >= max(C.WINDAGE_GRID):   # 0 is the physical lower bound, so only the upper edge matters
            out["notes"].append(f"optimum windage at a grid edge for {c}/{wv}/stokes={s}: {b.windage}")
    best = pd.DataFrame(best).reset_index(drop=True)
    best["stokes"] = best.stokes.astype(bool)
    top = best.med24.min()
    best["tie"] = best.med24 <= top * 1.02
    sel = best[best.tie].sort_values("med48").iloc[0]
    lost_choice = dict(current=sel.current, wind=sel.wind, stokes=bool(sel.stokes), windage=float(sel.windage),
                       med24_tuning=float(sel.med24), med48_tuning=float(sel.med48), n_tied=int(best.tie.sum()))
    out["choices"]["lost"] = lost_choice
    tables_v.append(best.assign(cls="lost", selected=(best.index == sel.name)))
    lines += ["## Lost-drogue class (OceanDrift, tuned windage)",
              f"Windows used for the comparison: {len(tl)} (windows that have both wind sources). Windage grid: {C.WINDAGE_GRID}.", "",
              "Best windage and tuning-set 24 h / 48 h median error (km) for every current / wind / Stokes variant:", "",
              "| current | wind | Stokes | best windage | median 24 h km | median 48 h km | within 2% of best |", "|---|---|---|---|---|---|---|"]
    for r in best.sort_values("med24").itertuples():
        lines.append(f"| {r.current} | {r.wind} | {'on' if r.stokes else 'off'} | {r.windage} | {r.med24:.2f} | {r.med48:.2f} | {'yes' if r.tie else ''} |")
    lines += ["", f"**Selected**: {lost_choice['current']}, wind {lost_choice['wind']}, Stokes {'on' if lost_choice['stokes'] else 'off'}, "
              f"windage {lost_choice['windage']} (24 h median {lost_choice['med24_tuning']:.2f} km; {lost_choice['n_tied']} variants were within 2%)."]
    if lost_choice["windage"] >= max(C.WINDAGE_GRID):
        out["notes"].append("SELECTED lost windage is at a grid edge")

    # ---------------------------------------------------------------- drogued class
    td = tw[tw.drogue_class == "drogued"]
    d = det[(det.method == "OD") & det.window_id.isin(td.window_id) & det.lead_h.isin([24, 48]) & (det.wind == "none")]
    p = d.groupby(["current", "lead_h"]).err_km.median().unstack("lead_h").reset_index()
    p.columns = ["current", "med24", "med48"]
    top = p.med24.min()
    p["tie"] = p.med24 <= top * 1.02
    sel = p[p.tie].sort_values("med48").iloc[0]
    drog_choice = dict(current=sel.current, med24_tuning=float(sel.med24), med48_tuning=float(sel.med48), n_tied=int(p.tie.sum()))
    out["choices"]["drogued"] = drog_choice
    tables_v.append(p.assign(cls="drogued", selected=(p.current == sel.current)))
    lines += ["", "## Drogued class (OceanDrift, no windage, no Stokes)", "",
              f"Windows: {len(td)}. Tuning-set median errors (km):", "",
              "| current | median 24 h km | median 48 h km | within 2% of best |", "|---|---|---|---|"]
    for r in p.sort_values("med24").itertuples():
        lines.append(f"| {r.current} | {r.med24:.2f} | {r.med48:.2f} | {'yes' if r.tie else ''} |")
    lines += ["", f"**Selected**: {drog_choice['current']} (24 h median {drog_choice['med24_tuning']:.2f} km; {drog_choice['n_tied']} variants within 2%)."]

    # ---------------------------------------------------------------- baseline
    lines += ["", "## Baseline (hand calculation), see work/baseline.md", ""]
    base = {}
    b1 = det[det.method == "B1"]
    for cls, tc, wsrc in (("lost", tw[tw.drogue_class == "lost"], lost_choice["wind"]), ("drogued", td, "none")):
        g = b1[b1.window_id.isin(tc.window_id) & (b1.wind == wsrc) & (b1.lead_h == 24)]
        tab = g.groupby("windage").err_km.median()
        if cls == "lost":
            wd = float(tab.idxmin())
            if wd >= max(C.WINDAGE_GRID):
                out["notes"].append(f"baseline windage at a grid edge: {wd}")
            lines.append("B1 windage grid, tuning-set 24 h median error (km), wind source " + wsrc + ": " +
                         ", ".join(f"{k}: {v:.2f}" for k, v in tab.items()))
        else:
            wd = 0.0
        sub = g[g.windage == wd].set_index("window_id").err_km
        # circle slopes: smallest radius that contains p of the tuning-set observations at 24 h, divided by 24 h
        a50 = float(np.quantile(sub.values, 0.5, method="higher") / 24.0)
        a95 = float(np.quantile(sub.values, 0.95, method="higher") / 24.0)
        # shift to the median offshore error of B1 and B0, B3 for reference
        ref = {m: med(det[(det.method == m) & det.window_id.isin(tc.window_id) & (det.lead_h == 24)].err_km.dropna()) for m in ("B0", "B3")}
        base[cls] = dict(wind=wsrc, windage=wd, a50_km_per_h=a50, a95_km_per_h=a95, med24_B1_tuning=med(sub.values),
                         med24_B0_tuning=ref["B0"], med24_B3_tuning=ref["B3"])
        lines.append(f"- {cls}: B1 windage {wd}, wind {wsrc}; circle radius slopes a_0.5 = {a50:.3f} km/h, a_0.95 = {a95:.3f} km/h; "
                     f"tuning 24 h median error B1 {med(sub.values):.2f} km, B0 {ref['B0']:.2f} km, B3 {ref['B3']:.2f} km")
    out["baseline"] = base
    # tuning-set model vs baseline
    lines += ["", "## Tuning-set skill (in-sample, optimistic)", ""]
    for cls, ch in (("lost", lost_choice), ("drogued", drog_choice)):
        tc = tw[tw.drogue_class == cls]
        if cls == "lost":
            m = det[(det.method == "OD") & (det.current == ch["current"]) & (det.wind == ch["wind"]) &
                    (det.stokes == ch["stokes"]) & (det.windage == ch["windage"]) & det.window_id.isin(tc.window_id)]
        else:
            m = det[(det.method == "OD") & (det.current == ch["current"]) & det.window_id.isin(tc.window_id)]
        bl = base[cls]
        bb = det[(det.method == "B1") & (det.wind == bl["wind"]) & (det.windage == bl["windage"]) & det.window_id.isin(tc.window_id)]
        for lead in C.LEADS_H:
            mm = med(m[m.lead_h == lead].err_km)
            bm = med(bb[bb.lead_h == lead].err_km)
            out.setdefault("tuning_skill", {}).setdefault(cls, {})[str(lead)] = dict(model=mm, B1=bm, improvement=1 - mm / bm)
            lines.append(f"- {cls} {lead} h: model median {mm:.2f} km, B1 {bm:.2f} km, improvement {100 * (1 - mm / bm):.1f}%")
    out["notes"] = sorted(set(out["notes"]))
    (C.RESULTS / "tuning_choices.json").write_text(json.dumps(out, indent=2, default=str))
    pd.concat(tables_w).to_csv(C.RESULTS / "tables_tuning_windage.csv", index=False)
    pd.concat(tables_v).to_csv(C.RESULTS / "tables_tuning_variants.csv", index=False)
    if out["notes"]:
        lines += ["", "## Notes", ""] + [f"- {n}" for n in out["notes"]]
    (C.WORK / "tuning.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    if out["notes"]:
        print("\nGRID EDGE WARNING:", out["notes"])
        sys.exit(2)


if __name__ == "__main__":
    main()

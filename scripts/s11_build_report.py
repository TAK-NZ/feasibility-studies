"""Build REPORT.md from report/REPORT.template.md and the CSV files in results/.

Every table in the report is rendered here from the CSV files, so the numbers cannot drift from the files.
Kilometre values are written as they are in the CSV (up to 5 significant digits); fractions are shown as percent
with one decimal. Usage: python scripts/s11_build_report.py   (writes REPORT.md)
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLASSES = ["all", "lost", "drogued"]
LEADS = [6, 12, 24, 48]


def g(x):
    if pd.isna(x):
        return ""
    x = float(x)
    if x.is_integer():
        return f"{int(x):,}"
    return f"{x:,g}"


def pct(x):
    return "" if pd.isna(x) else f"{100 * x:.1f}%"


def ci(lo, hi, f=g):
    return "" if pd.isna(lo) or pd.isna(hi) else f"{f(lo)} to {f(hi)}"


def table(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def t_step0():
    d = pd.read_csv(RES / "step0_counts.csv", keep_default_na=False)
    rows = []
    for _, r in d.iterrows():
        def n(v):
            return "" if v == "" else f"{int(float(v)):,}"
        rows.append([r.step, r.stratum, n(r.n_drifters), n(r.n_tracks), n(r.n_windows), n(r.n_drifter_days)])
    return table(["Step", "Stratum", "Drifters", "Tracks", "Windows", "Drifter-days"], rows)


def t_sizes():
    d = pd.read_csv(RES / "metrics_sample_sizes.csv")
    d = d[((d.set == "test") & (d.stratum_type == "zone")) | ((d.set == "tuning_in_sample") & (d.stratum_type == "headline"))]
    rows = [[r.set, r.stratum, r.drogue_class, g(r.n_windows), g(r.n_drifters), g(r.n_blocks), g(r.mean_windows_per_drifter),
             g(r.icc_24h_error), g(r.n_eff), "yes" if r.anecdotal else "no"] for _, r in d.iterrows()]
    return table(["Set", "Zone", "Class", "Windows", "Drifters", "96 h blocks", "Windows per drifter", "ICC (24 h error)",
                  "n_eff (approx.)", "Anecdotal"], rows)


def t_lead():
    d = pd.read_csv(RES / "metrics_test_by_lead.csv")
    rows = []
    for c in CLASSES:
        for lead in LEADS:
            r = d[(d.drogue_class == c) & (d.lead_h == lead)].iloc[0]
            rows.append([c, lead, g(r.n_windows), g(r.n_drifters),
                         f"{g(r.err_med_model_ens)} ({ci(r.err_med_model_ens_lo, r.err_med_model_ens_hi)})", g(r.err_p90_model_ens),
                         g(r.err_med_B0), g(r.err_med_B1), g(r.err_p90_B1), g(r.err_med_B3),
                         f"{pct(r.impr_B1_model_ens)} ({ci(r.impr_B1_model_ens_lo, r.impr_B1_model_ens_hi, pct)})",
                         f"{pct(r.impr_B0_model_ens)} ({ci(r.impr_B0_model_ens_lo, r.impr_B0_model_ens_hi, pct)})"])
    return table(["Class", "Lead h", "Windows", "Drifters", "Model median km (95% CI)", "Model p90 km", "B0 median km",
                  "B1 median km", "B1 p90 km", "B3 median km", "Reduction vs B1 (95% CI)", "Reduction vs B0 (95% CI)"], rows)


def t_lead_p90():
    d = pd.read_csv(RES / "metrics_test_by_lead.csv")
    rows = []
    for c in CLASSES:
        for lead in LEADS:
            r = d[(d.drogue_class == c) & (d.lead_h == lead)].iloc[0]
            rows.append([c, lead, f"{g(r.err_p90_model_ens)} ({ci(r.err_p90_model_ens_lo, r.err_p90_model_ens_hi)})",
                         g(r.err_p90_B0), g(r.err_p90_B1), g(r.err_p90_B3)])
    return table(["Class", "Lead h", "Model p90 km (95% CI)", "B0 p90 km", "B1 p90 km", "B3 p90 km"], rows)


def t_leeway():
    d = pd.read_csv(RES / "metrics_test_by_lead.csv")
    d = d[d.drogue_class == "lost"]
    rows = []
    for lead in LEADS:
        r = d[d.lead_h == lead].iloc[0]
        rows.append([lead, g(r.err_med_model_ens), g(r.err_med_leeway_native), g(r.err_med_leeway_offset),
                     g(r.err_p90_model_ens), g(r.err_p90_leeway_native), g(r.err_p90_leeway_offset),
                     pct(r.impr_B1_leeway_native), pct(r.impr_B1_leeway_offset)])
    return table(["Lead h", "OceanDrift median km", "Leeway, native keys median km", "Leeway, with offset median km",
                  "OceanDrift p90 km", "Leeway native p90 km", "Leeway offset p90 km", "Leeway native, reduction vs B1",
                  "Leeway offset, reduction vs B1"], rows)


def t_contours():
    d = pd.read_csv(RES / "metrics_test_contours_by_lead.csv")
    rows = []
    for c in CLASSES:
        for lead in LEADS:
            r = d[(d.drogue_class == c) & (d.lead_h == lead)].iloc[0]
            rows.append([c, lead,
                         f"{pct(r.hit50_model_ens)} ({ci(r.hit50_model_ens_lo, r.hit50_model_ens_hi, pct)})",
                         f"{pct(r.hit95_model_ens)} ({ci(r.hit95_model_ens_lo, r.hit95_model_ens_hi, pct)})",
                         g(r.area50_med_model_ens), g(r.area95_med_model_ens),
                         pct(r.hit50_B1_circle), pct(r.hit95_B1_circle), g(r.area95_med_B1_circle)])
    return table(["Class", "Lead h", "Model 50% hit (95% CI)", "Model 95% hit (95% CI)", "Model median 50% area km2",
                  "Model median 95% area km2", "Circle 50% hit", "Circle 95% hit", "Circle 95% area km2"], rows)


def t_reliab():
    d = pd.read_csv(RES / "reliability_offshore_test.csv")
    d = d[d.lead_h == 24]
    noms = sorted(d.nominal.unique())
    rows = []
    for nm in noms:
        row = [pct(nm)]
        for c in CLASSES:
            r = d[(d.drogue_class == c) & (d.nominal == nm)].iloc[0]
            row.append(f"{pct(r.hit_rate)} ({ci(r.lo, r.hi, pct)})")
        rows.append(row)
    return table(["Nominal contour level"] + [f"Observed share inside, {c} (95% CI)" for c in CLASSES], rows)


def t_strata():
    d = pd.read_csv(RES / "metrics_test_by_stratum.csv")
    d = d[(d.lead_h == 24) & (d.stratum_type != "headline")]
    rows = []
    for st in ["zone", "season", "wind_class", "ecmwf_grid"]:
        for s in d[d.stratum_type == st].stratum.unique():
            for c in CLASSES:
                m = d[(d.stratum_type == st) & (d.stratum == s) & (d.drogue_class == c)]
                if m.empty:
                    continue
                r = m.iloc[0]
                rows.append([st, s, c, g(r.n_windows), g(r.n_drifters), "yes" if r.anecdotal else "no",
                             g(r.err_med_model_ens), g(r.err_med_B1),
                             f"{pct(r.impr_B1_model_ens)} ({ci(r.impr_B1_model_ens_lo, r.impr_B1_model_ens_hi, pct)})",
                             pct(r.hit50_model_ens), pct(r.hit95_model_ens), g(r.area95_med_model_ens)])
    return table(["Stratum type", "Stratum", "Class", "Windows", "Drifters", "Anecdotal", "Model median km", "B1 median km",
                  "Reduction vs B1 (95% CI)", "50% hit", "95% hit", "Median 95% area km2"], rows)


def t_variants():
    d = pd.read_csv(RES / "metrics_test_variants.csv")
    rows = []
    for _, r in d.iterrows():
        rows.append([r.drogue_class, r.current, r.wind, "on" if r.stokes else "off", g(r.windage),
                     "chosen" if r.chosen_on_tuning else ("test-best (post hoc)" if "TEST" in r.role else ""),
                     g(r.tuning_med24), f"{g(r.med24_test)} ({ci(r.med24_lo, r.med24_hi)})",
                     f"{pct(r.impr24_vs_B1)} ({ci(r.impr24_lo, r.impr24_hi, pct)})",
                     f"{g(r.med24_minus_chosen)} ({ci(r.med24_minus_chosen_lo, r.med24_minus_chosen_hi)})"])
    return table(["Class", "Current", "Wind", "Stokes", "Windage", "Selection", "Tuning-set 24 h median km",
                  "Test-set 24 h median km (95% CI)", "Reduction vs B1 (95% CI)",
                  "Median minus chosen variant, km (95% CI)"], rows)


def t_pass():
    d = pd.read_csv(RES / "pass_criteria.csv")
    rows = []
    for _, r in d.iterrows():
        crit = r.criterion.split(" ")[0]
        lead = "" if pd.isna(r.lead_h) else g(r.lead_h)
        mm = g(r.model_median_km)
        bm = g(r.baseline_median_km)
        rows.append([crit, r["item"], lead, g(r.n_windows), g(r.n_drifters),
                     mm, bm, pct(r.value), ci(r.ci_lo_wider, r.ci_hi_wider, pct),
                     "" if pd.isna(r.threshold) else pct(r.threshold), r.verdict])
    return table(["Crit.", "Item", "Lead h", "Windows", "Drifters", "Model median km", "B1 median km",
                  "Value (reduction, or hit rate)", "95% CI (wider of drifter and block)", "Threshold", "Verdict"], rows)


def t_tuning_cal():
    d = pd.read_csv(RES / "metrics_tuning_calibration_by_lead.csv")
    rows = [[r.cls, r.setting, r.l1, g(r.sigma_b), r.lead_h, g(r.n_windows), pct(r.hit50), pct(r.hit95)] for _, r in d.iterrows()]
    return table(["Class", "Setting", "Level-1 key set", "sigma_b m/s", "Lead h", "Windows", "50% hit", "95% hit"], rows)


def t_tuning_skill():
    d = pd.read_csv(RES / "metrics_tuning_insample_by_lead.csv")
    rows = []
    for c in CLASSES:
        for lead in LEADS:
            r = d[(d.drogue_class == c) & (d.lead_h == lead)].iloc[0]
            rows.append([c, lead, g(r.n_windows), g(r.n_drifters), g(r.err_med_model_det), g(r.err_med_B1), g(r.err_med_B3),
                         pct(r.impr_B1_model_det)])
    return table(["Class", "Lead h", "Windows", "Drifters", "Model (deterministic) median km", "B1 median km",
                  "B3 median km", "Reduction vs B1"], rows)


TABLES = {k[2:]: v for k, v in globals().items() if k.startswith("t_") and callable(v)}


def main():
    tpl = (ROOT / "report" / "REPORT.template.md").read_text()

    def sub(m):
        return TABLES[m.group(1)]()
    out = re.sub(r"\{\{TABLE:(\w+)\}\}", sub, tpl)
    (ROOT / "REPORT.md").write_text(out)
    print("wrote REPORT.md", len(out), "chars")


if __name__ == "__main__":
    main()

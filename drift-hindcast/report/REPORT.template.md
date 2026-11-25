# OpenDrift drift hindcast for New Zealand SAR on free forcing data: final report

Spike folder: `drift-hindcast-spike`. Tables in this report are generated from the CSV files in `results/` by `scripts/s11_build_report.py` (kilometre values as in the CSV files, up to 5 significant digits; fractions shown as percent). Figures are in `figures/`. OpenDrift 1.14.12. This spike builds nothing for production. It produces evidence for a go / no-go decision.

## 1. Result in brief

**Recommendation: NO-GO for the free-data-only plan as specified. A bounded follow-up is justified before the idea is dropped (section 9).**

| Criterion (agreed with the user) | Result on the test set | Verdict |
|---|---|---|
| 1. Offshore: median 24 h separation error at least about 30% below the hand-calculation baseline | Pooled: model 10.451 km vs baseline 13.685 km, **23.6% lower** (95% CI 20.2% to 26.8%) | **FAIL** |
| 1, lost-drogue drifters only | 10.632 km vs 14.429 km, 26.3% lower (CI 21.9% to 29.7%) | BORDERLINE |
| 1, drogued drifters only | 10.2 km vs 12.636 km, 19.3% lower (CI 13.9% to 24.9%) | FAIL |
| 2. After tuning, the 95% contour holds at least about 80% of observed endpoints (24 h, offshore) | **92.0%** (CI 90.3% to 93.6%). The 50% contour holds 56.5% (CI 52.7% to 60.5%) | **PASS**, but only with a tuned layer outside OpenDrift (section 5.3) |
| 3. Coastal skill (not pass/fail) | 18 test windows from 8 drifters | reported, **anecdotal** |

- The model does not clearly beat the baseline offshore. The whole confidence interval for the pooled 24 h reduction lies below 30%.
- The model is better than the baseline at every lead (pooled reduction 5.4% at 6 h, 17.6% at 12 h, 23.6% at 24 h, 31.3% at 48 h), so the result is "useful but short of the target", not "no skill".
- Calibration was reached by adding a tuned persistent-velocity offset outside OpenDrift. OpenDrift's own uncertainty keys alone produce contours that hold about 0% of the endpoints.
- The hindcast used best-estimate analysis currents. CMEMS keeps no archived forecasts, so **live results will be somewhat worse than the numbers here**. The size of the loss was not measured.
- The tests used 35 cm GDP drifters, not people or boats, with the start position and time known exactly. Section 8 lists what this test cannot show.
- The intended use is live SAR support for RCCNZ-style users, so uncertainty must be shown to operators together with measured hit rates, not just a contour on a map (section 9).

![Separation error by lead](figures/fig1_error_vs_lead.png)

## 2. Data and method

### 2.1 Data

| Item | Source | What was used |
|---|---|---|
| Drifter tracks | NOAA AOML Global Drifter Program, ERDDAP `drifter_6hour_qc` (6-hourly, ends 2025-06-18) | 6-hourly positions, `drogue_lost_date`, `err_lat`/`err_lon`, in lon 158 to 187 E, lat -57 to -23, from 2022-12-01. The 1-hourly dataset ends 2022-10-31, before CMEMS GLO12 starts, so it was not usable. |
| Currents | Copernicus Marine GLOBAL_ANALYSISFORECAST_PHY_001_024 | `cmems_mod_glo_phy_anfc_merged-uv_PT1H-i` (hourly, 0.494 m, with tidal components `utide`/`vtide`, and Stokes drift `vsdx`/`vsdy`); `cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i` (6-hourly, no tides, 0.494 m and 15.81 m). Analysis fields only. |
| Winds | NOAA GFS 0.25 deg (AWS `noaa-gfs-bdp-pds`, archive starts 2021-01-01) and ECMWF IFS open data (AWS `ecmwf-forecasts`, CC-BY-4.0, archive starts 2023-01-18; 2023-01-17 returns 404) | 10 m wind in forecast mode: cycle at t0 minus 12 h, forecast hours 12 to 60 in 3 h steps (00Z and 12Z operational cycles only). ECMWF is 0.4 deg until 2024-01-31 and 0.25 deg from 2024-02-11. The ECMWF 12Z cycles of 2023-05-01 and 2023-05-02 are missing from the archive. |

- **No archived forecast currents.** CMEMS does not keep archived forecasts. The hindcast uses the best-estimate analysis, which has seen later observations than a live forecast would have. Live skill will be lower. The winds are the part that was done in forecast mode.
- **Stokes drift** comes from the merged-uv `vsdx`/`vsdy` fields (same wave model, MFWAM). The separate wave product was not downloaded or compared.
- The CMEMS system versions in the period 2023 to 2025 and whether drifter data are assimilated into GLO12 were not checked (see limitations).
- Data volumes: GDP 20 MB, CMEMS 6.2 GB (216 files), winds 350 MB. Everything is under `data/` (git-ignored). Details: `work/data_notes.md`.
- Forcing sanity check against the drifters' own velocity (test set, chosen currents): vector RMSE 0.22 m/s (lost) and 0.19 m/s (drogued), vector correlation 0.62 and 0.60 (`results/forcing_check.csv`).

### 2.2 Windows and the split

- **Time sampling.** 72 blocks of 96 h (60 random, 2 per calendar month from 2023-01 to 2025-06, plus 12 targeted at coastal drifters). Each block has two 00Z seed times, 48 h apart, so windows from one drifter in a block do not overlap. There is no forcing between the blocks.
- **Window.** A 48 h track of 9 consecutive 6-hourly observations that passes the filters (position error under 0.05 deg, speed under 2.5 m/s, not after the last good time). The simulation is seeded at the observed position at t0 and compared with the observed positions at 6, 12, 24 and 48 h.
- **Strata.** Offshore: the minimum distance to land over the window is more than 50 km. Coastal: under 10 km. Mid: in between. Drogued versus drogue-lost from `drogue_lost_date`. Windows straddling the loss date were excluded (387 candidate windows from 31 drifters).
- **Tuning and test sets.** Drifters were split once, at random, 40% tuning and 60% test. A drifter is in one set only (asserted in code). Everything tuned used the tuning set only (section 6). Every number in sections 4 and 5 is a test-set number unless it says otherwise.
- **Non-independence.** Windows from one drifter, and windows on the same date, are correlated. Confidence intervals come from a cluster bootstrap that resamples whole drifters (2,000 replicates, seed 42). The headline numbers are also bootstrapped with 96 h blocks as clusters, and the wider of the two intervals is quoted. The tables show an approximate effective sample size `n_eff`, estimated from the intra-drifter correlation of the 24 h error.

### 2.3 Decisions on the open points

| Open point | Decision | Why |
|---|---|---|
| Model for lost-drogue drifters | OpenDrift `OceanDrift` with the chosen hourly current (`cur_tide`: `uo + utide`), ECMWF 10 m wind as a windage fraction of 0.5% (`wind_drift_factor` 0.005), Stokes drift off. | A GDP drifter is a 35 cm float, so there is no calibrated Leeway class for it. Windage was tuned on the tuning set from a grid of 0 to 3%. Stokes drift on was worse on the tuning set (best Stokes-on variant 11.50 km, chosen variant 11.19 km, 24 h median). |
| Model for drogued drifters | `OceanDrift`, no windage, no Stokes, chosen current `cur_notide` (hourly `uo`, `vo`, 0.494 m). | The drogue holds the float in the water at about 15 m. The 15 m current variants were tested too (section 5.5). |
| Leeway | Run as a reference for the lost class only: object class 1 (person-in-water, unknown state), same current and wind, untuned. Not part of any criterion. | Leeway is what a SAR tool would use for a person. The drifter result says nothing about a person (section 8). |
| Hand-calculation baseline (B1) | Datum, plus one drift vector `V = U_cur + k x U_wind` taken once at the datum at t0 and held for 48 h, plus an expanding circle `R(t) = a x t`. `U_cur` is the same CMEMS current (0.494 m for lost, 15.81 m for drogued). `k` = 0.005 for lost, 0 for drogued, with ECMWF wind. The circle slopes were set on the tuning set to the 50% and 95% quantiles of the 24 h B1 error (lost: 0.582 and 1.406 km/h; drogued: 0.538 and 1.236 km/h). | The baseline gets the same tuning help as the model. It uses the model's own current field, so the test measures what particle tracking adds over a one-point snapshot. It does not test against a baseline that uses no ocean model. |
| Reference lines (not criteria) | B0: stay at the datum. B3: the drifter's own mean velocity over the 12 h before t0 carried forward. | B3 needs the drifter's own earlier track, which an operator does not have at t0. It shows how much plain persistence achieves. |
| Windows and independence | See 2.2. | |
| Combining GFS and ECMWF ("both") | Arithmetic mean of u and v, after ECMWF is bilinearly interpolated to the GFS grid. | Simplest option. It was tested as a variant. It was not the chosen variant. |
| CMEMS dataset ids | The two ids above, confirmed in the catalogue (`work/describe_*.json`). | |
| Model error definition | Distance from the observed position to the centroid of the 500-particle ensemble. The deterministic-particle and no-offset variants are shown as sensitivity rows. | Headline definition fixed before those variants were computed. |
| Ensemble and contours | 500 particles per window (300 in tuning), 900 s step, Runge-Kutta, 6-hourly output. Contours are highest-density regions of a Gaussian kernel density estimate with a 3 km minimum bandwidth. | The self-test in `drift_lib.py` checks hit fractions on Gaussian samples. |

Verdict thresholds were fixed in `scripts/config.py` before the test-set numbers were looked at. Criterion 1: PASS at 30% or more, BORDERLINE from 25% to under 30% (reading "about 30%"), FAIL below 25%. Criterion 2: PASS at 80% or more, BORDERLINE from 75%, FAIL below.

## 3. Sample sizes at each step

Funnel from the raw GDP query to the evaluation sets (`results/step0_counts.csv`; generated by `scripts/s02_windows_stepzero.py` before any forcing was downloaded). For steps 1 to 5 the "Windows" column holds six-hourly observation rows, not windows. A track is a run of at least 9 consecutive six-hourly observations.

{{TABLE:step0}}

Selected windows: 4,086 in total (offshore 3,940, mid 111, coastal 35), of which 1,570 are in the tuning set and 2,516 in the test set. The sample is large enough to interpret for offshore drifters. For coastal drifters it is anecdotal. Windows without ECMWF winds (seed time 2023-05-03) have no model result, so the offshore test windows that were evaluated number 2,415 (1,447 lost, 968 drogued) against 2,427 selected. Sample sizes behind every table, with the approximate effective sample size:

{{TABLE:sizes}}

The pooled offshore test sample is 2,415 windows, but only about 595 effectively independent ones from 75 drifters. The drogued sample is smaller still (n_eff about 222). "Anecdotal" marks groups with fewer than 30 windows or fewer than 8 drifters.

## 4. Results by lead time

### 4.1 Separation error (test set, offshore)

Median error with the 95% cluster-bootstrap interval for the model, and the median for the baselines. "Reduction" is 1 - median(model error) / median(baseline error) on the same windows. B3 uses information an operator does not have.

{{TABLE:lead}}

90th percentile error:

{{TABLE:lead_p90}}

### 4.2 Contours (test set, offshore)

A contour "hit" means the observed position at that lead is inside the contour. "Circle" is the baseline's expanding circle, calibrated on the tuning set at 24 h to the same 50% and 95% targets.

{{TABLE:contours}}

![Calibration](figures/fig2_calibration.png)

Reliability at 24 h: the share of observed endpoints inside the contour at each nominal level.

{{TABLE:reliab}}

- Up to about the 80% level the contours hold slightly more than nominal. At 90% and above they hold less. The 99% contour holds 96.7% pooled. The cloud is too concentrated in the middle and its tails are too light.
- Short leads are under-dispersed. The pooled 95% contour holds 80.2% at 6 h and 88.2% at 12 h, and 74.7% for drogued drifters at 6 h. The offset grows linearly from zero, while real short-lead errors do not.
- The nominal 95% is reached only by the lost class at 48 h (95.9%). The 24 h figure that passes criterion 2 (92.0%) is below the nominal 95%.
- At 24 h the model's median 95% contour area is 2,067 km2 against 3,575.4 km2 for the baseline circle (42% smaller), at hit rates of 92.0% and 94.2%. The coverage is not exactly matched, so this indicates a smaller search area at similar coverage and is not an equal-coverage comparison.

### 4.3 By stratum, 24 h (test set)

Season is austral (DJF = summer). Wind class is the mean GFS 10 m wind at the datum over the window (m/s). Strata are shown without any multiple-comparison correction, and many intervals overlap. Do not read a ranking into them.

{{TABLE:strata}}

![Strata](figures/fig4_strata.png)

- Reductions against the baseline vary from 16.1% (June to August) to 31.4% (December to February) across seasons. The intervals overlap in places, and no season stratum was set up or powered as a test.
- The drogued class in September to November has a 95% hit rate of 79.2% (101 windows, 17 drifters), below the 80% line. Among the other non-anecdotal strata the 95% hit rate is 87.3% or higher.
- Wind class shows no clear pattern in the reduction (21.6% to 24.0% pooled).
- The ECMWF 0.4 degree and 0.25 degree periods give similar pooled results (24.0% and 23.5%).

### 4.4 Coastal and mid zones

Coastal results are the `zone = coastal` rows of the table above. Across the whole study there are 35 coastal windows from 13 drifters (17 windows in the tuning set, 18 in the test set). The 18 test windows come from 8 drifters (15 lost-class windows from 6 drifters, 3 drogued windows from 2 drifters). **This is anecdotal and is not interpreted.** The pooled coastal reduction is 45.5% with an interval of -3.7% to 61.3%, which says nothing reliable. The 1/12 degree grid (about 8 km) cannot resolve harbours, estuaries or tidal straits, and a window that is coastal by distance may still be in open water. Particles that reach land are stranded by the landmask, and the coastal lost-class ensembles had 26% stranded particles on average. Mid-zone windows (10 to 50 km from land: 71 windows from 23 drifters) are also listed above and are too few to interpret.

### 4.5 Example windows

![Example windows](figures/fig3_example_windows.png)

## 5. Pass criteria evaluated

### 5.1 Table

Intervals are 95% percentile intervals; the quoted interval is the wider of the by-drifter and by-96 h-block intervals. "Threshold" is the PASS line in `scripts/config.py`. Machine-readable version: `results/pass_criteria.csv`; commentary: `results/pass_criteria.md`.

{{TABLE:pass}}

### 5.2 Criterion 1 (offshore skill against the hand calculation): FAIL

- Pooled 24 h median error: model 10.451 km, baseline B1 13.685 km, 23.6% lower (CI 20.2% to 26.8%). The whole interval is below 30%. The point estimate is also below the 25% lower edge of the BORDERLINE band, so even the lenient reading of "about 30%" fails on the point estimate (the upper end of the interval, 26.8%, is above 25%).
- The lost class (the one closer to the SAR question) is BORDERLINE at 26.3% (CI 21.9% to 29.7%); its interval includes 25% and excludes 30%. The drogued class fails at 19.3% (CI 13.9% to 24.9%).
- The three model definitions (ensemble centroid, a single deterministic particle, the centroid without the offset) agree within 0.3 percentage points, so the verdict does not depend on that choice.
- At 48 h the pooled reduction is 31.3% (CI 27.5% to 35.4%). The agreed lead is 24 h, so this is context.
- The baseline is weaker for drogued drifters. B1 is no better than staying at the datum there (B0 median 12.862 km against B1 12.636 km at 24 h). Against B0 (no calculation at all) the pooled reduction at 24 h is 28.2% (CI 24.8% to 31.6%).

### 5.3 Criterion 2 (calibration): PASS, with a caveat that matters

- 24 h, offshore, test set: the 95% contour holds 92.0% of observed endpoints (lost 93.8%, drogued 89.4%), and the 50% contour holds 56.5% (lost 61.2%, drogued 49.6%).
- The contour was calibrated with an additional tuned layer: a persistent velocity offset applied to each particle outside OpenDrift (sigma 0.125 m/s per component for lost, 0.1 m/s for drogued), chosen on the tuning set so that the 24 h 95% hit rate was at least 0.90 with the smallest area. The native keys that were tried (current and wind uncertainty, windage and current-scaling spread) did not produce enough spread in OpenDrift 1.14.12. With only the chosen native uncertainty keys, the 95% contour holds 0.0% pooled (lost 0.1%, drogued 0.0%) at 24 h. On the tuning set the best native-only setting held 22.8% (lost) and 38.4% (drogued).
- The offset was fitted on 40% of the drifters in this region and period and tested on the other 60%. Whether it transfers to other regions, seasons or years was not tested.
- Reliability is not uniform (section 4.2): the pass is specific to the 95% level at 24 h.

### 5.4 Criterion 3 (coastal): reported, anecdotal

See section 4.4. Not a pass/fail criterion.

### 5.5 Forcing variants

Every deterministic variant was run on the test set. The windage of each lost-class variant was tuned on the tuning set. "Chosen" marks the variant selected on the tuning set. The test-best row is selection on test data, so it is an optimistic ceiling. The last column is the paired difference to the chosen variant.

{{TABLE:variants}}

![Variants](figures/fig5_variants.png)

- Across the 18 lost-class variants the test-set 24 h median ranges from 10.506 to 11.303 km, and across the five drogued-class variants from 9.7976 to 10.263 km. For the 21 paired comparisons against the chosen variant, none has an interval that excludes zero.
- Even the best configuration chosen on the test set reaches only 27.2% (lost) and 22.5% (drogued). No forcing choice among the free options closes the gap to 30%.
- Tides, Stokes drift, GFS against ECMWF or their mean, and hourly against 6-hourly currents did not separate the variants. Stokes drift on made the lost class slightly worse (about 0.6 km at 24 h in the table).
- Drogued class: the 15 m current variants (the drogue depth) had the lowest test-set medians (22.2% and 22.5% reduction), but on the tuning set they were within the 2% tie band, and the pre-set tie rule (lower 48 h median) picked the surface no-tide current (19.3%). The rule was not changed after the test set was seen.
- Leeway class 1 (person in water, unknown state; untuned), lost class, same current and wind:

{{TABLE:leeway}}

At 24 h Leeway class 1 and the tuned OceanDrift configuration are indistinguishable (10.591 and 10.632 km). Leeway is a little better at 48 h (18.901 km against 20.162 km). No interval was computed for that difference. This is a float, not a person, so it is not evidence about Leeway skill for people.

## 6. What was tuned, and on which data

All tuning used the tuning set only (40% of drifters, offshore windows: lost 985 windows from 33 drifters, of which 975 have both wind sources; drogued 528 windows from 29 drifters). No test-set window or drifter was used for any of it, except the explicitly labelled post hoc test-best rows.

| Quantity | How chosen | Metric on the tuning set |
|---|---|---|
| Lost class: windage fraction | Grid 0, 0.005, 0.01, 0.015, 0.02, 0.03 per variant | 24 h median error |
| Lost class: current, wind source, Stokes on/off | Best of 18 variants. Ties within 2% resolved by the 48 h median. Chosen: `cur_tide`, ECMWF, Stokes off, windage 0.005 (9 variants within the tie band). | 24 h median error |
| Drogued class: current | Best of 5 variants, same tie rule. Chosen: `cur_notide` (5 within the tie band). | 24 h median error |
| Baseline: windage, wind source, circle slopes | Grid and quantiles of the baseline's own 24 h error. | 24 h median and quantiles |
| Uncertainty: OpenDrift keys | Lost: `drift:current_uncertainty` 0, `drift:wind_uncertainty` 3.0 m/s, no windage spread. Drogued: `drift:current_uncertainty` 0, no current-scaling spread. | 24 h 95% hit rate and area |
| Uncertainty: persistent velocity offset (outside OpenDrift) | Grid 0 to 0.3 m/s per component. Rule: smallest mean 95% contour area at 24 h among settings with a 24 h 95% hit rate of at least 0.90. Chosen: 0.125 m/s (lost), 0.1 m/s (drogued). The grid was refined once (0.125 added) after the first tuning pass, before any test-set run. | 24 h hit rate and area |

Fixed by hand, not tuned: 900 s time step, 6-hourly output, Runge-Kutta, 300 (tuning) and 500 (test) particles, the 3 km kernel floor, the grids themselves, `seed:jibe_probability` (default). Leeway was not tuned. Random seeds are fixed.

Tuning-set numbers are in-sample and are not evidence of skill. They are kept here for comparison with the test set.

{{TABLE:tuning_skill}}

Tuning-set hit rates of the chosen uncertainty setting (in-sample):

{{TABLE:tuning_cal}}

## 7. Limitations

- **Analysis fields, not forecasts.** The hindcast currents are the best estimate from the ocean analysis. CMEMS keeps no archived forecasts, so live skill will be worse by an unmeasured amount. The test-set winds are forecasts (12 to 60 h from the cycle 12 h before t0), which is the realistic case.
- **Possible optimistic bias from assimilation.** Whether GLO12 assimilates drifter data was not checked. If it does, the hindcast is partly not independent of the observations.
- **Float, not a SAR object.** GDP drifters are 35 cm floats with a drogue (or without one). Their windage, shape and response to waves are unlike a person, liferaft or small boat. The tuned windage of 0.5% does not transfer.
- **Perfect datum.** Every window starts from an exactly known position and time. A real datum has neither.
- **Region and period.** NZ box 160 to 185 E, 55 to 25 S, January 2023 to June 2025, sampled as 72 blocks of 96 h. This is a sample of the period, not the whole of it. Offshore windows dominate. There is no cyclone or storm-specific analysis. The highest wind class (mean wind 10 m/s or more) has 539 windows.
- **Effective sample size.** The pooled offshore test sample is 2,415 windows and about 595 effectively independent ones from 75 drifters. The drogued sample (about 222) is smaller and its intervals are wide.
- **Tuned layer.** The calibration relies on an offset layer outside OpenDrift, tuned on 40% of drifters in this region. Calibration is not uniform across levels or leads (section 4.2), and the 99% contour would not be reliable.
- **Baseline definition.** B1 is one reasonable "hand calculation" built from the model's own current field. It is not RCCNZ's actual method. A different baseline would change the percentages. For drogued drifters it is weaker than no calculation at all.
- **Strata.** Season, wind-class and ECMWF-grid results have no multiple-comparison correction, and overlapping intervals are common.
- **Stokes drift** comes from the merged-uv product, not from the separate wave product. Not compared.
- **ECMWF coverage.** ECMWF open data starts 2023-01-18 with a 0.4 degree grid until 2024-01-31. Two 12Z cycles in May 2023 are missing, which removed 22 lost-class windows (12 offshore test windows) from every ECMWF run.
- **Contours** are kernel-density estimates of finite ensembles. A 3 km floor on the bandwidth affects 6 h areas.
- **B3** and Leeway use information or settings that are not comparable with the main model, and are shown for context only.
- Two small defects in the generated figures: in `fig5_variants.png` the drogued-panel row labels extend into the left panel, and the legend covers the test-best marker. They do not change any number.

## 8. What this test cannot show

- **Live forecast skill.** It used analysis currents. How much worse a forecast is was not measured, and cannot be from CMEMS archives.
- **Skill for people, boats or liferafts.** No SAR object was tracked. The model, windage and Leeway settings were tested on 35 cm drifters only.
- **Skill when the datum is uncertain.** Start position and time were exact. Real datums carry errors of unknown size, which add to everything here.
- **Coastal, harbour, estuary and strait performance.** Only 18 coastal test windows from 8 drifters, and the 1/12 degree grid cannot resolve those places.
- **Whether higher-resolution forcing would reach the 30% target.** It was not tested. The data show only that the free 1/12 degree product did not.
- **Whether the calibration transfers.** The offset was tuned on one region and period. It was not tested on other regions, years or different drifter types.
- **Heavy-weather and rare-event performance.** No cyclone or storm-specific analysis was done.
- **Operational value.** Whether a 24% smaller median error, or a 42% smaller contour at similar coverage, changes a search plan or a search outcome was not tested, and the baseline is not the method any rescue coordination centre uses.
- **A ranking between forcing options.** None of the free options separated from the others at this sample size. The differences between them are within noise, so "no difference found" does not mean "no difference".

## 9. Recommendation and implications for the free-data plan

**No-go for the free-data-only plan as specified.** The agreed criterion 1 failed (23.6%, interval entirely below 30%), and the one criterion that passed (calibration) passed through an extra tuned layer, not through OpenDrift's own uncertainty settings. The test does not support presenting this as a validated operational drift service.

It is not a "no skill" result, which is why a bounded follow-up is worth considering:

- The model is better than the hand calculation at every lead, by 23.6% at 24 h and 31.3% at 48 h, and the 95% contour is about 42% smaller than the baseline circle at similar coverage.
- The lost class is borderline (26.3%, interval 21.9% to 29.7%).
- The pooled median is 0.87 km short of the 30% line (9.58 km needed).
- Free forcing choices (current dataset, tide, Stokes, wind source) do not move the result, so more free-data variants will not close the gap.
- B3, which uses the drifter's own earlier velocity, reaches 7.522 km at 24 h. That points to part of the error being a persistent drift that the 1/12 degree analysis misses. This is an inference from the data. No remedy was tested.

**Implications for the free-data plan.**

1. Free data (CMEMS currents, GFS and ECMWF winds) is enough to run the model and, after tuning, to produce contours with a measured hit rate at 24 h. It was not enough, in this test, to reach the agreed skill target against the hand calculation.
2. What would need to change to reach the target (all untested):
   - higher-resolution or regional current forcing, or a different or additional ocean model for eddy and boundary-current fields, which may not be free;
   - anchoring the cloud to an observed drift (for example a datum marker buoy or a drifter report), which addresses the persistent drift that B3 exposes;
   - a hindcast with archived forecast currents, which can only lower the numbers here, so it would be needed before any live claim.
3. If a prototype is built anyway for RCCNZ-style users, treat it as decision support under these conditions:
   - show the contour with its measured hit rate at that lead (for example, about 92% for the 95% contour at 24 h, lower at 6 and 12 h) and say that it was tuned on drifters, not on SAR objects;
   - do not present the contours as probabilities at levels above 90%;
   - state that the numbers come from analysis currents and live results will be worse;
   - re-run this validation for the SAR objects and region of concern, and for the live forecast configuration, before relying on it.
4. A decision between "stop" and "follow-up" belongs to the user. The cheapest follow-up that could change this recommendation is a test of one higher-resolution current product on the same windows, using the same scripts.

## 10. Reproduction and file index

`README.md` has one command sequence for setup, download, runs, metrics, plots and this report. Key files:

- Scripts: `scripts/s01` to `s11`; shared code `drift_lib.py`, `runner.py`, `config.py`.
- Tables: `results/metrics_*.csv`, `results/pass_criteria.csv`, `results/reliability_offshore_test.csv`, `results/step0_counts.csv`.
- Notes: `work/plan.md`, `work/step0_counts.md`, `work/data_notes.md`, `work/baseline.md`, `work/tuning.md`, `work/run_notes.md`.
- Licence: GPL-2.0 (OpenDrift is GPL-2.0-only). No code was copied from other TAK-NZ repositories.

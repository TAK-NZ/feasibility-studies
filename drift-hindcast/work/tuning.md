# Tuning record (part A: windage, variant choice, baseline)

Tuned on the **tuning set only**: windows with `split == tuning` and `zone == offshore` (more than 50 km from land). Tuning and test drifters are disjoint (checked in `s05_build_forcing.py`).

- lost: 985 offshore tuning windows from 33 drifters
- drogued: 528 offshore tuning windows from 29 drifters

## Lost-drogue class (OceanDrift, tuned windage)
Windows used for the comparison: 975 (windows that have both wind sources). Windage grid: [0.0, 0.005, 0.01, 0.015, 0.02, 0.03].

Best windage and tuning-set 24 h / 48 h median error (km) for every current / wind / Stokes variant:

| current | wind | Stokes | best windage | median 24 h km | median 48 h km | within 2% of best |
|---|---|---|---|---|---|---|
| cur_6h_surf | ecmwf | off | 0.005 | 11.13 | 21.13 | yes |
| cur_6h_surf | gfs | off | 0.005 | 11.17 | 21.03 | yes |
| cur_6h_surf | both | off | 0.005 | 11.18 | 21.02 | yes |
| cur_notide | ecmwf | off | 0.005 | 11.18 | 20.97 | yes |
| cur_tide | ecmwf | off | 0.005 | 11.19 | 20.96 | yes |
| cur_notide | both | off | 0.005 | 11.21 | 21.03 | yes |
| cur_tide | both | off | 0.005 | 11.21 | 21.00 | yes |
| cur_notide | gfs | off | 0.005 | 11.22 | 20.96 | yes |
| cur_tide | gfs | off | 0.005 | 11.26 | 20.97 | yes |
| cur_6h_surf | both | on | 0.0 | 11.50 | 21.59 |  |
| cur_6h_surf | ecmwf | on | 0.0 | 11.50 | 21.59 |  |
| cur_6h_surf | gfs | on | 0.0 | 11.50 | 21.59 |  |
| cur_notide | gfs | on | 0.0 | 11.57 | 21.68 |  |
| cur_notide | ecmwf | on | 0.0 | 11.57 | 21.68 |  |
| cur_notide | both | on | 0.0 | 11.57 | 21.68 |  |
| cur_tide | both | on | 0.0 | 11.57 | 21.68 |  |
| cur_tide | ecmwf | on | 0.0 | 11.57 | 21.68 |  |
| cur_tide | gfs | on | 0.0 | 11.57 | 21.68 |  |

**Selected**: cur_tide, wind ecmwf, Stokes off, windage 0.005 (24 h median 11.19 km; 9 variants were within 2%).

## Drogued class (OceanDrift, no windage, no Stokes)

Windows: 528. Tuning-set median errors (km):

| current | median 24 h km | median 48 h km | within 2% of best |
|---|---|---|---|
| cur_6h_15m_tide | 9.26 | 17.75 | yes |
| cur_6h_15m | 9.29 | 17.81 | yes |
| cur_notide | 9.38 | 17.35 | yes |
| cur_tide | 9.43 | 17.50 | yes |
| cur_6h_surf | 9.43 | 17.55 | yes |

**Selected**: cur_notide (24 h median 9.38 km; 5 variants within 2%).

## Baseline (hand calculation), see work/baseline.md

B1 windage grid, tuning-set 24 h median error (km), wind source ecmwf: 0.0: 14.02, 0.005: 13.97, 0.01: 14.60, 0.015: 15.86, 0.02: 17.64, 0.03: 22.20
- lost: B1 windage 0.005, wind ecmwf; circle radius slopes a_0.5 = 0.582 km/h, a_0.95 = 1.406 km/h; tuning 24 h median error B1 13.97 km, B0 15.83 km, B3 8.98 km
- drogued: B1 windage 0.0, wind none; circle radius slopes a_0.5 = 0.538 km/h, a_0.95 = 1.236 km/h; tuning 24 h median error B1 12.91 km, B0 12.25 km, B3 5.44 km

## Tuning-set skill (in-sample, optimistic)

- lost 6 h: model median 3.77 km, B1 3.91 km, improvement 3.6%
- lost 12 h: model median 6.44 km, B1 7.32 km, improvement 11.9%
- lost 24 h: model median 11.19 km, B1 13.97 km, improvement 19.9%
- lost 48 h: model median 20.96 km, B1 29.14 km, improvement 28.1%
- drogued 6 h: model median 3.34 km, B1 3.52 km, improvement 5.3%
- drogued 12 h: model median 5.22 km, B1 6.73 km, improvement 22.4%
- drogued 24 h: model median 9.38 km, B1 12.91 km, improvement 27.4%
- drogued 48 h: model median 17.35 km, B1 25.39 km, improvement 31.7%

## Part B: uncertainty tuning (tuning set only, offshore, 24 h)

Particles per drifter: 300. Level 1 = native OpenDrift keys; Level 2 = post hoc persistent velocity offset (sigma_b in [0.0, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3] m/s per component). Selection rule: smallest mean 95% contour area at 24 h among settings with a 24 h 95% hit rate of at least 0.9; if none reaches it, the highest hit rate. Areas use a fixed subsample of tuning windows (every k-th window of the sorted ids, about 300 per class; the count is `n_area` in results/tables_tuning_uncertainty.csv).

### lost: 18 Level-1 settings x 9 offsets, 975 windows
Selected (smallest area with hit95 >= target): Level 1 `cu0.0_wu3.0_rw0.0`, sigma_b = 0.125 m/s; tuning-set 24 h hit rate 50%: 0.581, 95%: 0.920, mean 95% area 2220 km2.
Native only (sigma_b = 0), best Level-1 setting `cu0.3_wu0.0_rw0.3`: hit 50%: 0.062, 95%: 0.228.

Hit rates at 24 h by sigma_b for the selected Level-1 setting:

| sigma_b | hit 50% | hit 95% | mean 95% area km2 |
|---|---|---|---|
| 0.0 | 0.000 | 0.000 | 0 |
| 0.03 | 0.053 | 0.197 | 126 |
| 0.05 | 0.143 | 0.462 | 350 |
| 0.075 | 0.273 | 0.692 | 792 |
| 0.1 | 0.440 | 0.849 | 1415 |
| 0.125 | 0.581 | 0.920 | 2220 |
| 0.15 | 0.673 | 0.970 | 3208 |
| 0.2 | 0.823 | 0.994 | 5731 |
| 0.3 | 0.956 | 1.000 | 12939 |

### drogued: 6 Level-1 settings x 9 offsets, 528 windows
Selected (smallest area with hit95 >= target): Level 1 `cu0.0_rc0.0`, sigma_b = 0.1 m/s; tuning-set 24 h hit rate 50%: 0.545, 95%: 0.907, mean 95% area 1430 km2.
Native only (sigma_b = 0), best Level-1 setting `cu0.3_rc0.2`: hit 50%: 0.110, 95%: 0.384.

Hit rates at 24 h by sigma_b for the selected Level-1 setting:

| sigma_b | hit 50% | hit 95% | mean 95% area km2 |
|---|---|---|---|
| 0.0 | 0.000 | 0.000 | 1 |
| 0.03 | 0.074 | 0.288 | 127 |
| 0.05 | 0.203 | 0.576 | 354 |
| 0.075 | 0.384 | 0.777 | 801 |
| 0.1 | 0.545 | 0.907 | 1430 |
| 0.125 | 0.667 | 0.953 | 2244 |
| 0.15 | 0.756 | 0.991 | 3243 |
| 0.2 | 0.877 | 0.994 | 5796 |
| 0.3 | 0.985 | 1.000 | 13085 |

## Part C: what was tuned, and on what data

Data for every tuned quantity below: the tuning split only (drifters assigned by `s02`, 40 percent of drifters), offshore windows (more than 50 km from land): lost class 985 windows from 33 drifters (975 of them have both wind sources), drogued class 528 windows from 29 drifters. Test-set windows and drifters were not used for any of it.

| quantity | how chosen | metric on the tuning set |
|---|---|---|
| lost class: windage fraction of the OceanDrift model | grid `config.WINDAGE_GRID`, per variant | 24 h median separation error |
| lost class: current, wind source, Stokes on/off | best of 18 variants; ties within 2 percent resolved by the 48 h median | 24 h median separation error |
| drogued class: current variant | best of 5; same tie rule | 24 h median separation error |
| baseline: windage k, wind source (same as the model), circle slopes a_0.5 and a_0.95 | grid / quantile of the 24 h baseline error | 24 h median, 24 h quantiles |
| uncertainty: Level-1 settings and sigma_b (post hoc offset), per class | rule above | 24 h 95 percent hit rate and 95 percent area |

Fixed by hand, not tuned: time step 900 s, output step 6 h, Runge-Kutta advection, 300 (tuning) and 500 (test) particles per drifter, the 3 km kernel floor of the contour estimator, the grids themselves. The sigma_b grid was refined once (0.125 m/s added) after the first tuning pass showed the 0.1 to 0.15 step was coarse (tuning data only, before any test-set run). Leeway was not tuned (object class 1).


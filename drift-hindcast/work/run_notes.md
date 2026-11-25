# Run notes: simulation step

All runs: OpenDrift 1.14.12, `opendrift.readers.reader_netCDF_CF_generic`, time step 900 s, output every 6 h, Runge-Kutta advection, no vertical motion, 48 h, fixed random seeds (`config.SIM_SEED`), `seed:ocean_only=False` (coastal drifters are seeded where they were observed; the roaring landmask strands particles that reach land), `seed_elements` at the observed position with radius 0.

## Windows
4,086 selected windows (`results/windows.csv`): tuning 1,570 (1,029 lost-drogue, 541 drogued), test 2,516 (1,524 lost-drogue, 992 drogued). By zone (this table fixes a number in the data step's hand-off, which said 149 mid windows; the file has 111): offshore 3,940, mid 111, coastal 35. Seed times: 143. A drifter is in one split only (asserted in `s05_build_forcing.py`).
Wind-speed class (GFS 10 m wind at the datum, mean of the 17 forecast steps), counts in the file: columns `wind_class`, `gfs_wind_ms`.

## Variants that ran
| group | variants | windows | OpenDrift runs |
|---|---|---|---|
| deterministic, lost class | 3 currents (`cur_tide`, `cur_notide`, `cur_6h_surf`) x 3 winds (gfs, ecmwf, both) x Stokes on/off, 6 windage values carried as separate elements in each run = 108 configurations | 2,553 with GFS, 2,531 with ECMWF and "both" | 2,562 (858 for each of gfs, ecmwf, both = 143 seed times x 6) |
| deterministic, drogued class | 5 currents (`cur_tide`, `cur_notide`, `cur_6h_surf`, `cur_6h_15m`, `cur_6h_15m_tide`), windage 0, Stokes off | 1,533 | 655 (131 seed times x 5) |
| baselines B0, B1, B3 | B1 for every windage and wind source | all 4,086 | none (closed form) |
| ensemble, tuning | lost 18 Level-1 settings, drogued 6 Level-1 settings, 300 particles per drifter, offshore tuning windows only; 9 post hoc offsets each | 985 lost (975 with ECMWF) + 528 drogued | 3,306 |
| ensemble, test | the chosen setting per class, 500 particles per drifter, all zones; OceanDrift for both classes, Leeway (object class 1, PIW unknown state) for the lost class | 1,512 lost + 992 drogued | 415 (OceanDrift and Leeway) |

Leeway ran as a reference with the chosen current and wind and the native uncertainty keys of the chosen setting (`drift:current_uncertainty`, `drift:wind_uncertainty`), once with no offset and once with the same post hoc offset. `seed:jibe_probability` was left at its default (0.04); it was not varied.

## Variants skipped or reduced, and why
- **ECMWF-only and "both" wind runs** have no runs for seed time 2023-05-03 (the ECMWF 12Z cycles of 2023-05-01 and 2023-05-02 are not in the archive). That removes 22 lost-class windows from every ECMWF/"both" configuration. The chosen lost-class model uses ECMWF, so those windows (12 offshore test windows, 10 offshore tuning windows, plus mid/coastal ones) have no model result in `results/ens_metrics.csv`; the baseline (also ECMWF) has no row for them either. Windows without ECMWF are flagged `has_ecmwf = False` in `results/windows.csv`.
- **Stokes drift on** exists only for the lost class (OceanDrift). Drogued windows run with Stokes off, as decided in the plan. Leeway does not use Stokes drift.
- **Separate CMEMS wave product** was not downloaded; Stokes drift comes from the merged-uv `vsdx`/`vsdy` fields (data step decision).
- **Level-1 uncertainty grid is run for the selected model only**, not for all 18 variants (the grid is 18 + 6 settings).
- **The tuning-set uncertainty selection** uses the 24 h hit rate pooled over windows. Area is evaluated on every k-th tuning window (about 300 per class, `n_area` in `results/tables_tuning_uncertainty.csv`) because contour areas on a grid cost about 30 ms each.
- **Mid-zone and coastal windows** were run, but they are not used for tuning. The 'coastal' test windows are anecdotal (15 lost and 3 drogued windows).
- **Drogued class baseline**: B1 is slightly worse than B0 (stay) on the tuning set (see `work/baseline.md`).

## Problems found and fixed in this step
- One GFS wind file (`data/winds/gfs/20230127T00.nc`) was a truncated file with only `u10`, and one (`20231208T00.nc`) had NaN latitude/longitude coordinates. Both came from the segmentation-fault runs of the data step and were not caught by its load test. Both were deleted and re-fetched; `s04_wind_download.py` now also checks both variables, the 17 steps and finite coordinates. All 865 forcing files and all 287 raw wind files were checked for monotonic, finite coordinates afterwards.
- `cur_6h_15m_tide` was first built with an outer join on slightly different float32 coordinates; it was caught by the smoke test and rebuilt on the hourly file's grid.
- The `hdr_level`/`hdr_area_km2` helpers from the plan's step 1 were replaced by `hdr_levels`/`kde_density`/`hdr_hits`/`hdr_areas` (leave-one-out densities, kernel covariance n^(-1/3) Cov + 3 km floor); the self-test checks the 95% and 50% hit fractions on Gaussian samples with n = 300 (0.946, 0.497) and n = 5000 (0.958, 0.507).

## Checks
- Forcing against drifter velocity (`results/forcing_check.csv`): current RMSE 0.19 to 0.23 m/s, vector correlation 0.55 to 0.65. No unit, sign or depth error.
- Benchmark (`results/benchmark.txt`): a run takes about 1.2 to 1.7 s; reader reuse gives identical results to fresh readers; the 6 h displacement of a drogued drifter in a zero-wind, zero-Stokes `cur_notide` run is 0.98 times the displacement predicted by the forcing velocity at the datum (correlation 0.80).
- Deactivation: 0.54 percent of deterministic lost-class elements are stranded at 48 h (all coastal), none were deactivated for missing data. In the test ensembles no offshore or mid window has a stranded or deactivated particle; coastal lost windows have 26 percent stranded particles on average.
- `s07_tune.py` and the tuning part of `s08` were re-run once and reproduce `results/tuning_choices.json` byte for byte.

## Raw outputs
- `data/runs/det/<seed>.parquet`, `data/runs/det_errors.parquet`: per window, method (OD, B0, B1, B3), current, wind, Stokes, windage, lead: lon, lat, error (km), status. 1,348,800 rows. Git-ignored.
- `data/runs/ens/tune_<seed>.npz`: particle displacements (km, east/north relative to the datum) at 6, 12, 24, 48 h for every tuning ensemble setting. `data/runs/ens/tune_hits.parquet`: per window, setting, offset, lead hit flags and areas. Git-ignored.
- `data/runs/ens/test_<seed>.npz` and `.parquet`: test-set particle displacements (OceanDrift, Leeway) and per-window contour rows. Git-ignored.
- Committed: `results/windows.csv`, `results/baseline_inputs.csv`, `results/forcing_check.csv`, `results/tuning_choices.json`, `results/tables_tuning_*.csv`, `results/ens_metrics.csv` (one row per window x method x lead for the test set: centroid error, 50/95 hit, 50/95 area, fraction of active particles), `results/benchmark.txt`.
- The deterministic test-set errors of the chosen configuration are in `data/runs/det_errors.parquet` (method `OD`, filter on the chosen current, wind, Stokes, windage).

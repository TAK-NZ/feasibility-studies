# Baseline: the hand calculation an operator would do

Implemented in `scripts/s06_run_deterministic.py` (`baseline_rows`) and `scripts/s07_tune.py` (constants), applied in `scripts/s08_run_ensemble.py` (`baseline_circle_rows`).
It uses only what an operator has at t0: the datum, a current and a wind at the datum, and a stored rule of thumb for the circle.

## Formula

For each window, with the datum = the observed drifter position at t0:

```
position(t) = datum + V * t          (geodesic on the WGS84 ellipsoid, V constant for the whole 48 h)
V           = U_cur + k * U_wind     (m/s; evaluated once, at the datum, at t0)
```

| term | lost-drogue windows | drogued windows |
|---|---|---|
| `U_cur` | hourly Eulerian current `uo`, `vo` at 0.494 m (`cur_notide`: no tide, no Stokes), bilinear in space at the datum, time t0 | 6-hourly current `uo`, `vo` at 15.81 m (`cur_6h_15m`, the drogue depth), same sampling |
| `U_wind` | forecast-mode 10 m wind valid at t0 (cycle t0 - 12 h, forecast hour 12) at the datum. The source is the one the primary model uses (chosen on the tuning set: ECMWF) | none (k = 0) |
| `k` (windage) | chosen from the grid `config.WINDAGE_GRID` = 0, 0.005, 0.01, 0.015, 0.02, 0.03 as the value with the lowest **tuning-set** 24 h median B1 error: **0.005** | 0 |

Expanding circle (the search-area part of the baseline):

```
R_p(t) = a_p * t            (km, t in hours)
a_p    = quantile_p( |observed position at 24 h - B1 position at 24 h| over the tuning windows ) / 24 h
```

with p = 0.5 and p = 0.95, per class, offshore tuning windows only. Constants chosen on the tuning set (`results/tuning_choices.json`):

| class | a_0.5 (km/h) | a_0.95 (km/h) | tuning windows |
|---|---|---|---|
| lost-drogue | 0.582 | 1.406 | 985 offshore tuning windows, 33 drifters |
| drogued | 0.538 | 1.236 | 528 offshore tuning windows, 29 drifters |

A circle "hit" at lead t means the observed position is within R_p(t) of the B1 position. The baseline gets the same calibration help as the model (a tuned windage and tuned circle growth rates), on the same tuning data.

## Reference lines (not part of any pass criterion)

- **B0**: stay at the datum.
- **B3**: the drifter's own mean velocity over the 12 h before t0 (`ve_pre`, `vn_pre`) carried forward. A real operator does not have this (it needs the drifter's own earlier track). It shows how much skill plain persistence has.

## Where B1 stands (tuning set, offshore, 24 h median error, in-sample for the tuned constants)

| class | B0 stay | B1 hand calculation | B3 persistence |
|---|---|---|---|
| lost-drogue | 15.83 km | 13.97 km | 8.98 km |
| drogued | 12.25 km | 12.91 km | 5.44 km |

For drogued windows B1 is slightly **worse** than leaving the datum where it is: a current at one point and one time does not hold for 24 to 48 h. The criterion in the brief is measured against B1, as agreed, so this makes B1 an easy baseline for the drogued class. B3 is much stronger than both B1 and the model (see `work/tuning.md`), so the report must not read "beats B1" as "beats every simple method".

## Choices to note

- The baseline's current is the 0.083 degree model field itself, so the comparison tests what the particle model adds over a single-point snapshot of the same data. It does not test the model against a baseline that uses no model at all (for example a chart of mean currents).
- Land (NaN) cells at the datum are filled with the mean of valid cells within two grid cells (`runner.sample_field`); the count is in `results/baseline_inputs.csv` (`cur_filled`).

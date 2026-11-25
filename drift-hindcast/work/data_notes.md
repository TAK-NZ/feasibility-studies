# Data notes: what was obtained

Downloaded 2026-11-25 by `scripts/s01`, `s03`, `s04`. Everything sits under `data/` (git-ignored, about 6.6 GB for forcing).
Time sampling: 72 blocks of 96 h (60 random, 2 per calendar month from 2023-01 to 2025-06, plus 12 coastal-targeted blocks),
144 seed times at 00Z (block start and block start + 48 h). Block starts run from 2023-01-21 to 2025-06-13. The days between blocks have no forcing.

## Drifters (NOAA GDP ERDDAP `drifter_6hour_qc`)
- Query: lon360 158..187, lat -57..-23, time 2022-12-01 to 2025-06-18T12Z (dataset end). 140,362 rows, 172 drifters, 20 MB.
- Buoy types in the box: SVPB 155, SVP 17.
- No position quality flag exists. The filter is `err_lat`/`err_lon` < 0.05 deg plus `time <= end_date`.
- `ve`/`vn` are missing (-999999) at the first row of each track.
- Drogue state: `drogue_lost_date` (empty = still attached, 1970-01-01 = uncertain; no candidate window had the 1970 value, so none was excluded for it).

## CMEMS (login via the stored Copernicus Marine account; credentials not read or printed)
| file type | dataset id (version) | files | per-file | total | time axis per file | levels |
|---|---|---|---|---|---|---|
| `uv1h_<block>.nc` | `cmems_mod_glo_phy_anfc_merged-uv_PT1H-i` (202211): `uo vo utide vtide vsdx vsdy` | 72 | 54 MB | 3.9 GB | 97 hourly steps | 0.494 m |
| `cur6h_s_<block>.nc` | `cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i` (202406): `uo vo` | 72 | 16 MB | 1.14 GB | 17 steps at 00/06/12/18Z | 0.494 m |
| `cur6h_15_<block>.nc` | same dataset | 72 | 16 MB | 1.14 GB | 17 steps | 15.81 m |

- Total 6.19 GB (the plan estimated about 4 GB; the 6-hourly files do not compress much). Box 158..187E, -57..-23, longitude stored as 0..360 up to 187 so the dateline needs no special handling.
- All 216 files: start and end times equal the block window, depth as expected, no file with all-NaN `uo`. Land is NaN (about 3 percent of cells).
- `data/cmems/manifest.csv` lists every file.
- Stokes drift comes from the merged-uv `vsdx`/`vsdy` for all variants. The separate wave product (`cmems_mod_glo_wav_anfc_0.083deg_PT3H-i`, `VSDX`/`VSDY`) was not downloaded, as decided in the plan, so the brief's wave product was substituted by the merged-uv Stokes fields from the same wave model (MFWAM). The wave product was not compared to them.
- These are best-estimate analysis fields (CMEMS keeps no archived forecasts). Not checked: GLO12 system changes during 2023 to 2025, and whether drifters are assimilated.

## Winds (forecast mode: cycle t0 - 12 h, forecast hours 12..60 step 3, 17 steps, 10 m U/V)
| model | files | range of seed times | grid | note |
|---|---|---|---|---|
| GFS 0.25 deg (AWS `noaa-gfs-bdp-pds`) | 144 of 144 | 2023-01-21 to 2025-06-15 | 0.25 deg | 247 MB. Archive starts 2021-01-01 (probed). |
| ECMWF IFS open data (AWS `ecmwf-forecasts`, CC-BY-4.0) | 143 of 144 | 2023-01-21 to 2025-06-15 | 0.4 deg (path `0p4-beta`) for 69 files up to 2024-01-31; 0.25 deg for 74 files from 2024-02-11 (path `0p25` for 3 files, `ifs/0p25` for 71) | 104 MB. Archive starts 2023-01-18 (probed: 2023-01-17 is a 404). |

- **Gap:** the ECMWF 12Z cycles of 2023-05-01 and 2023-05-02 are not in the archive (the `20230502` folder is empty). That affects seed time 2023-05-03T00Z only (block 20230503): 35 windows (13 offshore drogued, 22 offshore lost). The ECMWF-only and "both" wind variants have no data for them. Later steps must exclude these 35 windows from any ECMWF or "both" comparison, or run GFS-only variants on them and say so.
- Both models: 00Z/12Z `oper` (HRES) only. 06Z/18Z cycles were not used.
- ECMWF values were checked against GFS for three seed times (2023-01-21, 2024-02-28, 2025-06-15): domain-mean absolute wind-speed difference 0.89, 1.09 and 1.22 m/s; direction at (-40, 170) within 16 to 32 degrees. No unit or sign error.
- Not obtained: the optional GFS analysis-chain wind variant (`w_gfs_ana`, stretch item in the plan). Not downloaded because forecast-mode winds are what the study uses and nothing in this step needs it.
- `data/winds/manifest.csv` lists every file with its path variant and grid.

## Other notes
- Parallel threads crashed the first wind runs (eccodes and HDF5 are not thread-safe, segmentation fault). `s04` now uses processes and writes files through a temp name; truncated files from the crashed runs were detected and re-fetched.
- Nothing needing a login other than the Copernicus Marine account was used.

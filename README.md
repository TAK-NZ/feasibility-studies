# Drift hindcast spike

Goal: decide whether an OpenDrift-based ocean drift / search-area service for New Zealand SAR, built only on free forcing data, is worth building, by hindcasting against NOAA GDP drifter tracks.

## Setup

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/copernicusmarine login --check-credentials-valid   # run `copernicusmarine login` first if invalid
.venv/bin/python -c "import opendrift, opendrift.readers.reader_netCDF_CF_generic"
```

GRIB2 reading: cfgrib + eccodes, and herbie-data (with s3fs/boto3) for the AWS bucket noaa-gfs-bdp-pds.

## Data download
Run from the project folder, in this order. Every step skips files that already exist, so an interrupted run can be repeated.
Needs the one-off `.venv/bin/copernicusmarine login` (stored outside the project, never copied here).
```bash
export MPLCONFIGDIR=$PWD/data/mplconfig
.venv/bin/python scripts/drift_lib.py --selftest       # checks geodesy, land distance, contours, bootstrap
.venv/bin/python scripts/s01_gdp_download.py           # GDP drifter_6hour_qc tracks -> data/gdp/ (about 20 MB)
.venv/bin/python scripts/s02_windows_stepzero.py       # step zero counts, blocks, tuning/test split -> results/, work/step0_counts.md
.venv/bin/python scripts/s03_cmems_download.py         # CMEMS currents + Stokes for the chosen blocks -> data/cmems/ (about 6.2 GB, about 75 min)
.venv/bin/python scripts/s04_wind_download.py          # GFS and ECMWF forecast-mode 10 m winds -> data/winds/ (about 350 MB, about 15 min)
```
Step zero (`s02`) runs before any forcing is downloaded, and `s03`/`s04` read `results/blocks.csv` so they only fetch the 96 h blocks that are used.
What was obtained, and what was not, is recorded in `work/data_notes.md`.

Data sources: NOAA AOML Global Drifter Program (ERDDAP `drifter_6hour_qc`); Copernicus Marine GLOBAL_ANALYSISFORECAST_PHY_001_024
(`cmems_mod_glo_phy_anfc_merged-uv_PT1H-i`, `cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i`); NOAA GFS 0.25 degree (AWS `noaa-gfs-bdp-pds`);
ECMWF IFS open data (AWS `ecmwf-forecasts`, CC-BY-4.0).

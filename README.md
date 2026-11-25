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

"""Constants for the NZ drift hindcast spike. All tunable numbers live here."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
WORK = ROOT / "work"
FIGURES = ROOT / "figures"

# --- regions (lon in 0..360) ---
OBS_LON = (160.0, 185.0)      # drifter observation box
OBS_LAT = (-55.0, -25.0)
FORCING_LON = (158.0, 187.0)  # forcing subset box (a margin around the observation box)
FORCING_LAT = (-57.0, -23.0)

# --- GDP ERDDAP ---
GDP_URL = "https://erddap.aoml.noaa.gov/gdp/erddap/tabledap/drifter_6hour_qc.csv"
GDP_VARS = ["ID", "time", "latitude", "longitude", "lon360", "ve", "vn", "err_lat", "err_lon",
            "drogue_lost_date", "start_date", "end_date", "typedeath", "typebuoy"]
GDP_TIME_MIN = "2022-12-01T00:00:00Z"
GDP_TIME_MAX = "2025-06-18T12:00:00Z"   # end of dataset coverage
GDP_CSV = DATA / "gdp" / "gdp_6h_box.csv"

# --- study period (seed times) ---
SEED_START = "2023-01-21T00:00:00"   # ECMWF open-data archive on AWS starts 2023-01-18
SEED_END = "2025-06-14T00:00:00"     # GDP ends 2025-06-18T12Z, minus 96 h block
WINDOW_HOURS = 48
OBS_STEP_H = 6
N_OBS = WINDOW_HOURS // OBS_STEP_H + 1   # 9 observations per window
BLOCK_HOURS = 96
SEED_OFFSETS_H = (0, 48)
BLOCKS_PER_MONTH = 2
COASTAL_EXTRA_BLOCKS = 12
MIN_WINDOWS_PER_BLOCK = 15
RNG_SEED = 42

# --- QC ---
ERR_DEG_MAX = 0.05
MAX_SPEED_MS = 2.5
DROGUE_BUFFER_H = 24
OFFSHORE_KM = 50.0
COASTAL_KM = 10.0
TUNING_FRACTION = 0.4

LEADS_H = (6, 12, 24, 48)
ANECDOTAL_MIN_WINDOWS = 30
ANECDOTAL_MIN_DRIFTERS = 8
OFFSHORE_MIN_TOTAL_WINDOWS = 150

# --- CMEMS ---
CMEMS_UV1H = "cmems_mod_glo_phy_anfc_merged-uv_PT1H-i"
CMEMS_CUR6H = "cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i"
CMEMS_UV1H_VARS = ["uo", "vo", "utide", "vtide", "vsdx", "vsdy"]
CMEMS_DIR = DATA / "cmems"
CMEMS_MAX_GB = 8.0
CMEMS_WORKERS = 3

# --- winds ---
WIND_DIR = DATA / "winds"
WIND_CYCLE_OFFSET_H = 12            # use the cycle at t0 - 12 h
WIND_FHOURS = list(range(12, 61, 3))  # forecast hours 12..60, 17 steps
GFS_BASE = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
ECMWF_BASE = "https://ecmwf-forecasts.s3.amazonaws.com"
ECMWF_PREFIXES = ["ifs/0p25", "0p25", "ifs/0p4-beta", "0p4-beta"]
WIND_WORKERS = 4

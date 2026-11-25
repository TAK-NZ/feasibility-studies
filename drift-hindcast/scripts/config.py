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

# --- simulation (steps s05 to s08) ---
FORCING_DIR = DATA / "forcing"
RUNS_DIR = DATA / "runs"
TIME_STEP_S = 900
OUTPUT_STEP_H = 6
LEAD_IDX = {6: 1, 12: 2, 24: 4, 48: 8}          # index of each lead in the 6-hourly output
WINDAGE_GRID = [0.0, 0.005, 0.01, 0.015, 0.02, 0.03]   # extended by s07 if the optimum sits at an edge
LOST_CURRENTS = ["cur_tide", "cur_notide", "cur_6h_surf"]
DROGUED_CURRENTS = ["cur_tide", "cur_notide", "cur_6h_surf", "cur_6h_15m", "cur_6h_15m_tide"]
WINDS = ["gfs", "ecmwf", "both"]
STOKES = [True, False]
SIM_SEED = 1234                  # base random seed of all ensemble runs
N_PER_DRIFTER_TUNE = 300
N_PER_DRIFTER_TEST = 500
# Level 1 (native OpenDrift) uncertainty grid, D8 in work/plan.md
CUR_UNC_GRID = [0.0, 0.1, 0.3]                 # drift:current_uncertainty (m/s)
WIND_UNC_GRID = [0.0, 1.5, 3.0]                # drift:wind_uncertainty (m/s), lost class only
WINDAGE_REL_SIGMA_GRID = [0.0, 0.3]            # relative sigma of per-element windage, lost class only
CDF_REL_SIGMA_GRID = [0.0, 0.2]                # relative sigma of current_drift_factor, drogued class only
# Level 2 (post hoc persistent velocity offset, m/s per component)
SIGMA_B_GRID = [0.0, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3]
HIT_TARGET_TUNE = 0.90           # selection margin on the tuning set (pass line is 0.80)
N_AREA_SUBSAMPLE = 300           # tuning windows used for the area part of the selection rule
# --- metrics and pass criteria (s09, s10) ---
BOOT_B = 2000                     # cluster bootstrap replicates (clusters = drifters; sensitivity: clusters = 96 h blocks)
BOOT_SEED = 42
PASS_IMPROVEMENT = 0.30           # offshore 24 h median error must be at least this much below the hand-calculation baseline
BORDER_IMPROVEMENT = 0.25         # "about 30 percent": 25 to <30 is BORDERLINE
PASS_HIT95 = 0.80                 # 95 percent contour must hold at least this fraction of observed endpoints
BORDER_HIT95 = 0.75
RELIABILITY_PS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
FIG_DPI = 110                     # figures are at most 2000 x 1400 px

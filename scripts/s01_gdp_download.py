"""Download GDP drifter_6hour_qc tracks for the NZ box from NOAA AOML ERDDAP.

Usage: python scripts/s01_gdp_download.py [--force]
"""
import sys
import urllib.parse
import urllib.request
from io import StringIO

import numpy as np
import pandas as pd

import config as C


def build_url():
    # lon360 box wider than the observation box so edge tracks keep context
    cons = [
        f"lon360%3E={C.FORCING_LON[0]}", f"lon360%3C={C.FORCING_LON[1]}",
        f"latitude%3E={C.FORCING_LAT[0]}", f"latitude%3C={C.FORCING_LAT[1]}",
        f"time%3E={C.GDP_TIME_MIN}", f"time%3C={C.GDP_TIME_MAX}",
    ]
    order = urllib.parse.quote('orderBy("ID,time")', safe="()")
    return f"{C.GDP_URL}?{','.join(C.GDP_VARS)}&{'&'.join(cons)}&{order}"


def main():
    force = "--force" in sys.argv
    if C.GDP_CSV.exists() and not force:
        print(f"exists, skipping: {C.GDP_CSV} (use --force to refresh)")
    else:
        C.GDP_CSV.parent.mkdir(parents=True, exist_ok=True)
        url = build_url()
        print("GET", url)
        with urllib.request.urlopen(url, timeout=600) as r:
            raw = r.read().decode()
        lines = raw.split("\n")
        # line 0 = header, line 1 = units row: drop the units row
        C.GDP_CSV.write_text("\n".join([lines[0]] + lines[2:]))
    df = pd.read_csv(C.GDP_CSV, dtype={"ID": str}, low_memory=False)
    for v in ("ve", "vn"):
        df.loc[df[v] <= -999990, v] = np.nan
    print(f"rows={len(df)} drifters={df.ID.nunique()} "
          f"time {df.time.min()} .. {df.time.max()} size={C.GDP_CSV.stat().st_size / 1e6:.1f} MB")
    print("columns:", list(df.columns))


if __name__ == "__main__":
    main()

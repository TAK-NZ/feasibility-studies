"""Download CMEMS forcing for each 96 h block in results/blocks.csv.

Per block (start bs .. bs+96 h inclusive), box 158-187E 57-23S:
  uv1h_<bs>.nc     merged-uv hourly:  uo vo utide vtide vsdx vsdy  (0.494 m)
  cur6h_s_<bs>.nc  6-hourly currents: uo vo at 0.494 m
  cur6h_15_<bs>.nc 6-hourly currents: uo vo at 15.81 m
Needs a stored Copernicus Marine login (`copernicusmarine login`). Never reads or prints credentials.
Usage: python scripts/s03_cmems_download.py [--force] [--dry-run-only]
"""
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

import config as C

H = pd.Timedelta(hours=1)


def jobs():
    b = pd.read_csv(C.RESULTS / "blocks.csv", parse_dates=["block_start", "block_end"])
    out = []
    for _, r in b.iterrows():
        s, e = r.block_start, r.block_end
        common = dict(minimum_longitude=C.FORCING_LON[0], maximum_longitude=C.FORCING_LON[1],
                      minimum_latitude=C.FORCING_LAT[0], maximum_latitude=C.FORCING_LAT[1],
                      start_datetime=s.strftime("%Y-%m-%dT%H:%M:%S"), end_datetime=e.strftime("%Y-%m-%dT%H:%M:%S"))
        out.append(dict(common, kind="uv1h", block=r.block_id, dataset_id=C.CMEMS_UV1H,
                        variables=C.CMEMS_UV1H_VARS, minimum_depth=0.0, maximum_depth=1.0,
                        output_filename=f"uv1h_{r.block_id}.nc"))
        out.append(dict(common, kind="cur6h_s", block=r.block_id, dataset_id=C.CMEMS_CUR6H,
                        variables=["uo", "vo"], minimum_depth=0.4, maximum_depth=0.6,
                        output_filename=f"cur6h_s_{r.block_id}.nc"))
        out.append(dict(common, kind="cur6h_15", block=r.block_id, dataset_id=C.CMEMS_CUR6H,
                        variables=["uo", "vo"], minimum_depth=15.7, maximum_depth=15.9,
                        output_filename=f"cur6h_15_{r.block_id}.nc"))
    return out


def request(job, dry_run=False, tries=4):
    import copernicusmarine as cm
    kw = {k: v for k, v in job.items() if k not in ("kind", "block")}
    for i in range(tries):
        try:
            return cm.subset(**kw, output_directory=str(C.CMEMS_DIR), netcdf_compression_level=4,
                             skip_existing=True, disable_progress_bar=True, dry_run=dry_run)
        except Exception as e:   # transient network or server errors
            if i == tries - 1:
                raise
            print(f"retry {job['output_filename']} after {type(e).__name__}", flush=True)
            time.sleep(10 * (i + 1))


def run_one(job):
    t = time.time()
    request(job)
    return job["output_filename"], time.time() - t


def manifest(js):
    import xarray as xr
    rows = []
    for j in js:
        f = C.CMEMS_DIR / j["output_filename"]
        if not f.exists():
            rows.append(dict(file=f.name, kind=j["kind"], block=j["block"], present=False))
            continue
        with xr.open_dataset(f) as d:
            dep = d.depth.values if "depth" in d.coords else []
            rows.append(dict(file=f.name, kind=j["kind"], block=j["block"], present=True,
                             size_mb=round(f.stat().st_size / 1e6, 1), t_start=str(d.time.values[0])[:16],
                             t_end=str(d.time.values[-1])[:16], n_time=d.sizes["time"],
                             depth=",".join(f"{x:.3f}" for x in dep), lon_min=float(d.longitude.min()),
                             lon_max=float(d.longitude.max()), lat_min=float(d.latitude.min()),
                             lat_max=float(d.latitude.max()), variables=" ".join(d.data_vars),
                             uo_all_nan=bool(d.uo.isnull().all())))
    m = pd.DataFrame(rows)
    m.to_csv(C.CMEMS_DIR / "manifest.csv", index=False)
    return m


def main():
    force = "--force" in sys.argv
    C.CMEMS_DIR.mkdir(parents=True, exist_ok=True)
    js = jobs()
    todo = [j for j in js if not (C.CMEMS_DIR / j["output_filename"]).exists()]
    print(f"{len(js)} requests, {len(todo)} still to download")
    if todo:
        # dry run estimate on one block of each kind, scaled to the number of blocks
        est = {}
        for k in ("uv1h", "cur6h_s", "cur6h_15"):
            j = next(x for x in todo if x["kind"] == k)
            res = request(j, dry_run=True)
            est[k] = float(res.file_size or 0)
            print(f"dry run {k}: estimated file_size={res.file_size} (MB, uncompressed estimate)")
        n = {k: sum(1 for x in todo if x["kind"] == k) for k in est}
        tot_gb = sum(est[k] * n[k] for k in est) / 1024
        print(f"estimated uncompressed total {tot_gb:.1f} GB (compression level 4 reduces files about 6x)")
        # the server estimate is uncompressed: compare against the cap after the measured ~6x compression
        if tot_gb / 6 > C.CMEMS_MAX_GB and not force:
            print("estimate exceeds the cap, use --force", file=sys.stderr)
            sys.exit(3)
        if "--dry-run-only" in sys.argv:
            return
        t0 = time.time()
        with ProcessPoolExecutor(C.CMEMS_WORKERS) as ex:
            futs = [ex.submit(run_one, j) for j in todo]
            for i, f in enumerate(as_completed(futs), 1):
                name, dt = f.result()
                print(f"[{i}/{len(todo)}] {name} {dt:.0f}s (elapsed {time.time() - t0:.0f}s)", flush=True)
    m = manifest(js)
    print(m.groupby("kind").agg(files=("file", "size"), present=("present", "sum"), size_mb=("size_mb", "sum")))
    print("total size GB:", round(m.size_mb.sum() / 1000, 2))
    print("depth values:", m.groupby("kind").depth.unique().to_dict())
    print("all-NaN uo files:", int(m.uo_all_nan.sum()))


if __name__ == "__main__":
    main()

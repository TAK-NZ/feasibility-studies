"""Shared helpers: geodesy, distance to land, highest-density contours, cluster bootstrap.

Run `python scripts/drift_lib.py --selftest` to check them.
"""
import sys

import numpy as np
import pandas as pd
from pyproj import Geod

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import config as C  # noqa: E402

GEOD = Geod(ellps="WGS84")
R_EARTH_KM = 6371.0088


def geod_km(lon1, lat1, lon2, lat2):
    """Geodesic distance in km (WGS84), vectorised."""
    _, _, d = GEOD.inv(np.asarray(lon1, float), np.asarray(lat1, float),
                       np.asarray(lon2, float), np.asarray(lat2, float))
    return np.asarray(d) / 1000.0


def lon180(lon):
    return (np.asarray(lon, float) + 180.0) % 360.0 - 180.0


def local_km(lon, lat, lon0, lat0):
    """Local tangent-plane km (east, north) of points relative to (lon0, lat0)."""
    dlon = lon180(np.asarray(lon, float) - lon0)
    x = np.deg2rad(dlon) * R_EARTH_KM * np.cos(np.deg2rad(lat0))
    y = np.deg2rad(np.asarray(lat, float) - lat0) * R_EARTH_KM
    return x, y


# ---------------------------------------------------------------- land distance
class LandDistance:
    """Distance to land (km) from the roaring landmask, rasterised at `res` degrees."""

    def __init__(self, res=0.02, lon_range=(158.0, 187.0), lat_range=(-57.0, -23.0)):
        from roaring_landmask import RoaringLandmask
        from scipy.spatial import cKDTree
        self.res = res
        lons = np.arange(lon_range[0], lon_range[1] + res / 2, res)
        lats = np.arange(lat_range[0], lat_range[1] + res / 2, res)
        LON, LAT = np.meshgrid(lons, lats)
        mask = RoaringLandmask.new()
        land = mask.contains_many_par(lon180(LON).ravel(), LAT.ravel()).reshape(LON.shape)
        self.land_fraction = float(land.mean())
        # boundary = land cells with at least one ocean 4-neighbour
        pad = np.pad(land, 1, constant_values=False)
        ocean_nb = ~pad[:-2, 1:-1] | ~pad[2:, 1:-1] | ~pad[1:-1, :-2] | ~pad[1:-1, 2:]
        boundary = land & ocean_nb
        self.n_boundary = int(boundary.sum())
        self.tree = cKDTree(self._xyz(LON[boundary], LAT[boundary]))
        self._mask = mask

    @staticmethod
    def _xyz(lon, lat):
        lo, la = np.deg2rad(lon), np.deg2rad(lat)
        return R_EARTH_KM * np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])

    def km(self, lon, lat):
        """Chord distance (about equal to the surface distance for tens of km). 0 on land."""
        lon = np.asarray(lon, float)
        lat = np.asarray(lat, float)
        d, _ = self.tree.query(self._xyz(lon, lat))
        on_land = self._mask.contains_many_par(lon180(lon).ravel(), lat.ravel())
        d = np.where(on_land.reshape(d.shape), 0.0, d)
        return d


# ---------------------------------------------------------------- contours
def hdr_level(sample_xy, p, bw_min_km=3.0, grid_step_km=1.0):
    """Gaussian KDE of a 2-D sample (km). Returns (kde, level_p, area_km2 at level, pad).

    The p-contour is the density level that encloses fraction p of the sample's own KDE
    values (highest-density region). Bandwidth = max(Scott factor bandwidth, bw_min_km).
    """
    from scipy.stats import gaussian_kde
    xy = np.asarray(sample_xy, float)
    n = xy.shape[1]
    kde = gaussian_kde(xy)
    sd = xy.std(axis=1, ddof=1)
    scott = n ** (-1.0 / 6.0)
    bw = np.maximum(scott * sd, bw_min_km)
    kde.set_bandwidth(bw_method=float(np.mean(bw / np.maximum(sd, 1e-9))))
    dens = kde(xy)
    level = np.quantile(dens, 1.0 - p)
    return kde, level


def hdr_area_km2(kde, level, xy, grid_step_km=1.0, margin_km=30.0):
    xmin, ymin = xy.min(axis=1) - margin_km
    xmax, ymax = xy.max(axis=1) + margin_km
    gx, gy = np.meshgrid(np.arange(xmin, xmax, grid_step_km), np.arange(ymin, ymax, grid_step_km))
    dens = kde(np.vstack([gx.ravel(), gy.ravel()]))
    return float((dens >= level).sum() * grid_step_km ** 2)


# ---------------------------------------------------------------- bootstrap
def cluster_bootstrap(values, clusters, stat=np.median, B=2000, seed=C.RNG_SEED, alpha=0.05):
    """Resample clusters with replacement, return (point, lo, hi)."""
    values = np.asarray(values, float)
    clusters = np.asarray(clusters)
    uniq, inv = np.unique(clusters, return_inverse=True)
    groups = [values[inv == i] for i in range(len(uniq))]
    rng = np.random.default_rng(seed)
    out = np.empty(B)
    for b in range(B):
        pick = rng.integers(0, len(groups), len(groups))
        out[b] = stat(np.concatenate([groups[i] for i in pick]))
    return float(stat(values)), float(np.quantile(out, alpha / 2)), float(np.quantile(out, 1 - alpha / 2))


def _selftest():
    ok = True

    def check(name, cond, info=""):
        nonlocal ok
        print(("OK   " if cond else "FAIL ") + name, info)
        ok &= bool(cond)

    # (a) geodesic: one degree of latitude on the meridian is about 111.2 km at -40
    d = float(geod_km(175.0, -40.0, 175.0, -39.0))
    check("geod one degree latitude", abs(d - 111.0) < 1.0, f"{d:.2f} km")
    # (b) HDR contours on a Gaussian
    rng = np.random.default_rng(1)
    s = rng.normal(0, 20, size=(2, 5000))
    fresh = rng.normal(0, 20, size=(2, 20000))
    for p, tol in ((0.95, 0.02), (0.5, 0.03)):
        kde, lev = hdr_level(s, p)
        hit = float((kde(fresh) >= lev).mean())
        check(f"HDR {int(p * 100)} percent hit", abs(hit - p) < tol, f"{hit:.3f}")
    # (c) cluster bootstrap wider than naive on clustered data
    rng = np.random.default_rng(2)
    cl = np.repeat(np.arange(30), 40)
    vals = rng.normal(0, 1, 30)[cl] * 3 + rng.normal(0, 1, cl.size)
    _, lo, hi = cluster_bootstrap(vals, cl, B=500)
    naive = [np.median(rng.choice(vals, vals.size)) for _ in range(500)]
    nlo, nhi = np.quantile(naive, [0.025, 0.975])
    check("cluster CI wider than naive CI", (hi - lo) > (nhi - nlo), f"{hi - lo:.3f} vs {nhi - nlo:.3f}")
    # (d) land distance: Wellington airport area near land, mid-Tasman far from land
    ld = LandDistance()
    dd = ld.km([174.8, 170.0], [-41.3, -42.0])
    check("land distance near coast / open sea", dd[0] < 10 and dd[1] > 100, f"{dd.round(1)} km")
    return ok


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(0 if _selftest() else 1)

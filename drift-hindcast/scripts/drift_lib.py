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
BW_FLOOR_KM = 3.0


def _kernel_cov(sample_xy, bw_floor_km=BW_FLOOR_KM):
    """Gaussian kernel covariance: Scott's rule n^(-1/3) x sample covariance (2-D) + floor^2 x identity."""
    n = sample_xy.shape[1]
    S = np.cov(sample_xy) if n > 2 else np.zeros((2, 2))
    return n ** (-1.0 / 3.0) * S + bw_floor_km ** 2 * np.eye(2)


def _kern(dx, dy, Hinv):
    return np.exp(-0.5 * (Hinv[0, 0] * dx * dx + 2 * Hinv[0, 1] * dx * dy + Hinv[1, 1] * dy * dy))


def hdr_levels(sample_xy, ps=(0.5, 0.95), bw_floor_km=BW_FLOOR_KM):
    """Highest-density-region levels of a 2-D particle cloud (km).

    Gaussian KDE with covariance H = n^(-1/3) Cov + floor^2 I. The p-level is the density below which
    a fraction (1 - p) of the cloud's own leave-one-out KDE values lie, so the p-contour {density >= level}
    holds about p of the particles. Returns (H, {p: level}) with densities normalised by the number of
    kernels (the 1/(2 pi sqrt(det H)) factor is dropped, and the same is done in kde_density).
    """
    sample_xy = np.asarray(sample_xy, float)
    n = sample_xy.shape[1]
    H = _kernel_cov(sample_xy, bw_floor_km)
    Hinv = np.linalg.inv(H)
    dx = sample_xy[0][:, None] - sample_xy[0][None, :]
    dy = sample_xy[1][:, None] - sample_xy[1][None, :]
    K = _kern(dx, dy, Hinv)
    loo = (K.sum(1) - 1.0) / (n - 1)
    return H, {p: float(np.quantile(loo, 1.0 - p)) for p in ps}


def kde_density(sample_xy, pts_xy, H):
    """KDE value (same normalisation as hdr_levels) of the cloud at points pts_xy (2, m)."""
    sample_xy = np.asarray(sample_xy, float)
    pts_xy = np.asarray(pts_xy, float).reshape(2, -1)
    Hinv = np.linalg.inv(H)
    dx = pts_xy[0][:, None] - sample_xy[0][None, :]
    dy = pts_xy[1][:, None] - sample_xy[1][None, :]
    return _kern(dx, dy, Hinv).mean(1)


def hdr_hits(sample_xy, obs_xy, ps=(0.5, 0.95)):
    """Boolean hit per p: is the observed point inside the p highest-density contour of the cloud?"""
    H, lev = hdr_levels(sample_xy, ps)
    d = float(kde_density(sample_xy, np.asarray(obs_xy, float).reshape(2, 1), H)[0])
    return {p: bool(d >= lev[p]) for p in ps}


def hdr_areas(sample_xy, ps=(0.5, 0.95), max_grid=150, pad_sd=3.0):
    """Areas (km^2) of the p contours on a regular grid of at most max_grid x max_grid cells (step >= 1 km)."""
    sample_xy = np.asarray(sample_xy, float)
    H, lev = hdr_levels(sample_xy, ps)
    pad = pad_sd * np.sqrt(np.diag(H)) + 0.0
    lo = sample_xy.min(1) - pad
    hi = sample_xy.max(1) + pad
    step = max(1.0, float(np.max(hi - lo)) / max_grid)
    gx, gy = np.meshgrid(np.arange(lo[0], hi[0], step), np.arange(lo[1], hi[1], step))
    dens = kde_density(sample_xy, np.vstack([gx.ravel(), gy.ravel()]), H)
    return {p: float((dens >= lev[p]).sum() * step * step) for p in ps}


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
    # (b) HDR contours on a Gaussian: fresh draws must fall inside with the nominal fraction, for n = 300 and 5000
    rng = np.random.default_rng(1)
    for n in (300, 5000):
        hits = {0.5: [], 0.95: []}
        for _ in range(200 if n == 300 else 10):
            s = rng.normal(0, 20, size=(2, n))
            H, lev = hdr_levels(s)
            fresh = rng.normal(0, 20, size=(2, 400))
            dd = kde_density(s, fresh, H)
            for p in hits:
                hits[p].append((dd >= lev[p]).mean())
        for p, tol in ((0.95, 0.02), (0.5, 0.03)):
            h = float(np.mean(hits[p]))
            check(f"HDR {int(p * 100)} percent hit, n={n}", abs(h - p) < tol, f"{h:.3f}")
    a = hdr_areas(rng.normal(0, 20, size=(2, 2000)))
    check("HDR 95 percent area close to the Gaussian value", abs(a[0.95] / (np.pi * 5.991 * 400) - 1) < 0.1, f"{a[0.95]:.0f} km2")
    check("HDR of a point cloud does not fail", hdr_hits(np.zeros((2, 50)), [0.0, 0.0])[0.95])
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

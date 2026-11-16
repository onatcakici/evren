#!/usr/bin/env python3
"""
Evren spectrum: what a moving solar system body is made of, measured from
raw SPHEREx Level 2 images.
Code written with help from an AI assistant.

SPHEREx looks at the sky through a linear variable filter, so every pixel
sees its own wavelength. When a planet or asteroid lands on a pixel, that
exposure measures its brightness at that one wavelength. Exposures taken on
different days, on different detectors, add up to a spectrum.

For every exposure that contains the target this script:
  1. reads the image header and WCS with an HTTP Range request,
  2. places the target where JPL Horizons says it was at mid exposure,
  3. reads only the rows of pixels around the target's track,
  4. measures the light in a small aperture, then subtracts the light at the
     same spot in exposures from other days, when the target had moved on,
     so any star under it cancels out (the blink comparator idea again),
  5. looks up the wavelength and solid angle of that pixel in the SPHEREx
     calibration maps.

Frames come from saved archive answers (--rows, the Level 2 rows that
evren_build.py caches) and, optionally, from a scan of the archive's file
tree (--scan), which finds exposures where the target fell on other
detectors. The scan reads file headers directly, so it works even when the
archive's search service is down.

Example:
  python evren_spectrum.py --target 999 --name Pluto \\
      --rows cache/pluto2025_D1_89.npz --rows cache/pluto2026_D1_95.npz \\
      --bands D1,D4 --out docs/data/spectra/pluto
"""
import argparse
import concurrent.futures as cf
import hashlib
import io
import json
import math
import re
import sys
import time
import urllib.parse
import warnings
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from evren_build import get, read_rows, _read_header, _data_size  # noqa: E402

warnings.simplefilter("ignore", FITSFixedWarning)

S3 = "https://nasa-irsa-spherex.s3.amazonaws.com/"
IBE = "https://irsa.ipac.caltech.edu/ibe/data/spherex/"
CAL = {
    "wave": "dr1/spectral_wcs/cal-swcs-v5-2026-191/{d}/spectral_wcs_D{d}_spx_cal-swcs-v5-2026-191.fits",
    "omega": "dr1/solid_angle_pixel_map/cal-sapm-v3-2026-191/{d}/solid_angle_pixel_map_D{d}_spx_cal-sapm-v3-2026-191.fits",
}
BAND_RANGE = {1: (0.70, 1.15), 2: (1.05, 1.70), 3: (1.58, 2.50),
              4: (2.35, 3.90), 5: (3.75, 4.50), 6: (4.35, 5.10)}
CACHE = HERE / "cache"
NPIX = 2040
ARCSEC2_SR = (math.pi / 180 / 3600) ** 2
R_AP, R_IN, R_OUT = 2.0, 5.0, 9.0      # aperture and background annulus, pixels
AWAY_ARCSEC = 75.0                       # target counts as gone beyond this


# ---------------------------------------------------------------- archive

def s3_url(url):
    return S3 + url[len(IBE):] if url.startswith(IBE) else url


def band_url(url, band):
    d = int(str(band).strip().upper().lstrip("D"))
    return re.sub(r"/(\d)/level2_(.+?)D\d_spx",
                  lambda m: f"/{d}/level2_{m.group(2)}D{d}_spx", url)


def det_of(url):
    m = re.search(r"D(\d)_spx", url)
    return int(m.group(1)) if m else 0


def load_rows(path):
    z = np.load(path, allow_pickle=True)
    return json.loads(str(z["rows"]))


def header_pair(url, want=48 * 1024):
    """Primary header and IMAGE header of a SPHEREx L2 file, cached on disk."""
    fn = CACHE / "hdr" / (hashlib.md5(url.encode()).hexdigest() + ".json")
    if fn.exists():
        try:
            d = json.loads(fn.read_text())
            return fits.Header.fromstring(d["p"]), fits.Header.fromstring(d["h"]), int(d["s"])
        except Exception:                                        # noqa: BLE001
            pass
    prim, h, dstart = _header_pair(url, want)
    fn.parent.mkdir(parents=True, exist_ok=True)
    fn.write_text(json.dumps({"p": prim.tostring(), "h": h.tostring(), "s": dstart}))
    return prim, h, dstart


def _header_pair(url, want=48 * 1024):
    """Primary header and IMAGE header of a SPHEREx L2 file, by Range reads.

    Only header bytes are downloaded (see evren_build._read_header), so a file
    without an IMAGE extension stops at its end instead of looping.
    """
    buf, start, pos, prim = b"", 0, 0, None
    for _ in range(6):
        res = _read_header(url, pos, buf, start, chunk=want)
        if res is None:
            break
        hs, dstart, buf, start = res
        h = fits.Header.fromstring(hs)
        if prim is None:
            prim = h
        if str(h.get("EXTNAME", "")).strip().upper() == "IMAGE":
            return prim, h, dstart
        pos = dstart + _data_size(h)
    raise RuntimeError("IMAGE extension not found in " + url)


TIME_KEYS = ("MJD-AVG", "MJD-MID", "MJDMID", "MJD-OBS", "MJD-BEG", "MJD_OBS")
DATE_KEYS = ("DATE-AVG", "DATE-MID", "DATE-OBS", "DATE-BEG")


def mid_mjd(prim, h):
    for hh in (h, prim):
        for k in TIME_KEYS:
            if k in hh:
                try:
                    v = float(hh[k])
                    if k in ("MJD-OBS", "MJD-BEG", "MJD_OBS"):
                        exp = float(hh.get("EXPTIME", hh.get("TEXP", 0)) or 0)
                        v += 0.5 * exp / 86400.0
                    return v
                except (TypeError, ValueError):
                    pass
    for hh in (h, prim):
        for k in DATE_KEYS:
            if k in hh:
                s = str(hh[k]).strip().replace("Z", "")
                for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S"):
                    try:
                        dt = datetime.strptime(s, fmt)
                        return (dt - datetime(1858, 11, 17)).total_seconds() / 86400.0
                    except ValueError:
                        pass
    return None


def s3_list(prefix, delimiter=None):
    """Keys (or common prefixes when delimiter is set) under prefix in the bucket."""
    out, token = [], None
    ns = "{http://s3.amazonaws.com/doc/2006-03-01/}"
    while True:
        q = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if delimiter:
            q["delimiter"] = delimiter
        if token:
            q["continuation-token"] = token
        root = ET.fromstring(get(S3 + "?" + urllib.parse.urlencode(q), timeout=60))
        if delimiter:
            out += [el.text for el in root.iter(ns + "Prefix") if el.text != prefix]
        else:
            out += [el.text for el in root.iter(ns + "Key")]
        trunc = root.find(ns + "IsTruncated")
        if trunc is None or trunc.text != "true":
            return out
        token = root.find(ns + "NextContinuationToken").text


# --------------------------------------------------------------- ephemeris

def mjd_to_iso(mjd):
    return (datetime(1858, 11, 17) + timedelta(days=mjd)).strftime("%Y-%m-%d %H:%M")


def horizons(target, t0, t1, step="20m"):
    """Geocentric astrometric RA, Dec (deg) of target between MJD t0 and t1."""
    p = {
        "format": "text", "COMMAND": f"'{target}'", "OBJ_DATA": "'NO'",
        "MAKE_EPHEM": "'YES'", "EPHEM_TYPE": "'OBSERVER'", "CENTER": "'500@399'",
        "START_TIME": f"'{mjd_to_iso(t0)}'", "STOP_TIME": f"'{mjd_to_iso(t1)}'",
        "STEP_SIZE": f"'{step}'", "QUANTITIES": "'1'", "ANG_FORMAT": "'DEG'",
        "CSV_FORMAT": "'YES'",
    }
    txt = get("https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(p),
              timeout=180).decode("utf-8", "replace")
    if "$$SOE" not in txt:
        raise RuntimeError("Horizons gave no ephemeris:\n" + txt[:1500])
    t, ra, de = [], [], []
    for line in txt.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines():
        parts = [s.strip() for s in line.split(",")]
        nums = []
        for s in parts[1:]:
            try:
                nums.append(float(s))
            except ValueError:
                pass
        if len(nums) < 2:
            continue
        dt = datetime.strptime(parts[0], "%Y-%b-%d %H:%M")
        t.append((dt - datetime(1858, 11, 17)).total_seconds() / 86400.0)
        ra.append(nums[0])
        de.append(nums[1])
    return np.array(t), np.unwrap(np.radians(ra)), np.radians(de)


class Track:
    """Interpolated target positions, built from Horizons in time clusters."""

    def __init__(self, target, times):
        times = np.sort(np.asarray(times, float))
        groups, cur = [], [times[0]]
        for x in times[1:]:
            if x - cur[-1] > 2.0:
                groups.append(cur)
                cur = [x]
            else:
                cur.append(x)
        groups.append(cur)
        self.parts = []
        for g in groups:
            t, ra, de = horizons(target, g[0] - 0.05, g[-1] + 0.05)
            self.parts.append((t, ra, de))

    def __call__(self, mjd):
        mjd = np.atleast_1d(np.asarray(mjd, float))
        ra = np.full(mjd.shape, np.nan)
        de = np.full(mjd.shape, np.nan)
        for t, r, d in self.parts:
            sel = (mjd >= t[0]) & (mjd <= t[-1])
            ra[sel] = np.degrees(np.interp(mjd[sel], t, r)) % 360.0
            de[sel] = np.degrees(np.interp(mjd[sel], t, d))
        return ra, de


def sep_arcsec(ra1, de1, ra2, de2):
    r1, d1, r2, d2 = map(np.radians, (ra1, de1, ra2, de2))
    s = np.sin((d2 - d1) / 2) ** 2 + np.cos(d1) * np.cos(d2) * np.sin((r2 - r1) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(s, 0, 1)))) * 3600.0


# ------------------------------------------------------------ calibration

def cal_map(kind, d):
    out = CACHE / "cal" / f"{kind}_D{d}.npy"
    if out.exists():
        return np.load(out, mmap_mode="r")
    out.parent.mkdir(parents=True, exist_ok=True)
    raw = get(S3 + CAL[kind].format(d=d), timeout=600)
    planes = []
    with fits.open(io.BytesIO(raw)) as hl:
        for hdu in hl:
            if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                cols = {c.upper(): c for c in hdu.columns.names}
                if {"X", "Y", "VALUES"} <= set(cols):
                    row = hdu.data[0]
                    X = np.asarray(row[cols["X"]], float).ravel()
                    Y = np.asarray(row[cols["Y"]], float).ravel()
                    V = np.asarray(row[cols["VALUES"]], float)
                    if X.min() >= 1 and X.max() > NPIX - 2:      # FITS pixels are 1 based
                        X = X - 1
                    if Y.min() >= 1 and Y.max() > NPIX - 2:
                        Y = Y - 1
                    print(f"  {kind} table D{d}: X {X[0]:.0f}..{X[-1]:.0f}, Y {Y[0]:.0f}..{Y[-1]:.0f}, "
                          f"values {V.shape}", flush=True)
                    # FITS -TAB order: VALUES[y][x][m], m = 0 wavelength, 1 bandwidth
                    for m in range(V.shape[-1]):
                        G = V[..., m]
                        xs = np.arange(NPIX, dtype=float)
                        tmp = np.array([np.interp(xs, X, G[j, :]) for j in range(len(Y))])
                        full = np.empty((NPIX, NPIX), np.float32)
                        for i in range(NPIX):
                            full[:, i] = np.interp(xs, Y, tmp[:, i])
                        planes.append((f"{hdu.name}.VALUES[{m}]", full))
                continue
            a = getattr(hdu, "data", None)
            if a is None or not hasattr(a, "ndim"):
                continue
            try:
                a = np.asarray(a, dtype=np.float32)
            except (TypeError, ValueError):
                continue
            if a.ndim == 2 and a.shape == (NPIX, NPIX):
                planes.append((hdu.name, a))
            elif a.ndim == 3 and a.shape[1:] == (NPIX, NPIX):
                planes += [(f"{hdu.name}[{i}]", a[i]) for i in range(a.shape[0])]
    if not planes:
        raise RuntimeError(f"no {NPIX}x{NPIX} plane in {kind} map for D{d}")
    if kind == "wave":
        lo, hi = BAND_RANGE[d]
        best = None
        for name, a in planes:
            m = float(np.nanmedian(a))
            for scale in (1.0, 1e-3, 1e6):          # microns, nm, meters
                if lo <= m * scale <= hi:
                    score = (0 if "WAVE" in name.upper() else 1)
                    if best is None or score < best[0]:
                        best = (score, name, a * scale)
        if best is None:
            raise RuntimeError(f"no wavelength plane fits D{d}: " +
                               ", ".join(f"{n}={np.nanmedian(a):.4g}" for n, a in planes))
        _, name, data = best
    else:
        name, data = planes[0]
        if float(np.nanmedian(data)) > 1e-3:          # arcsec^2 per pixel
            data = data * ARCSEC2_SR
    print(f"  {kind} map D{d}: plane {name}, median {np.nanmedian(data):.5g}", flush=True)
    np.save(out, data.astype(np.float32))
    return np.load(out, mmap_mode="r")


# ------------------------------------------------------------- photometry

_SUB = (np.arange(5) + 0.5) / 5.0 - 0.5


def aperture(img, x, y, r_ap=R_AP, r_in=R_IN, r_out=R_OUT):
    """Background subtracted sum (image units) in a circle at (x, y).

    Returns (sum, error, ok). x and y are 0 based pixel coordinates in img.
    """
    H, W = img.shape
    R = int(math.ceil(r_out)) + 1
    ix, iy = int(round(x)), int(round(y))
    x0, x1, y0, y1 = ix - R, ix + R + 1, iy - R, iy + R + 1
    if x0 < 0 or y0 < 0 or x1 > W or y1 > H:
        return np.nan, np.nan, False
    sub = img[y0:y1, x0:x1].astype(np.float64)
    yy, xx = np.mgrid[y0:y1, x0:x1].astype(np.float64)
    d = np.hypot(xx - x, yy - y)
    ann = (d >= r_in) & (d <= r_out) & np.isfinite(sub)
    if ann.sum() < 20:
        return np.nan, np.nan, False
    vals = sub[ann]
    for _ in range(4):
        med = np.median(vals)
        mad = 1.4826 * np.median(np.abs(vals - med)) + 1e-30
        keep = np.abs(vals - med) < 3.0 * mad
        if keep.all():
            break
        vals = vals[keep]
    bg, sig = float(np.median(vals)), float(1.4826 * np.median(np.abs(vals - np.median(vals))))
    near = d <= r_ap + 1.0
    if not np.isfinite(sub[near]).all():
        return np.nan, np.nan, False
    w = np.zeros_like(sub)
    for oy in _SUB:
        for ox in _SUB:
            w += (np.hypot(xx + ox - x, yy + oy - y) <= r_ap)
    w /= len(_SUB) ** 2
    area = float(w.sum())
    total = float(np.nansum(w * (sub - bg)))
    err = sig * math.sqrt(area * (1.0 + area / len(vals)))
    return total, err, True


def planck_nu(lam_um, T=5772.0):
    """Planck B_nu, arbitrary units. Stands in for sunlight."""
    x = 14387.77 / (lam_um * T)
    return (1.0 / lam_um ** 3) / np.expm1(x)


def runs(rows_needed):
    """Group sorted row numbers into contiguous (r0, r1) runs."""
    out = []
    for r in sorted(rows_needed):
        if out and r == out[-1][1] + 1:
            out[-1][1] = r
        else:
            out.append([r, r])
    return out


# ------------------------------------------------------------------ frames

class Frame:
    def __init__(self, url, obs_id, t=None):
        self.url, self.obs_id, self.t = url, obs_id, t
        self.det = det_of(url)
        self.prim = self.h = self.dstart = self.wcs = None

    def load(self):
        self.prim, self.h, self.dstart = header_pair(self.url)
        self.wcs = WCS(self.h)
        c = self.wcs.all_pix2world([[(NPIX - 1) / 2, (NPIX - 1) / 2]], 0)[0]
        self.center = (float(c[0]), float(c[1]))
        if self.t is None:
            self.t = mid_mjd(self.prim, self.h)
        return self

    def pix(self, ra, dec):
        ra, dec = np.atleast_1d(np.asarray(ra, float)), np.atleast_1d(np.asarray(dec, float))
        x = np.full(ra.shape, np.nan)
        y = np.full(ra.shape, np.nan)
        near = np.isfinite(ra) & (sep_arcsec(ra, dec, *self.center) < 4.0 * 3600)
        if near.any():
            xx, yy = self.wcs.all_world2pix(ra[near], dec[near], 0, quiet=True)
            x[near], y[near] = xx, yy
        return x, y


def frames_from_rows(paths, bands):
    out, seen = [], set()
    for p in paths:
        for r in load_rows(p):
            t = 0.5 * (float(r["t_min"]) + float(r["t_max"]))
            for b in bands:
                url = s3_url(band_url(r["access_url"], b))
                if url in seen:
                    continue
                seen.add(url)
                out.append(Frame(url, r["obs_id"], t))
    return out


def scan_frames(weeks, dets, track_fn, collection="dr1", workers=24, step_one_first=True):
    """Find exposures whose detector `d` contains the target, by reading headers.

    Reads every step 1 header of each exposure sequence first; only sequences
    whose step 1 frame lands within ~500 pixels of the target get all of their
    steps read. Returns a list of loaded Frame objects that contain the target.
    """
    found = []
    folders = []
    for wk in weeks:
        if "_" in wk:
            folders.append(wk)
        else:
            top = s3_list(f"{collection}/level2/", delimiter="/")
            folders += [p.rstrip("/").split("/")[-1] for p in top
                        if p.rstrip("/").split("/")[-1].startswith(wk)]
    for wk in folders:
        base = f"{collection}/level2/{wk}/"
        try:
            vers = s3_list(base, delimiter="/")
        except Exception as e:                                   # noqa: BLE001
            print(f"  scan {wk}: cannot list ({e})", flush=True)
            continue
        for ver in vers:
            for d in dets:
                keys = [k for k in s3_list(f"{ver}{d}/") if k.endswith(".fits")]
                if not keys:
                    continue
                seqs = {}
                for k in keys:
                    m = re.search(r"level2_(\w+?_\w+?_\d+)_(\d+)D\d_spx", k)
                    if m:
                        seqs.setdefault(m.group(1), []).append((int(m.group(2)), k))
                firsts = [sorted(v)[0][1] for v in seqs.values()] if step_one_first else keys
                print(f"  scan {ver}{d}/: {len(keys)} files, {len(seqs)} sequences", flush=True)

                def probe(key):
                    try:
                        f = Frame(S3 + key, re.search(r"level2_(.+?)D\d_spx", key).group(1)).load()
                    except Exception as e:                       # noqa: BLE001
                        if not getattr(probe, "warned", False):
                            probe.warned = True
                            print(f"    header read failed for {key}: {e}", flush=True)
                        return None
                    if f.t is None:
                        if not getattr(probe, "warned_t", False):
                            probe.warned_t = True
                            keys = [k for k in list(f.prim.keys()) + list(f.h.keys())
                                    if any(s in k for s in ("MJD", "DATE", "TIME", "EXP"))]
                            print(f"    no time found in {key}; keys: {keys}", flush=True)
                        return None
                    ra, de = track_fn(f.t)
                    if not np.isfinite(ra[0]):
                        return None
                    x, y = f.pix(ra, de)
                    return f, float(x[0]), float(y[0])

                hits = []
                with cf.ThreadPoolExecutor(workers) as ex:
                    for res in ex.map(probe, firsts):
                        if res is None:
                            continue
                        f, x, y = res
                        if -520 <= x <= NPIX + 520 and -520 <= y <= NPIX + 520:
                            hits.append(f)
                if not step_one_first:
                    found += [f for f in hits if 12 <= f.pix(*track_fn(f.t))[0][0] <= NPIX - 13]
                    continue
                todo = []
                for f in hits:
                    seq = re.search(r"level2_(\w+?_\w+?_\d+)_\d+D", f.url).group(1)
                    todo += [k for _, k in sorted(seqs[seq])]
                with cf.ThreadPoolExecutor(workers) as ex:
                    for res in ex.map(probe, todo):
                        if res is None:
                            continue
                        f, x, y = res
                        if 12 <= x <= NPIX - 13 and 12 <= y <= NPIX - 13:
                            found.append(f)
                print(f"    {len(hits)} sequences near the target, {len(found)} frames so far", flush=True)
    return found


# ----------------------------------------------------------------- measure

def measure(frames, track, workers=12):
    """Forced photometry with same spot, other day subtraction.

    Every frame is measured at the target position of every other frame of
    the same detector, so frames taken after the target moved on give the
    static sky under it.
    """
    t = np.array([f.t for f in frames])
    ra, de = track(t)
    ok = np.isfinite(ra)
    frames = [f for f, k in zip(frames, ok) if k]
    t, ra, de = t[ok], ra[ok], de[ok]
    n = len(frames)
    for d in sorted({f.det for f in frames}):
        cal_map("wave", d)
        cal_map("omega", d)
    A = np.full((n, n), np.nan)
    E = np.full((n, n), np.nan)
    L = np.full((n, n), np.nan)
    own = [None] * n

    by_det = {}
    for i, f in enumerate(frames):
        by_det.setdefault(f.det, []).append(i)

    def work(g):
        f = frames[g]
        members = by_det[f.det]
        wave = cal_map("wave", f.det)
        xs, ys = f.pix(ra[members], de[members])
        inside = (xs >= 12) & (xs <= NPIX - 13) & (ys >= 12) & (ys <= NPIX - 13)
        need = set()
        for x, y, k in zip(xs, ys, inside):
            if k:
                need.update(range(max(0, int(y) - 11), min(NPIX, int(y) + 13)))
        if not need:
            return g, [], None
        strips = {}
        for r0, r1 in runs(need):
            strips[(r0, r1)] = read_rows(f.url, f.h, f.dstart, r0, r1)
        res = []
        for j, x, y, k in zip(members, xs, ys, inside):
            if not k:
                continue
            for (r0, r1), img in strips.items():
                if r0 <= y - 10 and y + 11 <= r1:
                    s, e, good = aperture(img, x, y - r0)
                    if good:
                        lam = float(wave[int(round(y)), int(round(x))])
                        res.append((j, s, e, lam))
                    break
        stamp = None
        j0 = members.index(g) if g in members else None
        if j0 is not None and inside[j0]:
            y0 = ys[j0]
            for (r0, r1), img in strips.items():
                if r0 <= y0 - 10 and y0 + 11 <= r1:
                    cy, cx = int(round(y0)) - r0, int(round(xs[j0]))
                    stamp = np.array(img[cy - 10:cy + 11, cx - 10:cx + 11])
                    break
        return g, res, stamp

    with cf.ThreadPoolExecutor(workers) as ex:
        futs = [ex.submit(work, g) for g in range(n)]
        for k, fu in enumerate(cf.as_completed(futs), 1):
            try:
                g, res, stamp = fu.result()
            except Exception as e:                               # noqa: BLE001
                print(f"  frame failed: {e}", flush=True)
                continue
            for j, s, e, lam in res:
                A[g, j], E[g, j], L[g, j] = s, e, lam
            own[g] = stamp
            if k % 20 == 0 or k == n:
                print(f"  measured {k}/{n} frames", flush=True)

    pts = []
    for f_i, f in enumerate(frames):
        if not np.isfinite(A[f_i, f_i]):
            continue
        same = np.array([frames[g].det == f.det for g in range(n)])
        away = sep_arcsec(ra, de, ra[f_i], de[f_i]) > AWAY_ARCSEC
        cand = np.where(same & away & np.isfinite(A[:, f_i]))[0]
        if len(cand) == 0:
            continue
        lam_f = L[f_i, f_i]
        cand = cand[np.argsort(np.abs(L[cand, f_i] - lam_f))][:7]
        ref = A[cand, f_i]
        sky = float(np.median(ref))
        sky_err = 1.2533 * float(np.std(ref)) / math.sqrt(len(ref)) if len(ref) > 1 else float(E[cand[0], f_i])
        xf, yf = f.pix(ra[f_i], de[f_i])
        xi = int(round(float(np.clip(xf[0], 0, NPIX - 1))))
        yi = int(round(float(np.clip(yf[0], 0, NPIX - 1))))
        omega = float(cal_map("omega", f.det)[yi, xi])
        unit = 1e12 if "MJY" in str(f.h.get("BUNIT", "MJy/sr")).upper() else 1.0
        flux = (A[f_i, f_i] - sky) * omega * unit
        err = math.hypot(E[f_i, f_i], sky_err) * omega * unit
        pts.append(dict(w=round(lam_f, 5), f=flux, e=err, sky=sky * omega * unit,
                        det=f.det, obs=f.obs_id, t=round(f.t, 5), nref=int(len(cand)),
                        ra=round(float(ra[f_i]), 6), dec=round(float(de[f_i]), 6)))
    stamps = {frames[i].obs_id + f"_D{frames[i].det}": own[i] for i in range(n) if own[i] is not None}
    return pts, stamps


# ------------------------------------------------------------------ output

def reflectance(pts, norm_lo=0.92, norm_hi=1.08):
    for p in pts:
        p["r"] = p["f"] / planck_nu(p["w"])
        p["re"] = p["e"] / planck_nu(p["w"])
    ref = [p["r"] for p in pts if norm_lo <= p["w"] <= norm_hi and p["f"] > 0]
    if len(ref) < 3:
        ref = sorted(p["r"] for p in pts if p["f"] > 0)[-5:]
    scale = float(np.median(ref)) if ref else 1.0
    for p in pts:
        p["r"] = round(float(p["r"] / scale), 5)
        p["re"] = round(float(p["re"] / scale), 5)
        p["f"] = round(float(p["f"]), 3)
        p["e"] = round(float(p["e"]), 3)
        p["sky"] = round(float(p["sky"]), 3)
    return pts


def clean(pts, k=4.0):
    """Drop single exposure outliers (cosmic rays, hot pixels, bad subtractions).

    A point is dropped if it is far outside what Pluto like bodies can reflect,
    or more than k sigma away from the median of its spectral neighbors.
    """
    if not pts:
        return pts, 0
    w = np.array([p["w"] for p in pts])
    r = np.array([p["r"] for p in pts])
    e = np.array([p["re"] for p in pts])
    keep = (r < 1.6) & (r > -0.3)
    idx = np.arange(len(pts))
    for i in idx:
        if not keep[i]:
            continue
        nb = keep & (np.abs(w - w[i]) < 0.04 * w[i]) & (idx != i)
        if nb.sum() >= 2:
            med = float(np.median(r[nb]))
            mad = 1.4826 * float(np.median(np.abs(r[nb] - med)))
            if abs(r[i] - med) > k * math.hypot(e[i], max(mad, 0.03)):
                keep[i] = False
    return [p for p, ok in zip(pts, keep) if ok], int((~keep).sum())


def binned(pts, R=40.0):
    if not pts:
        return []
    w = np.array([p["w"] for p in pts])
    r = np.array([p["r"] for p in pts])
    e = np.array([max(p["re"], 1e-6) for p in pts])
    edges = [w.min() * 0.999]
    while edges[-1] < w.max():
        edges.append(edges[-1] * (1 + 1 / R))
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (w >= a) & (w < b)
        if sel.sum() == 0:
            continue
        rr, ee = r[sel], e[sel]
        med = float(np.median(rr))
        if sel.sum() >= 3:
            sig = 1.2533 * float(np.std(rr)) / math.sqrt(sel.sum())
        else:
            sig = float(np.sqrt(np.sum(ee ** 2))) / sel.sum()
        out.append(dict(w=round(float(np.median(w[sel])), 4), r=round(med, 4),
                        e=round(float(sig), 4), n=int(sel.sum())))
    return out


def save_stamp(arr, path, scale=8):
    from PIL import Image
    a = np.asarray(arr, float)
    fin = a[np.isfinite(a)]
    if fin.size < 10:
        return False
    med = np.median(fin)
    sig = 1.4826 * np.median(np.abs(fin - med)) + 1e-12
    lo, hi = med - 2.5 * sig, med + 30 * sig
    v = np.clip((np.nan_to_num(a, nan=med) - lo) / (hi - lo), 0, 1) ** 0.5
    img = Image.fromarray((v * 255).astype(np.uint8), "L")
    img = img.resize((a.shape[1] * scale, a.shape[0] * scale), Image.NEAREST)
    img.save(path, optimize=True)
    return True


def pick_stamps(pts, stamps, wanted, out_dir):
    chosen = []
    for target_w in wanted:
        best = None
        for p in pts:
            key = p["obs"] + f"_D{p['det']}"
            if key not in stamps or abs(p["w"] - target_w) > 0.04 * target_w:
                continue
            # cleanest spot first: least starlight under the target, then closest wavelength
            score = (abs(p["sky"]) / max(p["e"], 1e-9), abs(p["w"] - target_w))
            if best is None or score < best[0]:
                best = (score, p, key)
        if best is None:
            continue
        _, p, key = best
        name = f"stamp_{p['w']:.3f}um.png"
        if save_stamp(stamps[key], out_dir / name):
            chosen.append(dict(file=name, w=p["w"], obs=p["obs"], det=p["det"],
                               snr=round(p["f"] / p["e"], 1) if p["e"] > 0 else None))
    return chosen


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", required=True, help="JPL Horizons id, e.g. 999 for Pluto")
    ap.add_argument("--name", required=True)
    ap.add_argument("--rows", action="append", default=[], help="cached Level 2 rows (.npz from evren_build)")
    ap.add_argument("--bands", default="D1,D4", help="detectors to read for the cached rows")
    ap.add_argument("--scan", default="", help="week folders to scan, e.g. 2025W40_1A,2025W40_2B")
    ap.add_argument("--scan-dets", default="2,3", help="detectors to scan")
    ap.add_argument("--scan-collection", default="dr1")
    ap.add_argument("--stamps", default="0.95,2.6,3.3,3.75", help="wavelengths for image stamps, microns")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--out", required=True, help="output folder")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    bands = [b.strip().upper() for b in args.bands.split(",") if b.strip()]

    frames = frames_from_rows(args.rows, bands) if args.rows else []
    print(f"{len(frames)} frames from cached rows ({', '.join(bands)})", flush=True)
    times = [f.t for f in frames]

    track = None
    if args.scan:
        weeks = [w.strip() for w in args.scan.split(",") if w.strip()]
        days = []
        for wk in weeks:
            m = re.search(r"(\d{4})W(\d{2})", wk)
            if m:
                d0 = datetime.strptime(f"{m.group(1)}-W{int(m.group(2)):02d}-1", "%G-W%V-%u")
                days += [(d0 - datetime(1858, 11, 17)).days + k for k in range(-1, 9)]
        track = Track(args.target, sorted(set(times + days)))
        dets = [int(x) for x in args.scan_dets.split(",") if x.strip()]
        more = []
        by_col = {}
        for wk in weeks:                         # tokens may read "qr3:2026W17"
            col, _, w = wk.rpartition(":")
            by_col.setdefault(col or args.scan_collection, []).append(w)
        for col, wks in by_col.items():
            more += scan_frames(wks, dets, track, collection=col, workers=24)
        pair = {1: 4, 2: 5, 3: 6, 4: 1, 5: 2, 6: 3}
        more += [Frame(band_url(f.url, pair[f.det]), f.obs_id, f.t) for f in list(more) if f.det in pair]
        have = {f.url for f in frames}
        frames += [f for f in more if f.url not in have]
        print(f"{len(frames)} frames after the scan", flush=True)
    if not frames:
        sys.exit("no frames")

    with cf.ThreadPoolExecutor(args.workers) as ex:
        loaded = []
        for f in ex.map(lambda fr: fr if fr.wcs is not None else fr.load(), frames):
            loaded.append(f)
    frames = [f for f in loaded if f.t is not None]
    if track is None:
        track = Track(args.target, [f.t for f in frames])
    print(f"headers read for {len(frames)} frames in {time.time() - t_start:.0f}s; "
          f"BUNIT {frames[0].h.get('BUNIT')}", flush=True)

    pts, stamps = measure(frames, track, workers=args.workers)
    pts.sort(key=lambda p: p["w"])
    pts = reflectance(pts)
    pts, dropped = clean(pts)
    print(f"dropped {dropped} outlier exposures", flush=True)
    wanted = [float(x) for x in args.stamps.split(",") if x.strip()]
    chosen = pick_stamps(pts, stamps, wanted, out)
    res = dict(target=args.name, horizons_id=args.target, n_points=len(pts),
               detectors=sorted({p["det"] for p in pts}),
               t_first=min(p["t"] for p in pts) if pts else None,
               t_last=max(p["t"] for p in pts) if pts else None,
               aperture_px=R_AP, sunlight="Planck 5772 K", points=pts,
               binned=binned(pts), stamps=chosen)
    (out / "spectrum.json").write_text(json.dumps(res, separators=(",", ":")))
    print(f"{len(pts)} spectrum points from detectors {res['detectors']} in "
          f"{time.time() - t_start:.0f}s -> {out / 'spectrum.json'}", flush=True)
    for b in res["binned"]:
        print(f"  {b['w']:.3f} um  R={b['r']:.3f} +- {b['e']:.3f}  (n={b['n']})")


if __name__ == "__main__":
    main()

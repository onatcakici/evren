#!/usr/bin/env python3
"""
Evren: a blink comparator for NASA SPHEREx.
Code written with help from an AI assistant.

For one patch of sky this script:
  1. finds every SPHEREx Level 2 spectral image that covers it (IRSA SIA2),
  2. reads only the pixels it needs from each image (HTTP Range requests),
  3. aligns all frames on one north up grid using each image's own WCS,
  4. builds a median "static sky" and subtracts it from every frame,
  5. finds what changed, links detections into tracklets,
  6. checks tracklets against known small bodies (IMCCE SkyBoT, optionally JPL SB Ident),
  7. writes PNG frames and a meta.json for the web viewer.

Example:
  python evren_build.py --id ecl160 --name "Ecliptic field in Leo" --ra 161.54 --dec 7.82 --band D1
"""
import argparse
import concurrent.futures as cf
import json
import math
import os
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import warnings
import xml.etree.ElementTree as ET

import certifi
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
from PIL import Image
from scipy.ndimage import median_filter

warnings.simplefilter("ignore", FITSFixedWarning)
warnings.simplefilter("ignore", RuntimeWarning)
warnings.simplefilter("ignore", UserWarning)

CTX = ssl.create_default_context(cafile=certifi.where())
UA = {"User-Agent": "Evren-SPHEREx-blink/1.0 (NASA Space Apps 2026)"}
SIA_URL = "https://irsa.ipac.caltech.edu/SIA"
COLLECTIONS = ("spherex_dr1", "spherex_qr3")
SBIDENT_URL = "https://ssd-api.jpl.nasa.gov/sb_ident.api"
SKYBOT_URL = "https://ssp.imcce.fr/webservices/skybot/api/conesearch.php"
PIX = 6.15          # SPHEREx pixel scale, arcsec
BLOCK = 2880        # FITS block size


# ----------------------------------------------------------------- network

def get(url, headers=None, timeout=90, tries=6):
    h = dict(UA)
    if headers:
        h.update(headers)
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=h)
            with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(1.0 + 1.5 * i)
                continue
            raise
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            last = e
            time.sleep(1.0 + 1.5 * i)
    raise last


def votable_rows(xml_bytes):
    root = ET.fromstring(xml_bytes)
    tag = lambda el: el.tag.split("}")[-1]
    fields = [el.get("name") for el in root.iter() if tag(el) == "FIELD"]
    rows = []
    for tr in (el for el in root.iter() if tag(el) == "TR"):
        rows.append(dict(zip(fields, [td.text for td in tr if tag(td) == "TD"])))
    return rows


def find_images(ra, dec, band):
    out, seen = [], set()
    for coll in COLLECTIONS:
        url = (f"{SIA_URL}?COLLECTION={coll}"
               f"&POS=circle+{ra:.6f}+{dec:.6f}+0.002&MAXREC=10000")
        for r in votable_rows(get(url)):
            if r.get("energy_bandpassname") != f"SPHEREx-{band}":
                continue
            if r.get("dataproduct_subtype") not in (None, "science"):
                continue
            if r.get("obs_id") in seen:
                continue
            seen.add(r.get("obs_id"))
            r["collection"] = coll
            out.append(r)
    out.sort(key=lambda r: float(r["t_min"]))
    return out


# --------------------------------------------------------------- FITS bits

def _header_at(buf, pos):
    cards, p = [], pos
    while p + BLOCK <= len(buf):
        block = buf[p:p + BLOCK]
        p += BLOCK
        for i in range(36):
            card = block[i * 80:(i + 1) * 80]
            cards.append(card)
            if card[:8] == b"END     ":
                return b"".join(cards).decode("ascii", "replace"), p
    return None


def _data_size(h):
    naxis = h.get("NAXIS", 0)
    if naxis == 0:
        return 0
    n = 1
    for k in range(1, naxis + 1):
        n *= h[f"NAXIS{k}"]
    size = (abs(h["BITPIX"]) // 8) * h.get("GCOUNT", 1) * (h.get("PCOUNT", 0) + n)
    return int(math.ceil(size / BLOCK) * BLOCK)


def image_header(url):
    """Read headers with Range requests until the IMAGE extension header."""
    want = 96 * 1024
    buf = get(url, {"Range": f"bytes=0-{want - 1}"})
    pos = 0
    for _ in range(6):
        res = _header_at(buf, pos)
        while res is None:
            want *= 2
            buf = get(url, {"Range": f"bytes=0-{want - 1}"})
            res = _header_at(buf, pos)
        hstr, dstart = res
        h = fits.Header.fromstring(hstr)
        if str(h.get("EXTNAME", "")).strip().upper() == "IMAGE":
            return h, dstart
        pos = dstart + _data_size(h)
    raise RuntimeError("IMAGE extension not found")


def read_rows(url, h, dstart, r0, r1):
    nx = h["NAXIS1"]
    bpp = abs(h["BITPIX"]) // 8
    a = dstart + r0 * nx * bpp
    b = dstart + (r1 + 1) * nx * bpp - 1
    raw = get(url, {"Range": f"bytes={a}-{b}"})
    expected = b - a + 1
    if len(raw) != expected:          # server ignored Range and sent everything
        raw = raw[a:b + 1]
    dtype = {-32: ">f4", -64: ">f8", 16: ">i2", 32: ">i4"}[h["BITPIX"]]
    arr = np.frombuffer(raw, dtype=dtype).reshape(r1 - r0 + 1, nx).astype(np.float32)
    bscale, bzero = float(h.get("BSCALE", 1.0)), float(h.get("BZERO", 0.0))
    if bscale != 1.0 or bzero != 0.0:
        arr = arr * bscale + bzero
    return arr


# ------------------------------------------------------------- geometry

def out_wcs(ra, dec, n):
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra, dec]
    w.wcs.crpix = [(n + 1) / 2.0, (n + 1) / 2.0]
    w.wcs.cdelt = [-PIX / 3600.0, PIX / 3600.0]
    return w


def bilinear(img, x, y):
    hgt, wid = img.shape
    good = np.isfinite(x) & np.isfinite(y)
    x = np.where(good, x, -10.0)
    y = np.where(good, y, -10.0)
    x0 = np.floor(x).astype(np.int64)
    y0 = np.floor(y).astype(np.int64)
    fx = (x - x0).astype(np.float32)
    fy = (y - y0).astype(np.float32)
    ok = good & (x0 >= 0) & (y0 >= 0) & (x0 < wid - 1) & (y0 < hgt - 1)
    x0 = np.clip(x0, 0, wid - 2)
    y0 = np.clip(y0, 0, hgt - 2)
    v = (img[y0, x0] * (1 - fx) * (1 - fy) + img[y0, x0 + 1] * fx * (1 - fy)
         + img[y0 + 1, x0] * (1 - fx) * fy + img[y0 + 1, x0 + 1] * fx * fy)
    v = v.astype(np.float32)
    v[~ok] = np.nan
    return v


def fetch_frame(row, n, world):
    """Return an n x n frame on the output grid (row 0 = north), or None."""
    url = row["access_url"]
    h, dstart = image_header(url)
    w = WCS(h)
    nx, ny = h["NAXIS1"], h["NAXIS2"]
    xin, yin = w.all_world2pix(world[0], world[1], 0, quiet=True)
    inside = np.isfinite(xin) & np.isfinite(yin) & (xin > -1) & (yin > -1) & (xin < nx) & (yin < ny)
    if inside.mean() < 0.5:
        return None
    r0 = max(int(np.floor(yin[inside].min())) - 2, 0)
    r1 = min(int(np.ceil(yin[inside].max())) + 2, ny - 1)
    c0 = max(int(np.floor(xin[inside].min())) - 2, 0)
    c1 = min(int(np.ceil(xin[inside].max())) + 2, nx - 1)
    sub = read_rows(url, h, dstart, r0, r1)[:, c0:c1 + 1]
    img = bilinear(sub, xin - c0, yin - r0).reshape(n, n)
    return img[::-1].copy()   # flip so row 0 is north (top)


# ------------------------------------------------------------- detection

def robust_sigma(a):
    a = a[np.isfinite(a)]
    if a.size < 20:
        return np.nan
    return 1.4826 * np.median(np.abs(a - np.median(a)))


def detect(snr, frame, resid, thr, static_snr=None, med=None):
    """Local maxima above thr. Skips spots where the static sky already has a source:
    stars change brightness from one SPHEREx wavelength to the next, movers do not sit on stars."""
    n = snr.shape[0]
    s = np.where(np.isfinite(snr), snr, -99.0)
    pad = np.pad(s, 1, constant_values=-99.0)
    neigh = np.max([pad[1 + dy:1 + dy + n, 1 + dx:1 + dx + n]
                    for dy in (-1, 0, 1) for dx in (-1, 0, 1) if (dy, dx) != (0, 0)], axis=0)
    ys, xs = np.nonzero((s > thr) & (s >= neigh))
    found = []
    oy, ox = np.mgrid[-1:2, -1:2]
    for y, x in zip(ys, xs):
        if x < 3 or y < 3 or x >= n - 3 or y >= n - 3:
            continue
        if not np.isfinite(frame[y - 2:y + 3, x - 2:x + 3]).all():
            continue
        flux = float(np.nansum(resid[y - 1:y + 2, x - 1:x + 2]))
        if static_snr is not None and np.nanmax(static_snr[y - 1:y + 2, x - 1:x + 2]) > 3.0:
            # a star sits here: keep it only if the new light clearly outshines what
            # wavelength to wavelength changes of that star could explain
            sflux = float(np.nansum(np.clip(med[y - 1:y + 2, x - 1:x + 2], 0, None)))
            if flux < 0.5 * sflux:
                continue
        win = np.clip(s[y - 1:y + 2, x - 1:x + 2], 0, None)
        tot = win.sum()
        found.append((float(x + (win * ox).sum() / tot), float(y + (win * oy).sum() / tot), float(s[y, x]), flux))
    return found


def link_tracklets(dets, mjd, max_dt_h=3.0, max_rate=150.0, tol_px=1.3, min_frames=3):
    """Union detections from different frames that are consistent with one moving source.

    dets: (frame, x, y, snr, flux). A real asteroid shows up in every exposure of a
    visit (a few minutes apart) at nearly the same spot with similar brightness;
    cosmic rays, hot pixels and glints show up once.
    """
    parent = list(range(len(dets)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    dets = sorted(dets, key=lambda d: mjd[d[0]])
    for a in range(len(dets)):
        fa, xa, ya, _, fla = dets[a]
        for b in range(a + 1, len(dets)):
            fb, xb, yb, _, flb = dets[b]
            dt_h = (mjd[fb] - mjd[fa]) * 24.0
            if dt_h > max_dt_h:
                break
            if fa == fb or dt_h == 0:
                continue
            if fla <= 0 or flb <= 0 or max(fla, flb) / min(fla, flb) > 4.0:
                continue
            d = math.hypot(xb - xa, yb - ya)
            if d <= tol_px + max_rate * dt_h / PIX:
                parent[root(a)] = root(b)
    groups = {}
    for i in range(len(dets)):
        groups.setdefault(root(i), []).append(i)

    tracks = []
    for members in groups.values():
        frames = sorted({dets[i][0] for i in members})
        if len(frames) < min_frames:
            continue
        # keep the strongest detection per frame
        best = {}
        for i in members:
            f = dets[i][0]
            if f not in best or dets[i][3] > dets[best[f]][3]:
                best[f] = i
        pts = [dets[best[f]] for f in frames]
        t = np.array([mjd[p[0]] for p in pts]) * 24.0
        x = np.array([p[1] for p in pts])
        y = np.array([p[2] for p in pts])
        t0 = t.mean()
        if np.ptp(t) > 0:
            A = np.vstack([t - t0, np.ones_like(t)]).T
            vx, x0 = np.linalg.lstsq(A, x, rcond=None)[0]
            vy, y0 = np.linalg.lstsq(A, y, rcond=None)[0]
        else:
            vx = vy = 0.0
            x0, y0 = x.mean(), y.mean()
        rms = float(np.sqrt(np.mean((x - (x0 + vx * (t - t0))) ** 2 + (y - (y0 + vy * (t - t0))) ** 2)))
        combined = math.sqrt(sum(p[3] ** 2 for p in pts))
        if rms > 0.8 or max(p[3] for p in pts) < 6.0 or combined < 11.0:
            continue
        rate = float(math.hypot(vx, vy) * PIX)   # arcsec per hour
        tracks.append({
            "frames": [int(p[0]) for p in pts],
            "x": [round(float(v), 2) for v in x],
            "y": [round(float(v), 2) for v in y],
            "snr": [round(float(p[3]), 1) for p in pts],
            "rate_arcsec_h": round(rate, 1),
            "rms_px": round(rms, 2),
            "mjd_mid": float(t0 / 24.0),
        })
    return tracks


# ------------------------------------------------------- known objects

def sexa_ra(deg):
    h = (deg % 360.0) / 15.0
    hh = int(h)
    mm = int((h - hh) * 60)
    ss = ((h - hh) * 60 - mm) * 60
    return f"{hh:02d}-{mm:02d}-{ss:05.2f}"


def sexa_dec(deg):
    sign = "M" if deg < 0 else ""
    d = abs(deg)
    dd = int(d)
    mm = int((d - dd) * 60)
    ss = ((d - dd) * 60 - mm) * 60
    return f"{sign}{dd:02d}-{mm:02d}-{ss:04.1f}"


def parse_sexa(s, hours):
    s = s.replace(":", " ").replace("'", " ").replace('"', " ").strip()
    neg = s.startswith("-")
    parts = [float(p) for p in s.replace("+", "").replace("-", "").split()]
    while len(parts) < 3:
        parts.append(0.0)
    v = parts[0] + parts[1] / 60 + parts[2] / 3600
    v = -v if neg else v
    return v * 15.0 if hours else v


def known_skybot(ra, dec, radius_deg, mjd, vlim=21.0):
    """Known asteroids and comets in the field at time mjd (IMCCE SkyBoT, geocentric)."""
    q = {"-ep": f"{mjd + 2400000.5:.6f}", "-ra": f"{ra:.6f}", "-dec": f"{dec:.6f}",
         "-rd": f"{radius_deg:.4f}", "-mime": "text", "-output": "all", "-loc": "500",
         "-filter": "120", "-objFilter": "111", "-refsys": "EQJ2000", "-from": "Evren"}
    txt = get(SKYBOT_URL + "?" + urllib.parse.urlencode(q), timeout=120).decode("utf-8", "replace")
    names, out = None, []
    for line in txt.splitlines():
        if line.startswith("# Num"):
            names = [c.strip().lower() for c in line[1:].split("|")]
            continue
        if line.startswith("#") or not line.strip() or names is None:
            continue
        cols = [c.strip() for c in line.split("|")]
        rec = dict(zip(names, cols))
        try:
            v = float(rec.get("mv", "nan"))
            if not v <= vlim:
                continue
            out.append({
                "name": rec.get("name", "?"),
                "ra": parse_sexa(rec["ra(h)"], True),
                "dec": parse_sexa(rec["de(deg)"], False),
                "vmag": v,
                "cls": rec.get("class", ""),
                "rate_ra": float(rec.get("dra(arcsec/h)", 0) or 0),
                "rate_dec": float(rec.get("ddec(arcsec/h)", 0) or 0),
            })
        except (KeyError, ValueError):
            continue
    return out


def known_objects(ra, dec, half_deg, mjd, vlim=20.0):
    t = Time(mjd, format="mjd", scale="utc").strftime("%Y-%m-%d_%H:%M:%S")
    q = {
        "mpc-code": "500", "obs-time": t,
        "fov-ra-center": sexa_ra(ra), "fov-dec-center": sexa_dec(dec),
        "fov-ra-hwidth": f"{half_deg:.4f}", "fov-dec-hwidth": f"{half_deg:.4f}",
        "two-pass": "true", "suppress-first-pass": "true",
        "mag-required": "true", "req-elem": "false", "vmag-lim": f"{vlim}",
    }
    url = SBIDENT_URL + "?" + urllib.parse.urlencode(q, safe=":-_")
    d = json.loads(get(url, timeout=120))
    fields = d.get("fields_second") or d.get("fields_first") or []
    data = d.get("data_second_pass") or d.get("data_first_pass") or []

    def col(*keys):
        for i, f in enumerate(fields):
            fl = f.lower()
            if all(k in fl for k in keys):
                return i
        return None

    i_name, i_ra, i_dec = col("name"), col("ra", "hh"), col("dec", "dd")
    i_v = col("magnitude") if col("magnitude") is not None else col("v")
    i_rra, i_rdec = col("ra", "rate"), col("dec", "rate")
    out = []
    for row in data:
        try:
            out.append({
                "name": row[i_name].strip(),
                "ra": parse_sexa(row[i_ra], True),
                "dec": parse_sexa(row[i_dec], False),
                "vmag": float(row[i_v]) if i_v is not None and row[i_v] not in (None, "", "n.a.") else None,
                "rate_ra": float(row[i_rra]) if i_rra is not None else 0.0,
                "rate_dec": float(row[i_rdec]) if i_rdec is not None else 0.0,
            })
        except (TypeError, ValueError, IndexError):
            continue
    return out, fields


# ---------------------------------------------------------------- output

def stretch(img, lo, hi, soft=None):
    soft = max(soft if soft else (hi - lo) * 0.04, 1e-9)
    v = np.arcsinh((img - lo) / soft) / np.arcsinh((hi - lo) / soft)
    v = np.where(np.isfinite(v), v, 0.0)
    return (np.clip(v, 0, 1) * 255).astype(np.uint8)


def resid_png(snr, lim=6.0):
    v = np.where(np.isfinite(snr), snr / lim * 127.5 + 127.5, 20.0)
    return np.clip(v, 0, 255).astype(np.uint8)


def build(args):
    ra, dec, band = args.ra, args.dec, args.band
    n = int(round(args.size * 3600 / PIX)) | 1
    ow = out_wcs(ra, dec, n)
    jj, ii = np.mgrid[0:n, 0:n]
    wra, wdec = ow.all_pix2world(ii.ravel(), jj.ravel(), 0)

    cache_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache")
    os.makedirs(cache_dir, exist_ok=True)
    cpath = os.path.join(cache_dir, f"{args.id}_{band}_{n}.npz")
    if os.path.exists(cpath) and not args.refresh:
        z = np.load(cpath)
        frames = list(zip(json.loads(str(z["rows"])), list(z["cube"])))
        print(f"loaded {len(frames)} cached frames", flush=True)
        return analyse(args, frames, n, ow)

    rows = find_images(ra, dec, band)
    print(f"{len(rows)} SPHEREx {band} images cover RA {ra:.4f} Dec {dec:.4f}", flush=True)
    if args.max_frames and len(rows) > args.max_frames:
        keep = np.linspace(0, len(rows) - 1, args.max_frames).round().astype(int)
        rows = [rows[k] for k in sorted(set(keep))]

    frames = []
    t0 = time.time()
    with cf.ThreadPoolExecutor(args.workers) as ex:
        futs = {ex.submit(fetch_frame, r, n, (wra, wdec)): r for r in rows}
        for f in cf.as_completed(futs):
            r = futs[f]
            try:
                img = f.result()
            except Exception as e:     # keep going if one image fails
                print(f"  skipped {r.get('obs_id')}: {e}", flush=True)
                continue
            if img is None or np.isfinite(img).mean() < 0.6:
                continue
            frames.append((r, img))
            print(f"  {len(frames)}/{len(rows)} frames  {time.time() - t0:.0f}s", flush=True)
    frames.sort(key=lambda fr: float(fr[0]["t_min"]))
    if len(frames) < 3:
        sys.exit("Not enough usable frames here.")
    np.savez_compressed(cpath, cube=np.stack([im for _, im in frames]),
                        rows=json.dumps([r for r, _ in frames]))
    return analyse(args, frames, n, ow)


def analyse(args, frames, n, ow):
    ra, dec, band = args.ra, args.dec, args.band

    cube = np.stack([im for _, im in frames]).astype(np.float32)
    for k in range(len(cube)):       # remove zodiacal light gradients and stray light
        c = cube[k]
        fill = np.where(np.isfinite(c), c, np.nanmedian(c))
        cube[k] = c - median_filter(fill, size=15, mode="nearest")
    med = np.nanmedian(cube, axis=0)
    bright = np.isfinite(med) & (med > np.nanpercentile(med, 97))
    scale = np.array([np.nanmedian(c[bright] / med[bright]) for c in cube])
    scale = np.where(np.isfinite(scale) & (scale > 0.05), scale, 1.0)
    cube /= scale[:, None, None]
    med = np.nanmedian(cube, axis=0)
    resid = cube - med[None]

    faint = np.isfinite(med) & (med < np.nanpercentile(med, 70))
    static_snr = med / robust_sigma(med[faint])
    noise_disp = float(np.nanmedian([robust_sigma(c[faint]) for c in cube]))
    mjd = [(float(r["t_min"]) + float(r["t_max"])) / 2.0 for r, _ in frames]
    dets, snrs, noisy = [], [], []
    for k in range(len(cube)):
        s_bg = robust_sigma(resid[k][faint])
        sig = np.sqrt(s_bg ** 2 + (args.alpha * np.clip(med, 0, None)) ** 2)
        snr = resid[k] / sig
        snrs.append(snr)
        fin = snr[np.isfinite(snr)]
        if fin.size == 0 or np.mean(np.abs(fin) > 5) > 0.015:
            noisy.append(k)          # stray light edges or streaks: show it, do not search it
            continue
        for x, y, s, fl in detect(snr, cube[k], resid[k], args.thr, static_snr, med):
            dets.append((k, x, y, s, fl))
    print(f"{len(noisy)} noisy frames skipped for detection", flush=True)
    print(f"{len(dets)} single frame detections above {args.thr} sigma", flush=True)
    tracks = link_tracklets(dets, mjd)
    print(f"{len(tracks)} tracklets (seen in 3 or more frames)", flush=True)

    # known small bodies, one JPL query per visit (frames within 1 hour)
    visits, cur = [], []
    for k in range(len(mjd)):
        if cur and (mjd[k] - mjd[cur[0]]) * 24 > 1.0:
            visits.append(cur)
            cur = []
        cur.append(k)
    if cur:
        visits.append(cur)
    known = {}
    half = args.size / 2.0 + 0.02
    kpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache", f"known_{args.id}_{band}.json")
    kcache = {}
    if os.path.exists(kpath) and not args.refresh:
        with open(kpath) as fh:
            kcache = json.load(fh)
    if not args.no_jpl:
        for v in visits:
            tq = float(np.mean([mjd[k] for k in v]))
            key = f"{tq:.5f}"
            try:
                if key in kcache:
                    objs = kcache[key]
                elif args.jpl:
                    objs, _ = known_objects(ra, dec, half, tq, args.vlim)
                else:
                    objs = known_skybot(ra, dec, half * 1.42, tq, args.vlim)
                kcache[key] = objs
            except Exception as e:
                print(f"  known object query failed for MJD {tq:.4f}: {e}", flush=True)
                continue
            for o in objs:
                xq, yq = ow.all_world2pix([[o["ra"], o["dec"]]], 0)[0]
                for k in v:
                    dt_h = (mjd[k] - tq) * 24.0
                    x = xq - o["rate_ra"] * dt_h / PIX
                    y_wcs = yq + o["rate_dec"] * dt_h / PIX
                    row = (n - 1) - y_wcs
                    if -2 <= x <= n + 1 and -2 <= row <= n + 1:
                        ent = known.setdefault(o["name"], {"name": o["name"], "vmag": o["vmag"],
                                                           "cls": o.get("cls", ""), "pos": {}})
                        ent["pos"][str(k)] = [round(float(x), 2), round(float(row), 2)]
        print(f"{len(known)} known small bodies cross this field during the visits", flush=True)
        with open(kpath, "w") as fh:
            json.dump(kcache, fh)

    # forced check: is there signal where each known object should be?
    for ent in known.values():
        ent["snr"] = {}
        for k, (x, y) in ent["pos"].items():
            xi, yi = int(round(x)), int(round(y))
            if 1 <= xi < n - 1 and 1 <= yi < n - 1 and int(k) not in noisy:
                w = snrs[int(k)][yi - 1:yi + 2, xi - 1:xi + 2]
                if np.isfinite(w).any():
                    ent["snr"][k] = round(float(np.nanmax(w)), 1)
        cap = 300.0 if (ent["vmag"] or 99) > 16 else 1e9   # a faint asteroid that bright is a coincident artifact
        vals = [v for v in ent["snr"].values() if v < cap]
        hits = sum(1 for v in vals if v >= 4.0)
        ent["recovered"] = bool(vals) and len(vals) == len(ent["snr"]) and (
            (len(vals) == 1 and vals[0] >= 8.0) or (len(vals) >= 2 and hits >= max(2, math.ceil(0.5 * len(vals)))))

    for tr in tracks:
        tr["known"] = None
        best = 99.0
        for ent in known.values():
            dsum, cnt = 0.0, 0
            for f, x, y in zip(tr["frames"], tr["x"], tr["y"]):
                p = ent["pos"].get(str(f))
                if p:
                    dsum += math.hypot(p[0] - x, p[1] - y)
                    cnt += 1
            if cnt and dsum / cnt < best:
                best = dsum / cnt
                if best <= args.match_px:
                    tr["known"] = ent["name"]
        tr["kind"] = "stationary" if tr["rate_arcsec_h"] < 2.0 and (max(mjd[f] for f in tr["frames"]) - min(mjd[f] for f in tr["frames"])) * 24 > 0.5 else "mover"
    # the same unmatched spot on different days does not move: a variable source, not a mover
    for a in tracks:
        if a["known"]:
            continue
        xa, ya = int(round(a["x"][0])), int(round(a["y"][0]))
        if 0 <= xa < n and 0 <= ya < n and np.nanmax(static_snr[max(ya - 1, 0):ya + 2, max(xa - 1, 0):xa + 2]) > 5.0:
            a["kind"] = "variable"
        for b in tracks:
            if b is a or b["known"]:
                continue
            if abs(b["mjd_mid"] - a["mjd_mid"]) > 0.5 and math.hypot(a["x"][0] - b["x"][0], a["y"][0] - b["y"][0]) < 1.5:
                a["kind"] = b["kind"] = "variable"
    tracks.sort(key=lambda t: (t["known"] is not None, -max(t["snr"])))
    for i, tr in enumerate(tracks):
        tr["id"] = i + 1

    out = os.path.join(args.out, args.id)
    os.makedirs(out, exist_ok=True)
    lo = -2.0 * noise_disp
    hi = max(float(np.nanpercentile(med, 99.8)), 25.0 * noise_disp)
    soft = 2.5 * noise_disp
    Image.fromarray(stretch(med, lo, hi, soft)).save(os.path.join(out, "median.png"))
    meta_frames = []
    gif = []
    for k, ((r, _), c) in enumerate(zip(frames, cube)):
        fn, rn = f"f{k:03d}.png", f"r{k:03d}.png"
        im8 = stretch(c, lo, hi, soft)
        Image.fromarray(im8).save(os.path.join(out, fn))
        Image.fromarray(resid_png(snrs[k])).save(os.path.join(out, rn))
        gif.append(Image.fromarray(im8).resize((n * 4, n * 4), Image.NEAREST))
        meta_frames.append({
            "i": k, "file": fn, "resid": rn,
            "t": Time(mjd[k], format="mjd").isot[:19] + "Z",
            "mjd": round(mjd[k], 6),
            "obs_id": r.get("obs_id"), "collection": r.get("collection"),
            "um": [round(float(r["em_min"]) * 1e6, 3), round(float(r["em_max"]) * 1e6, 3)],
            "cover": round(float(np.isfinite(c).mean()), 3),
            "noisy": k in noisy,
        })
    if gif:
        gif[0].save(os.path.join(out, "blink.gif"), save_all=True, append_images=gif[1:], duration=350, loop=0)

    meta = {
        "id": args.id, "name": args.name, "note": args.note,
        "ra": ra, "dec": dec, "band": band, "size_px": n, "pix_arcsec": PIX,
        "built": Time.now().isot[:19] + "Z",
        "frames": meta_frames, "tracklets": tracks, "known": list(known.values()),
    }
    with open(os.path.join(out, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)

    idx_path = os.path.join(args.out, "targets.json")
    idx = []
    if os.path.exists(idx_path):
        with open(idx_path) as fh:
            idx = [t for t in json.load(fh) if t["id"] != args.id]
    idx.append({"id": args.id, "name": args.name, "note": args.note, "band": band,
                "frames": len(meta_frames), "tracklets": len(tracks),
                "unknown": sum(1 for t in tracks if not t["known"] and t["kind"] == "mover")})
    with open(idx_path, "w") as fh:
        json.dump(idx, fh, indent=1)
    unk = [t for t in tracks if not t["known"] and t["kind"] == "mover"]
    print(f"done: {len(meta_frames)} frames, {len(tracks)} tracklets, {len(unk)} unmatched movers -> {out}")


def main():
    p = argparse.ArgumentParser(description="SPHEREx blink comparator builder")
    p.add_argument("--id", required=True)
    p.add_argument("--name", default="")
    p.add_argument("--note", default="")
    p.add_argument("--ra", type=float, required=True)
    p.add_argument("--dec", type=float, required=True)
    p.add_argument("--band", default="D1", choices=["D1", "D2", "D3", "D4", "D5", "D6"])
    p.add_argument("--size", type=float, default=0.2, help="field width in degrees")
    p.add_argument("--max-frames", type=int, default=90)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--thr", type=float, default=4.5)
    p.add_argument("--alpha", type=float, default=0.12)
    p.add_argument("--vlim", type=float, default=20.0)
    p.add_argument("--match-px", type=float, default=3.0)
    p.add_argument("--no-jpl", action="store_true", help="skip the known object check")
    p.add_argument("--jpl", action="store_true", help="use JPL SB Ident instead of IMCCE SkyBoT")
    p.add_argument("--refresh", action="store_true", help="download again even if cached")
    p.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "docs", "data"))
    build(p.parse_args())


if __name__ == "__main__":
    main()

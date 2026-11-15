#!/usr/bin/env python3
"""List SPHEREx frames that contain a moving solar system object.

Uses JPL Horizons for the object's path and IRSA SIA2 to find SPHEREx images
along it, then keeps the frames whose footprint holds the object at the moment
the frame was taken. Examples:
  python find_object_frames.py --target 999 --start 2025-05-01 --stop 2026-10-01        (Pluto)
  python find_object_frames.py --target "DES=C/2025 N1;CAP;" --start 2025-08-01 --stop 2025-08-25 --step 1
"""
import argparse
import concurrent.futures as cf
import json
import sys
import urllib.parse

import numpy as np
from astropy.time import Time

import evren_build as eb


def ephemeris(target, start, stop, step="6 h"):
    q = {"format": "text", "COMMAND": f"'{target}'", "OBJ_DATA": "'NO'", "MAKE_EPHEM": "'YES'",
         "EPHEM_TYPE": "'OBSERVER'", "CENTER": "'500@399'", "START_TIME": f"'{start}'",
         "STOP_TIME": f"'{stop}'", "STEP_SIZE": f"'{step}'", "QUANTITIES": "'1,9'", "ANG_FORMAT": "'DEG'"}
    txt = eb.get("https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(q), timeout=60, tries=3).decode()
    if "$$SOE" not in txt:
        sys.exit("Horizons did not return an ephemeris:\n" + txt[:800])
    body = txt.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    t, ra, de, mag = [], [], [], []
    for line in body:
        p = line.split()
        try:
            t.append(Time.strptime(p[0] + " " + p[1], "%Y-%b-%d %H:%M").mjd)
            ra.append(float(p[2]))
            de.append(float(p[3]))
            mag.append(float(p[4]) if len(p) > 4 and p[4] != "n.a." else np.nan)
        except (ValueError, IndexError):
            continue
    return np.array(t), np.array(ra), np.array(de), np.array(mag)


def inside(poly, ra, dec):
    """Point in polygon on a small patch of sky (handles RA wrap)."""
    xs = np.array(poly[0::2]); ys = np.array(poly[1::2])
    xs = (xs - ra + 180.0) % 360.0 - 180.0
    xs *= np.cos(np.radians(dec))
    ys = ys - dec
    hit = False
    j = len(xs) - 1
    for i in range(len(xs)):
        if (ys[i] > 0) != (ys[j] > 0) and 0 < (xs[j] - xs[i]) * (0 - ys[i]) / (ys[j] - ys[i]) + xs[i]:
            hit = not hit
        j = i
    return hit


def search(t, ra, de, half, coll):
    url = f"{eb.SIA_URL}?COLLECTION={coll}&POS=circle+{ra:.5f}+{de:.5f}+{half['r']:.3f}&MAXREC=50000"
    if half.get("t"):
        url += f"&TIME={t - half['t']:.4f}+{t + half['t']:.4f}"
    try:
        return eb.votable_rows(eb.get(url, timeout=60, tries=3))
    except Exception as e:
        print("archive search failed:", e, file=sys.stderr)
        return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, help="JPL Horizons id, for example 999 for Pluto")
    ap.add_argument("--start", required=True)
    ap.add_argument("--stop", required=True)
    ap.add_argument("--band", default="D1")
    ap.add_argument("--step", type=float, default=20.0, help="days between archive searches")
    ap.add_argument("--json", help="also write the matching frames here")
    a = ap.parse_args()
    T, RA, DE, MAG = ephemeris(a.target, a.start, a.stop)
    # sample the track every 0.4 degrees and search a 0.3 degree circle around each point
    pts = [0]
    for i in range(1, len(T)):
        d = np.hypot((RA[i] - RA[pts[-1]]) * np.cos(np.radians(DE[i])), DE[i] - DE[pts[-1]])
        if d > 0.4:
            pts.append(i)
    if pts[-1] != len(T) - 1:
        pts.append(len(T) - 1)
    jobs = []
    with cf.ThreadPoolExecutor(6) as ex:
        for i in pts:
            for coll in eb.COLLECTIONS:
                jobs.append(ex.submit(search, T[i], RA[i], DE[i], {"r": 0.3, "t": None}, coll))
        rows = [r for f in jobs for r in f.result()]
    rows = [r for r in rows if T[0] <= float(r["t_min"]) <= T[-1]]
    found = {}
    for r in rows:
        if r.get("energy_bandpassname") != f"SPHEREx-{a.band}" or r["obs_id"] in found:
            continue
        tm = (float(r["t_min"]) + float(r["t_max"])) / 2
        ra, de = float(np.interp(tm, T, RA)), float(np.interp(tm, T, DE))
        reg = (r.get("s_region") or "").replace("POLYGON", "").replace("ICRS", "").split()
        try:
            poly = [float(v) for v in reg]
        except ValueError:
            continue
        if len(poly) >= 6 and inside(poly, ra, de):
            found[r["obs_id"]] = {"t": Time(tm, format="mjd").isot[:16], "mjd": tm, "obs_id": r["obs_id"],
                                  "ra": round(ra, 5), "dec": round(de, 5),
                                  "vmag": None if np.isnan(np.interp(tm, T, MAG)) else round(float(np.interp(tm, T, MAG)), 1)}
    out = sorted(found.values(), key=lambda f: f["mjd"])
    for f in out:
        print(f["t"], f["obs_id"], f"object at RA {f['ra']:.4f} Dec {f['dec']:.4f}", f"V {f['vmag']}")
    print(f"{len(out)} frames")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()

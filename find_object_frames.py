#!/usr/bin/env python3
"""List SPHEREx frames that contain a moving solar system object.

Uses JPL Horizons for the object's path and IRSA SIA2 to find SPHEREx images
along it. Example (Pluto is Horizons id 999):
  python find_object_frames.py --target 999 --start 2025-05-01 --stop 2026-10-01
"""
import argparse
import urllib.parse

import numpy as np
from astropy.time import Time

import evren_build as eb


def ephemeris(target, start, stop):
    q = {"format": "text", "COMMAND": f"'{target}'", "OBJ_DATA": "'NO'", "MAKE_EPHEM": "'YES'",
         "EPHEM_TYPE": "'OBSERVER'", "CENTER": "'500@399'", "START_TIME": f"'{start}'",
         "STOP_TIME": f"'{stop}'", "STEP_SIZE": "'6 h'", "QUANTITIES": "'1'", "ANG_FORMAT": "'DEG'"}
    txt = eb.get("https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(q)).decode()
    body = txt.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    t, ra, de = [], [], []
    for line in body:
        p = line.split()
        t.append(Time.strptime(p[0] + " " + p[1], "%Y-%b-%d %H:%M").mjd)
        ra.append(float(p[2]))
        de.append(float(p[3]))
    return np.array(t), np.array(ra), np.array(de)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, help="JPL Horizons id, for example 999 for Pluto")
    ap.add_argument("--start", required=True)
    ap.add_argument("--stop", required=True)
    ap.add_argument("--band", default="D1")
    ap.add_argument("--step", type=float, default=8.0, help="days between archive searches")
    a = ap.parse_args()
    T, RA, DE = ephemeris(a.target, a.start, a.stop)
    found = {}
    for t in np.arange(T[0] + a.step / 2, T[-1] - a.step / 2, a.step):
        ra, de = np.interp(t, T, RA), np.interp(t, T, DE)
        for coll in eb.COLLECTIONS:
            url = (f"{eb.SIA_URL}?COLLECTION={coll}&POS=circle+{ra:.5f}+{de:.5f}+0.05"
                   f"&TIME={t - a.step / 2:.3f}+{t + a.step / 2:.3f}&MAXREC=5000")
            try:
                rows = eb.votable_rows(eb.get(url))
            except Exception as e:
                print("archive search failed:", e)
                continue
            for r in rows:
                if r.get("energy_bandpassname") == f"SPHEREx-{a.band}":
                    found[r["obs_id"]] = r
    for r in sorted(found.values(), key=lambda r: float(r["t_min"])):
        tm = (float(r["t_min"]) + float(r["t_max"])) / 2
        print(Time(tm, format="mjd").isot[:16], r["obs_id"],
              f"object at RA {np.interp(tm, T, RA):.4f} Dec {np.interp(tm, T, DE):.4f}")
    print(f"{len(found)} frames")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Render the 30 second demo video for the Space Apps submission (1280x720, H.264)."""
import json
import os
import subprocess
import sys

from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "docs", "data")
W, H, FPS = 1280, 720, 30
BG = (11, 13, 16)
TXT = (232, 230, 227)
MUTED = (154, 160, 168)
ACCENT = (255, 90, 31)
KNOWN = (124, 199, 255)
URL = sys.argv[1] if len(sys.argv) > 1 else ""


def font(size, bold=False):
    for path, idx in (("/System/Library/Fonts/Helvetica.ttc", 1 if bold else 0),
                      ("/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf", 0)):
        try:
            return ImageFont.truetype(path, size, index=idx)
        except OSError:
            continue
    return ImageFont.load_default()


F_HUGE, F_BIG, F_MED, F_SM = font(96, True), font(44, True), font(30), font(22)


def load(field):
    with open(os.path.join(DATA, field, "meta.json")) as fh:
        return json.load(fh)


def panel(field, meta, k, kind, size=600, mark=None, path=None):
    src = meta["frames"][k]["resid" if kind == "resid" else "file"]
    im = Image.open(os.path.join(DATA, field, src)).convert("RGB")
    n = meta["size_px"]
    im = im.resize((size, size), Image.LANCZOS)
    d = ImageDraw.Draw(im)
    s = size / n
    if path:
        pts = [((x + 0.5) * s, (y + 0.5) * s) for x, y in path]
        for a, b in zip(pts, pts[1:]):
            d.line([a, b], fill=(90, 140, 180), width=2)
    for (x, y, label, color) in (mark or []):
        cx, cy, r = (x + 0.5) * s, (y + 0.5) * s, 4.0 * s
        d.ellipse([cx - r, cy - r, cx + r, cy + r], outline=color, width=3)
        tw = d.textlength(label, font=F_SM)
        lx = cx + r + 8 if cx + r + 8 + tw < size - 6 else cx - r - 8 - tw
        ly = max(cy - r - 8, 6)
        d.text((lx, ly), label, font=F_SM, fill=color)
    return im


def text_block(d, x, y, lines, gap=12):
    for f, color, line in lines:
        d.text((x, y), line, font=f, fill=color)
        y += f.size + gap
    return y


def frame_base():
    return Image.new("RGB", (W, H), BG)


def fade(img, a):
    if a >= 1:
        return img
    return Image.blend(frame_base(), img, max(a, 0.0))


def main():
    pluto = load("pluto2025")
    leo = load("ecl160")
    pent = next(e for e in pluto["known"] if "luto" in e["name"])
    pk = sorted(int(k) for k in pent["pos"])
    ppath = [pent["pos"][str(k)] for k in pk]

    out = os.path.join(HERE, "demo.mp4")
    ff = subprocess.Popen(["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
                           "-s", f"{W}x{H}", "-r", str(FPS), "-i", "-", "-c:v", "libx264",
                           "-preset", "veryfast", "-crf", "20", "-pix_fmt", "yuv420p", out],
                          stdin=subprocess.PIPE)
    count = 0

    def emit(img, seconds=None, frames=None):
        nonlocal count
        nfr = frames if frames is not None else int(round(seconds * FPS))
        raw = img.tobytes()
        for _ in range(nfr):
            ff.stdin.write(raw)
        count += nfr

    # 1. title
    t = frame_base()
    d = ImageDraw.Draw(t)
    d.ellipse([110, 262, 140, 292], fill=ACCENT)
    d.text((160, 220), "Evren", font=F_HUGE, fill=TXT)
    d.text((164, 340), "A blink comparator for NASA SPHEREx", font=F_MED, fill=MUTED)
    d.text((164, 600), "Team Evren, Onat Cakici", font=F_SM, fill=MUTED)
    d.text((164, 632), "NASA Space Apps Challenge 2026: Planet X and SPHEREx", font=F_SM, fill=MUTED)
    for i in range(15):
        emit(fade(t, (i + 1) / 15), frames=1)
    emit(t, 2.5)

    # 2. the idea
    t = frame_base()
    d = ImageDraw.Draw(t)
    text_block(d, 110, 230, [
        (F_BIG, TXT, "1930: Clyde Tombaugh found Pluto"),
        (F_BIG, TXT, "by blinking two photos of the same stars."),
        (F_MED, MUTED, ""),
        (F_MED, ACCENT, "2026: we ran his method on NASA's SPHEREx images."),
    ])
    emit(t, 3.0)

    # 3a. the blink itself: first and last visit, flipping like Tombaugh's plates
    good = sorted(int(k) for k, v in pent["snr"].items() if v >= 8)
    fa, fb = good[0], good[-1]
    for i in range(12):
        k = fa if i % 2 == 0 else fb
        img = frame_base()
        p = pent["pos"][str(k)]
        img.paste(panel("pluto2025", pluto, k, "frame", 600, mark=[(p[0], p[1], "Pluto", KNOWN)]), (80, 60))
        d = ImageDraw.Draw(img)
        tag = ("A  " if k == fa else "B  ") + pluto["frames"][k]["t"][:10]
        d.rectangle([92, 72, 92 + 16 + d.textlength(tag, font=F_SM), 108], fill=(0, 0, 0))
        d.text((100, 77), tag, font=F_SM, fill=TXT if k == fa else ACCENT)
        text_block(d, 730, 220, [(F_BIG, TXT, "Blink two dates."), (F_BIG, TXT, "One dot jumps.")])
        text_block(d, 730, 360, [(F_SM, MUTED, line) for line in (
            "October 6 and October 12, 2025.",
            "Exactly how Tombaugh saw Pluto",
            "move on his glass plates in 1930.")])
        emit(img, 0.5)

    # 3. Pluto blinking across the stars
    def pluto_scene(kind, caption, sub, loops, per):
        for _ in range(loops):
            for k in pk:
                img = frame_base()
                p = pent["pos"][str(k)]
                pan = panel("pluto2025", pluto, k, kind, 600,
                            mark=[(p[0], p[1], "Pluto", KNOWN)], path=ppath)
                img.paste(pan, (80, 60))
                d = ImageDraw.Draw(img)
                d.text((80, 672), pluto["frames"][k]["t"][:16].replace("T", "  ") + " UTC", font=F_SM, fill=MUTED)
                text_block(d, 730, 220, [(F_BIG, TXT, caption[0]), (F_BIG, TXT, caption[1])])
                text_block(d, 730, 360, [(F_SM, MUTED, line) for line in sub])
                emit(img, frames=per)

    pluto_scene("frame", ("Stars stay put.", "One dot drifts."),
                ["That dot is Pluto, October 6 to 12, 2025.",
                 "Evren found it on its own in four visits",
                 "and matched it to Pluto's orbit within",
                 "about two arcseconds."], loops=1, per=7)
    pluto_scene("resid", ("Subtract the sky", "that never changes."),
                ["What is left is only what changed.",
                 "A cosmic ray hits once. A real mover",
                 "shows up in every exposure of a visit."], loops=1, per=7)

    # 4. known asteroids recovered in the Leo field
    for name, color in (("Lexell", KNOWN), ("Rosa-Brusin", KNOWN), ("Irvine", KNOWN)):
        ent = next((e for e in leo["known"] if e["name"] == name), None)
        if not ent or not ent.get("snr"):
            continue
        k = int(max(ent["snr"].items(), key=lambda kv: kv[1])[0])
        p = ent["pos"][str(k)]
        img = frame_base()
        img.paste(panel("ecl160", leo, k, "frame", 600, mark=[(p[0], p[1], name, color)]), (80, 60))
        d = ImageDraw.Draw(img)
        d.text((80, 672), leo["frames"][k]["t"][:16].replace("T", "  ") + " UTC, ecliptic field in Leo", font=F_SM, fill=MUTED)
        text_block(d, 730, 220, [(F_BIG, TXT, "Checked against"), (F_BIG, TXT, "every known asteroid.")])
        text_block(d, 730, 360, [(F_SM, MUTED, line) for line in (
            f"{name}, V {ent['vmag']}, is right where",
            "its orbit says. Anything with no match",
            "becomes a candidate worth a second look.")])
        emit(img, 1.5)

    # 5. end card
    t = frame_base()
    d = ImageDraw.Draw(t)
    d.ellipse([110, 222, 140, 252], fill=ACCENT)
    d.text((160, 180), "Evren", font=F_HUGE, fill=TXT)
    text_block(d, 164, 310, [
        (F_MED, MUTED, "Every SPHEREx frame of a patch of sky, aligned,"),
        (F_MED, MUTED, "differenced and checked against every known asteroid."),
    ])
    if URL:
        d.text((164, 470), URL, font=F_MED, fill=ACCENT)
    d.text((164, 632), "Data: NASA SPHEREx via NASA IPAC IRSA. Known objects: IMCCE SkyBoT.", font=F_SM, fill=MUTED)
    remaining = max(30 * FPS - count, int(2.5 * FPS))
    emit(t, frames=remaining)

    ff.stdin.close()
    ff.wait()
    print(f"wrote {out}: {count / FPS:.1f} s")


if __name__ == "__main__":
    main()

# Evren: a blink comparator for NASA SPHEREx

In 1930 Clyde Tombaugh found Pluto by flipping between two photographs of the same stars taken six days apart. Stars stay put. Anything that jumps is moving.

Evren does the same with SPHEREx, the NASA space telescope that maps the whole sky in 102 infrared colors every six months. Pick a patch of sky and Evren pulls every SPHEREx image of it, lines them up, subtracts the sky that never changes and links what is left into moving objects. Every object is checked against all known asteroids and comets. Anyone can flip through the frames in a browser.

Built for the 2026 NASA Space Apps Challenge, **Planet X and SPHEREx**, by team Evren.

## What it found

* **Pluto, found again.** In SPHEREx images from October 2025 and April 2026, Evren detected a moving source in six separate visits and matched it to Pluto's orbit within about two arcseconds, with no hint about where to look. Up to 318 sigma.
* **Known asteroids recovered.** In an ecliptic field in Leo it picked up Lexell (V 17.2), Irvine (V 18.3) and Rosa-Brusin (V 18.6) right where their orbits put them. SPHEREx sees main belt asteroids down to about V 18.5 in a single exposure.
* **Candidates.** Anything that moves and matches no known object is listed for a second look, with image stamps from its own exposures and from other dates.

## How it works

1. **Find.** Query the NASA IPAC Infrared Science Archive (IRSA) SIA2 service for every SPHEREx Level 2 spectral image that covers the position (Data Release 1 and Quick Release 3).
2. **Read only what is needed.** Each image is about 70 MB. Evren reads the FITS headers and then only the rows it needs with HTTP Range requests, about 1 MB per image.
3. **Align.** Every frame is resampled onto one north up grid using its own TAN SIP world coordinate system.
4. **Subtract the static sky.** A local background removes zodiacal light and stray light. The median of all frames is the static sky. Each frame minus the static sky leaves only change.
5. **Detect and link.** Peaks above 4.5 sigma are found in each difference image. Spots on top of stars are kept only if the new light clearly outshines what the star's color changes between SPHEREx wavelengths could explain. Detections in several exposures of one visit that agree in position and brightness become a tracklet. A cosmic ray hits once, an asteroid stays.
6. **Check against known objects.** For every visit, IMCCE SkyBoT lists the known asteroids and comets in the field. Tracklets that match are labeled. Evren also measures the signal at every predicted position, which shows which known objects SPHEREx actually recorded.
7. **Show it.** The pipeline writes PNG frames and a JSON file per field. The viewer in `docs/` is plain HTML and JavaScript with no build step.

## Run it yourself

```
python3 -m venv .venv
.venv/bin/pip install numpy astropy scipy pillow certifi
.venv/bin/python evren_build.py --id myfield --name "My field" --ra 161.54 --dec 7.82 --band D1 --size 0.2
.venv/bin/python -m http.server 8765 --directory docs
```

Then open http://localhost:8765. Each run adds the field to `docs/data/targets.json`.

To find the SPHEREx frames that contain a moving solar system object (this is how the Pluto fields were chosen):

```
.venv/bin/python find_object_frames.py --target 999 --start 2025-05-01 --stop 2026-10-01
```

## Data and credits

* SPHEREx Level 2 spectral images: NASA, JPL Caltech, served by the NASA IPAC Infrared Science Archive (IRSA).
* Known asteroids and comets: IMCCE SkyBoT.
* Ephemerides used to locate Pluto: JPL Horizons.
* Built with Python, NumPy, Astropy, SciPy and Pillow.
* Code written with help from an AI assistant.

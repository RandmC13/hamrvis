"""Every scenario's bundle in a real browser: no page errors, every lane drawn; and what the page computes and prints
for a write equals what the bundle holds. Skipped without node and Playwright."""
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import unittest
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCEN = os.path.join(HERE, "scenarios")
sys.path.insert(0, ROOT)
from hamrvis import reduce  # noqa: E402


def have_playwright():
    return shutil.which("node") and subprocess.run(["node", "-e", "require('playwright')"], cwd=HERE,
                                                   capture_output=True).returncode == 0


@unittest.skipUnless(have_playwright() or os.environ.get("HAMRVIS_REQUIRE_VIEWER_TEST"), "needs node and Playwright")
class ViewerSmoke(unittest.TestCase):
    def test_every_scenario_draws(self):
        for name in sorted(os.listdir(SCEN)):
            with self.subTest(name), tempfile.TemporaryDirectory() as tmp:
                out = os.path.join(tmp, "bundle")
                reduce.reduce([os.path.join(SCEN, name)], out, workers=1)
                runs = [r["id"] for r in json.load(open(os.path.join(out, "manifest.json")))["runs"]]
                with socket.socket() as s:
                    s.bind(("127.0.0.1", 0))
                    port = s.getsockname()[1]
                view = subprocess.Popen([sys.executable, "-m", "hamrvis", "view", out, "-p", str(port), "--no-browser"],
                                        cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                try:
                    url = f"http://127.0.0.1:{port}/"
                    for _ in range(100):
                        try:
                            urllib.request.urlopen(url)
                            break
                        except OSError:
                            time.sleep(0.1)
                    r = subprocess.run(["node", "viewer_smoke.js", url] + runs, cwd=HERE, capture_output=True, text=True)
                    self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
                finally:
                    view.kill()
                    view.wait()


if __name__ == "__main__":
    unittest.main()


# ------------------------------------------------------------------------------------------------ viewer fidelity
import math  # noqa: E402
import re  # noqa: E402
from collections import Counter  # noqa: E402

import numpy as np  # noqa: E402


def serve(bundle):
    """Start `hamrvis view` on a free port; returns (process, url)."""
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        port = s.getsockname()[1]
    view = subprocess.Popen([sys.executable, "-m", "hamrvis", "view", bundle, "-p", str(port), "--no-browser"],
                            cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    url = f"http://127.0.0.1:{port}/"
    for _ in range(100):
        try:
            urllib.request.urlopen(url)
            break
        except OSError:
            time.sleep(0.1)
    return view, url


def num(text):
    """First number in a printed value (VAMPIRE-style minus signs included)."""
    return float(re.search(r"[-−+]?\d+(\.\d+)?", text).group(0).replace("−", "-"))


class Expected:
    """What the viewer should show for a write run, worked out from the bundle files alone (not from app.js)."""

    def __init__(self, bundle, e):
        self.e, self.p = e, e["hamr"]
        d = os.path.join(bundle, e["id"])
        nG, nF = len(e["grains"]["id"]), e["grain_series"]["frames"]
        g = np.fromfile(os.path.join(d, "grains.i16"), np.int16).reshape(nF, nG, 4)[-1].astype(np.float64) / 32767
        self.mz, self.ml = g[:, 2] * g[:, 3], g[:, 3]                    # last frame: m_z |m| and |m|
        self.gm = e.get("grain_map")
        if self.gm:
            self.gmap = np.fromfile(os.path.join(d, "grain_map.i16"), np.int16).reshape(self.gm["ny"], self.gm["nx"])
            grid = (e.get("snapshots") or {}).get("grid", {})
            self.a = self.gm.get("atom_pitch_nm") or grid.get("atom_pitch_nm") or self.gm["pitch_nm"]
            self.f = round(self.gm["pitch_nm"] / self.a)

    def centre(self, i):
        """Mean position of the atom columns k = i f ... i f + f - 1 of pixel i (the reducer bins k = rint(x / a) // f)."""
        return self.a * (i * self.f + (self.f - 1) / 2)

    def cells(self, lag):
        p, e = self.p, self.e
        Tc = e["Tc"]
        half = p["sigma_y"] * math.sqrt(2 * math.log((p["Tmax"] - p["Tmin"]) / (Tc - p["Tmin"]))) if Tc and p["Tmin"] < Tc < p["Tmax"] else 0
        hw, out = half or p["TW"] / 2, []
        for t in range(p["n_tracks"]):
            for b in range(p["bits_per_track"]):
                yc, x0 = p["pad"] + (t + .5) * p["TW"], b * p["BL"] - lag
                x1 = x0 + p["BL"]
                if self.gm:
                    ys = [iy for iy in range(self.gm["ny"]) if abs(self.centre(iy) - yc) <= hw]
                    xs = [ix for ix in range(self.gm["nx"]) if x0 <= self.centre(ix) < x1]
                    parts = Counter(int(k) for k in self.gmap[np.ix_(ys, xs)].ravel() if k >= 0)
                else:
                    gx, gy, at = e["grains"]["x"], e["grains"]["y"], e["grains"]["atoms"]
                    parts = Counter({k: at[k] for k in range(len(gx)) if x0 <= gx[k] < x1 and abs(gy[k] - yc) <= hw})
                out.append([t, b, int(np.sign(p["sequence"][t * p["bits_per_track"] + b] * p["field_uz"])), dict(parts)])
        return out

    def bits_ok(self, cells):
        ok = n = 0
        for _, _, want, parts in cells:
            if want and parts:
                n += 1
                ok += np.sign(sum(w * self.mz[k] for k, w in parts.items()) / sum(parts.values())) == want
        return ok, n

    def grains_ok(self, lag):
        p, ok, n = self.p, 0, 0
        for a, x in zip(self.e["analysis"]["grains"], self.e["grains"]["x"]):
            b = math.floor((x + lag) / p["BL"])
            k = a["track"] * p["bits_per_track"] + b
            want = int(np.sign((p["sequence"][k] if 0 <= b < p["bits_per_track"] and k < len(p["sequence"]) else 0) * p["field_uz"]))
            if a["category"] != "unheated" and want:
                n += 1
                ok += a["category"] != "undecided" and np.sign(a["final"]) == want
        return ok, n


@unittest.skipUnless(have_playwright() or os.environ.get("HAMRVIS_REQUIRE_VIEWER_TEST"), "needs node and Playwright")
class ViewerFidelity(unittest.TestCase):
    """What the page computes and prints for a write equals what the bundle holds (P5 layer 2)."""

    def check_bundle(self, bundle):
        M = json.load(open(os.path.join(bundle, "manifest.json"), encoding="utf-8"))
        runs = M.get("runs") or [dict(e, id=k) for k, e in M["stages"].items()]
        writes = [e for e in runs if e.get("analysis") and e.get("grain_series")]
        view, url = serve(bundle)
        try:
            r = subprocess.run(["node", "viewer_fidelity.js", url] + [e["id"] for e in writes], cwd=HERE, capture_output=True, text=True)
        finally:
            view.kill()
            view.wait()
        self.assertEqual(r.returncode, 0, r.stderr)
        got = json.loads(r.stdout)
        for e in writes:
            X = Expected(bundle, e)
            for mode, lag in (("inferred", e["freeze_lag_nm"]), ("none", 0.0)):
                with self.subTest(bundle=os.path.basename(bundle), run=e["id"], lag=mode):
                    v = got[e["id"]][mode]
                    self.assertEqual(v["lag"], lag)
                    cells = X.cells(lag)
                    self.assertEqual([[t, b, w, {int(k): c for k, c in parts}] for t, b, w, parts in v["cells"]], cells)
                    facts = dict(v["facts"])
                    ok, n = X.grains_ok(lag)
                    self.assertEqual(facts["melted grains matching their bit at end"], f"{ok} / {n} ({100 * ok / n:.0f} %)")
                    bok, bn = X.bits_ok(cells)
                    if bn:
                        self.assertEqual(facts["bits matching intended at end"], f"{bok} / {bn}")
                    if X.gm:
                        for x, y, g in v["hover"]:
                            ix, iy = int(np.rint(x / X.a)) // X.f, int(np.rint(y / X.a)) // X.f
                            self.assertEqual(g, int(X.gmap[iy, ix]), (x, y))
                    for g, rows in v["inspector"].items():
                        a, rows = e["analysis"]["grains"][int(g)], dict(rows)
                        self.assertAlmostEqual(num(rows["mz"]), X.mz[int(g)], delta=5e-4 + 1e-12)
                        self.assertAlmostEqual(num(rows["|m|"]), X.ml[int(g)], delta=5e-4 + 1e-12)
                        self.assertAlmostEqual(num(rows["final mz"]), a["final"], delta=5e-4 + 1e-12)
                        self.assertAlmostEqual(num(rows["head Bz then"]), a["B_dec"] * X.p["field_uz"], delta=5e-3 + 1e-12)
                        self.assertAlmostEqual(num(rows["Bz at grain then"]), a["Bz_dec"], delta=5e-3 + 1e-12)
                        self.assertEqual("(outside the field box)" in rows["Bz at grain then"], a["Bz_dec"] == 0 and a["B_dec"] != 0)
        return len(writes)

    def test_scenarios(self):
        n = 0
        with tempfile.TemporaryDirectory() as tmp:
            for name in ("chain", "two-tracks", "rewrite", "write-only"):
                out = os.path.join(tmp, name)
                reduce.reduce([os.path.join(SCEN, name)], out, workers=1)
                n += self.check_bundle(out)
        self.assertGreaterEqual(n, 5)

    def test_real_bundles(self):
        found = [os.path.join(ROOT, "bundles", b) for b in ("h3", "2x11") if os.path.exists(os.path.join(ROOT, "bundles", b, "manifest.json"))]
        if not found:
            self.skipTest("real-run bundles in bundles/ (git-ignored) not present")
        for b in found:
            self.check_bundle(b)

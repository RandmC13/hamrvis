"""What the viewer prints and plots equals what the bundle holds: the clock, film tooltips, bits lane and stray field
tooltips, and the chart series (scrubber, life story, space-time, score against lag, read-back waveform), at several
frames of every run. Expected values come from the bundle files and the independent head model in
hamrvis/verify/reference.py, never from app.js. Printed numbers must agree within the digits shown, chart values to
float rounding. Skipped without node and Playwright."""
import gzip
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # also when run as tests.<module>
from test_viewer_smoke import HERE, ROOT, SCEN, Expected, have_playwright, reduce, serve
from hamrvis.verify.reference import Head

CAT = {"full-field": "full field", "reversal": "field reversal", "zero-field": "field off", "unheated": "not melted",
       "undecided": "undecided", "empty": "no magnetic atoms"}
NUM = re.compile(r"[-−+]?\d[\d,]*(?:\.(\d+))?")


def nums(text):
    """[(value, half a unit of the last digit shown)] of the numbers in a printed value."""
    out = []
    for m in NUM.finditer(text):
        v = float(m.group(0).replace("−", "-").replace(",", ""))
        out.append((v, 0.5 * 10.0 ** -(len(m.group(1)) if m.group(1) else 0)))
    return out


def printed(text, value, unit=1.0):
    """True if the first number of text equals value (in units of unit) to the digits shown."""
    (v, tol), = nums(text)[:1] or [(math.nan, 0)]
    return abs(v - value / unit) <= tol + 1e-9 * abs(value / unit)


def printed_time(text, t, dt, mc):
    """A time as the viewer prints it: 'MC step N', 'x ps' (2 decimals, trailing zeros dropped) or 'x ns' (3)."""
    if mc:
        return nums(text)[0][0] == round(t / dt)
    v = nums(text)[0][0]
    if text.strip().endswith("ps"):
        return abs(v - t * 1e12) <= 0.005 + 1e-9
    return abs(v - t * 1e9) <= 0.0005 + 1e-9


class Run:
    """One run of a bundle, read from its files."""

    def __init__(self, bundle, e):
        self.e, self.d = e, os.path.join(bundle, e["id"])
        self.dt = e["dt"]
        # version 1 bundles hold no monte_carlo flag: their relax stage is Monte Carlo when its integrator is (bundle-format.md)
        self.mc = bool(e["monte_carlo"]) if "monte_carlo" in e else e["id"] == "relax" and "monte-carlo" in e["input"].get("sim:integrator", "")
        nG = len(e["grains"]["id"])
        gs = e.get("grain_series")
        if gs:
            f32 = gs.get("dtype") == "f32"
            raw = np.fromfile(os.path.join(self.d, "grains.f32" if f32 else "grains.i16"), np.float32 if f32 else np.int16)
            self.g = raw.reshape(gs["frames"], nG, 4).astype(np.float64)
            self.gs = 1.0 if f32 else 32767.0
            self.time = np.fromfile(os.path.join(self.d, "grain_lead.f64")).reshape(gs["frames"], -1)[:, 0]
        else:
            self.g = None
            snaps = (e.get("snapshots") or {}).get("frames") or []
            self.time = np.array([s["time"] for s in snaps]) if len(snaps) > 1 else np.zeros(1)
        self.steps = [int(k) for k in np.floor(self.time / self.dt + 0.5)]
        self.nF = len(self.time)
        self.head = Head(e["hamr"]) if e.get("hamr") else None
        self.gm = e.get("grain_map")
        if self.gm:
            self.gmap = np.fromfile(os.path.join(self.d, "grain_map.i16"), np.int16).reshape(self.gm["ny"], self.gm["nx"])
        sn = e.get("snapshots")
        if sn:
            G = self.grid = sn["grid"]
            self.snaps = [(s["index"], math.floor(s["time"] / self.dt + 0.5)) for s in sn["frames"]]
            self.count = np.fromfile(os.path.join(self.d, "top_count.u8"), np.uint8).reshape(G["ny"], G["nx"])
            self.minm = np.fromfile(os.path.join(self.d, "min_m.u8"), np.uint8).reshape(G["ny"], G["nx"])
        if e.get("dipole"):
            P = e["dipole"]
            self.B = np.fromfile(os.path.join(self.d, "dipole.f32"), np.float32).reshape(3, P["nz"], P["ny"], P["nx"]).astype(np.float64)
        cols = e["output"]["columns"]
        f64 = e["output"].get("dtype") == "f64"
        self.out = np.fromfile(os.path.join(self.d, "output.f64" if f64 else "output.f32"), np.float64 if f64 else np.float32)
        self.out = self.out.astype(np.float64).reshape(-1, len(cols))
        self.cols = cols

    def mz(self, f, g):
        return self.g[f, g, 2] / self.gs * self.g[f, g, 3] / self.gs

    def ml(self, f, g):
        return self.g[f, g, 3] / self.gs

    def grain_at(self, x, y):
        if not self.gm:   # no grain map: the first of VAMPIRE's Voronoi outlines that holds the point (even-odd rule)
            for g, o in enumerate(self.e["grains"].get("outline") or []):
                inside = False
                for (x0, y0), (x1, y1) in zip(o, o[1:] + o[:1]):
                    if (y0 > y) != (y1 > y) and x < x0 + (y - y0) * (x1 - x0) / (y1 - y0):
                        inside = not inside
                if o and inside:
                    return g
            return -1
        a = self.gm.get("atom_pitch_nm") or self.gm["pitch_nm"]
        f = round(self.gm["pitch_nm"] / a)
        ix, iy = int(np.rint(x / a)) // f, int(np.rint(y / a)) // f
        return int(self.gmap[iy, ix]) if 0 <= ix < self.gm["nx"] and 0 <= iy < self.gm["ny"] else -1

    def column(self, f, x, y):
        """(stored column m_z or nan, atoms, lowest |m|) at the snapshot nearest frame f."""
        G = self.grid
        idx = min(self.snaps, key=lambda s: abs(s[1] - self.steps[f]))[0]
        tag = self.e["snapshots"].get("dtype", "i8")
        dtype, sc, empty = {"i8": (np.int8, 127, -128), "i16": (np.int16, 32767, -32768), "f32": (np.float32, 1, None)}[tag]
        a = np.frombuffer(gzip.open(os.path.join(self.d, "top", f"{idx:08d}.{tag}.gz")).read(), dtype).reshape(3, G["ny"], G["nx"])
        p = G["atom_pitch_nm"]
        k = round(G["pitch_nm"] / p)
        ix, iy = int(np.rint(x / p)) // k, int(np.rint(y / p)) // k
        v = a[2, iy, ix]
        v = math.nan if (empty is not None and v == empty) or v != v else v / sc
        return v, int(self.count[iy, ix]), self.minm[iy, ix] / 255

    def cells(self, lag):
        """Bit cells as tests/test_viewer_smoke.py Expected.cells works them out from the bundle, cached."""
        if getattr(self, "_cells", (None,))[0] != lag:
            self._cells = (lag, Expected(os.path.dirname(self.d), self.e).cells(lag))
        return self._cells[1]

    def scores(self, lag):
        """Per grain (bit, intended, correct) at a freeze lag: bit floor((x + lag) / BL) of its track."""
        e, p = self.e, self.e["hamr"]
        seq, n, res = p["sequence"], p["bits_per_track"], []
        for a, x in zip(e["analysis"]["grains"], e["grains"]["x"]):
            b = math.floor((x + lag) / p["BL"])
            k = a["track"] * n + b
            want = int(np.sign((seq[k] if 0 <= b < n and k < len(seq) else 0) * p["field_uz"]))
            ok = a["category"] not in ("unheated", "undecided", "empty") and want != 0 and np.sign(a["final"]) == want
            res.append((b, want, bool(ok)))
        return res


def melted(a):
    return a["category"] not in ("unheated", "empty")


@unittest.skipUnless(have_playwright() or os.environ.get("HAMRVIS_REQUIRE_VIEWER_TEST"), "needs node and Playwright")
class ViewerReadouts(unittest.TestCase):

    def check_bundle(self, bundle):
        M = json.load(open(os.path.join(bundle, "manifest.json"), encoding="utf-8"))
        runs = M.get("runs") or [dict(e, id=k) for k, e in M["stages"].items()]
        view, url = serve(bundle)
        try:
            r = subprocess.run(["node", "viewer_readouts.js", url] + [e["id"] for e in runs], cwd=HERE, capture_output=True, text=True)
        finally:
            view.kill()
            view.wait()
        self.assertEqual(r.returncode, 0, r.stderr)
        got = json.loads(r.stdout)
        byid = {e["id"]: e for e in runs}
        for e in runs:
            if not M.get("runs"):   # version 1: read compares with write
                e.setdefault("ref", "write" if e["id"] == "read" and "write" in byid else None)
            R = Run(bundle, e)
            v = got[e["id"]]
            bad = []
            self.check_frames(R, byid, v, bad)
            self.check_charts(R, v, bad)
            with self.subTest(bundle=os.path.basename(bundle), run=e["id"]):
                self.assertEqual(bad[:12], [], f"{len(bad)} readouts differ")
        return len(runs)

    # ---------------------------------------------------------------- printed values
    def check_frames(self, R, byid, v, bad):
        e, h = R.e, R.head
        lag = v["lag"] or 0.0
        sc = R.scores(lag) if e.get("analysis") and h else None
        for fr in v["frames"]:
            f, step = fr["f"], R.steps[fr["f"]]
            clock = dict(fr["clock"])
            if R.nF > 1:
                a, b = clock["t"].split(" / ")
                if not (printed_time(a, R.time[f], R.dt, R.mc) and printed_time(b, R.time[-1], R.dt, R.mc)):
                    bad.append(("clock t", f, clock["t"], R.time[f]))
            if h:
                hx, _, on, _, _ = h.at(step)
                if (clock["laser x"] == "off") == on or (on and not printed(clock["laser x"], hx)):
                    bad.append(("clock laser x", f, clock["laser x"], hx, on))
                bt = h.bit(step)
                want_bit = f"{bt[1] + 1} / {e['hamr']['bits_per_track']}" if bt else "–"
                if clock["bit"] != want_bit:
                    bad.append(("clock bit", f, clock["bit"], want_bit))
                Bz = h.field(step) * e["hamr"]["field_uz"]
                if not (clock["Bz head"] == "0 T" if Bz == 0 else printed(clock["Bz head"], Bz)):
                    bad.append(("clock Bz", f, clock["Bz head"], Bz))
            for view, tips in fr["tips"].items():
                for x, y, rows in tips:
                    self.check_tip(R, sc, lag, f, view, x, y, rows, bad)
            for x, rows in fr.get("bits", []):
                self.check_bits(R, sc, lag, f, x, dict(rows), bad)
            for x, y, rows in fr.get("bmap", []):
                rows, P = dict(rows), e["dipole"]
                cx, cy = (P.get("cell_nm_xyz") or [P["cell_nm"]] * 3)[:2]
                ix, iy, z = math.floor(x / cx), math.floor(y / cy), fr["wave"]["z"]
                for c, k in enumerate(("Bx", "By", "Bz")):
                    if not printed(rows[k], R.B[c, z, iy, ix] * 1e3):
                        bad.append(("stray field tip", k, x, y, rows[k], R.B[c, z, iy, ix] * 1e3))
            if fr.get("wave"):
                W, P = fr["wave"], e["dipole"]
                cy = (P.get("cell_nm_xyz") or [P["cell_nm"]] * 3)[1]
                rows_in = [iy for iy in range(P["ny"] - 1) if abs((iy + .5) * cy - W["ysec"]) <= W["rw"] / 2]
                want = R.B[2, W["z"], rows_in, :P["nx"] - 1].sum(0) / max(len(rows_in), 1)
                if len(W["v"]) != P["nx"] - 1 or np.abs(np.array(W["v"]) - want).max() > 1e-12:
                    bad.append(("read-back waveform", f))

    def check_tip(self, R, sc, lag, f, view, x, y, rows, bad):
        e, h, got = R.e, R.head, dict(rows)
        want = {}
        xs = nums(got.get("x, y", ""))
        if len(xs) != 2 or abs(xs[0][0] - x) > xs[0][1] + 1e-9 or abs(xs[1][0] - y) > xs[1][1] + 1e-9:
            bad.append(("tip x, y", view, x, y, got.get("x, y")))
        if view == "atoms":
            mzc, n, low = R.column(f, x, y)
            want["column mz"] = ("–" if mzc != mzc else mzc)
            want["atoms in column"] = n
            want["lowest column |m|"] = low
        if h:
            want["laser T now"] = h.T(x, y, R.steps[f])
            want["peak laser T"] = h.peak_T(x, y)
        g = R.grain_at(x, y)
        if view == "cells" and sc is not None:
            p = e["hamr"]
            t = min(p["n_tracks"] - 1, max(0, math.floor((y - p["pad"]) / p["TW"])))
            b = math.floor((x + lag) / p["BL"])
            if 0 <= b < p["bits_per_track"]:
                _, _, w, parts = R.cells(lag)[t * p["bits_per_track"] + b]
                m = sum(c * R.mz(f, k) for k, c in parts.items()) / sum(parts.values()) if parts else 0.0
                want[f"bit {b + 1} mean mz"] = (m, "" if not w else "✓" if np.sign(m) == w and parts else "✗")
        if g >= 0:
            want["grain"] = e["grains"]["id"][g]
            if R.g is not None:
                want["mz"], want["|m|"] = R.mz(f, g), R.ml(f, g)
            if sc is not None:
                a, (_, w, ok) = e["analysis"]["grains"][g], sc[g]
                lab = CAT[a["category"]]
                want["outcome"] = lab if not melted(a) or a["category"] == "undecided" or not w else ("✓ " if ok else "✗ ") + lab
            elif e.get("analysis"):
                want["outcome"] = None          # analysis without a head model (not expected); compared loosely
        missing = set(want) ^ (set(got) - {"x, y"})
        if missing:
            bad.append(("tip rows", view, f, x, y, sorted(missing)))
            return
        for k, w in want.items():
            s = got[k]
            if isinstance(w, tuple):
                ok = printed(s, w[0]) and s.endswith(w[1])
            elif isinstance(w, str):
                ok = s == w
            elif w is None:
                ok = True
            elif k == "grain":
                ok = s == str(w)
            else:
                ok = printed(s, w)
            if not ok:
                bad.append(("tip", view, f, x, y, k, s, w))

    def check_bits(self, R, sc, lag, f, x, got, bad):
        e, p = R.e, R.e["hamr"]
        b = math.floor((x + lag) / p["BL"])
        if not 0 <= b < p["bits_per_track"]:
            if got:
                bad.append(("bits tip outside the bits", x, got))
            return
        an = e["analysis"]["grains"]
        gs = [g for g, a in enumerate(an) if melted(a) and sc[g][1] and a["track"] == 0 and sc[g][0] == b]
        w = int(np.sign(p["sequence"][b] * p["field_uz"]))
        written = f == R.nF - 1 or any(an[g].get("decision_step") is not None and an[g]["decision_step"] <= R.steps[f] for g in gs)
        now = sum(np.sign(R.mz(f, g)) == sc[g][1] for g in gs)
        end = sum(sc[g][2] for g in gs)
        want = {"bit": f"{b + 1} / {p['bits_per_track']}", "intended": "up (+z)" if w > 0 else "down (−z)" if w < 0 else "no field"}
        at = next((k for k in got if k.startswith("✓ at ") and k != "✓ at end"), None)
        lo, hi = (n for n, _ in nums(got.get("x", "")))
        checks = [got.get("bit") == want["bit"], got.get("intended") == want["intended"],
                  abs(lo - (b * p["BL"] - lag)) <= 0.05 + 1e-9, abs(hi - ((b + 1) * p["BL"] - lag)) <= 0.05 + 1e-9,
                  ("written" in got) == (not written), at is not None and got[at] == f"{now} / {len(gs)}",
                  got.get("✓ at end") == f"{end} / {len(gs)}"]
        if not all(checks):
            bad.append(("bits tip", f, x, got, want, now, end, len(gs), checks))

    # ---------------------------------------------------------------- chart series
    def check_charts(self, R, v, bad):
        e, h = R.e, R.head
        S = v["scrub"]
        if S:
            cols = R.cols
            ci = lambda n: next((i for i, c in enumerate(cols) if c.endswith(n)), -1)
            iz, il, tc = ci("M_norm_z"), ci("M_norm_l"), cols.index("Real_time") if "Real_time" in cols else -1
            rows, t1 = len(R.out), (R.time[-1] or 1)
            tOf = (lambda r: R.out[r, tc]) if tc >= 0 else (lambda r: r / rows * t1)
            val = {"M_l": lambda r: R.out[r, il], "M_z": lambda r: R.out[r, iz] * (R.out[r, il] if il >= 0 else 1)}
            for s in S["series"]:
                fn = val[s["key"]]
                if S["mode"] == "rows":
                    step = max(1, math.floor(rows / (S["w"] * 2)))
                    want = [[r, tOf(r), fn(r)] for r in range(0, rows, step)]
                    if len(want) != len(s["pts"]) or any(a[0] != b[0] or abs(a[1] - b[1]) > 1e-24 or abs(a[2] - b[2]) > 1e-12 for a, b in zip(want, s["pts"])):
                        bad.append(("scrubber", s["key"], "rows"))
                else:
                    n = math.ceil(S["w"]) + 1
                    lo, hi = np.full(n, np.inf), np.full(n, -np.inf)
                    for r in range(rows):
                        k = max(0, min(n - 1, math.floor(tOf(r) / t1 * S["w"])))
                        lo[k], hi[k] = min(lo[k], fn(r)), max(hi[k], fn(r))
                    want = [[k, lo[k], hi[k]] for k in range(n) if hi[k] >= lo[k]]
                    if len(want) != len(s["pts"]) or any(a[0] != b[0] or abs(a[1] - b[1]) > 1e-12 or abs(a[2] - b[2]) > 1e-12 for a, b in zip(want, s["pts"])):
                        bad.append(("scrubber", s["key"], "min/max per pixel"))
        for st in v["story"]:
            g, gx, gy = st["g"], e["grains"]["x"][st["g"]], e["grains"]["y"][st["g"]]
            f0, f1 = 0, R.nF - 1
            if h and not st["all"]:
                best, bf = 0, 0
                for f in range(0, R.nF, 2):
                    T = h.T(gx, gy, R.steps[f])
                    if T > best:
                        best, bf = T, f
                p = e["hamr"]
                half = math.floor(3 * p["sigma_x"] / p["v"] * 1e-9 / ((R.time[1] - R.time[0]) or 1e-12) + 0.5)
                f0, f1 = max(0, bf - half), min(R.nF - 1, bf + 2 * half)
            pst = max(1, (f1 - f0) // 300)
            if (st["f0"], st["f1"], st["pst"]) != (f0, f1, pst):
                bad.append(("life story window", g, st["all"], (st["f0"], st["f1"], st["pst"]), (f0, f1, pst)))
                continue
            fs = range(f0, f1 + 1, pst)
            fn = {"T": lambda f: h.T(gx, gy, R.steps[f]), "m": lambda f: R.ml(f, g), "mz": lambda f: R.mz(f, g),
                  "Bz": lambda f: h.field(R.steps[f]) * e["hamr"]["field_uz"] if h.in_box(gx, gy, R.steps[f]) else 0.0}
            keys = (["T"] if h else []) + ["m", "mz"] + (["Bz"] if h else [])
            if [P["key"] for P in st["panels"]] != keys:
                bad.append(("life story panels", g, [P["key"] for P in st["panels"]]))
                continue
            for P in st["panels"]:
                if [q[0] for q in P["pts"]] != list(fs) or any(abs(q[1] - fn[P["key"]](q[0])) > 1e-9 * max(1, abs(q[1])) for q in P["pts"]):
                    bad.append(("life story", g, st["all"], P["key"]))
        if v["st"]:
            an = e["analysis"]["grains"] if e.get("analysis") else None
            rows = [g for g in range(len(e["grains"]["id"])) if not an or melted(an[g])]
            rows.sort(key=lambda g: (an[g]["track"] if an else 0, e["grains"]["x"][g]))
            if v["st"]["rows"] != rows:
                bad.append(("space-time rows", v["st"]["rows"][:8], rows[:8]))
            elif h:
                if v["st"]["pass"] != [h.pass_step(e["grains"]["x"][g], an[g]["track"] if an else 0) for g in rows]:
                    bad.append(("space-time laser passes",))
                if an and v["st"]["dec"] != [an[g].get("decision_step") for g in rows]:
                    bad.append(("space-time decision ticks",))
        if v["lagscan"]:
            p, inf = e["hamr"], e.get("freeze_lag_nm") or 0
            xmax = max(3 * p["BL"], 2 * inf, (v["lag"] or 0) * 1.2, 1)
            an = e["analysis"]["grains"]
            for i, (l, ok, n, out) in enumerate(v["lagscan"]):
                sc = R.scores(i / 150 * xmax)
                m = [g for g, a in enumerate(an) if melted(a)]
                w = [g for g in m if sc[g][1]]
                want = (sum(sc[g][2] for g in w), len(w), len(m) - len(w))
                if abs(l - i / 150 * xmax) > 1e-12 or (ok, n, out) != want:
                    bad.append(("score against lag", l, (ok, n, out), want))
                    break

    # ---------------------------------------------------------------- bundles
    def test_scenarios(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name in ("chain", "two-tracks", "rewrite", "continuous", "llg-relax"):
                out = os.path.join(tmp, name)
                reduce.reduce([os.path.join(SCEN, name)], out, workers=1)
                self.check_bundle(out)

    def test_real_bundles(self):
        found = [os.path.join(ROOT, "bundles", b) for b in ("h3", "h3-v1", "2x11", "h3-scale-window", "h3-scale-full") if os.path.exists(os.path.join(ROOT, "bundles", b, "manifest.json"))]
        if not found:
            self.skipTest("real-run bundles in bundles/ (git-ignored) not present")
        for b in found:
            self.check_bundle(b)


if __name__ == "__main__":
    unittest.main()

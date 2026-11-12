"""The pixels the viewer draws hold the data under them: sampled canvas pixels of the film (grains, atom
columns, one layer, in every colour mode), the x-z section and the stray field map must have the colour the colour map
gives the bundle's value at that film position, and a selected grain's outline (the pixels that change when it is selected) must lie on its grain-map boundary.
This catches drawing at the wrong place (flipped y, an offset of a pixel or more, the wrong grid), the wrong layer and
the wrong frame. The film position of a screen pixel is worked out here from the lane size and the x window alone.

Tolerances: one screen pixel of placement (the browser's rounding where a data pixel edge falls inside a screen pixel),
one colour map entry (1/127.5 in m_z, rounding of values on an entry edge). Left out: the section line, the head
(isotherms, field box and its label) while the laser is on, the snapshot time label, the reader band. Skipped without
node and Playwright."""
import gzip
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # also when run as tests.<module>
from test_viewer_readouts import Run, melted
from test_viewer_smoke import HERE, ROOT, SCEN, have_playwright, reduce, serve

OFFSETS = ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1))     # screen pixels around the sampled one
VOID, MID = -1, -2                                       # expected colour: no data, or the neutral colour


def idx_div(v):
    return np.clip(np.floor((np.asarray(v, float) + 1) * 127.5 + 0.5), 0, 255).astype(int)


def idx_seq(v):
    return np.clip(np.floor(np.asarray(v, float) * 255 + 0.5), 0, 255).astype(int)


class Film:
    """Screen pixel -> film position for a lane on the film axis: x from the x window, y from the film height."""

    def __init__(self, cap, L):
        self.W, self.H, self.Lx, self.Ly = cap["W"], cap["H"], L[0], L[1]
        x0, x1 = cap["win"]["x0"], cap["win"]["x1"]
        self.x0, self.sx = x0, cap["W"] / (x1 - x0)
        full = self.Ly * self.sx <= self.H + 1
        self.sy = self.H / self.Ly if full else self.sx
        self.y0 = 0.0 if full else max(0.0, min(self.Ly - self.H / self.sy, cap["y0"]))

    def xy(self, q, r):
        return self.x0 + (q + 0.5) / self.sx, self.y0 + (self.H - (r + 0.5)) / self.sy

    def row(self, y):
        return self.H - (y - self.y0) * self.sy


class Expect:
    """Expected colour (a colour map and entry, VOID or MID) of film positions, per view and colour mode."""

    def __init__(self, R, ref, luts):
        self.R, self.ref, self.luts = R, ref, {k: np.array(v, int).reshape(-1, 3) for k, v in luts.items() if k not in ("void", "mid")}
        self.void, self.mid = np.array(luts["void"]), np.array(luts["mid"])
        e = R.e
        self.lag = (e.get("freeze_lag_nm") or 0) if e.get("hamr") else None
        W = R if e.get("analysis") and e.get("hamr") else ref
        self.want = None
        if W is not None and W.e.get("analysis") and W.e.get("hamr"):
            sc = W.scores(W.e.get("freeze_lag_nm") or 0)
            self.want = np.array([s[1] if melted(a) else 0 for a, s in zip(W.e["analysis"]["grains"], sc)])

    def grains(self, x, y):
        """Grain index at film positions (-1: none): the grain map, else the last Voronoi outline drawn over the point."""
        R = self.R
        if R.gm:
            a = R.gm.get("atom_pitch_nm") or R.gm["pitch_nm"]
            f = round(R.gm["pitch_nm"] / a)
            ix, iy = np.rint(x / a).astype(int) // f, np.rint(y / a).astype(int) // f
            ok = (ix >= 0) & (ix < R.gm["nx"]) & (iy >= 0) & (iy < R.gm["ny"])
            return np.where(ok, R.gmap[np.clip(iy, 0, R.gm["ny"] - 1), np.clip(ix, 0, R.gm["nx"] - 1)], -1)
        g = np.full(x.shape, -1)
        for k, o in enumerate(R.e["grains"].get("outline") or []):
            if not o:
                continue
            o = np.array(o)
            inside = np.zeros(x.shape, bool)
            for (xa, ya), (xb, yb) in zip(o, np.roll(o, -1, axis=0)):
                with np.errstate(divide="ignore", invalid="ignore"):
                    inside ^= ((ya > y) != (yb > y)) & (x < xa + (y - ya) * (xb - xa) / (yb - ya))
            g[inside] = k
        return g

    def snapshot(self, f, kind):
        R = self.R
        idx = min(R.snaps, key=lambda s: abs(s[1] - R.steps[f]))[0]
        tag = R.e["snapshots"].get("dtype", "i8")
        dtype, sc, empty = {"i8": (np.int8, 127, -128), "i16": (np.int16, 32767, -32768), "f32": (np.float32, 1, None)}[tag]
        a = np.frombuffer(gzip.open(os.path.join(R.d, kind, f"{idx:08d}.{tag}.gz")).read(), dtype).astype(np.float64)
        if empty is not None:
            a[a == empty] = np.nan
        return a / sc

    def cols(self, x, y, pitch):
        """Column indices (and validity) of film positions on a grid of the given pitch (atom grids start at -a/2)."""
        G = self.R.grid
        a = G["atom_pitch_nm"]
        k = round(pitch / a)
        return np.rint(x / a).astype(int) // k, np.rint(y / a).astype(int) // k, a, k

    def colour(self, cap, x, y):
        """(map name or None, entry or VOID/MID) arrays for film positions."""
        R, e, view, mode, f = self.R, self.R.e, cap["view"], cap.get("colour"), cap["f"]
        if view == "grains":
            g = self.grains(x, y)
            gg = np.clip(g, 0, None)
            if mode == "mz":
                v = np.array([R.mz(f, k) for k in range(len(e["grains"]["id"]))])
                return "m", np.where(g >= 0, idx_div(v[gg]), VOID)
            if mode == "agree":
                v = np.array([R.mz(f, k) for k in range(len(e["grains"]["id"]))])
                w = self.want[gg]
                return "agree", np.where(g < 0, VOID, np.where(w != 0, idx_div(v[gg] * w), MID))
            if mode == "peak":
                p = e["hamr"]
                pk = np.array([a["peakT"] for a in e["analysis"]["grains"]])
                return "s", np.where(g >= 0, idx_seq((pk[gg] - p["Tmin"]) / (p["Tmax"] - p["Tmin"])), VOID)
            if mode == "minm":
                low = np.minimum(1.0, R.g[:, :, 3].min(0) / R.gs)
                return "s", np.where(g >= 0, idx_seq(1 - low[gg]), VOID)
        G = R.grid
        if view == "atoms":
            ix, iy, a, k = self.cols(x, y, G["pitch_nm"])
            ok = (ix >= 0) & (ix < G["nx"]) & (iy >= 0) & (iy < G["ny"])
            ix, iy = np.clip(ix, 0, G["nx"] - 1), np.clip(iy, 0, G["ny"] - 1)
            filled = ok & (R.count[iy, ix] > 0)
            if mode == "mz":
                v = self.snapshot(f, "top").reshape(3, G["ny"], G["nx"])[2, iy, ix]
                return "m", np.where(ok & (v == v), idx_div(np.nan_to_num(v)), VOID)
            if mode == "peak":
                p, c = e["hamr"], lambda i: a * (i * k + (k - 1) / 2)
                T = np.array([R.head.peak_T(cx, cy) for cx, cy in zip(c(ix), c(iy))])
                return "s", np.where(filled, idx_seq((T - p["Tmin"]) / (p["Tmax"] - p["Tmin"])), VOID)
            if mode == "minm":
                return "s", np.where(filled, idx_seq(1 - R.minm[iy, ix] / 255), VOID)
        if view in ("layer", "xz"):
            lp = G.get("layer_pitch_nm") or 2 * G["pitch_nm"]
            vol = self.snapshot(f, "vol").reshape(3, G["nz"], G["vy"], G["vx"])[2]
            ix, iy, a, k = self.cols(x, y, lp)
            if view == "xz":   # x along the film, y here is the layer, the row is the one holding the section line
                iz, iy = y, np.full(x.shape, math.floor((cap["ysec"] + a / 2) / lp))
            else:
                iz = np.full(x.shape, cap["layer"])
            ok = (ix >= 0) & (ix < G["vx"]) & (iy >= 0) & (iy < G["vy"]) & (iz >= 0) & (iz < G["nz"])
            v = vol[np.clip(iz, 0, G["nz"] - 1), np.clip(iy, 0, G["vy"] - 1), np.clip(ix, 0, G["vx"] - 1)]
            return "m", np.where(ok & (v == v), idx_div(np.nan_to_num(v)), VOID)
        raise AssertionError(f"no expectation for {view}/{mode}")

    def matches(self, px, lut, ent):
        """True where a screen colour is the expected one, within one colour map entry."""
        good = np.zeros(len(px), bool)
        good |= (ent == VOID) & (np.abs(px - self.void).max(1) <= 1)
        good |= (ent == MID) & (np.abs(px - self.mid).max(1) <= 1)
        if lut:
            L = self.luts[lut]
            for dk in (-1, 0, 1):
                k = np.clip(ent + dk, 0, 255)
                good |= (ent >= 0) & (np.abs(px - L[k]).max(1) <= 1)
        return good


def isotherm_reach(e):
    """Largest distance from the laser centre (in sigma) of the isotherms the viewer draws."""
    p, Tc = e["hamr"], e.get("Tc")
    Ts = ([Tc] if Tc and p["Tmin"] < Tc < p["Tmax"] else []) + [T for T in (600, 500, 400) if p["Tmin"] < T < p["Tmax"] and not (Tc and abs(T - Tc) <= 15)]
    return max((math.sqrt(2 * math.log((p["Tmax"] - p["Tmin"]) / (T - p["Tmin"]))) for T in Ts), default=0)


@unittest.skipUnless(have_playwright() or os.environ.get("HAMRVIS_REQUIRE_VIEWER_TEST"), "needs node and Playwright")
class ViewerPixels(unittest.TestCase):

    def check_bundle(self, bundle):
        M = json.load(open(os.path.join(bundle, "manifest.json"), encoding="utf-8"))
        runs = M.get("runs") or [dict(e, id=k, ref="write" if k == "read" and "write" in M["stages"] else None) for k, e in M["stages"].items()]
        view, url = serve(bundle)
        try:
            r = subprocess.run(["node", "viewer_pixels.js", url] + [e["id"] for e in runs], cwd=HERE, capture_output=True, text=True)
        finally:
            view.kill()
            view.wait()
        self.assertEqual(r.returncode, 0, r.stderr)
        got = json.loads(r.stdout)
        Rs = {e["id"]: Run(bundle, e) for e in runs}
        n = 0
        for e in runs:
            R, v = Rs[e["id"]], got[e["id"]]
            X = Expect(R, Rs.get(e.get("ref")), v["luts"])
            for cap in v["captures"]:
                with self.subTest(bundle=os.path.basename(bundle), run=e["id"], lane=cap["lane"], view=cap["view"], colour=cap.get("colour"), frame=cap["f"]):
                    checked, bad = self.check_capture(R, X, cap)
                    self.assertGreater(checked, 1000, "too few pixels left to check")
                    self.assertLessEqual(bad, 0, f"{bad} of {checked} pixels do not hold the data under them")
                    n += checked
            if v.get("outline"):
                with self.subTest(bundle=os.path.basename(bundle), run=e["id"], check="outline"):
                    self.check_outline(R, v["outline"])
        return n

    def samples(self, cap):
        st, W, H = cap["stride"], cap["W"], cap["H"]
        q, r = np.meshgrid(np.arange(0, W, st), np.arange(0, H, st))
        px = np.array(cap["px"], int).reshape(-1, 3)
        return q.ravel(), r.ravel(), px

    def check_capture(self, R, X, cap):
        q, r, px = self.samples(cap)
        keep = np.ones(len(q), bool)
        if cap["lane"] == "bmap":
            return self.check_field(R, X, cap, q, r, px)
        if cap["lane"] == "xz":
            G = R.grid
            F = Film(cap, R.e["size_nm"])
            x = F.xy(q, r)[0]
            iz = np.floor((cap["H"] - (r + 0.5)) / cap["H"] * G["nz"]).astype(int)
            ok = np.zeros(len(q), bool)
            for dq, dr in OFFSETS:
                xx = F.xy(q + dq, r)[0]
                izz = np.floor((cap["H"] - (r + dr + 0.5)) / cap["H"] * G["nz"]).astype(int)
                lut, ent = X.colour(cap, xx, izz)
                ok |= X.matches(px, lut, ent)
            return len(q), int((~ok).sum())
        F = Film(cap, R.e["size_nm"])
        keep &= np.abs(r + 0.5 - F.row(cap["ysec"])) > 3                     # section line
        if R.e.get("dipole"):                                                  # reader band
            keep &= np.abs(r + 0.5 - F.row(cap["ysec"])) > cap["rw"] / 2 * F.sy + 3
        if cap["view"] != "grains":                                            # snapshot time label, bottom left
            keep &= ~((q < 260) & (r > cap["H"] - 30))
        if R.head:
            hx, _, on, _, _ = R.head.at(R.steps[cap["f"]])
            if on:
                p, reach = R.e["hamr"], isotherm_reach(R.e)
                x = F.xy(q, r)[0]
                near_laser = np.abs(x - hx) <= reach * p["sigma_x"] + 1
                near_box = np.abs(x - (hx - p["NPS"])) <= p["box_x"] / 2 + 1
                keep &= ~(near_laser | near_box)
        q, r, px = q[keep], r[keep], px[keep]
        ok = np.zeros(len(q), bool)
        ents = []
        for dq, dr in OFFSETS:
            x, y = F.xy(q + dq, r + dr)
            lut, ent = X.colour(cap, x, y)
            ents.append(ent)
            ok |= X.matches(px, lut, ent)
        if not R.gm and cap["view"] == "grains":   # Voronoi polygons are drawn smoothed: their edges blend two grains
            edge = np.any([en != ents[0] for en in ents[1:]], axis=0)
            ok |= edge
        return len(q), int((~ok).sum())

    def check_field(self, R, X, cap, q, r, px):
        """Stray B_z map: each drawn cell is one image pixel in the colour of its value, and the browser draws the image
        smoothed, so a screen pixel holds the bilinear blend of the four nearest cell colours (texel centres at the cell
        centres). Compared within 3 per colour channel plus 1/8 of the four cells' colour range (the browser's
        interpolation weights are coarse; a cell drawn half a cell or more out of place differs by about half the range); pixels near grain boundaries (darkened over the map), the
        reader band and the outermost half cell are left out."""
        P = R.e["dipole"]
        cx, cy = (P.get("cell_nm_xyz") or [P["cell_nm"]] * 3)[:2]
        Bz = R.B[2, cap["z"]]
        a = np.sort(np.abs(Bz[np.isfinite(Bz)]).ravel())
        sc = a[int(len(a) * 0.995)] or 1e-3
        C = X.luts["b"][idx_div(Bz[:P["ny"] - 1, :P["nx"] - 1] / sc)].astype(float)      # the drawn cells
        F = Film(cap, R.e["size_nm"])
        x, y = F.xy(q, r)
        u, v = x / cx - 0.5, y / cy - 0.5
        keep = (u >= 0) & (u <= P["nx"] - 2) & (v >= 0) & (v <= P["ny"] - 2)
        keep &= np.abs(r + 0.5 - F.row(cap["ysec"])) > cap["rw"] / 2 * F.sy + 3
        if R.gm:
            for dq in (-2, 0, 2):
                for dr in (-2, 0, 2):
                    keep &= X.grains(*F.xy(q + dq, r + dr)) >= 0
        u, v, px = u[keep], v[keep], px[keep]
        i0, j0 = np.minimum(np.floor(u).astype(int), P["nx"] - 3), np.minimum(np.floor(v).astype(int), P["ny"] - 3)
        tu, tv = (u - i0)[:, None], (v - j0)[:, None]
        want = (C[j0, i0] * (1 - tu) * (1 - tv) + C[j0, i0 + 1] * tu * (1 - tv) + C[j0 + 1, i0] * (1 - tu) * tv
                + C[j0 + 1, i0 + 1] * tu * tv)
        quad = np.stack([C[j0, i0], C[j0, i0 + 1], C[j0 + 1, i0], C[j0 + 1, i0 + 1]])
        tol = 3 + (quad.max(0) - quad.min(0)).max(1) / 8
        bad = np.abs(px - want).max(1) > tol
        return len(px), int(bad.sum())

    def check_outline(self, R, o):
        """Every pixel that changes when a grain is selected lies on its grain-map boundary (within one grain-map pixel,
        as the outline is drawn 1.2 pixels wide), and every boundary pixel is drawn."""
        g, gm = o["g"], R.gmap
        ny, nx = gm.shape
        mine = gm == g
        edge = mine & ~(np.roll(mine, 1, 0) & np.roll(mine, -1, 0) & np.roll(mine, 1, 1) & np.roll(mine, -1, 1))
        F = Film(dict(o, W=o["W"], H=o["H"]), R.e["size_nm"])
        pts = np.array(o["pts"], int).reshape(-1, 2)
        self.assertGreater(len(pts), 20, "no outline drawn")
        a = R.gm.get("atom_pitch_nm") or R.gm["pitch_nm"]
        f = round(R.gm["pitch_nm"] / a)
        x, y = F.xy(pts[:, 0], pts[:, 1])
        ix, iy = np.rint(x / a).astype(int) // f, np.rint(y / a).astype(int) // f
        near = np.zeros(len(pts), bool)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                near |= edge[(iy + dy) % ny, (ix + dx) % nx]
        self.assertEqual(int((~near).sum()), 0, f"{int((~near).sum())} of {len(pts)} outline pixels are not on grain {g}'s boundary")
        ink = np.zeros((o["H"], o["W"]), bool)
        ink[pts[:, 1], pts[:, 0]] = True
        ey, ex = np.nonzero(edge)
        cq = np.floor((a * (ex * f + (f - 1) / 2) - F.x0) * F.sx).astype(int)
        cr = np.floor(F.row(a * (ey * f + (f - 1) / 2))).astype(int)
        vis = (cq >= 1) & (cq < o["W"] - 1) & (cr >= 1) & (cr < o["H"] - 1) & (np.abs(cr + 0.5 - F.row(o["ysec"])) > 3)
        if o.get("dipole"):
            vis &= np.abs(cr + 0.5 - F.row(o["ysec"])) > o["rw"] / 2 * F.sy + 3
        drawn = np.zeros(len(cq), bool)
        for dq in (-1, 0, 1):
            for dr in (-1, 0, 1):
                drawn |= ink[np.clip(cr + dr, 0, o["H"] - 1), np.clip(cq + dq, 0, o["W"] - 1)]
        missing = int((vis & ~drawn).sum())
        self.assertEqual(missing, 0, f"{missing} of {int(vis.sum())} boundary pixels of grain {g} have no outline")

    def test_scenarios(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name in ("chain", "two-tracks", "continuous"):
                out = os.path.join(tmp, name)
                reduce.reduce([os.path.join(SCEN, name)], out, workers=1)
                self.check_bundle(out)

    def test_real_bundles(self):
        found = [os.path.join(ROOT, "bundles", b) for b in ("h3", "2x11", "h3-scale-window", "h3-scale-full") if os.path.exists(os.path.join(ROOT, "bundles", b, "manifest.json"))]
        if not found:
            self.skipTest("real-run bundles in bundles/ (git-ignored) not present")
        for b in found:
            self.check_bundle(b)


if __name__ == "__main__":
    unittest.main()

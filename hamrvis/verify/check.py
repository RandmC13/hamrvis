"""hamrvis verify: every array and manifest number of a bundle against the run directories it was made from, read by
hamrvis/verify/reference.py, which shares no code with the rest of hamrvis. The only allowed differences are the
stored encodings and roundings:
    output table          exact (float64); joined jobs: every stored row is a row of one job, in time order
    grain list            ids and listed atom counts exact
    grain centres         periodic mean of the grain's magnetic atoms, stored to 1e-3 nm; listed centre without any
    magnetic atoms        exact; atom counts per column (top_count) and per layer cell (vol_count), capped at 255, exact
    grid                  atom pitch, columns, layers, layer cells from the atoms exactly as docs/bundle-format.md says
    grain map             the grain with most magnetic atoms in each column, -1 without any, exact
    grain series          int16: within half a step (0.5/32767); float32: float32 of grain.txt's value
    snapshot columns and  int8: 0.5/127, int16: 0.5/32767, float32: float32 rounding; empty pixels marked empty;
    layers                window snapshots are rebuilt from the complete snapshot before them
    lowest |m|            within one step of 1/255 (hamrvis keeps the column |m| in float32)
    dipole cells          float32 of dipole-field.txt's fields, atoms per cell exact
    head parameters       input values in VAMPIRE's units

Run on a bundle next to the run directories it was made from (their paths are the manifest's segments):
    hamrvis verify <bundle> [--run ID ...] [--every N] [--json FILE]
--every N compares every N-th snapshot (lowest |m| is then not checked). The result is stored in the manifest
("verify") and shown by the viewer under Run details. Exit status 0 if everything agrees, 1 if anything differs,
2 if run directories are missing.
"""
import argparse
import gzip
import json
import math
import os
import sys
import time

import numpy as np

from .. import __version__                               # the only thing taken from the rest of hamrvis
from . import reference as ref

STEP = {"i8": 127, "i16": 32767, "f32": None}
DTYPE = {"i8": np.int8, "i16": np.int16, "f32": np.float32}


class Results:
    """Rows (run, quantity, values compared, largest difference, tolerance, ok, note)."""

    def __init__(self, out=None):
        self.rows, self.out = [], out

    def add(self, run, what, n, diff, tol, ok, note=""):
        row = (run, what, int(n), float(diff), tol, bool(ok), note)
        self.rows.append(row)
        if self.out:
            self.out(f"  {'ok  ' if ok else 'FAIL'} {what:<34} {int(n):>12,d} values  max diff {diff:.3g}  (tolerance {tol}){'  ' + note if note else ''}")

    def failed(self):
        return [r for r in self.rows if not r[5]]


def compare(res, run, what, got, want, tol, note=""):
    """Element-wise |got - want| <= tol (NaN must match NaN)."""
    got, want = np.asarray(got, float), np.asarray(want, float)
    if got.shape != want.shape:
        res.add(run, what, 0, math.inf, tol, False, f"shape {got.shape} != {want.shape}")
        return
    nan = np.isnan(want)
    ok_nan = bool((np.isnan(got) == nan).all())
    d = np.abs(got[~nan] - want[~nan])
    diff = float(d.max()) if d.size else 0.0
    res.add(run, what, d.size, diff, tol, ok_nan and diff <= tol + 1e-12, note if ok_nan else "empty pixels differ " + note)


def check_run(bundle, e, res, every=1, log=None):
    """Every check of one run (manifest entry e) of the bundle."""
    rid, d, src, segs = e["id"], os.path.join(bundle, e["id"]), e["source"], e["segments"]
    say = log or (lambda m: None)

    # output table
    cols, parts = None, []
    for seg in segs:
        c, r = ref.output_table(seg)
        cols = cols or c
        parts.append(r)
    stored = np.fromfile(os.path.join(d, "output.f64")).reshape(-1, len(cols))
    if e["output"]["columns"] != cols:
        res.add(rid, "output columns", len(cols), math.inf, "exact", False)
    elif len(parts) == 1:
        compare(res, rid, "output table", stored, parts[0], 0)
    else:
        rows = {tuple(r) for p in parts for r in p}
        t = stored[:, cols.index("Real_time")]
        res.add(rid, "output table (joined jobs)", stored.size, 0, "rows of the jobs, in time order",
                all(tuple(r) in rows for r in stored) and bool((np.diff(t) > 0).all()))

    # head parameters
    if "hamr" in e:
        keys, p = ref.parse_input(ref.input_path(src)), e["hamr"]
        got = [p["v"], p["BL"], p["dt"], p["Lx"], p["Hmax"], p["TW"], p["NPS"]]
        want = [ref.value(keys, "hamr:head-speed", "velocity", 30e10) / 1e10, ref.value(keys, "hamr:bit-size", "length", 0) / 10,
                ref.value(keys, "sim:time-step", "time", 1e-15), ref.value(keys, "dimensions:system-size-x", "length", 100) / 10,
                ref.value(keys, "sim:maximum-applied-field-strength", "field", 1), ref.value(keys, "hamr:track-size", "length", 0) / 10,
                ref.value(keys, "hamr:NPS", "length", 0) / 10]
        compare(res, rid, "head parameters", np.array(got) / np.array(want), np.ones(len(got)), 1e-12, "(relative)")

    G = e["grains"]
    nG = len(G["id"])
    listed = ref.grain_list(src) if os.path.exists(os.path.join(src, "grain-coordinates.txt")) else None
    if listed is not None:
        res.add(rid, "grain list (ids, listed atoms)", 2 * nG, 0, "exact",
                G["id"] == listed[:, 0].astype(int).tolist() and G["atoms"] == listed[:, 3].astype(int).tolist())

    # grain series
    gs = e.get("grain_series")
    if gs and len(segs) == 1:
        say("  reading grain.txt")
        want = ref.grain_magnetisation(src, nG)
        k = gs.get("every", 1)
        if k > 1:
            keep = list(range(0, len(want), k))
            want = want[keep + ([len(want) - 1] if keep[-1] != len(want) - 1 else [])]
        if gs.get("dtype") == "f32":
            got = np.fromfile(os.path.join(d, "grains.f32"), np.float32).reshape(-1, nG, 4)
            compare(res, rid, "grain series (float32)", got, want.astype(np.float32), 0)
        else:
            got = np.fromfile(os.path.join(d, "grains.i16"), np.int16).reshape(-1, nG, 4) / 32767
            compare(res, rid, "grain series (int16)", got, np.clip(want, -1, 1), 0.5 / 32767)

    # atoms: grid, counts, grain map, centres, snapshots
    sn = e.get("snapshots")
    coords = os.path.exists(os.path.join(src, "atoms-coords.meta"))
    if coords and (sn or "grain_map" in e) and not e.get("grain_map", {}).get("from_run"):
        say("  reading the atom coordinates")
        typ, xyz = ref.coordinates(src)
        gid = ref.grain_ids(src)
        mats = ref.magnetic_materials(src)
        res.add(rid, "magnetic materials", 1, 0, "exact", e["magnetic_types"] == mats)
        mag = np.isin(typ, mats) if mats is not None else np.ones(len(typ), bool)
        xyz_m = xyz[mag]
        xs = np.unique(np.round(xyz_m[:, 0], 3))
        dx = np.diff(xs)
        a = float(dx[dx > 2e-3].min())
        zs = np.unique(np.round(xyz_m[:, 2], 3))
        dz = np.diff(zs)
        pz = float(dz[dz > 2e-3].min()) if len(dz) else 1.0
        z0 = float(xyz_m[:, 2].min())
        Lx, Ly = (v * 10 for v in e["size_nm"][:2])
        nx0, ny0 = int(round(Lx / a)), int(round(Ly / a))
        nz = int(round((float(xyz_m[:, 2].max()) - z0) / pz)) + 1
        grid = sn["grid"] if sn else dict(e["grain_map"], pitch_z_nm=None)
        tf = int(round(grid["pitch_nm"] * 10 / a))
        nx, ny = -(-nx0 // tf), -(-ny0 // tf)
        got = [grid["atom_pitch_nm"] * 10, grid["nx"], grid["ny"]]
        want = [a, nx, ny]
        if sn:
            lf = int(round(grid["layer_pitch_nm"] * 10 / a))
            vx, vy = -(-nx0 // lf), -(-ny0 // lf)
            got += [grid["pitch_z_nm"] * 10, grid["z0_nm"] * 10, grid["nz"], grid["vx"], grid["vy"]]
            want += [pz, z0, nz, vx, vy]
        compare(res, rid, "grid (pitch, columns, layers)", got, want, 1e-6,
                "" if tf == max(1, math.ceil(math.sqrt(nx0 * ny0 / 4e6))) else "(top map coarsened by an option)")
        sites = (nx0, ny0)
        ix = np.rint(xyz_m[:, 0] / a).astype(np.int64) % nx0 // tf
        iy = np.rint(xyz_m[:, 1] / a).astype(np.int64) % ny0 // tf
        col = iy * nx + ix
        count = np.bincount(col, minlength=nx * ny)
        if sn:
            got = np.fromfile(os.path.join(d, "top_count.u8"), np.uint8)
            compare(res, rid, "atoms per column", got, np.minimum(count, 255), 0)
            jx = np.rint(xyz_m[:, 0] / a).astype(np.int64) % nx0 // lf
            jy = np.rint(xyz_m[:, 1] / a).astype(np.int64) % ny0 // lf
            jz = np.clip(np.rint((xyz_m[:, 2] - z0) / pz).astype(np.int64), 0, nz - 1)
            vcount = np.bincount((jz * vy + jy) * vx + jx, minlength=nz * vy * vx)
            got = np.fromfile(os.path.join(d, "vol_count.u8"), np.uint8)
            compare(res, rid, "atoms per layer cell", got, np.minimum(vcount, 255), 0)
        if gid is not None and listed is not None:
            g_m = gid[mag]
            index = {int(g): i for i, g in enumerate(G["id"])}
            gi = np.array([index.get(int(g), -1) for g in np.unique(g_m)])[np.searchsorted(np.unique(g_m), g_m)]
            # grain map: most magnetic atoms of a listed grain in the column (atoms of unlisted ids count as none)
            pair, n = np.unique(col * (nG + 1) + np.where(gi < 0, nG, gi), return_counts=True)
            c, g = pair // (nG + 1), pair % (nG + 1)
            best = np.full(nx * ny, -1)
            most = np.zeros(nx * ny, np.int64)
            for ci, gg, nn in zip(c, g, n):
                if nn > most[ci]:
                    most[ci], best[ci] = nn, (-1 if gg == nG else gg)
            got = np.fromfile(os.path.join(d, "grain_map.i16"), np.int16)
            compare(res, rid, "grain map", got, best, 0)
            res.add(rid, "mixed columns", 1, abs(e["grain_map"]["mixed_columns"] - int(len(pair) - len(np.unique(c)))), "exact",
                    e["grain_map"]["mixed_columns"] == int(len(pair) - len(np.unique(c))))
            # centres: periodic mean of each grain's magnetic atoms; grains with none keep their listed centre
            cx = np.zeros(nG); cy = np.zeros(nG); n_at = np.zeros(nG, int)
            th_x, th_y = 2 * np.pi * xyz_m[:, 0] / Lx, 2 * np.pi * xyz_m[:, 1] / Ly
            ok = gi >= 0
            for arr, vals in ((cx, np.cos(th_x)), (cy, np.cos(th_y))):
                arr += np.bincount(gi[ok], vals[ok], minlength=nG)
            sx = np.bincount(gi[ok], np.sin(th_x)[ok], minlength=nG)
            sy = np.bincount(gi[ok], np.sin(th_y)[ok], minlength=nG)
            n_at = np.bincount(gi[ok], minlength=nG)
            mx = (np.arctan2(sx, cx) % (2 * np.pi)) * Lx / (2 * np.pi) / 10
            my = (np.arctan2(sy, cy) % (2 * np.pi)) * Ly / (2 * np.pi) / 10
            mx = np.where(n_at > 0, mx, listed[:, 1] / 10)
            my = np.where(n_at > 0, my, listed[:, 2] / 10)
            compare(res, rid, "grain centres (nm)", np.r_[G["x"], G["y"]], np.r_[mx, my], 5e-4 + 1e-9, "(stored to 1e-3 nm)")
            compare(res, rid, "magnetic atoms per grain", G["magnetic_atoms"], n_at, 0)
        if sn:
            tag = sn.get("dtype", "i8")
            step = STEP[tag]
            listed_snaps = {}
            for seg in segs:                     # a later job's snapshot replaces an earlier one with the same number
                for idx, t, w in ref.snapshot_list(seg):
                    listed_snaps[idx] = (seg, t, w)
            frames = sn["frames"]
            res.add(rid, "snapshot list (numbers, times)", len(frames), 0, "exact",
                    all(f["index"] in listed_snaps and listed_snaps[f["index"]][1] == f["time"] for f in frames))
            left = sn.get("skipped") or []
            res.add(rid, "snapshots left out", len(frames) + len(left), len(left), "none", not left,
                    f"left out by the reduce: {left[:5]}{' ...' if len(left) > 5 else ''}" if left else "")
            check = {f["index"] for f in frames[::every]} | {frames[-1]["index"]}
            tol = 0.5 / step if step else 1e-7
            worst_t = worst_v = 0.0
            n_t = n_v = 0
            bad = []
            low = np.full(nx * ny, np.inf)
            want_snap = {f["index"] for f in frames} if every == 1 else check
            start = time.time()
            atoms = (xyz_m, mag)
            for seg in dict.fromkeys(listed_snaps[f["index"]][0] for f in frames):
                if seg != src and os.path.exists(os.path.join(seg, "atoms-coords.meta")):
                    # a job continued on another number of processes lists its atoms in another order: its own list
                    t2, x2 = ref.coordinates(seg)
                    m2 = np.isin(t2, mats) if mats is not None else np.ones(len(t2), bool)
                    atoms = (x2[m2], m2)
                xyz_m, mag = atoms
                for idx, s in ref.spin_states(seg, want_snap):
                    if listed_snaps[idx][0] != seg or idx not in want_snap:
                        continue
                    s_m = s[mag]
                    top = ref.column_means(xyz_m, s_m, a, nx, ny, tf, sites)
                    if every == 1:
                        with np.errstate(invalid="ignore"):
                            low = np.fmin(low, np.sqrt((top ** 2).sum(0)).ravel())
                    if idx not in check:
                        continue
                    for kind, w in (("top", top), ("vol", ref.layer_means(xyz_m, s_m, a, pz, z0, nz, vx, vy, lf, sites) if sn["volume"] else None)):
                        if w is None:
                            continue
                        q = np.frombuffer(gzip.open(os.path.join(d, kind, f"{idx:08d}.{tag}.gz")).read(), DTYPE[tag]).reshape(w.shape).astype(float)
                        if step:
                            q[q == -step - 1] = np.nan
                            q /= step
                        nan = np.isnan(w)
                        diff = float(np.abs(q[~nan] - w[~nan]).max()) if (~nan).any() else 0.0
                        if not (np.isnan(q) == nan).all() or diff > tol + 1e-12:
                            bad.append((kind, idx, diff))
                        if kind == "top":
                            worst_t, n_t = max(worst_t, diff), n_t + int((~nan).sum())
                        else:
                            worst_v, n_v = max(worst_v, diff), n_v + int((~nan).sum())
                    say(f"  snapshot {idx} ({len(bad)} differing so far, {time.time() - start:.0f} s)")
            res.add(rid, f"snapshot columns ({len(check)} snapshots)", n_t, worst_t, round(tol, 9), not any(b[0] == "top" for b in bad),
                    f"differ: {[b[1] for b in bad if b[0] == 'top'][:5]}" if bad else "")
            if sn["volume"]:
                res.add(rid, f"snapshot layers ({len(check)} snapshots)", n_v, worst_v, round(tol, 9), not any(b[0] == "vol" for b in bad),
                        f"differ: {[b[1] for b in bad if b[0] == 'vol'][:5]}" if bad else "")
            if every == 1:
                low = np.where(count > 0, np.rint(np.clip(low, 0, 1) * 255), 0)
                got = np.fromfile(os.path.join(d, "min_m.u8"), np.uint8)
                compare(res, rid, "lowest column |m| (x 255)", got, low, 1, "(one step: float32 in hamrvis)")
    elif e.get("grain_map", {}).get("from_run"):
        res.add(rid, "grain map", 1, 0, "copied", True, f"shared from run {e['grain_map']['from_run']} (no coordinates in this run)")

    # dipole cells
    if e.get("dipole"):
        P = e["dipole"]
        cells = ref.dipole_cells(segs[-1])
        n = (P["nx"], P["ny"], P["nz"])
        ok = len(cells) == n[0] * n[1] * n[2] and (cells[:, 0] == np.arange(len(cells))).all()
        if not ok:
            res.add(rid, "dipole cells", len(cells), math.inf, "exact", False, "cell count or order")
        else:
            B = np.fromfile(os.path.join(d, "dipole.f32"), np.float32).reshape(3, n[2], n[1], n[0])
            want = cells[:, 6:9].reshape(n[0], n[1], n[2], 3).transpose(3, 2, 1, 0).astype(np.float32)
            compare(res, rid, "dipole field (float32)", B, want, 0)
            A = np.fromfile(os.path.join(d, "dipole_atoms.u16"), np.uint16).reshape(n[2], n[1], n[0])
            compare(res, rid, "atoms per dipole cell", A, np.minimum(cells[:, 1], 65535).reshape(n).transpose(2, 1, 0), 0)


def add_arguments(ap):
    ap.add_argument("bundle", help="bundle directory made by 'reduce'")
    ap.add_argument("--run", nargs="+", help="run ids (default: all)")
    ap.add_argument("--every", type=int, default=1, help="compare every N-th snapshot (default 1: all)")
    ap.add_argument("--json", help="write every check's result to this file")


def verify(bundle, runs=None, every=1, json_path=None):
    """Check the bundle, store the result in its manifest and return the exit status (0 ok, 1 differences, 2 no raw files)."""
    mp = os.path.join(bundle, "manifest.json")
    if not os.path.isfile(mp):
        print(f"{bundle} is not a hamrvis bundle: it has no manifest.json")
        return 2
    man = json.load(open(mp, encoding="utf-8"))
    todo = [e for e in man["runs"] if not runs or e["id"] in runs]
    if runs and len(todo) != len(set(runs)):
        print(f"no run {', '.join(sorted(set(runs) - {e['id'] for e in todo}))} in this bundle (it has {', '.join(e['id'] for e in man['runs'])})")
        return 2
    gone = [s for e in todo for s in e["segments"] if not os.path.isdir(s)]
    if gone:
        print("run director" + ("ies" if len(gone) > 1 else "y") + " not found:\n  " + "\n  ".join(gone) +
              "\nverify reads the raw run files the bundle was made from, at the paths recorded when it was reduced; "
              "run it where they are (a bundle copied away from them cannot be checked)")
        return 2
    res = Results(print)
    start = time.time()
    for e in todo:
        print(f"{e['id']} ({e['type']}): {e['source']}", flush=True)
        check_run(bundle, e, res, every, lambda m: print(m, flush=True))
    bad = res.failed()
    print(f"{len(res.rows)} checks, {len(bad)} failed, {time.time() - start:.0f} s")
    if json_path:
        with open(json_path, "w") as f:
            json.dump([dict(zip(("run", "quantity", "values", "max_diff", "tolerance", "ok", "note"), r)) for r in res.rows], f, indent=1)
    man["verify"] = dict(hamrvis_version=__version__, date=time.strftime("%Y-%m-%d %H:%M:%S"), checks=len(res.rows), failed=len(bad),
                         runs=None if not runs else [e["id"] for e in todo], every=every)
    tmp = mp + ".tmp"                                   # manifest only, replaced in one step
    with open(tmp, "w") as f:
        json.dump(man, f, separators=(",", ":"))
    os.replace(tmp, mp)
    return 1 if bad else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="Check a hamrvis bundle against the VAMPIRE runs it was made from.")
    add_arguments(ap)
    a = ap.parse_args(argv)
    if a.every < 1:
        raise SystemExit("--every must be 1 or more")
    return verify(a.bundle, a.run, a.every, a.json)


if __name__ == "__main__":
    sys.exit(main())

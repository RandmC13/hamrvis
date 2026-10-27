"""Reduce VAMPIRE HAMR run directories (relax / write / read) to a compact bundle the viewer loads.

Bundle layout (see docs/bundle-format.md):
    manifest.json                     run parameters, stage list, array index
    <stage>/output.f32                'output' file, [rows, cols]
    <stage>/grains.i16                grain.txt per-grain block, [frames, grains, 4] x 1/32767
    <stage>/grain_lead.f64            grain.txt leading columns (time, temperature), [frames, n_lead]
    <stage>/top/NNNNNNNN.i8.gz        snapshot NNNNNNNN: column-averaged magnetisation, [3, ny, nx] x 1/127
    <stage>/vol/NNNNNNNN.i8.gz        snapshot NNNNNNNN: per-layer magnetisation, [3, nz, ny/2, nx/2] x 1/127
    <stage>/top_count.u8, vol_count.u8  magnetic atoms per pixel (0 = no magnetic atom: grain boundary)
    <stage>/min_m.u8                  lowest column |m| any snapshot reached, x 1/255
    <stage>/grain_map.i16             grain index of the magnetic atoms in each column, [ny, nx] (-1: none)
    read/dipole.f32                   hierarchical dipole cell field, [3, nz, ny, nx] (T)
    read/dipole_atoms.u16             atoms per dipole cell (0 = field probe cell)
Empty pixels hold -128 in the int8 arrays. Everything is little-endian.
"""
import gzip
import json
import math
import multiprocessing
import os
import shutil
import time

import numpy as np

from . import vampire as vp

FORMAT_VERSION = 1

# write analysis thresholds (see docs/analysis.md)
DECISION = 0.05     # m_z|m| below this in the final direction = not yet decided
HEATED = 0.1        # |m| below this at some point = the grain was melted
SETTLED = 0.3       # |final m_z|m|| below this = still disordered at the end
FULL_FIELD = 0.4    # |B| / Hpeak at the decision above this = full field, below = reversal


def log(msg):
    print(f"[hamrvis {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ------------------------------------------------------------------------------------------------ helpers
def magnetic_types(run_dir):
    """Material indices (0-based, as in the coordinate files) that are magnetic, from the material file."""
    inp = vp.Input(run_dir)
    mat = vp.read_keys(os.path.join(run_dir, inp.str("material:file")))
    n = max(int(k.split("[")[1].split("]")[0]) for k in mat if k.startswith("material["))
    nonmag = {int(k.split("[")[1].split("]")[0]) for k in mat if k.endswith(":non-magnetic")}
    return [i - 1 for i in range(1, n + 1) if i not in nonmag]


def _pitch(values, tol=1e-3):
    """Smallest spacing between distinct coordinate values (A); values closer than 2 tol count as equal."""
    u = np.unique(np.round(np.asarray(values) / tol) * tol)
    d = np.diff(u)
    d = d[d > 2 * tol]
    return float(d.min()) if len(d) else 1.0


class Grid:
    """Pixel grids: top-down columns at the atomic in-plane pitch, and a per-layer volume at twice that pitch
    (one atomic site per layer per pixel can be empty in non-cubic lattices; 2x2 sites never are inside a grain)."""

    def __init__(self, Lx_nm, Ly_nm, xyz):
        self.p = _pitch(xyz[:, 0])                                # A
        self.pz = _pitch(xyz[:, 2])
        self.z0 = float(xyz[:, 2].min())
        self.nx = int(round(Lx_nm * 10 / self.p))
        self.ny = int(round(Ly_nm * 10 / self.p))
        self.nz = int(round((xyz[:, 2].max() - self.z0) / self.pz)) + 1
        self.vx, self.vy = (self.nx + 1) // 2, (self.ny + 1) // 2

    def indices(self, xyz):
        ix = np.rint(xyz[:, 0] / self.p).astype(np.int64) % self.nx
        iy = np.rint(xyz[:, 1] / self.p).astype(np.int64) % self.ny
        iz = np.clip(np.rint((xyz[:, 2] - self.z0) / self.pz).astype(np.int64), 0, self.nz - 1)
        return iy * self.nx + ix, (iz * self.vy + iy // 2) * self.vx + ix // 2

    def info(self):
        return dict(pitch_nm=self.p / 10, pitch_z_nm=self.pz / 10, z0_nm=self.z0 / 10, nx=self.nx, ny=self.ny,
                    nz=self.nz, vx=self.vx, vy=self.vy)


def _q8(sums, count):
    """int8 mean direction components (x 1/127), -128 where a pixel has no atoms."""
    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.rint(127 * sums / count)
    return np.where(count > 0, np.clip(q, -127, 127), -128).astype(np.int8)


def _write_gz(path, arr):
    with gzip.open(path, "wb", compresslevel=6) as f:
        f.write(np.ascontiguousarray(arr).tobytes())


# ------------------------------------------------------------------------------------------------ snapshots
_W = {}   # per-worker state, inherited through fork


def _bin_snapshot(job):
    idx, files, out_dir = job
    g, nG = _W["grid"], _W["n_grains"]
    top = np.zeros((3, g.ny * g.nx))
    vol = np.zeros((3, g.nz * g.vy * g.vx)) if _W["volume"] else None
    gs = np.zeros((3, nG))
    for f, (mask, kt, kv, gi) in zip(files, _W["bins"]):
        s = vp.read_spins(f)
        if len(s) != len(mask):
            raise IOError(f"{f}: {len(s)} spins, but its coordinate file has {len(mask)} atoms")
        s = s[mask]
        for c in range(3):
            top[c] += np.bincount(kt, s[:, c], minlength=top.shape[1])
            gs[c] += np.bincount(gi, s[:, c], minlength=nG)
            if vol is not None:
                vol[c] += np.bincount(kv, s[:, c], minlength=vol.shape[1])
    qt = _q8(top, _W["ct"])
    _write_gz(os.path.join(out_dir, "top", f"{idx:08d}.i8.gz"), qt.reshape(3, g.ny, g.nx))
    if vol is not None:
        _write_gz(os.path.join(out_dir, "vol", f"{idx:08d}.i8.gz"), _q8(vol, _W["cv"]).reshape(3, g.nz, g.vy, g.vx))
    m = np.sqrt((top ** 2).sum(0)) / np.maximum(_W["ct"], 1)
    # largest rounding error of the stored columns against the full-precision means (checks the int8 encoding)
    with np.errstate(invalid="ignore", divide="ignore"):
        q_err = np.nanmax(np.abs(qt[:, _W["ct"] > 0] / 127 - top[:, _W["ct"] > 0] / _W["ct"][_W["ct"] > 0]))
    return idx, m.astype(np.float32), gs, float(q_err)


def reduce_atoms(run_dir, out_dir, entry, grain_steps=None, grain_m=None, volume=True, every=1, workers=None):
    """From the atom coordinates: exact top-down grain map, grain centres; then bin every spin snapshot onto the
    grids and check each snapshot's per-grain mean spin against grain.txt at the same time step."""
    cfiles = vp.coord_files(run_dir)
    if not cfiles:
        return
    Lx, Ly = entry["size_nm"][:2]
    mtypes, ids = entry["magnetic_types"], entry["grains"]["id"]
    nG = len(ids)
    gindex = np.full(max(ids) + 2, nG, np.int64)          # grain ID -> grain index; nG collects other IDs
    gindex[ids] = np.arange(nG)
    log(f"  reading {len(cfiles)} coordinate files")
    coords = [vp.read_coords(f) for f in cfiles]
    allmag = np.concatenate([xyz[np.isin(t, mtypes)] for t, xyz, _ in coords])
    grid = Grid(Lx, Ly, allmag)
    del allmag
    bins, ct, cv = [], np.zeros(grid.ny * grid.nx), np.zeros(grid.nz * grid.vy * grid.vx)
    pairs = []
    cx, sx, cy, sy, na = (np.zeros(nG + 1) for _ in range(5))
    for t, xyz, g in coords:
        mask = np.isin(t, mtypes)
        kt, kv = grid.indices(xyz[mask])
        gi = gindex[np.clip(g[mask], 0, len(gindex) - 1)] if g is not None else np.full(mask.sum(), nG)
        ct += np.bincount(kt, minlength=ct.size)
        cv += np.bincount(kv, minlength=cv.size)
        ax, ay = 2 * np.pi * xyz[mask, 0] / (Lx * 10), 2 * np.pi * xyz[mask, 1] / (Ly * 10)
        for acc, v in ((cx, np.cos(ax)), (sx, np.sin(ax)), (cy, np.cos(ay)), (sy, np.sin(ay)), (na, 1.0)):
            acc += np.bincount(gi, np.broadcast_to(v, gi.shape), minlength=nG + 1)
        pairs.append(kt * (nG + 1) + gi)                    # (column, grain) key
        bins.append((mask, kt, kv, gi))
    del coords
    files_bytes = sum(os.path.getsize(f) for f in cfiles)
    entry["atoms_source"] = dict(coordinate_files=len(cfiles), bytes=files_bytes)

    # grain map: the grain of the magnetic atoms in each column (columns hold one grain in a columnar film; the
    # majority grain otherwise); -1 = no magnetic atom
    key, cnt = np.unique(np.concatenate(pairs), return_counts=True)
    del pairs
    col, gr = key // (nG + 1), key % (nG + 1)
    order = np.lexsort((-cnt, col))                       # per column, most atoms first
    col, gr = col[order], gr[order]
    first = np.r_[True, col[1:] != col[:-1]]
    gm = np.full(grid.ny * grid.nx, -1, np.int64)
    gm[col[first]] = np.where(gr[first] == nG, -1, gr[first])
    mixed = int(len(col) - first.sum())
    gm.astype(np.int16).tofile(os.path.join(out_dir, "grain_map.i16"))
    entry["grain_map"] = dict(nx=grid.nx, ny=grid.ny, pitch_nm=grid.p / 10, mixed_columns=mixed)
    # periodic-aware centres (circular means) and magnetic atom counts replace grain-coordinates.txt values,
    # which include non-magnetic atoms assigned to a grain
    ok = na[:nG] > 0
    gx = (np.arctan2(sx[:nG], cx[:nG]) % (2 * np.pi)) * Lx / (2 * np.pi)
    gy = (np.arctan2(sy[:nG], cy[:nG]) % (2 * np.pi)) * Ly / (2 * np.pi)
    G = entry["grains"]
    G["x"] = np.round(np.where(ok, gx, G["x"]), 3).tolist()
    G["y"] = np.round(np.where(ok, gy, G["y"]), 3).tolist()
    G["magnetic_atoms"] = na[:nG].astype(int).tolist()
    G.pop("outline", None)

    snaps = vp.snapshots(run_dir)
    if not snaps:
        return
    snaps = snaps[::every] if snaps[-1] in snaps[::every] else snaps[::every] + [snaps[-1]]
    log(f"  {len(snaps)} snapshots")
    for d in ("top", "vol") if volume else ("top",):
        os.makedirs(os.path.join(out_dir, d), exist_ok=True)
    np.minimum(ct, 255).astype(np.uint8).tofile(os.path.join(out_dir, "top_count.u8"))
    np.minimum(cv, 255).astype(np.uint8).tofile(os.path.join(out_dir, "vol_count.u8"))
    _W.update(grid=grid, bins=bins, ct=ct, cv=cv, volume=volume, n_grains=nG)
    jobs = [(i, files, out_dir) for i, _, files in snaps]
    step_of = {i: int(round(t / entry["dt"])) for i, t, _ in snaps}
    workers = workers or min(len(jobs), max(1, (os.cpu_count() or 2) - 1))
    min_m = np.full(ct.size, np.inf, np.float32)
    checks, q_max, t0 = [], 0.0, time.time()
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        for n, (idx, m, gs, q_err) in enumerate(pool.imap_unordered(_bin_snapshot, jobs), 1):
            np.minimum(min_m, m, out=min_m)
            q_max = max(q_max, q_err)
            if grain_m is not None:
                # per-grain means against grain.txt at the same step (they agree exactly then), else one step
                # earlier: VAMPIRE writes snapshot k at step k * rate + 1, one step after the grain row of that output
                for lag in (0, 1):
                    f = np.nonzero(grain_steps == step_of[idx] - lag)[0]
                    if len(f):
                        ref = grain_m[f[0], :, :3] * grain_m[f[0], :, 3:4]       # direction x |m|
                        mine = (gs / np.maximum(na[:nG], 1)).T
                        checks.append((idx, lag, float(np.abs(mine - ref)[ok].max())))
                        break
            if n % max(1, len(jobs) // 20) == 0 or n == len(jobs):
                el = time.time() - t0
                log(f"  snapshot {n}/{len(jobs)} ({el:.0f} s, ~{el / n * (len(jobs) - n):.0f} s left)")
    _W.clear()
    min_m[ct == 0] = 0
    np.rint(np.clip(min_m, 0, 1) * 255).astype(np.uint8).tofile(os.path.join(out_dir, "min_m.u8"))
    entry["snapshots"] = dict(grid=grid.info(), volume=volume, frames=[dict(index=i, time=t) for i, t, _ in snaps])
    same = [c[2] for c in checks if c[1] == 0]
    prev = [c[2] for c in checks if c[1] == 1]
    entry["checks"] = dict(column_rounding_max=q_max,
                           grain_txt_same_step=dict(n=len(same), max_diff=max(same, default=None)),
                           grain_txt_one_step_earlier=dict(n=len(prev), max_diff=max(prev, default=None),
                                                           median_diff=float(np.median(prev)) if prev else None))
    log(f"  checks: stored columns within {q_max:.4f} of the exact column means; per-grain means vs grain.txt: "
        f"{len(same)} snapshots at the same step (max diff {max(same, default=0):.1e}), {len(prev)} one step "
        f"earlier (max {max(prev, default=0):.1e}, median {np.median(prev) if prev else 0:.1e}: one step of motion)")


# ------------------------------------------------------------------------------------------------ stage tables
def reduce_tables(run_dir, out_dir, entry):
    cols, out = vp.read_output(run_dir)
    out.astype(np.float32).tofile(os.path.join(out_dir, "output.f32"))
    entry["output"] = dict(columns=cols, rows=len(out))
    gc = vp.read_grain_coordinates(run_dir)
    outlines = vp.grain_outlines(gc, vp.read_grain_shapes(run_dir))
    entry["grains"] = dict(
        id=gc[:, 0].astype(int).tolist(), x=np.round(gc[:, 1], 3).tolist(), y=np.round(gc[:, 2], 3).tolist(),
        atoms=gc[:, 3].astype(int).tolist(),
        outline=[None if o is None else np.round(o, 3).tolist() for o in outlines])
    if os.path.exists(os.path.join(run_dir, "grain.txt")):
        log("  parsing grain.txt")
        n_lead, lead, m = vp.read_grain_file(run_dir, len(gc))
        np.rint(np.clip(m, -1, 1) * 32767).astype(np.int16).tofile(os.path.join(out_dir, "grains.i16"))
        lead.tofile(os.path.join(out_dir, "grain_lead.f64"))
        names = ["time", "temperature", "extra"][:n_lead]
        entry["grain_series"] = dict(frames=len(m), lead=names)
        return lead, m
    return None, None


def reduce_dipole(run_dir, out_dir, entry):
    """Grid dipole-field.txt by cell ID: macrocells are numbered x-major, z fastest, on a grid of
    (L / cell size + 1) cells per axis. Its x, y, z columns are atom centroids, not cell centres."""
    d = vp.read_dipole_cells(run_dir)
    cell = vp.Input(run_dir).num("cells:macro-cell-size", 2.0)              # nm
    n = [int(round(L / cell)) + 1 for L in entry["size_nm"]]                # nx, ny, nz
    if n[0] * n[1] * n[2] != len(d) or not np.array_equal(d[:, 0], np.arange(len(d))):
        raise ValueError(f"dipole-field.txt has {len(d)} cells, expected {n[0]} x {n[1]} x {n[2]} in ID order")
    # file order is (x, y, z) with z fastest: reshape to [nx, ny, nz], store as [nz, ny, nx]
    B = d[:, 6:9].reshape(n[0], n[1], n[2], 3).transpose(3, 2, 1, 0).astype(np.float32)
    A = np.minimum(d[:, 1], 65535).reshape(n[0], n[1], n[2]).transpose(2, 1, 0).astype(np.uint16)
    np.ascontiguousarray(B).tofile(os.path.join(out_dir, "dipole.f32"))
    np.ascontiguousarray(A).tofile(os.path.join(out_dir, "dipole_atoms.u16"))
    entry["dipole"] = dict(nx=n[0], ny=n[1], nz=n[2], cell_nm=cell)
    # shortcut: atomistic_dipole_field.txt (per-atom read-back field) is not reduced yet; add it when the
    # cell field's 1 nm resolution is not enough for read-back analysis


# ------------------------------------------------------------------------------------------------ write analysis
def write_analysis(h, entry, lead, m, Tc=None):
    """Per-grain outcome of the write: intended polarity, final state, decision moment and its circumstances."""
    g = entry["grains"]
    gx, gy = np.array(g["x"]), np.array(g["y"])
    steps = np.rint(lead[:, 0] / h.dt).astype(np.int64)
    mz = m[:, :, 2] * m[:, :, 3]
    T = h.temperature(gx[None, :], gy[None, :], steps[:, None])
    B = h.field(steps) * h.field_uz                                  # z field inside the box
    final = mz[-1]
    sgn = np.sign(final)
    below = mz * sgn < DECISION
    heated = (m[:, :, 3].min(0) < HEATED) & below.any(0)
    # decision frame: the frame after the last one below the threshold in the final direction
    dec = np.where(below.any(0), len(steps) - np.argmax(below[::-1], axis=0), 0)
    dec = np.minimum(dec, len(steps) - 1)

    if Tc is None:   # effective freezing temperature: median centre temperature at the decision of heated grains
        sel = heated & (np.abs(final) >= SETTLED)
        Tc = float(np.median(T[dec[sel], np.nonzero(sel)[0]])) if sel.any() else None
        entry["Tc_source"] = "median grain-centre laser temperature at the decision moment of heated grains"
    else:
        entry["Tc_source"] = "given on the command line"
    entry["Tc"] = Tc
    lag = h.freeze_lag(Tc) if Tc and h.Tmin < Tc < h.Tmax else 0.0
    entry["freeze_lag_nm"] = lag

    res = []
    for i in range(len(gx)):
        t = int(h.track_of(gy[i]))
        b = int(math.floor((gx[i] + lag) / h.BL))
        seq = h.sequence
        pol = seq[t * h.bits_per_track + b] if 0 <= b < h.bits_per_track and t * h.bits_per_track + b < len(seq) else 0
        want = int(np.sign(pol * h.field_uz))                         # the z direction the field asked for
        r = dict(track=t, bit=b, intended=want, final=round(float(final[i]), 4),
                 peakT=round(float(T[:, i].max()), 1), minM=round(float(m[:, i, 3].min()), 4))
        if not heated[i]:
            r["category"] = "unheated"
        elif abs(final[i]) < SETTLED:
            r["category"] = "undecided"
        else:
            d = int(dec[i])
            hx, hy, _ = h.head(steps[d])
            bz = float(B[d]) if h.in_field_box(gx[i], gy[i], steps[d]) else 0.0
            bd = float(h.field(steps[d]))
            r["category"] = ("zero-field" if bd == 0 else "reversal" if abs(bd) < FULL_FIELD * h.Hpeak
                             else "full-field")
            r.update(decision_frame=d, decision_step=int(steps[d]), Bz_dec=round(bz, 3), B_dec=round(bd, 3),
                     T_dec=round(float(T[d, i]), 1), m_dec=round(float(m[d, i, 3]), 4),
                     lag_dec=round(float(hx - gx[i]), 3), dy=round(float(gy[i] - hy), 3))
        r["correct"] = bool(r["category"] not in ("unheated", "undecided") and want != 0 and sgn[i] == want)
        res.append(r)
    entry["analysis"] = dict(grains=res, thresholds=dict(decision=DECISION, heated=HEATED, settled=SETTLED,
                                                         full_field=FULL_FIELD))
    cats = {}
    for r in res:
        cats[r["category"]] = cats.get(r["category"], 0) + 1
    written = [r for r in res if r["intended"] != 0 and r["category"] != "unheated"]
    log(f"  Tc {Tc:.0f} K, lag {lag:.2f} nm; categories {cats}; correct {sum(r['correct'] for r in res)}"
        f" of {len(written)} heated grains inside the bit windows")


def share_grain_map(out, manifest):
    """Stages without coordinate output (e.g. a read that wrote no config) take the grain map and centres of a stage
    of the same structure (same size and grains); only then fall back to the approximate Voronoi outlines."""
    st = manifest["stages"]
    donor = next((s for s in ("write", "relax", "read") if s in st and "grain_map" in st[s]), None)
    for s, e in st.items():
        if "grain_map" in e or donor is None:
            continue
        d = st[donor]
        if d["size_nm"] == e["size_nm"] and d["grains"]["id"] == e["grains"]["id"]:
            shutil.copy(os.path.join(out, donor, "grain_map.i16"), os.path.join(out, s, "grain_map.i16"))
            e["grain_map"] = dict(d["grain_map"], from_stage=donor)
            for k in ("x", "y", "magnetic_atoms"):
                e["grains"][k] = d["grains"][k]
            e["grains"].pop("outline", None)
    with open(os.path.join(out, "manifest.json"), "w") as f:
        json.dump(manifest, f, separators=(",", ":"))


# ------------------------------------------------------------------------------------------------ driver
def find_runs(paths):
    """Run directories (folders with an 'input' file) among paths and their immediate subfolders."""
    runs = []
    for p in paths:
        if os.path.exists(os.path.join(p, "input")):
            runs.append(p)
        else:
            runs += sorted(os.path.join(p, d) for d in os.listdir(p) if os.path.exists(os.path.join(p, d, "input")))
    return runs


def reduce(paths, out, volume=True, every=1, Tc=None, workers=None, name=None):
    runs = find_runs(paths)
    if not runs:
        raise SystemExit(f"no run directories (folders with an 'input' file) in {paths}")
    os.makedirs(out, exist_ok=True)
    manifest = dict(format="hamrvis-bundle", version=FORMAT_VERSION, created=time.strftime("%Y-%m-%d %H:%M:%S"),
                    name=name or os.path.basename(os.path.commonpath([os.path.abspath(r) for r in runs]) if len(runs) > 1
                                                 else os.path.abspath(runs[0])), stages={})
    for run in runs:
        stage = vp.stage_of(run)
        if stage in manifest["stages"]:
            raise SystemExit(f"two {stage} runs ({manifest['stages'][stage]['source']} and {run}): reduce them "
                             f"into separate bundles")
        log(f"{stage}: {run}")
        sd = os.path.join(out, stage)
        os.makedirs(sd, exist_ok=True)
        inp = vp.Input(run)
        e = dict(source=os.path.abspath(run), input={k: v for k, (v, u) in inp.keys.items()},
                 input_units={k: u for k, (v, u) in inp.keys.items() if u},
                 dt=inp.num("sim:time-step", 1e-15, "s"),
                 size_nm=[inp.num(f"dimensions:system-size-{a}") for a in "xyz"],
                 magnetic_types=magnetic_types(run))
        for f in ("log", "stdout", "time.txt"):
            if os.path.exists(os.path.join(run, f)):
                e.setdefault("files", []).append(f)
        meta = os.path.join(run, "atoms-coords.meta")
        if os.path.exists(meta):
            e["n_atoms"] = int(vp._meta_files(meta)[0].get("Number of atoms", "0").split()[0])
        e["source_files"] = {f: [os.path.getsize(os.path.join(run, f)), int(os.path.getmtime(os.path.join(run, f)))]
                             for f in ("input", "output", "grain.txt", "dipole-field.txt") if os.path.exists(os.path.join(run, f))}
        lead, m = reduce_tables(run, sd, e)
        gsteps = np.rint(lead[:, 0] / e["dt"]).astype(np.int64) if lead is not None else None
        reduce_atoms(run, sd, e, gsteps, m, volume, every, workers)
        if stage == "write":
            h = vp.Hamr(run)
            e["hamr"] = h.params()
            if lead is not None:
                write_analysis(h, e, lead, m, Tc)
        if stage == "read" and os.path.exists(os.path.join(run, "dipole-field.txt")):
            log("  gridding dipole-field.txt")
            reduce_dipole(run, sd, e)
        manifest["stages"][stage] = e
        with open(os.path.join(out, "manifest.json"), "w") as f:   # after each stage, so a partial bundle opens
            json.dump(manifest, f, separators=(",", ":"))
    share_grain_map(out, manifest)
    size = sum(os.path.getsize(os.path.join(d, f)) for d, _, fs in os.walk(out) for f in fs)
    log(f"bundle {out}: {size / 1e6:.1f} MB")

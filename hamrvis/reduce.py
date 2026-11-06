"""Reduce VAMPIRE runs (any mix of HAMR writes, Monte Carlo runs, dipole reads, ...) to a compact bundle the viewer loads.

Bundle layout (see docs/bundle-format.md):
    manifest.json                     run list (in run order), parameters, array index
    <run>/output.f64                  'output' file, [rows, cols]
    <run>/grains.i16                grain.txt per-grain block, [frames, grains, 4] x 1/32767
    <run>/grain_lead.f64            grain.txt leading columns (time, temperature), [frames, n_lead]
    <run>/grains_overview.i16       long series only: per time bin mz min, mz max, |m| min, |m| max, [bins, grains, 4] x 1/32767
    <run>/top/NNNNNNNN.i8.gz        snapshot NNNNNNNN: column-averaged magnetisation, [3, ny, nx] x 1/127
    <run>/vol/NNNNNNNN.i8.gz        snapshot NNNNNNNN: per-layer magnetisation, [3, nz, ny/2, nx/2] x 1/127
                                      (.i16.gz x 1/32767 or .f32.gz with --precision 16 / 32)
    <run>/top_count.u8, vol_count.u8  magnetic atoms per pixel (0 = no magnetic atom: grain boundary)
    <run>/min_m.u8                  lowest column |m| any snapshot reached, x 1/255
    <run>/grain_map.i16             grain index of the magnetic atoms in each column, [ny, nx] (-1: none)
    <run>/dipole.f32                  dipole cell field, [3, nz, ny, nx] (T)
    <run>/dipole_atoms.u16             atoms per dipole cell (0 = field probe cell)
Empty pixels hold -128 (int8), -32768 (int16) or NaN (float32); grains.f32 replaces grains.i16 with
--grain-precision 32. Everything is little-endian.
"""
import contextlib
import gzip
import json
import math
import multiprocessing
import os
import re
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool

import numpy as np

from . import __version__
from . import runs as rn
from . import vampire as vp

FORMAT_VERSION = 2
WORK = ".reduce-work"   # inside the bundle while reducing: per-atom indices and progress; removed at the end

# write analysis thresholds (see docs/analysis.md)
DECISION = 0.05     # m_z|m| below this in the final direction = not yet decided
HEATED = 0.1        # |m| below this at some point = the grain was melted
SETTLED = 0.3       # |final m_z|m|| below this = still disordered at the end
FULL_FIELD = 0.4    # |B| / Hpeak at the decision above this = full field, below = reversal
THRESHOLDS = dict(decision=DECISION, heated=HEATED, settled=SETTLED, full_field=FULL_FIELD)


def log(msg):
    print(f"[hamrvis {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ------------------------------------------------------------------------------------------------ helpers
def magnetic_types(run_dir):
    """Material indices (0-based, as in the coordinate files) that are magnetic, from the material file, or None
    (every material counts as magnetic) when the run's material file is not there."""
    inp = vp.Input(run_dir)
    name = inp.str("material:file")
    path = os.path.join(run_dir, name) if name else None
    if not path or not os.path.isfile(path):
        log(f"  WARNING: material file {name!r} not in {run_dir}: every material is treated as magnetic "
            f"(non-magnetic atoms would then count in the atom views)" if name else
            f"  WARNING: the input names no material:file: every material is treated as magnetic")
        return None
    mat = vp.read_keys(path)
    idx = [int(k.split("[")[1].split("]")[0]) for k in mat if k.startswith("material[")]
    if not idx:
        log(f"  WARNING: {name} defines no material[n] keys: every material is treated as magnetic")
        return None
    nonmag = {int(k.split("[")[1].split("]")[0]) for k in mat if k.startswith("material[") and k.endswith(":non-magnetic")}
    return [i - 1 for i in range(1, max(idx) + 1) if i not in nonmag]


TOL = 1e-3         # coordinates closer than 2 TOL (A) count as equal


def axis_stats(xyz, stats=None):
    """Add atoms to the extent of a run's lattice, which Grid needs: distinct x and z values and the z range."""
    ux, uz = np.unique(np.round(xyz[:, 0] / TOL) * TOL), np.unique(np.round(xyz[:, 2] / TOL) * TOL)
    if stats is None:
        return dict(ux=ux, uz=uz, z0=float(xyz[:, 2].min()), z1=float(xyz[:, 2].max()))
    return dict(ux=np.union1d(stats["ux"], ux), uz=np.union1d(stats["uz"], uz),
                z0=min(stats["z0"], float(xyz[:, 2].min())), z1=max(stats["z1"], float(xyz[:, 2].max())))


def _pitch(u):
    """Smallest spacing between distinct coordinate values (A)."""
    d = np.diff(u)
    d = d[d > 2 * TOL]
    return float(d.min()) if len(d) else 1.0


TOP_MAX = 4e6      # top-down pixels per snapshot before the column map is coarsened
LAYER_MAX = 2.5e6  # per-layer cells per snapshot before the layer grid is coarsened


class Grid:
    """Pixel grids on the atomic in-plane pitch p (smallest spacing between atom x coordinates):
    top-down columns of tf x tf sites (tf = 1, i.e. single atom columns, unless that exceeds TOP_MAX pixels) and a
    per-layer grid of lf x lf sites per atomic layer (lf = 2 by default: one site per layer can be empty in
    non-cubic lattices, 2 x 2 sites never are inside a grain; larger when LAYER_MAX would be exceeded)."""

    def __init__(self, Lx_nm, Ly_nm, stats, layer_factor=None, top_max=TOP_MAX, layer_max=LAYER_MAX):
        self.p = _pitch(stats["ux"])                              # A
        self.pz = _pitch(stats["uz"])
        self.z0 = stats["z0"]
        self.nx0 = int(round(Lx_nm * 10 / self.p))
        self.ny0 = int(round(Ly_nm * 10 / self.p))
        self.nz = int(round((stats["z1"] - self.z0) / self.pz)) + 1
        self.tf = max(1, math.ceil(math.sqrt(self.nx0 * self.ny0 / top_max)))
        self.lf = layer_factor or max(2, math.ceil(math.sqrt(self.nx0 * self.ny0 * self.nz / layer_max)))
        self.nx, self.ny = -(-self.nx0 // self.tf), -(-self.ny0 // self.tf)
        self.vx, self.vy = -(-self.nx0 // self.lf), -(-self.ny0 // self.lf)

    def indices(self, xyz):
        ix = np.rint(xyz[:, 0] / self.p).astype(np.int64) % self.nx0
        iy = np.rint(xyz[:, 1] / self.p).astype(np.int64) % self.ny0
        iz = np.clip(np.rint((xyz[:, 2] - self.z0) / self.pz).astype(np.int64), 0, self.nz - 1)
        return (iy // self.tf) * self.nx + ix // self.tf, (iz * self.vy + iy // self.lf) * self.vx + ix // self.lf

    def info(self):
        return dict(pitch_nm=self.p * self.tf / 10, layer_pitch_nm=self.p * self.lf / 10, atom_pitch_nm=self.p / 10,
                    pitch_z_nm=self.pz / 10, z0_nm=self.z0 / 10, nx=self.nx, ny=self.ny, nz=self.nz, vx=self.vx,
                    vy=self.vy)


# snapshot grid encodings (--precision): bits -> (file tag, dtype, scale; None = stored as float32)
ENCODINGS = {8: ("i8", np.int8, 127), 16: ("i16", np.int16, 32767), 32: ("f32", np.float32, None)}


def _quantise(sums, count, bits=8):
    """Mean direction components: int8 x 1/127 (-128 where a pixel has no atoms), int16 x 1/32767 (-32768) or
    float32 (NaN)."""
    _, dtype, sc = ENCODINGS[bits]
    with np.errstate(invalid="ignore", divide="ignore"):
        if sc is None:
            return np.where(count > 0, sums / count, np.nan).astype(np.float32)
        q = np.rint(sc * sums / count)
    return np.where(count > 0, np.clip(q, -sc, sc), -sc - 1).astype(dtype)


def _decode(q, bits=8):
    sc = ENCODINGS[bits][2]
    return q.astype(np.float64) if sc is None else q / sc


def _write_gz(path, arr):
    with gzip.open(path, "wb", compresslevel=6) as f:
        f.write(np.ascontiguousarray(arr).tobytes())


# ------------------------------------------------------------------------------------------------ snapshots
_W = {}   # per-worker state, set by _init_worker


def _init_worker(state):
    """Worker set-up (works with spawn, so on every platform): the per-atom bin indices stay in files mapped read-only,
    so all workers share one copy in the page cache instead of each holding its own."""
    _W.update(state)
    _W["bins"] = [tuple(np.load(os.path.join(state["work"], f"bin{k}.{n}.npy"), mmap_mode="r") for n in BIN_ARRAYS)
                  for k in range(state["n_files"])]
    _W["ct"], _W["cv"] = (np.load(os.path.join(state["work"], f"{n}.npy")) for n in ("ct", "cv"))


BIN_ARRAYS = ("mask", "kt", "kv", "gi")


def _bin_snapshot(job):
    try:
        return _bin_snapshot_checked(job)
    except (IOError, ValueError) as err:          # truncated copy, missing file
        return job[0], None, None, str(err)


def _new_sums():
    g, nG = _W["grid"], _W["n_grains"]
    return (np.zeros((3, g.ny * g.nx)), np.zeros((3, g.nz * g.vy * g.vx)) if _W["volume"] else None,
            np.zeros((3, nG + 1)))                 # last grain bin: atoms of no listed grain


def _add_spins(sums, s, kt, kv, gi):
    """Add spins of magnetic atoms to the column, volume and grain sums."""
    top, vol, gs = sums
    for c in range(3):
        top[c] += np.bincount(kt, s[:, c], minlength=top.shape[1])
        gs[c] += np.bincount(gi, s[:, c], minlength=gs.shape[1])
        if vol is not None:
            vol[c] += np.bincount(kv, s[:, c], minlength=vol.shape[1])


def _bin_snapshot_checked(job):
    idx, files, out_dir = job
    sums = _new_sums()
    for f, (mask, kt, kv, gi) in zip(files, _W["bins"]):
        a, m = 0, 0                                # atoms and magnetic atoms done in this file
        for s in vp.spin_chunks(f):
            sel = np.asarray(mask[a:a + len(s)])
            if len(sel) != len(s):
                raise IOError(f"{f}: more spins than the {len(mask)} atoms of its coordinate file")
            s, c = s[sel], int(sel.sum())
            _add_spins(sums, s, kt[m:m + c], kv[m:m + c], gi[m:m + c])
            a, m = a + len(sel), m + c
        if a != len(mask):
            raise IOError(f"{f}: {a} spins, but its coordinate file has {len(mask)} atoms")
    return _write_snapshot(idx, out_dir, *sums)


def _window_series(snaps, out_dir, keep):
    """Results of the snapshots in keep, in order, for a run with window snapshots (hamr:active-window): the
    sums of a complete snapshot are binned in full, and a window snapshot changes them by the binned
    difference of the spins it holds (new minus old), so only the atoms near the head are touched."""
    state, sums = None, None
    mask_index = [np.cumsum(np.asarray(mask)).astype(np.int64) - 1 for mask, *_ in _W["bins"]]
    for idx, t, files, window in snaps:
        if not window:
            state = []
            sums = _new_sums()
            for f, (mask, kt, kv, gi) in zip(files, _W["bins"]):
                s = vp.read_spins(f)
                if len(s) != len(mask):
                    raise IOError(f"{f}: {len(s)} spins, but its coordinate file has {len(mask)} atoms")
                state.append(s)
                _add_spins(sums, s[np.asarray(mask)], kt, kv, gi)
        elif state is None:
            raise IOError(f"snapshot {idx} is a window snapshot with no complete snapshot before it")
        else:
            top, vol, gs = sums
            nG = _W["n_grains"]
            for s, f, (mask, kt, kv, gi), midx in zip(state, files, _W["bins"], mask_index):
                n, pos, new = vp.read_window_spins(f)
                if n != len(mask):
                    raise IOError(f"{f}: {n} atoms, but its coordinate file has {len(mask)}")
                d, sel = new - s[pos], np.asarray(mask[pos])
                k = midx[pos[sel]]
                for c in range(3):
                    top[c] += np.bincount(kt[k], d[sel, c], minlength=top.shape[1])
                    gs[c] += np.bincount(gi[k], d[sel, c], minlength=nG + 1)
                    if vol is not None:
                        vol[c] += np.bincount(kv[k], d[sel, c], minlength=vol.shape[1])
                s[pos] = new
        if idx in keep:
            yield _write_snapshot(idx, out_dir, *sums)


def _write_snapshot(idx, out_dir, top, vol, gs):
    g = _W["grid"]
    bits = _W["bits"]
    tag = ENCODINGS[bits][0]
    qt = _quantise(top, _W["ct"], bits)
    _write_gz(os.path.join(out_dir, "top", f"{idx:08d}.{tag}.gz"), qt.reshape(3, g.ny, g.nx))
    if vol is not None:
        _write_gz(os.path.join(out_dir, "vol", f"{idx:08d}.{tag}.gz"),
                  _quantise(vol, _W["cv"], bits).reshape(3, g.nz, g.vy, g.vx))
    m = np.sqrt((top ** 2).sum(0)) / np.maximum(_W["ct"], 1)
    # largest rounding error of the stored columns against the full-precision means (checks the int8 encoding)
    with np.errstate(invalid="ignore", divide="ignore"):
        q_err = np.nanmax(np.abs(_decode(qt[:, _W["ct"] > 0], bits) - top[:, _W["ct"] > 0] / _W["ct"][_W["ct"] > 0]))
    return idx, m.astype(np.float32), gs, float(q_err)


def available_memory():
    """Bytes this job may still use: free memory, or what is left of its cgroup (Slurm) limit if that is lower."""
    best = None
    try:
        with open("/proc/meminfo") as f:
            best = int(re.search(r"MemAvailable:\s*(\d+)", f.read()).group(1)) * 1024
    except (OSError, AttributeError):
        try:
            best = os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
        except (ValueError, OSError, AttributeError):
            pass
    slurm = os.environ.get("SLURM_MEM_PER_NODE") or (os.environ.get("SLURM_MEM_PER_CPU") and
                                                     int(os.environ["SLURM_MEM_PER_CPU"]) * _cpus())
    if slurm:                                      # MB; a tenth is left for the main process
        best = 0.9 * int(slurm) * 1e6 if best is None else min(best, 0.9 * int(slurm) * 1e6)
    try:
        with open("/sys/fs/cgroup/memory.max") as f, open("/sys/fs/cgroup/memory.current") as g:
            limit = f.read().strip()
            if limit != "max":
                left = int(limit) - int(g.read())
                best = left if best is None else min(best, left)
    except (OSError, ValueError):
        pass
    return best


def _cpus():
    try:
        return len(os.sched_getaffinity(0))        # the CPUs of this job, not of the node (Linux only)
    except AttributeError:
        return os.cpu_count() or 2


def worker_memory(grid, max_atoms_per_file, volume, bits=8):
    """Peak bytes of one snapshot worker: a piece of a spin file (twice during the mask copy), the sums and quantised grids."""
    atoms = min(max_atoms_per_file, vp.CHUNK_ATOMS)
    cells = grid.nx * grid.ny + (grid.nz * grid.vy * grid.vx if volume else 0)
    return 2 * 24 * atoms + cells * (24 + 3 * bits // 8) + 100e6


def worker_count(jobs, each, requested=None, mem=None):
    """(workers, why): the -j value, else one per CPU (one left free) as far as the memory budget allows."""
    if requested:
        return requested, "as requested"
    n = min(jobs, max(1, _cpus() - 1))
    budget = mem if mem else available_memory()
    if budget and n * each > budget:
        n = max(1, int(budget // each))
        return n, f"limited by memory: {each / 1e9:.1f} GB each, budget {budget / 1e9:.1f} GB"
    return n, "one per CPU"


def _resume_steps(dirs):
    """Step at which each segment after the first took over (None for the first)."""
    return [None] + [vp.checkpoint_resume_step(d) for d in dirs[1:]]


def run_snapshots(dirs, dt):
    """Snapshots of a run over its segments: each segment's own up to the step the next one resumed from (the later
    job wins where they overlap)."""
    res, stops = [], _resume_steps(dirs)[1:] + [None]
    for d, stop in zip(dirs, stops):
        res += [s for s in vp.snapshots(d) if stop is None or s[1] / dt <= stop + 0.5]
    idx = [s[0] for s in res]
    if len(set(idx)) != len(idx):
        raise SystemExit(f"{dirs}: snapshot numbers repeat across the jobs of this run")
    return res


def _join_rows(parts, times, dirs, dt):
    """Rows of each segment's table up to the step the next segment resumed from, concatenated."""
    stops = _resume_steps(dirs)[1:] + [None]
    keep = [a if stop is None or t is None else a[t / dt <= stop + 0.5] for a, t, stop in zip(parts, times, stops)]
    return np.concatenate(keep) if len(keep) > 1 else keep[0]


def _index_type(n):
    return np.int32 if n < 2 ** 31 else np.int64


def _merge_counts(keys, counts, new):
    """Add the keys in new to sorted unique keys with counts."""
    k, c = np.unique(new, return_counts=True)
    if keys is None:
        return k, c
    merged, inv = np.unique(np.concatenate([keys, k]), return_inverse=True)
    return merged, np.bincount(inv, np.concatenate([counts, c]), minlength=len(merged)).astype(np.int64)


def _coordinate_pass(cfiles, entry, out_dir, work, opts, layer_factor):
    """Read the coordinate files a piece at a time: lattice extent, then per-atom bin indices (written to work/, one
    set per data file), atom counts per pixel, grain map, grain centres. Returns what the snapshot pass needs."""
    Lx, Ly = entry["size_nm"][:2]
    mtypes, ids = entry["magnetic_types"], entry["grains"]["id"]
    nG = len(ids)
    gindex = np.full(max(ids, default=-1) + 2, nG, np.int64)          # grain ID -> grain index; nG collects other IDs
    gindex[ids] = np.arange(nG)
    magnetic = (lambda t: np.ones(len(t), bool)) if mtypes is None else (lambda t: np.isin(t, mtypes))
    log(f"  reading {len(cfiles)} coordinate files")
    stats, have_ids, atoms_per_file, magnetic_per_file = None, True, [], []
    for f in cfiles:
        n, nm = 0, 0
        for t, xyz, g in vp.coord_chunks(f):
            n += len(t)
            have_ids = have_ids and g is not None
            m = xyz[magnetic(t)]
            nm += len(m)
            if len(m):
                stats = axis_stats(m, stats)
        atoms_per_file.append(n)
        magnetic_per_file.append(nm)
    if stats is None:
        raise SystemExit(f"{cfiles[0]}: no magnetic atoms in the coordinate files")
    grid = Grid(Lx, Ly, stats, layer_factor, opts.get("top_max", TOP_MAX), opts.get("layer_max", LAYER_MAX))
    if not have_ids:
        log("  WARNING: the coordinate files hold no grain IDs (older text format): grain shapes stay VAMPIRE's "
            "approximate Voronoi outlines and snapshots are not checked against grain.txt")
    ct, cv = np.zeros(grid.ny * grid.nx), np.zeros(grid.nz * grid.vy * grid.vx)
    cx, sx, cy, sy, na = (np.zeros(nG + 1) for _ in range(5))
    keys, counts = None, None
    tt, tv = _index_type(ct.size), _index_type(cv.size)
    for k, f in enumerate(cfiles):
        out = {}
        for n, dtype, count in (("mask", bool, atoms_per_file[k]), ("kt", tt, magnetic_per_file[k]),
                                ("kv", tv, magnetic_per_file[k]), ("gi", np.int32, magnetic_per_file[k])):
            out[n] = open(os.path.join(work, f"bin{k}.{n}.npy"), "wb")   # raw arrays after a header, for np.load's mmap
            np.lib.format.write_array_header_1_0(out[n], dict(descr=np.lib.format.dtype_to_descr(np.dtype(dtype)),
                                                              fortran_order=False, shape=(count,)))
        for t, xyz, g in vp.coord_chunks(f):
            mask = magnetic(t)
            kt, kv = grid.indices(xyz[mask])
            gi = gindex[np.clip(g[mask], 0, len(gindex) - 1)] if have_ids else np.full(mask.sum(), nG)
            ct += np.bincount(kt, minlength=ct.size)
            cv += np.bincount(kv, minlength=cv.size)
            ax, ay = 2 * np.pi * xyz[mask, 0] / (Lx * 10), 2 * np.pi * xyz[mask, 1] / (Ly * 10)
            for acc, v in ((cx, np.cos(ax)), (sx, np.sin(ax)), (cy, np.cos(ay)), (sy, np.sin(ay)), (na, 1.0)):
                acc += np.bincount(gi, np.broadcast_to(v, gi.shape), minlength=nG + 1)
            if have_ids and nG:
                keys, counts = _merge_counts(keys, counts, kt.astype(np.int64) * (nG + 1) + gi)   # (column, grain)
            for n, a in (("mask", mask), ("kt", kt.astype(tt)), ("kv", kv.astype(tv)), ("gi", gi.astype(np.int32))):
                a.tofile(out[n])
        for fh in out.values():
            fh.close()
    entry["atoms_source"] = dict(coordinate_files=len(cfiles), bytes=sum(os.path.getsize(f) for f in cfiles))

    # grain map: the grain of the magnetic atoms in each column (columns hold one grain in a columnar film; the
    # majority grain otherwise); -1 = no magnetic atom
    if have_ids and nG:
        col, gr = keys // (nG + 1), keys % (nG + 1)
        order = np.lexsort((-counts, col))                # per column, most atoms first
        col, gr = col[order], gr[order]
        first = np.r_[True, col[1:] != col[:-1]]
        gm = np.full(grid.ny * grid.nx, -1, np.int64)
        gm[col[first]] = np.where(gr[first] == nG, -1, gr[first])
        gm.astype(np.int16).tofile(os.path.join(out_dir, "grain_map.i16"))
        entry["grain_map"] = dict(nx=grid.nx, ny=grid.ny, pitch_nm=grid.p * grid.tf / 10,
                                  mixed_columns=int(len(col) - first.sum()))
    # periodic-aware centres (circular means) and magnetic atom counts replace grain-coordinates.txt values,
    # which include non-magnetic atoms assigned to a grain
    ok = na[:nG] > 0
    if not have_ids:
        ok[:] = False
    gx = (np.arctan2(sx[:nG], cx[:nG]) % (2 * np.pi)) * Lx / (2 * np.pi)
    gy = (np.arctan2(sy[:nG], cy[:nG]) % (2 * np.pi)) * Ly / (2 * np.pi)
    G = entry["grains"]
    G["x"] = np.round(np.where(ok, gx, G["x"]), 3).tolist()
    G["y"] = np.round(np.where(ok, gy, G["y"]), 3).tolist()
    if have_ids:
        G["magnetic_atoms"] = na[:nG].astype(int).tolist()
        G.pop("outline", None)
    np.save(os.path.join(work, "ct.npy"), ct)
    np.save(os.path.join(work, "cv.npy"), cv)
    return dict(grid=grid.__dict__, have_ids=have_ids, ok=ok.tolist(), na=na[:nG].tolist(), max_atoms=max(atoms_per_file))


def _save_json(path, obj):
    with open(path + ".tmp", "w") as f:
        json.dump(obj, f)
    os.replace(path + ".tmp", path)


SAVE_EVERY = 60     # seconds between saves of a run's snapshot progress (for resuming an interrupted reduce)


def _save_progress(path, done, min_m, q_max, checks, bad):
    with open(path + ".tmp", "wb") as f:
        np.savez(f, done=np.array(sorted(done), np.int64), min_m=min_m, q_max=q_max,
                 checks=np.array(checks, np.float64).reshape(-1, 3), bad=np.array(sorted(bad), np.int64))
    os.replace(path + ".tmp", path)


def _until_worker_lost(results):
    """A worker killed from outside, nearly always for lack of memory, stops the reduce with advice, not a traceback."""
    try:
        yield from results
    except BrokenProcessPool:
        raise SystemExit("A snapshot worker was killed, most likely because memory ran out. Run the same command to carry "
                         "on, with fewer workers (-j) or a lower memory budget (--mem).") from None


def reduce_atoms(dirs, out_dir, entry, series=None, volume=True, every=1, workers=None, layer_factor=None, opts=None,
                 work=None):
    """From the atom coordinates: exact top-down grain map, grain centres; then bin every spin snapshot onto the
    grids and check each snapshot's per-grain mean spin against grain.txt at the same time step. work: folder for
    the per-atom indices and the progress, which an interrupted reduce carries on from."""
    run_dir = dirs[0]                                     # the jobs of a run share its atoms
    cfiles = vp.coord_files(run_dir)
    if not cfiles:
        return
    opts = opts or {}
    state_file = os.path.join(work, "atoms.json")
    if os.path.exists(state_file):
        A = json.load(open(state_file))
        entry.clear()
        entry.update(A["entry"])
        log("  coordinates already read by the interrupted reduce")
    else:
        A = _coordinate_pass(cfiles, entry, out_dir, work, opts, layer_factor)
        A["entry"] = entry
        _save_json(state_file, A)
    grid = Grid.__new__(Grid)
    grid.__dict__.update(A["grid"])
    have_ids, ok, na = A["have_ids"], np.array(A["ok"], bool), np.array(A["na"])
    nG = len(entry["grains"]["id"])
    ct, cv = np.load(os.path.join(work, "ct.npy")), np.load(os.path.join(work, "cv.npy"))
    if not opts.get("snapshots", True):
        log("  spin snapshots left out (--no-snapshots)")
        return
    snaps = all_snaps = run_snapshots(dirs, entry["dt"])
    if not snaps:
        return
    snaps = snaps[::every] if snaps[-1] in snaps[::every] else snaps[::every] + [snaps[-1]]
    log(f"  {len(snaps)} snapshots")
    for d in ("top", "vol") if volume else ("top",):
        os.makedirs(os.path.join(out_dir, d), exist_ok=True)
    np.minimum(ct, 255).astype(np.uint8).tofile(os.path.join(out_dir, "top_count.u8"))
    np.minimum(cv, 255).astype(np.uint8).tofile(os.path.join(out_dir, "vol_count.u8"))
    bits = opts.get("precision", 8)
    state = dict(grid=grid, volume=volume, n_grains=nG, bits=bits, work=work, n_files=len(cfiles))
    step_of = {i: int(round(t / entry["dt"])) for i, t, _, _ in snaps}
    min_m = np.full(ct.size, np.inf, np.float32)
    checks, q_max, t0, bad, done = [], 0.0, time.time(), set(), set()
    progress = os.path.join(work, "progress.npz")
    if os.path.exists(progress):
        with np.load(progress) as z:
            done, bad = set(z["done"].tolist()), set(z["bad"].tolist())
            min_m, q_max, checks = z["min_m"], float(z["q_max"]), [tuple(c) for c in z["checks"].tolist()]
        log(f"  carrying on: {len(done)} of {len(snaps)} snapshots done by the interrupted reduce")
    jobs = [(i, files, out_dir) for i, _, files, _ in snaps if i not in done]
    windowed = any(w for *_, w in all_snaps)
    if windowed:
        log("  window snapshots (hamr:active-window): read in order, each updating the previous one")
        _init_worker(state)
        pool = contextlib.nullcontext()
    else:
        each = worker_memory(grid, A["max_atoms"], volume, bits)
        n, why = worker_count(len(jobs) or 1, each, workers, opts.get("mem"))
        log(f"  {n} snapshot workers ({why})")
        pool = ProcessPoolExecutor(n, multiprocessing.get_context("spawn"), _init_worker, (state,))
    saved = time.time()
    with pool as P:
        results = _window_series(all_snaps, out_dir, {j[0] for j in jobs}) if windowed else \
            (f.result() for f in as_completed([P.submit(_bin_snapshot, j) for j in jobs]))
        for n, (idx, m, gs, q_err) in enumerate(_until_worker_lost(results), 1):
            if m is None:
                log(f"  WARNING: snapshot {idx} skipped: {q_err}")
                bad.add(idx)
                done.add(idx)
                continue
            done.add(idx)
            np.minimum(min_m, m, out=min_m)
            q_max = max(q_max, q_err)
            if series is not None and have_ids:
                # per-grain means against grain.txt at the same step (they agree exactly then), else one step
                # earlier: VAMPIRE writes snapshot k at step k * rate + 1, one step after the grain row of that output
                for lag in (0, 1):
                    frame = series["frames"].get(step_of[idx] - lag)
                    if frame is not None:
                        ref = frame[:, :3] * frame[:, 3:4]                       # direction x |m|
                        mine = (gs[:, :nG] / np.maximum(na, 1)).T
                        checks.append((idx, lag, float(np.abs(mine - ref)[ok].max())))
                        break
            if n % max(1, len(jobs) // 20) == 0 or n == len(jobs):
                el = time.time() - t0
                log(f"  snapshot {n}/{len(jobs)} ({el:.0f} s, ~{el / n * (len(jobs) - n):.0f} s left)")
            if time.time() - saved > SAVE_EVERY:
                _save_progress(progress, done, min_m, q_max, checks, bad)
                saved = time.time()
    _W.clear()
    min_m[ct == 0] = 0
    np.rint(np.clip(min_m, 0, 1) * 255).astype(np.uint8).tofile(os.path.join(out_dir, "min_m.u8"))
    entry["snapshots"] = dict(grid=grid.info(), volume=volume, skipped=sorted(bad),
                              frames=[dict(index=i, time=t) for i, t, *_ in snaps if i not in bad])
    if bits != 8:
        entry["snapshots"]["dtype"] = ENCODINGS[bits][0]
    same = [c[2] for c in checks if c[1] == 0]
    prev = [c[2] for c in checks if c[1] == 1]
    entry["checks"] = dict(column_rounding_max=q_max,
                           grain_txt_same_step=dict(n=len(same), max_diff=max(same, default=None)),
                           grain_txt_one_step_earlier=dict(n=len(prev), max_diff=max(prev, default=None),
                                                           median_diff=float(np.median(prev)) if prev else None))
    log(f"  checks: stored columns within {q_max:.4f} of the exact column means; per-grain means vs grain.txt: "
        f"{len(same)} snapshots at the same step (max diff {max(same, default=0):.1e}), {len(prev)} one step "
        f"earlier (max {max(prev, default=0):.1e}, median {np.median(prev) if prev else 0:.1e}: one step of motion)")


# ------------------------------------------------------------------------------------------------ tables
def reduce_tables(dirs, out_dir, entry, opts=None, check_steps=frozenset()):
    """'output' and grain.txt of a run (its segments joined), grain list and outlines. Returns (grain.txt leading columns,
    {acc: GrainAccumulator, frames: float32 frames at check_steps}), or (None, None) without a grain series."""
    run_dir, dt = dirs[0], entry["dt"]
    parts = [vp.read_output(d) for d in dirs]
    cols = parts[0][0]
    if any(c != cols for c, _ in parts):
        raise SystemExit(f"{dirs}: the jobs of this run have different output columns")
    ti = cols.index("Real_time") if "Real_time" in cols else None
    out = _join_rows([o for _, o in parts], [None if ti is None else o[:, ti] for _, o in parts], dirs, dt)
    out.astype(np.float64).tofile(os.path.join(out_dir, "output.f64"))   # float32 loses steps beyond 2^24
    entry["output"] = dict(columns=cols, rows=len(out), dtype="f64")
    if not os.path.exists(os.path.join(run_dir, "grain-coordinates.txt")):
        # VAMPIRE writes it only for more than one grain: a continuous film (its grain.txt, if any, holds the
        # whole film, which the output table has already)
        log("  no grain-coordinates.txt (continuous film or not granular): no grain views")
        entry["grains"] = dict(id=[], x=[], y=[], atoms=[], outline=[])
        entry["notes"] = entry.get("notes", []) + ["no grain-coordinates.txt: VAMPIRE writes it only for granular "
                                                   "systems, so this run has atom views only"]
        return None, None
    gc = vp.read_grain_coordinates(run_dir)
    outlines = vp.grain_outlines(gc, vp.read_grain_shapes(run_dir))
    entry["grains"] = dict(
        id=gc[:, 0].astype(int).tolist(), x=np.round(gc[:, 1], 3).tolist(), y=np.round(gc[:, 2], 3).tolist(),
        atoms=gc[:, 3].astype(int).tolist(),
        outline=[None if o is None else np.round(o, 3).tolist() for o in outlines])
    if not all(os.path.exists(os.path.join(d, "grain.txt")) for d in dirs):
        return None, None
    log("  parsing grain.txt")
    opts = opts or {}
    try:
        return _stream_grains(dirs, out_dir, entry, len(gc), dt, opts, check_steps)
    except ValueError as err:
        for f in ("grains.i16", "grains.f32"):
            with contextlib.suppress(OSError):
                os.remove(os.path.join(out_dir, f + ".tmp"))
        log(f"  WARNING: grain.txt not used, so no grain magnetisation views: {err}")
        entry["notes"] = entry.get("notes", []) + [f"grain.txt not used: {err}"]
        return None, None


class GrainAccumulator:
    """What the write analysis needs from the grain series, gathered frame by frame so the series itself never has
    to be in memory (float32, as the analysis has always used). A grain's decision depends on the sign of its final
    m_z, unknown until the last frame, so both signs are followed: last frame below the threshold in that direction,
    and |m| on the frame after it."""

    def __init__(self, n_grains, decision):
        self.decision = decision
        self.n = 0
        self.final = np.zeros(n_grains, np.float32)
        self.last_m = np.zeros(n_grains, np.float32)
        self.min_m = np.full(n_grains, np.inf, np.float32)
        self.last = np.full((2, n_grains), -1, np.int64)            # [0]: direction +, [1]: direction -
        self.after = np.zeros((2, n_grains), np.float32)

    def add(self, m):
        mz, ml = m[:, 2] * m[:, 3], m[:, 3]
        np.minimum(self.min_m, ml, out=self.min_m)
        for i, v in enumerate((mz, -mz)):
            arrived = self.last[i] == self.n - 1
            self.after[i][arrived] = ml[arrived]
            self.last[i][v < self.decision] = self.n
        self.final, self.last_m = mz, ml
        self.n += 1

    def outcome(self):
        """(final mz, lowest |m|, any frame below the threshold, decision frame, |m| at that frame), per grain."""
        final = self.final
        d = np.where(final < 0, 1, 0)
        i = np.arange(len(final))
        last = np.where(final == 0, self.n - 1, self.last[d, i])    # direction 0: m_z x 0 is below every threshold
        last = np.where(np.isnan(final), -1, last)
        dec = np.minimum(np.where(last >= 0, last + 1, 0), self.n - 1)
        m_dec = np.where(last + 1 >= self.n, self.last_m, self.after[d, i])
        return final, self.min_m, last >= 0, dec, m_dec

    ARRAYS = ("final", "last_m", "min_m", "last", "after")


def _stream_grains(dirs, out_dir, entry, n_grains, dt, opts, check_steps):
    """Write the grain series from the runs' grain.txt, one frame in memory at a time. Returns (lead columns [frames,
    2], {accumulator, frames at check_steps}), or raises ValueError for an unusable grain.txt."""
    every = opts.get("grain_every", 1)
    f32 = opts.get("grain_precision", 16) == 32
    name = os.path.join(out_dir, "grains.f32" if f32 else "grains.i16")
    th = dict(THRESHOLDS, **(opts.get("thresholds") or {}))
    acc = GrainAccumulator(n_grains, th["decision"])
    stops = _resume_steps(dirs)[1:] + [None]
    lead, frames, count, held = [], {}, 0, None

    def keep(row_lead, m):
        lead.append(row_lead)
        # stored from grain.txt's values (float64): quantising a float32 copy could miss the half-step bound
        (m.astype(np.float32) if f32 else np.rint(np.clip(m, -1, 1) * 32767).astype(np.int16)).tofile(out)
        m32 = m.astype(np.float32)                              # the analysis and checks work on float32, as before
        acc.add(m32)
        step = int(np.rint(row_lead[0] / dt))
        if step in check_steps and step not in frames:
            frames[step] = m32

    with open(name + ".tmp", "wb") as out:
        for d, stop in zip(dirs, stops):
            names, rows = vp.grain_rows(d, n_grains)
            ti, si, ci = (names.index(k) if k in names else None for k in ("time", "step", "temperature"))
            if ti is None and si is None:   # neither grain:real-time nor grain:time-steps: frames every output-rate steps
                if len(dirs) > 1:
                    raise SystemExit(f"{dirs}: grain.txt needs grain:real-time or grain:time-steps to join the jobs")
                rate = vp.Input(d).num("grain:output-rate", 1, "int")
                entry["grain_time_note"] = "grain.txt has no time column: frame times assumed from grain:output-rate"
            for i, v in enumerate(rows):
                t = v[ti] if ti is not None else v[si] * dt if si is not None else (i + 1) * rate * dt
                if stop is not None and t / dt > stop + 0.5:
                    continue
                row_lead = (t, v[ci] if ci is not None else np.nan)
                m = v[len(names):].reshape(n_grains, 4)
                if count % every == 0:      # thinned frames: the last one is kept too (below); the analysis uses the kept frames
                    keep(row_lead, m)
                    held = None
                else:
                    held = (row_lead, m)
                count += 1
        if held:
            keep(*held)
    os.replace(name + ".tmp", name)
    lead = np.array(lead, np.float64).reshape(-1, 2)
    lead.tofile(os.path.join(out_dir, "grain_lead.f64"))
    entry["grain_series"] = dict(frames=len(lead), lead=["time", "temperature"], columns=names)
    if f32:
        entry["grain_series"]["dtype"] = "f32"
    if every > 1:
        entry["grain_series"]["every"] = every
    if len(lead) * n_grains * 8 * (2 if f32 else 1) > GRAIN_WHOLE_MAX:
        write_grain_overview(name, out_dir, entry, n_grains, f32)
    return lead, dict(acc=acc, frames=frames)


GRAIN_WHOLE_MAX = 512e6    # above this size (bytes) the viewer reads the grain series in blocks, not all at once
OVERVIEW_BINS = 1024       # time bins of the overview of such a series


def write_grain_overview(name, out_dir, entry, n_grains, f32):
    """grains_overview.i16 [bins, grains, 4]: lowest and highest m_z = mz x |m|, lowest and highest |m| in each time bin of
    a long grain series, x 32767, so the viewer's whole-run views need not read it all and miss no short event."""
    frames = entry["grain_series"]["frames"]
    dtype = np.float32 if f32 else np.int16
    scale = 1.0 if f32 else 1.0 / 32767
    bins = min(OVERVIEW_BINS, frames)
    row_bytes = n_grains * 4 * np.dtype(dtype).itemsize
    with open(os.path.join(out_dir, "grains_overview.i16"), "wb") as f:
        for b in range(bins):
            f0, f1 = b * frames // bins, (b + 1) * frames // bins
            part = np.fromfile(name, dtype, (f1 - f0) * n_grains * 4, offset=f0 * row_bytes).reshape(f1 - f0, n_grains, 4)
            mz, ml = part[:, :, 2] * scale * part[:, :, 3] * scale, part[:, :, 3] * scale
            row = np.stack([mz.min(0), mz.max(0), ml.min(0), ml.max(0)], axis=1)
            np.rint(np.clip(row, -1, 1) * 32767).astype(np.int16).tofile(f)
    entry["grain_series"]["overview"] = dict(bins=bins)
    log(f"  grain series {frames * row_bytes / 1e9:.1f} GB: the viewer reads it in blocks; overview of {bins} time bins written")


def reduce_dipole(run_dir, out_dir, entry):
    """Grid dipole-field.txt by cell ID: macrocells are numbered x-major, z fastest, on a grid of
    (L / cell size + 1) cells per axis. Its x, y, z columns are atom centroids, not cell centres."""
    d = vp.read_dipole_cells(run_dir)
    inp = vp.Input(run_dir)
    cell = []
    for a in "xyz":   # cells:macro-cell-size sets all three axes, -x/-y/-z one; the last line wins (src/cells)
        keys = [k for k in inp.keys if k in ("cells:macro-cell-size", f"cells:macro-cell-size-{a}")]
        cell.append(inp.num(keys[-1]) if keys else 1.0)                      # nm; VAMPIRE's default is 10 A
    n = [int(round(L / c)) + 1 for L, c in zip(entry["size_nm"], cell)]     # nx, ny, nz
    if n[0] * n[1] * n[2] != len(d) or not np.array_equal(d[:, 0], np.arange(len(d))):
        raise ValueError(f"dipole-field.txt has {len(d)} cells, expected {n[0]} x {n[1]} x {n[2]} in ID order")
    # file order is (x, y, z) with z fastest: reshape to [nx, ny, nz], store as [nz, ny, nx]
    B = d[:, 6:9].reshape(n[0], n[1], n[2], 3).transpose(3, 2, 1, 0).astype(np.float32)
    A = np.minimum(d[:, 1], 65535).reshape(n[0], n[1], n[2]).transpose(2, 1, 0).astype(np.uint16)
    np.ascontiguousarray(B).tofile(os.path.join(out_dir, "dipole.f32"))
    np.ascontiguousarray(A).tofile(os.path.join(out_dir, "dipole_atoms.u16"))
    entry["dipole"] = dict(nx=n[0], ny=n[1], nz=n[2], cell_nm=cell[0])
    if len(set(cell)) > 1:
        entry["dipole"]["cell_nm_xyz"] = cell
    # shortcut: atomistic_dipole_field.txt (per-atom read-back field) is not reduced yet; add it when the
    # cell field's 1 nm resolution is not enough for read-back analysis


# ------------------------------------------------------------------------------------------------ write analysis
def write_analysis(h, entry, lead, series, Tc=None, thresholds=None):
    """Per-grain outcome of the write: intended polarity, final state, decision moment and its circumstances."""
    th = dict(THRESHOLDS, **(thresholds or {}))
    DECISION, HEATED, SETTLED, FULL_FIELD = (th[k] for k in ("decision", "heated", "settled", "full_field"))
    g = entry["grains"]
    gx, gy = np.array(g["x"]), np.array(g["y"])
    steps = np.rint(lead[:, 0] / h.dt).astype(np.int64)
    B = h.field(steps) * h.field_uz                                  # z field inside the box
    final, minM, any_below, dec, m_dec = series["acc"].outcome()
    sgn = np.sign(final)
    heated = (minM < HEATED) & any_below
    # laser temperature at the grain centres: at the decision frame (the frame after the last one below the threshold
    # in the final direction) and the highest over the run, a block of frames at a time
    Tdec = h.temperature(gx, gy, steps[dec])
    peak = np.full(len(gx), -np.inf)
    for a in range(0, len(steps), 4096):
        np.maximum(peak, h.temperature(gx[None, :], gy[None, :], steps[a:a + 4096, None]).max(0), out=peak)

    if Tc is None:   # effective freezing temperature: median centre temperature at the decision of heated grains
        sel = heated & (np.abs(final) >= SETTLED)
        Td = Tdec[sel]
        Td = Td[Td > h.Tmin + 0.05 * (h.Tmax - h.Tmin)]     # decisions under the laser, not thermal flips far from it
        Tc = float(np.median(Td)) if len(Td) else None
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
                 peakT=round(float(peak[i]), 1), minM=round(float(minM[i]), 4))
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
                     T_dec=round(float(Tdec[i]), 1), m_dec=round(float(m_dec[i]), 4),
                     lag_dec=round(float(hx - gx[i]), 3), dy=round(float(gy[i] - hy), 3))
        r["correct"] = bool(r["category"] not in ("unheated", "undecided") and want != 0 and sgn[i] == want)
        res.append(r)
    entry["analysis"] = dict(grains=res, thresholds=dict(decision=DECISION, heated=HEATED, settled=SETTLED,
                                                         full_field=FULL_FIELD))
    cats = {}
    for r in res:
        cats[r["category"]] = cats.get(r["category"], 0) + 1
    written = [r for r in res if r["intended"] != 0 and r["category"] != "unheated"]
    log(f"  Tc {Tc or 0:.0f} K, lag {lag:.2f} nm; categories {cats}; correct {sum(r['correct'] for r in res)}"
        f" of {len(written)} heated grains inside the bit windows")


def share_grain_map(out, manifest):
    """Runs without coordinate output (e.g. a read that wrote no config) take the grain map and centres of the
    nearest run of the same structure (same size and grains); only then fall back to the approximate Voronoi
    outlines."""
    runs = manifest["runs"]
    for i, e in enumerate(runs):
        if "grain_map" in e:
            continue
        near = sorted(range(len(runs)), key=lambda j: (abs(j - i), j > i))
        donor = next((runs[j] for j in near if "grain_map" in runs[j] and "from_run" not in runs[j]["grain_map"]
                      and runs[j]["size_nm"] == e["size_nm"] and runs[j]["grains"]["id"] == e["grains"]["id"]), None)
        if donor is None:
            continue
        shutil.copy(os.path.join(out, donor["id"], "grain_map.i16"), os.path.join(out, e["id"], "grain_map.i16"))
        e["grain_map"] = dict(donor["grain_map"], from_run=donor["id"])
        for k in ("x", "y", "magnetic_atoms"):
            e["grains"][k] = donor["grains"][k]
        e["grains"].pop("outline", None)
    _write_manifest(out, manifest)


def _write_manifest(out, manifest):
    tmp = os.path.join(out, "manifest.json.tmp")
    with open(tmp, "w") as f:
        json.dump(manifest, f, separators=(",", ":"))
    os.replace(tmp, os.path.join(out, "manifest.json"))


def bundle_runs(man):
    """[(run id, entry)] of a manifest, in run order; version 1 bundles had the fixed stages relax, write, read."""
    if "runs" in man:
        return [(e["id"], e) for e in man["runs"]]
    return [(s, man["stages"][s]) for s in ("relax", "write", "read") if s in man["stages"]]


def _prepare_out(out, request):
    """Make out ready for a bundle; True if an interrupted reduce of the same runs with the same options carries on.
    Otherwise an existing bundle's run folders and manifest are removed (so no files of an earlier reduction stay
    behind); any other non-empty directory is refused."""
    state = os.path.join(out, WORK, "state.json")
    if os.path.exists(state) and json.load(open(state)) == request:
        return True
    if os.path.isdir(out) and os.listdir(out):
        mp = os.path.join(out, "manifest.json")
        work = os.path.join(out, WORK)
        if not os.path.exists(mp) and not os.path.isdir(work):
            raise SystemExit(f"{out} is not empty and holds no hamrvis bundle: choose a new or empty directory")
        if os.path.isdir(work):                           # an unfinished reduce's folder holds only what it wrote
            for name in os.listdir(out):
                shutil.rmtree(os.path.join(out, name)) if os.path.isdir(os.path.join(out, name)) else os.remove(os.path.join(out, name))
        elif os.path.exists(mp):
            for rid, _ in bundle_runs(json.load(open(mp))):
                if re.fullmatch(r"[A-Za-z0-9._-]+", rid) and rid not in (".", ".."):   # never leave the bundle
                    shutil.rmtree(os.path.join(out, rid), ignore_errors=True)
            os.remove(mp)
    os.makedirs(os.path.join(out, WORK), exist_ok=True)
    _save_json(state, request)
    return False


# ------------------------------------------------------------------------------------------------ driver
def _tables_stage(r, sd, e, opts, work):
    """reduce_tables, or its saved results if an interrupted reduce got past it (grain.txt may take long to parse)."""
    state = os.path.join(work, "tables.json")
    if os.path.exists(state):
        s = json.load(open(state))
        e.clear()
        e.update(s["entry"])
        log("  grain.txt already parsed by the interrupted reduce")
        if not s["series"]:
            return None, None
        with np.load(os.path.join(work, "tables.npz")) as z:
            acc = GrainAccumulator(len(e["grains"]["id"]), s["decision"])
            acc.n = s["n"]
            for k in acc.ARRAYS:
                setattr(acc, k, z[k])
            frames = dict(zip(z["frame_steps"].tolist(), z["frames"]))
            return z["lead"], dict(acc=acc, frames=frames)
    check_steps = set()
    if opts.get("snapshots", True) and vp.coord_files(r.dir):
        for _, t, *_ in run_snapshots(r.dirs, e["dt"]):
            step = int(round(t / e["dt"]))
            check_steps.update((step, step - 1))
    lead, series = reduce_tables(r.dirs, sd, e, opts, check_steps)
    if series:
        acc = series["acc"]
        steps = sorted(series["frames"])
        frames = np.array([series["frames"][k] for k in steps], np.float32).reshape(len(steps), len(e["grains"]["id"]), 4)
        with open(os.path.join(work, "tables.tmp"), "wb") as f:
            np.savez(f, lead=lead, frame_steps=np.array(steps, np.int64), frames=frames,
                     **{k: getattr(acc, k) for k in acc.ARRAYS})
        os.replace(os.path.join(work, "tables.tmp"), os.path.join(work, "tables.npz"))
    _save_json(state, dict(entry=e, series=bool(series), n=series["acc"].n if series else 0,
                           decision=series["acc"].decision if series else 0))
    return lead, series


def reduce(paths, out, volume=True, every=1, Tc=None, workers=None, name=None, layer_factor=None, order=None,
           opts=None):
    """opts (all optional, defaults give today's bundle): precision 8|16|32 (snapshot grids), grain_precision
    16|32, grain_every, snapshots (False: leave the spin snapshots out), top_max, layer_max, thresholds
    {decision, heated, settled, full_field} of the write analysis, mem (bytes: memory budget for the snapshot
    workers; default what is free). Run again with the same arguments after an interruption to carry on."""
    opts = opts or {}
    runs = rn.find(paths, order)
    if not runs:
        raise SystemExit(f"no run directories (folders with an 'input' file) in {paths}")
    dirs = [os.path.abspath(d) for r in runs for d in r.dirs]
    # anything that changes the bundle; the same request again carries on from where an interrupted reduce stopped
    request = dict(version=FORMAT_VERSION, hamrvis_version=__version__, dirs=dirs, every=every, volume=volume, Tc=Tc, name=name, layer_factor=layer_factor,
                   opts={k: v for k, v in opts.items() if k not in ("workers", "mem")},
                   inputs=[[os.path.getsize(f), int(os.path.getmtime(f))] for f in
                           (os.path.join(d, vp.input_name(d)) for d in dirs)])
    resuming = _prepare_out(out, request)
    mp = os.path.join(out, "manifest.json")
    if resuming and os.path.exists(mp):
        manifest = json.load(open(mp))
    else:
        manifest = dict(format="hamrvis-bundle", version=FORMAT_VERSION, hamrvis_version=__version__, created=time.strftime("%Y-%m-%d %H:%M:%S"),
                        name=name or os.path.basename(os.path.commonpath(dirs) if len(dirs) > 1 else dirs[0]), runs=[])
    if resuming:
        log(f"carrying on the interrupted reduce in {out}: {len(manifest['runs'])} of {len(runs)} runs were complete")
    finished = {e["id"] for e in manifest["runs"]}
    last_hamr = None
    for r in runs:
        if r.id in finished:
            if r.type == "hamr":
                last_hamr = r.id
            continue
        log(f"{r.id} ({r.type}): " + ", ".join(r.dirs))
        sd = os.path.join(out, r.id)
        os.makedirs(sd, exist_ok=True)
        work = os.path.join(out, WORK, r.id)
        os.makedirs(work, exist_ok=True)
        run = r.dir
        inp = vp.Input(run)
        e = dict(id=r.id, type=r.type, title=r.title, source=os.path.abspath(run),
                 segments=[os.path.abspath(d) for d in r.dirs],
                 input={k: v for k, (v, u) in inp.keys.items()},
                 input_units={k: u for k, (v, u) in inp.keys.items() if u},
                 dt=inp.num("sim:time-step", 1e-15, "time"), monte_carlo=rn.is_monte_carlo(inp),
                 size_nm=inp.system_size(),
                 magnetic_types=magnetic_types(run))
        if r.type != "hamr" and last_hamr:
            e["ref"] = last_hamr                       # the write whose bits this run is compared with
        for f in ("log", "stdout", "time.txt"):
            if os.path.exists(os.path.join(run, f)):
                e.setdefault("files", []).append(f)
        meta = os.path.join(run, "atoms-coords.meta")
        if os.path.exists(meta):
            e["n_atoms"] = int(vp._meta_files(meta)[0].get("Number of atoms", "0").split()[0])
        e["source_files"] = {f: [os.path.getsize(os.path.join(run, f)), int(os.path.getmtime(os.path.join(run, f)))]
                             for f in (inp.name, vp.output_name(run), "grain.txt", "dipole-field.txt")
                             if f and os.path.exists(os.path.join(run, f))}
        lead, series = _tables_stage(r, sd, e, opts, work)
        reduce_atoms(r.dirs, sd, e, series, volume, every, workers, layer_factor, opts, work)
        if r.type == "hamr":
            h = vp.Hamr(run, r.dirs)
            e["hamr"] = h.params()
            if lead is not None:
                write_analysis(h, e, lead, series, Tc, opts.get("thresholds"))
            last_hamr = r.id
        last = r.dirs[-1]                              # field maps are written at the end of a run
        if os.path.exists(os.path.join(last, "dipole-field.txt")):
            log("  gridding dipole-field.txt")
            reduce_dipole(last, sd, e)
        manifest["runs"].append(e)
        _write_manifest(out, manifest)                 # after each run, so a partial bundle opens
        shutil.rmtree(work)
    share_grain_map(out, manifest)
    shutil.rmtree(os.path.join(out, WORK))
    size = sum(os.path.getsize(os.path.join(d, f)) for d, _, fs in os.walk(out) for f in fs)
    log(f"bundle {out}: {size / 1e6:.1f} MB")

"""Readers for VAMPIRE output files and the HAMR head model.

Lengths are in nm, times in s, fields in T, temperatures in K unless a name says otherwise.
Everything here only reads files; the bundle layout lives in reduce.py.
"""
import glob
import math
import os
import re

import numpy as np

# unit tag -> factor to nm, s, m/s or T (VAMPIRE's defaults are Angstrom, s, T)
_UNITS = {"nm": 1.0, "A": 0.1, "m": 1e9, "m/s": 1.0, "s": 1.0, "ns": 1e-9, "ps": 1e-12, "fs": 1e-15, "T": 1.0,
          "mT": 1e-3, "muB": 1.0, "K": 1.0}


def read_keys(path):
    """{key: (value, unit or None)} for every 'key = value !unit' line of a VAMPIRE input or material file.
    Keys without a value (flags) map to ("", None). Later lines win, as in VAMPIRE."""
    keys = {}
    for line in open(path):
        line = line.split("#")[0].strip()
        if not line:
            continue
        k, _, v = (s.strip() for s in line.partition("="))
        v, _, unit = (s.strip() for s in v.partition("!"))
        keys[k.replace(" ", "")] = (v, unit or None)
    return keys


class Input:
    """A run's input file with typed access."""

    def __init__(self, run_dir):
        self.keys = read_keys(os.path.join(run_dir, "input"))

    def has(self, key):
        return key in self.keys

    def str(self, key, default=None):
        return self.keys[key][0] if key in self.keys else default

    def num(self, key, default=None, default_unit="A"):
        """Value converted to nm / s / m/s / T; default_unit is VAMPIRE's unit for the key when no '!unit' is given."""
        if key not in self.keys:
            return default
        v, unit = self.keys[key]
        return float(v) * _UNITS[unit or default_unit]


# ----------------------------------------------------------------------------------------------- run classification
def stage_of(run_dir):
    """'relax', 'write' or 'read' from the input: HAMR program = write, dipole output = read, anything else = relax."""
    inp = Input(run_dir)
    if inp.str("sim:program") == "hamr-simulation":
        return "write"
    if inp.has("dipole:output-dipole-field") or inp.has("dipole:output-atomistic-dipole-field"):
        return "read"
    return "relax"


# ----------------------------------------------------------------------------------------------- text tables
def read_output(run_dir):
    """(column names, float64 [rows, cols]) of the 'output' file; needs output:column-headers."""
    lines = [l for l in open(os.path.join(run_dir, "output")) if not l.startswith("#") and l.strip()]
    head = lines[0].split()
    if re.fullmatch(r"[-+.\deE]+", head[0]):
        raise ValueError(f"{run_dir}/output has no header row: add output:column-headers to the input")
    data = np.loadtxt(lines[1:], ndmin=2)
    return head, data


def read_grain_coordinates(run_dir):
    """float64 [grains, 5]: id, x (nm), y (nm), atom count, moment (muB)."""
    g = np.loadtxt(os.path.join(run_dir, "grain-coordinates.txt"), comments="#", ndmin=2)
    g[:, 1:3] /= 10
    return g


def read_grain_shapes(run_dir):
    """List of [n, 2] outlines (nm) of every Voronoi cell, in file order (see grain_outlines)."""
    path = os.path.join(run_dir, "grain_shapes.txt")
    if not os.path.exists(path):
        return []
    return [np.array(l.split()[1:], float).reshape(-1, 2) / 10 for l in open(path) if l.strip()]


def _contains(poly, x, y):
    px, py = poly[:, 0], poly[:, 1]
    qx, qy = np.roll(px, -1), np.roll(py, -1)
    cross = ((py > y) != (qy > y)) & (x < (qx - px) * (y - py) / np.where(qy == py, 1e-12, qy - py) + px)
    return bool(np.count_nonzero(cross) % 2)


def grain_outlines(coords, shapes):
    """One outline per grain. grain_shapes.txt lists every Voronoi cell made before a geometry cut, so its order is
    not the grain order once cells lose all their atoms: each grain takes the cell containing its centre (or the
    nearest cell centre). Periodic grains whose centre lies outside every cell also fall back to the nearest."""
    if not shapes:
        return [None] * len(coords)
    cen = np.array([s.mean(axis=0) for s in shapes])
    res = []
    for x, y in coords[:, 1:3]:
        inside = [s for s in shapes if _contains(s, x, y)]
        res.append(inside[0] if len(inside) == 1 else shapes[int(np.argmin(np.hypot(*(cen - (x, y)).T)))])
    return res


# grain.txt columns per grain:* keyword (src/vio/match_grain_list.cpp), written in the order of the input file
_GRAIN_LEAD = {"time-steps": ["step"], "real-time": ["time"], "temperature": ["temperature"],
               "electron-temperature": ["temperature"], "phonon-temperature": ["phonon_temperature"],
               "applied-field-strength": ["field"], "applied-field-unit-vector": ["field_ux", "field_uy", "field_uz"]}


def read_grain_file(run_dir, n_grains):
    """Per-frame grain data from grain.txt.

    The columns follow the input's grain:* keywords in order: single columns (step, time, temperature, field) and
    then blocks per grain, of which grain:magnetisation gives 4 per grain (mx, my, mz direction, |m|).
    Returns ({name: float64 [frames]} of the columns before the magnetisation block, float32 [frames, grains, 4])."""
    keys = [k.split(":", 1)[1] for k in Input(run_dir).keys if k.startswith("grain:")]
    names = []
    for k in keys:
        if k == "magnetisation":
            break
        if k in _GRAIN_LEAD:
            names += _GRAIN_LEAD[k]
        elif k != "output-rate":
            raise ValueError(f"grain:{k} comes before grain:magnetisation in the input; its column count is not "
                             f"known to hamrvis (add it to _GRAIN_LEAD in hamrvis/vampire.py)")
    else:
        raise ValueError("the input has no grain:magnetisation, so grain.txt holds no grain magnetisation")
    path = os.path.join(run_dir, "grain.txt")
    cols = list(range(len(names) + 4 * n_grains))
    d = np.loadtxt(path, comments="#", usecols=cols, ndmin=2, dtype=np.float64)
    lead = {n: d[:, i] for i, n in enumerate(names)}
    m = d[:, len(names):].reshape(len(d), n_grains, 4).astype(np.float32)
    return lead, m


def read_dipole_cells(run_dir):
    """dipole-field.txt: float64 [cells, 9] = id, atoms in cell, x, y, z (A), moment, Bx, By, Bz (T)."""
    return np.loadtxt(os.path.join(run_dir, "dipole-field.txt"), ndmin=2)


# ----------------------------------------------------------------------------------------------- config snapshots
def _meta_files(meta_path):
    """(header dict, data file names) of a config .meta file; data files are listed one per line after
    'Number of ... files: N' (v5+ format) or default to the .meta name with .data (single-file output)."""
    text = open(meta_path).read()
    head = dict(re.findall(r"^([A-Za-z][A-Za-z ]*?):\s*(.*)$", text, re.M))
    m = re.search(r"Number of (?:spin|coord) files:\s*(\d+)", text)
    d = os.path.dirname(meta_path)
    if m:
        n = int(m.group(1))
        after = text[m.end():].split()
        names = [w for w in after if w.endswith(".data")][:n]
    else:
        names = [os.path.basename(meta_path)[:-4] + "data"]
    return head, [os.path.join(d, f) for f in names]


def _is_binary(path):
    with open(path, "rb") as f:
        b = f.read(64)
    return b"\0" in b


def coord_files(run_dir):
    """Data files of the atom coordinates (in snapshot file order), or [] if the run wrote none."""
    meta = os.path.join(run_dir, "atoms-coords.meta")
    if os.path.exists(meta):
        return _meta_files(meta)[1]
    old = os.path.join(run_dir, "atoms-coords.data")
    return [old] if os.path.exists(old) else []


def read_coords(path):
    """(type int32 [n], xyz float64 [n, 3] in A, grain int32 [n] or None) of one coordinate data file."""
    if _is_binary(path):
        with open(path, "rb") as f:
            n = int(np.fromfile(f, np.uint64, 1)[0])
            typ = np.fromfile(f, np.int32, n)
            np.fromfile(f, np.int32, n)                     # category
            xyz = np.fromfile(f, np.float64, 3 * n).reshape(n, 3)
            grain = np.fromfile(f, np.int32, n)
        return typ, xyz, (grain if len(grain) == n else None)
    d = np.loadtxt(path, skiprows=1, ndmin=2)
    return d[:, 0].astype(np.int32), d[:, 2:5], (d[:, 5].astype(np.int32) if d.shape[1] > 5 else None)


def read_spins(path):
    """float64 [n, 3] spin directions of one spin data file."""
    if _is_binary(path):
        with open(path, "rb") as f:
            n = int(np.fromfile(f, np.uint64, 1)[0])
            s = np.fromfile(f, np.float64, 3 * n)
        if len(s) != 3 * n:
            raise IOError(f"{path}: truncated ({len(s) // 3} of {n} spins), copy caught mid-transfer?")
        return s.reshape(n, 3)
    s = np.loadtxt(path, skiprows=1, ndmin=2)
    return s


def snapshots(run_dir):
    """[(index, time s, [data files])] of the spin snapshots, in order. Snapshot files are found by probing
    spins-00000000.meta, ...01.meta, ... (directory listings of copied folders have been seen to miss files),
    falling back to a glob for runs that do not start at 0 (restarts)."""
    res, i = [], 0
    metas = []
    while os.path.exists(p := os.path.join(run_dir, f"spins-{i:08d}.meta")):
        metas.append(p)
        i += 1
    if not metas:
        metas = sorted(glob.glob(os.path.join(run_dir, "spins-[0-9]*.meta")))
    for p in metas:
        head, files = _meta_files(p)
        idx = int(re.search(r"spins-(\d+)", os.path.basename(p)).group(1))
        res.append((idx, float(head["Time"].split()[0]), files))
    return res


# ----------------------------------------------------------------------------------------------- HAMR head model
class Hamr:
    """Head position, laser temperature and field box of a hamr-simulation run, as in src/hamr/hamr_continuous.cpp.
    Steps count from the start of the run's sim::time (equilibration included)."""

    def __init__(self, run_dir):
        inp = Input(run_dir)
        n = inp.num
        self.dt = n("sim:time-step", 1e-15, "s")
        self.eq = int(n("sim:equilibration-time-steps", 0, "s"))
        self.Tmin = n("sim:minimum-temperature", 0.0, "K")
        self.Tmax = n("sim:maximum-temperature", 0.0, "K")
        # the head field ramps by (Hmax - Hmin) in each bit (src/hamr/field_profiles.cpp), so its plateau is
        # Hmax - Hmin, with VAMPIRE's defaults Hmax = +1 T, Hmin = -1 T (i.e. 2 T) unless the input sets them
        self.Hmin = n("sim:minimum-applied-field-strength", -1.0, "T")
        self.Hmax = n("sim:maximum-applied-field-strength", 1.0, "T")
        self.Hpeak = self.Hmax - self.Hmin
        uv = [float(c) for c in inp.str("sim:applied-field-unit-vector", "0,0,1").split(",")]
        self.field_uz = uv[2] / math.sqrt(sum(c * c for c in uv))
        self.rise = int(round(n("hamr:field-rise-time", 1e-12, "s") / self.dt))
        self.fall = int(round(n("hamr:field-fall-time", 1e-12, "s") / self.dt))
        self.Lx, self.Ly, self.Lz = (n(f"dimensions:system-size-{a}") for a in "xyz")
        self.v = n("hamr:head-speed", 30.0, "m/s")                      # m/s = nm/ns
        self.BL = n("hamr:bit-size", 0.0)
        self.TW = n("hamr:track-size", 0.0)
        self.pad = n("hamr:track-padding", 0.0)
        self.NPS = n("hamr:NPS", n("hamr:NFT-to-pole-spacing", 0.0))
        self.fwhm_x, self.fwhm_y = n("hamr:laser-FWHM-x", 20.0), n("hamr:laser-FWHM-y", 20.0)
        self.box_x = min(n("hamr:head-field-x", 20.0), self.Lx)
        self.box_y = min(n("hamr:head-field-y", 20.0), self.Ly)
        self.sigma_x = self.fwhm_x / math.sqrt(8 * math.log(2))
        self.sigma_y = self.fwhm_y / math.sqrt(8 * math.log(2))
        # tracks, bits and sequence as src/hamr/initialize.cpp and bit_sequence.cpp decide them (lengths in A there)
        self.bits_per_track = int(math.floor((self.Lx - 0.1) / self.BL))
        max_tracks = int(math.floor((self.Ly - 0.1 - 2 * self.pad) / self.TW))
        n_bits = int(n("hamr:number-of-bits", 0, "s"))
        try_tracks = int(round(n_bits / self.bits_per_track)) if self.bits_per_track else 0
        self.n_tracks = try_tracks if try_tracks < max_tracks else max_tracks
        if inp.str("hamr:bit-sequence-type", "user-defined") == "single-tone-predefined" or not inp.has("hamr:bit-sequence"):
            seq = [-1 if i % 2 == 0 else 1 for i in range(n_bits)]
        else:
            seq = [int(v) for v in re.findall(r"-?\d+", inp.str("hamr:bit-sequence"))]
            seq = seq[:n_bits] if n_bits < len(seq) else seq
        total = self.n_tracks * self.bits_per_track
        self.sequence = (seq + [0] * total)[:total]                    # trimmed, or padded with field-free bits
        vdt = self.v * 1e9 * self.dt                                     # nm per step
        self.step_nm = vdt
        self.nb = int(round(self.BL / vdt))
        self.ne = int(round(self.BL / 2 / vdt)) + int(round(self.NPS / vdt))
        self.nt = self.bits_per_track * self.nb + 2 * self.ne
        self.laser_steps = self.n_tracks * self.nt

    def params(self):
        """Everything the viewer needs to rebuild the head, as plain numbers."""
        keys = ("dt", "eq", "Tmin", "Tmax", "Hmin", "Hmax", "Hpeak", "field_uz", "rise", "fall", "Lx", "Ly", "Lz", "v", "BL", "TW", "pad", "NPS", "fwhm_x",
                "fwhm_y", "box_x", "box_y", "sigma_x", "sigma_y", "step_nm", "nb", "ne", "bits_per_track", "nt",
                "n_tracks", "laser_steps")
        p = {k: getattr(self, k) for k in keys}
        p["sequence"] = self.sequence
        return p

    def head(self, steps):
        """Laser centre x, y (nm) and laser-on flag at the given steps."""
        k = np.asarray(steps, np.int64) - self.eq
        on = (k >= 1) & (k <= self.laser_steps)
        kk = np.clip(k, 1, None)
        track, j = (kk - 1) // self.nt, (kk - 1) % self.nt + 1
        return -self.BL / 2 - self.NPS + (j + 1) * self.step_nm, self.TW * (track + 0.5) + self.pad, on

    def temperature(self, x, y, steps):
        hx, hy, on = self.head(steps)
        g = np.exp(-(x - hx) ** 2 / (2 * self.sigma_x ** 2) - (y - hy) ** 2 / (2 * self.sigma_y ** 2))
        return self.Tmin + np.where(on, (self.Tmax - self.Tmin) * g, 0.0)

    def in_field_box(self, x, y, steps):
        hx, hy, on = self.head(steps)
        return on & (np.abs(x - (hx - self.NPS)) <= self.box_x / 2) & (np.abs(y - hy) <= self.box_y / 2)

    def field(self, steps):
        """Signed head field B_applied (T) at the given steps, as VAMPIRE's output column reports it: plateau Hpeak,
        linear ramps of rise / fall steps inside each bit, sign of the bit. The z field is field() * field_uz."""
        k = np.asarray(steps, np.int64) - self.eq
        on = (k >= 1) & (k <= self.laser_steps)
        j = (np.clip(k, 1, None) - 1) % self.nt                       # 0-based step within the track
        b, i = np.divmod(j - self.ne, self.nb)                        # bit within track, 0-based step within bit
        inbit = on & (j >= self.ne) & (b < self.bits_per_track)
        mag = np.minimum((i + 1) / self.rise, 1.0)
        mag = np.where(i >= self.nb - self.fall, (self.nb - 1 - i) / self.fall, mag) * self.Hpeak
        track = (np.clip(k, 1, None) - 1) // self.nt
        seq = np.asarray(self.sequence[: self.n_tracks * self.bits_per_track] or [0])
        pol = seq[np.clip(track * self.bits_per_track + b, 0, len(seq) - 1)]
        return np.where(inbit, mag * pol, 0.0)

    def bit_windows(self):
        """[(first step, last step, polarity)] of every bit, in write order."""
        res = []
        for b, pol in enumerate(self.sequence[: self.n_tracks * self.bits_per_track]):
            t, i = divmod(b, self.bits_per_track)
            s0 = self.eq + t * self.nt + self.ne + i * self.nb + 1
            res.append((s0, s0 + self.nb - 1, pol))
        return res

    def track_of(self, y):
        return np.clip(np.floor((np.asarray(y) - self.pad) / self.TW), 0, self.n_tracks - 1).astype(int)

    def freeze_lag(self, Tc):
        """Distance (nm) behind the laser centre where the track centre line cools through Tc."""
        return self.sigma_x * math.sqrt(2 * math.log((self.Tmax - self.Tmin) / (Tc - self.Tmin)))

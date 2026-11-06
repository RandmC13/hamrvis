"""Readers for VAMPIRE output files and the HAMR head model.

Lengths are in nm, times in s, fields in T, temperatures in K unless a name says otherwise.
Everything here only reads files; the bundle layout lives in reduce.py.
"""
import glob
import itertools
import math
import os
import re

import numpy as np

# VAMPIRE's unit table (src/utility/units.cpp): unit -> (factor to VAMPIRE's internal unit, kind). Internal units are
# A, s, T, A/s, J, J/T; a value without a unit is in the internal unit. VAMPIRE stops on an unknown unit or one of
# the wrong kind for the key, so a run that ran only used these.
_UNITS = {
    "A": (1.0, "length"), "m": (1e10, "length"), "cm": (1e8, "length"), "mm": (1e7, "length"), "um": (1e4, "length"),
    "nm": (1e1, "length"), "pm": (1e-3, "length"), "in": (2.54e8, "length"), "feet": (3.12e9, "length"),
    "J": (1.0, "energy"), "eV": (1.602176634e-19, "energy"), "meV": (1.602176634e-22, "energy"),
    "erg": (1e-7, "energy"), "Ryd": (2.17987208e-18, "energy"), "mRyd": (2.17987208e-21, "energy"),
    "Ht": (4.35974417e-18, "energy"), "zJ": (1e-21, "energy"),
    "J/T": (1.0, "moment"), "Amm": (1.0, "moment"), "erg/G": (1e-3, "moment"), "abAcmcm": (1e-3, "moment"),
    "muB": (9.2740100783e-24, "moment"), "eV/T": (1.602176634e-19, "moment"), "erg/Oe": (1e-3, "moment"),
    "A/m": (1e-30, "magnetisation"), "emu/cc": (1e-27, "magnetisation"), "J/T/AAA": (1.0, "magnetisation"),
    "J/atom": (1.0, "anisotropy"), "J/mmm": (1e-30, "anisotropy"), "erg/cc": (1e-31, "anisotropy"),
    "T": (1.0, "field"), "mT": (1e-3, "field"), "uT": (1e-6, "field"), "Oe": (1e-4, "field"), "kOe": (1e-1, "field"),
    "s": (1.0, "time"), "ms": (1e-3, "time"), "us": (1e-6, "time"), "ns": (1e-9, "time"), "ps": (1e-12, "time"),
    "fs": (1e-15, "time"), "as": (1e-18, "time"), "zs": (1e-21, "time"),
    "Ohm-m": (1.0, "resistivity"),
    "V": (1.0, "potential"), "mV": (1e-3, "potential"), "uV": (1e-6, "potential"), "nV": (1e-9, "potential"),
    "Hz": (1.0, "frequency"), "kHz": (1e3, "frequency"), "MHz": (1e6, "frequency"), "GHz": (1e9, "frequency"),
    "THz": (1e12, "frequency"),
    "J/m": (1e-10, "exchange"), "erg/cm": (1e-15, "exchange"), "erg/cm2": (1e-23, "mm_energy"),
    "A/s": (1.0, "velocity"), "m/s": (1e10, "velocity"), "cm/s": (1e8, "velocity"), "mm/s": (1e7, "velocity"),
    "um/s": (1e4, "velocity"), "nm/s": (1e1, "velocity"), "pm/s": (1e-3, "velocity"), "in/s": (2.54e8, "velocity"),
    "feet/s": (3.12e9, "velocity"), "A/ms": (1e3, "velocity"), "A/us": (1e6, "velocity"), "A/ns": (1e9, "velocity"),
    "A/ps": (1e12, "velocity"), "nm/ns": (1e8, "velocity"),
}
# hamrvis works in nm, s, T and m/s: internal VAMPIRE value / factor. "none" (temperatures in K, plain numbers)
# takes no unit; "int" keys are read as integers by VAMPIRE, which ignores any unit on them.
_OUT = {"length": 10.0, "time": 1.0, "field": 1.0, "velocity": 1e10, "none": 1.0}


def read_keys(path):
    """{key: (value, unit or None)} for every 'key = value !unit' line of a VAMPIRE input or material file.
    Keys without a value (flags) map to ("", None). Later lines win, as in VAMPIRE, and the keys are in the order
    of the line that last set them. Whitespace inside keys is ignored (VAMPIRE removes all of it)."""
    keys = {}
    for line in open(path, errors="replace"):
        line = line.split("#")[0].strip()
        if not line:
            continue
        k, _, v = (s.strip() for s in line.partition("="))
        v, _, unit = (s.strip() for s in v.partition("!"))
        k = "".join(k.split())
        keys.pop(k, None)
        keys[k] = (v, "".join(unit.split()) or None)
    return keys


# names VAMPIRE gives its main files unless run with --input-file / --output-file; set by the command line flags
NAMES = {"input": None, "output": None}


def _log_head(run_dir, n=40):
    try:
        with open(os.path.join(run_dir, "log"), errors="replace") as f:
            return [f.readline() for _ in range(n)]
    except OSError:
        return []


def input_name(run_dir):
    """File name of a run's main input: --input-file if given and present, else 'input', else the name the log says VAMPIRE
    opened ('Opening main input file "..."', written when VAMPIRE starts), or None."""
    if NAMES["input"] and os.path.isfile(os.path.join(run_dir, NAMES["input"])):   # else found as without the flag
        return NAMES["input"]
    if os.path.isfile(os.path.join(run_dir, "input")):
        return "input"
    for line in _log_head(run_dir):
        m = re.search(r'Opening main input file "(.+)"', line)
        if m and os.path.isfile(os.path.join(run_dir, m.group(1))):
            return m.group(1)
    return None


def output_name(run_dir):
    """File name of a run's 'output' table: --output-file if given and present, else 'output', else the one file of the run
    directory whose second line is VAMPIRE's '# Output file for vampire simulation' (VAMPIRE does not log the
    name of the file given with --output-file; a job continued from a checkpoint appends to it without the header,
    so there the name must be given), or None."""
    if NAMES["output"] and os.path.isfile(os.path.join(run_dir, NAMES["output"])):   # else found as without the flag
        return NAMES["output"]
    if os.path.isfile(os.path.join(run_dir, "output")):
        return "output"
    found = []
    for f in sorted(os.listdir(run_dir)):
        p = os.path.join(run_dir, f)
        if f.endswith((".data", ".meta", ".chk")) or not os.path.isfile(p):
            continue
        try:
            with open(p, "rb") as fh:
                head = fh.read(400).split(b"\n")
        except OSError:
            continue
        if len(head) > 1 and head[1].startswith(b"# Output file for vampire simulation"):
            found.append(f)
    if len(found) > 1:
        raise SystemExit(f"{run_dir}: several VAMPIRE output tables ({', '.join(found)}): choose one with --output-file")
    return found[0] if found else None


class Input:
    """A run's input file with typed access."""

    def __init__(self, run_dir):
        name = input_name(run_dir)
        if name is None:
            raise SystemExit(f"{run_dir}: no VAMPIRE input file (looked for 'input' and the name in the log; "
                             f"use --input-file NAME)")
        self.name = name
        self.keys = read_keys(os.path.join(run_dir, name))

    def has(self, key):
        return key in self.keys

    def str(self, key, default=None):
        return self.keys[key][0] if key in self.keys else default

    def num(self, key, default=None, kind="length"):
        """Value of a key VAMPIRE reads as kind ("length", "time", "field", "velocity", "none" or "int"),
        converted as VAMPIRE converts it and returned in nm / s / T / m/s (default is returned as given)."""
        if key not in self.keys:
            return default
        v, unit = self.keys[key]
        if kind == "int":
            return int(float(v))
        factor, got = _UNITS.get(unit, (None, None)) if unit else (1.0, kind)
        if factor is None:
            raise ValueError(f"{key} = {v} !{unit}: '{unit}' is not one of VAMPIRE's units (src/utility/units.cpp)")
        if got != kind:
            raise ValueError(f"{key} = {v} !{unit}: VAMPIRE reads {key} as a {kind}, but '{unit}' is a {got} unit")
        return float(v) * factor / _OUT[kind]

    def system_size(self):
        """[x, y, z] in nm: dimensions:system-size sets all three, -x/-y/-z one, the last line wins (src/create/interface.cpp);
        10 nm each unless set (mp::default_system)."""
        size = [10.0, 10.0, 10.0]
        for k in self.keys:
            if k == "dimensions:system-size":
                size = [self.num(k)] * 3
            elif k in ("dimensions:system-size-x", "dimensions:system-size-y", "dimensions:system-size-z"):
                size["xyz".index(k[-1])] = self.num(k)
        return size


def vampire_integers(text):
    """Integer list as VAMPIRE's vin::integers_from_string reads it: comma separated, each field read with
    'stream >> int', which takes a leading [+-]digits and gives 0 when the field does not start with one. So
    '{1,0,-1}' (braces are not VAMPIRE syntax) reads as 0, 0, -1, as VAMPIRE writes it."""
    return [int(m.group()) if (m := re.match(r"[+-]?\d+", f)) else 0 for f in "".join(text.split()).split(",")]


# ----------------------------------------------------------------------------------------------- text tables
def read_output(run_dir):
    """(column names, float64 [rows, cols]) of the 'output' file; needs output:column-headers."""
    name = output_name(run_dir)
    if name is None:
        raise ValueError(f"{run_dir}: no VAMPIRE output table (looked for 'output' and for VAMPIRE's header in the "
                         f"other files; use --output-file NAME)")
    lines = [l for l in open(os.path.join(run_dir, name)) if not l.startswith("#") and l.strip()]
    if not lines:
        raise ValueError(f"{run_dir}/{name} has no rows: the run stopped before its first output")
    head = lines[0].split()
    if re.fullmatch(r"[-+.\deE]+", head[0]):
        raise ValueError(f"{run_dir}/{name} has no header row: add output:column-headers to the input")
    rows = lines[1:]
    if len(rows) > 1 and len(rows[-1].split()) < len(rows[0].split()):
        print(f"WARNING: {run_dir}/{name}: last row incomplete (job stopped while writing it): left out")
        rows = rows[:-1]
    data = np.loadtxt(rows, ndmin=2)
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
    lo, hi = np.array([s.min(axis=0) for s in shapes]), np.array([s.max(axis=0) for s in shapes])
    res = []
    for x, y in coords[:, 1:3]:
        # only cells whose bounding box holds the point can contain it
        near = np.nonzero((lo[:, 0] <= x) & (x <= hi[:, 0]) & (lo[:, 1] <= y) & (y <= hi[:, 1]))[0]
        inside = [shapes[i] for i in near if _contains(shapes[i], x, y)]
        res.append(inside[0] if len(inside) == 1 else shapes[int(np.argmin(np.hypot(*(cen - (x, y)).T)))])
    return res


# grain.txt columns per grain:* keyword (src/vio/match_grain_list.cpp), written in the order of the input file
_GRAIN_LEAD = {"time-steps": ["step"], "real-time": ["time"], "temperature": ["temperature"],
               "electron-temperature": ["temperature"], "phonon-temperature": ["phonon_temperature"],
               "applied-field-strength": ["field"], "applied-field-unit-vector": ["field_ux", "field_uy", "field_uz"]}


def grain_columns(run_dir, n_grains):
    """(names of the single columns before the magnetisation block, columns needed) of grain.txt, checked against its
    first row; ValueError if grain.txt holds no usable grain magnetisation. The columns follow the input's grain:*
    keywords in order: single columns (step, time, temperature, field) and then blocks per grain, of which
    grain:magnetisation gives 4 per grain (mx, my, mz direction, |m|)."""
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
    need = len(names) + 4 * n_grains
    with open(os.path.join(run_dir, "grain.txt"), errors="replace") as f:
        first = next((l for l in f if l.strip() and not l.startswith("#")), "")
    if len(first.split()) < need:
        raise ValueError(f"grain.txt rows have {len(first.split())} columns, but the input's grain keywords and the "
                         f"{n_grains} grains of grain-coordinates.txt need at least {need}")
    return names, need


def grain_rows(run_dir, n_grains):
    """(names of the columns before the magnetisation block, iterator of float64 rows of grain.txt, each cut to the
    columns needed), one row in memory at a time: grain.txt can be tens of GB. A short or damaged last row (a job
    stopped while writing it) is left out; one anywhere else raises ValueError."""
    names, need = grain_columns(run_dir, n_grains)
    path = os.path.join(run_dir, "grain.txt")

    def parse(lines):
        return np.loadtxt(lines, usecols=range(need), ndmin=2, dtype=np.float64)

    def rows():
        with open(path, errors="replace") as f:
            while True:
                batch = list(itertools.islice(f, 16))
                if not batch:
                    return
                batch = [l for l in batch if l.strip() and not l.startswith("#")]
                if not batch:
                    continue
                try:
                    yield from parse(batch)
                except ValueError:                        # find the row, which may only be a last one left half-written
                    for k, line in enumerate(batch):
                        try:
                            v = parse([line])
                        except ValueError:
                            if any(l.strip() for l in batch[k + 1:]) or any(l.strip() for l in f):
                                raise ValueError("grain.txt has a short or damaged row (a copy caught mid-transfer?)") from None
                            print(f"WARNING: {path}: last row incomplete (job stopped while writing it): left out")
                            return
                        yield v[0]
    return names, rows()


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


CHUNK_ATOMS = 1 << 24      # atoms read at a time, so a single data file of any size fits in memory


def _text_chunks(path, cols=None, chunk=CHUNK_ATOMS):
    """Rows of a text config data file (first line: number of rows) in pieces, checked to be complete."""
    with open(path) as f:
        head = f.readline().split()
        got = 0
        try:
            n = int(head[0])
            while got < n:
                d = np.loadtxt(f, ndmin=2, usecols=range(cols) if cols else None, max_rows=chunk)
                if not len(d):
                    break
                got += len(d)
                yield d
            if got == n and f.read(4096).strip():
                got = n + 1
        except (ValueError, IndexError):
            raise IOError(f"{path}: truncated or damaged text file (unreadable rows), copy caught mid-transfer?") from None
    if got != n:
        raise IOError(f"{path}: truncated or damaged text file ({got} of {n} rows), copy caught mid-transfer?")


def coord_chunks(path, chunk=CHUNK_ATOMS):
    """Pieces (type int32 [m], xyz float64 [m, 3] in A, grain int32 [m] or None) of one coordinate data file."""
    if _is_binary(path):
        n = int(np.fromfile(path, np.uint64, 1)[0])            # n, type[n], category[n], xyz[3n], grain[n] (newer files)
        size = os.path.getsize(path)
        if size < 8 + 32 * n:
            raise IOError(f"{path}: truncated (fewer than its {n} atoms), copy caught mid-transfer?")
        grains = size >= 8 + 36 * n
        for a in range(0, n, chunk):
            m = min(chunk, n - a)
            yield (np.fromfile(path, np.int32, m, offset=8 + 4 * a),
                   np.fromfile(path, np.float64, 3 * m, offset=8 + 8 * n + 24 * a).reshape(m, 3),
                   np.fromfile(path, np.int32, m, offset=8 + 32 * n + 4 * a) if grains else None)
    else:
        for d in _text_chunks(path, chunk=chunk):
            yield d[:, 0].astype(np.int32), d[:, 2:5], (d[:, 5].astype(np.int32) if d.shape[1] > 5 else None)


def spin_chunks(path, chunk=CHUNK_ATOMS):
    """Pieces (float64 [m, 3]) of one spin data file, in file order."""
    if _is_binary(path):
        n = int(np.fromfile(path, np.uint64, 1)[0])
        if os.path.getsize(path) < 8 + 24 * n:
            got = (os.path.getsize(path) - 8) // 24
            raise IOError(f"{path}: truncated ({got} of {n} spins), copy caught mid-transfer?")
        for a in range(0, n, chunk):
            m = min(chunk, n - a)
            yield np.fromfile(path, np.float64, 3 * m, offset=8 + 24 * a).reshape(m, 3)
    else:
        yield from _text_chunks(path, 3, chunk)


def read_spins(path):
    """float64 [n, 3] spin directions of one spin data file."""
    parts = list(spin_chunks(path))
    return np.concatenate(parts) if parts else np.zeros((0, 3))


def read_window_spins(path):
    """(atoms in the file, positions int64 [m], float64 [m, 3] spins) of one window snapshot data file
    (hamr:active-window): the spins of the atoms that changed since the previous snapshot, stored as
    (first position, count) ranges of positions in the coordinate file's order."""
    with open(path, "rb") as f:
        n, r = (int(v) for v in np.fromfile(f, np.uint64, 2))
        ranges = np.fromfile(f, np.uint64, 2 * r).astype(np.int64).reshape(r, 2)
        m = int(ranges[:, 1].sum())
        s = np.fromfile(f, np.float64, 3 * m)
    if len(s) != 3 * m:
        raise IOError(f"{path}: truncated ({len(s) // 3} of {m} spins), copy caught mid-transfer?")
    pos = np.concatenate([np.arange(a, a + c) for a, c in ranges]) if r else np.zeros(0, np.int64)
    if r and pos.max() >= n:
        raise IOError(f"{path}: positions beyond its {n} atoms")
    return n, pos, s.reshape(m, 3)


def spin_series(snaps):
    """Yield (index, time, [float64 [n, 3] spins per data file]) for snapshots(run_dir), with every window
    snapshot applied to the spins of the snapshots before it. The arrays are updated in place: copy to keep."""
    state = None
    for idx, t, files, window in snaps:
        if not window:
            state = [read_spins(f) for f in files]
        elif state is None:
            raise IOError(f"snapshot {idx} is a window snapshot with no complete snapshot before it")
        else:
            for s, f in zip(state, files):
                n, pos, new = read_window_spins(f)
                if n != len(s):
                    raise IOError(f"{f}: {n} atoms, but the snapshot before it has {len(s)}")
                s[pos] = new
        yield idx, t, state


def snapshots(run_dir):
    """[(index, time s, [data files], window)] of the spin snapshots, in order; window is True for a window
    snapshot (hamr:active-window), which holds only the atoms that changed since the previous snapshot (see
    spin_series). Snapshot files are found by probing spins-00000000.meta, ...01.meta, ... (directory listings
    of copied folders have been seen to miss files), falling back to a glob for runs that do not start at 0
    (restarts)."""
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
        res.append((idx, float(head["Time"].split()[0]), files, head.get("Spin data", "full").strip() == "window"))
    return res


# ----------------------------------------------------------------------------------------------- HAMR head model
def _bits_written_per_track(run_dirs):
    """Bits VAMPIRE wrote on the first track, from the "Written bit N" lines of the logs of a run's segments (in
    order) before the first end of track; None without a whole track. Distinct N are counted, since a run resumed
    from a checkpoint replays its loop and logs the bits before the checkpoint again."""
    seen = set()
    for d in run_dirs:
        try:
            with open(os.path.join(d, "log"), errors="replace") as f:
                for line in f:
                    m = re.search(r"Written bit (\d+)", line)
                    if m:
                        seen.add(int(m.group(1)))
                    elif "Head position at end of track" in line:
                        return len(seen)
        except OSError:
            pass
    return None


def log_start(run_dir):
    """Wall-clock time the run started, from the first time stamp of its log ("dd-mm-yyyy [hh:mm:ss]"), as a
    sortable string, or None."""
    try:
        with open(os.path.join(run_dir, "log"), errors="replace") as f:
            m = re.match(r"(\d\d)-(\d\d)-(\d{4}) \[(\d\d:\d\d:\d\d)\]", f.readline())
    except OSError:
        return None
    return f"{m.group(3)}-{m.group(2)}-{m.group(1)} {m.group(4)}" if m else None


def checkpoint_resume_step(run_dir):
    """Step at which a run continued from a checkpoint (sim:load-checkpoint = continue), from its log, or None."""
    if Input(run_dir).str("sim:load-checkpoint") != "continue":
        return None
    try:
        with open(os.path.join(run_dir, "log"), errors="replace") as f:
            for line in f:
                m = re.search(r"Checkpoint file loaded at sim::time (\d+)", line)
                if m:
                    return int(m.group(1))
    except OSError:
        pass
    return None


class Hamr:
    """Head position, laser temperature and field box of a hamr-simulation run, as in src/hamr/hamr_continuous.cpp.
    Steps count from the start of the run's sim::time (equilibration included)."""

    def __init__(self, run_dir, segments=None):
        """segments: all directories of a run continued over several jobs (default: run_dir alone)."""
        inp = Input(run_dir)
        n = inp.num
        self.dt = n("sim:time-step", 1e-15, "time")
        self.eq = n("sim:equilibration-time-steps", 0, "int")
        self.Tmin = n("sim:minimum-temperature", 0.0, "none")
        self.Tmax = n("sim:maximum-temperature", 300.0, "none")              # defaults as src/simulate/sim.cpp
        # the head field ramps by (Hmax - Hmin) in each bit (src/hamr/field_profiles.cpp), so its plateau is
        # Hmax - Hmin, with VAMPIRE's defaults Hmax = +1 T, Hmin = -1 T (i.e. 2 T) unless the input sets them
        self.Hmin = n("sim:minimum-applied-field-strength", -1.0, "field")
        self.Hmax = n("sim:maximum-applied-field-strength", 1.0, "field")
        self.Hpeak = self.Hmax - self.Hmin
        uv = [float(c) for c in inp.str("sim:applied-field-unit-vector", "0,0,1").split(",")]
        self.field_uz = uv[2] / math.sqrt(sum(c * c for c in uv))
        self.rise = int(round(n("hamr:field-rise-time", 1e-12, "time") / self.dt))
        self.fall = int(round(n("hamr:field-fall-time", 1e-12, "time") / self.dt))
        self.Lx, self.Ly, self.Lz = inp.system_size()
        self.v = n("hamr:head-speed", 30.0, "velocity")                      # m/s = nm/ns
        self.BL = n("hamr:bit-size", 0.0)
        self.TW = n("hamr:track-size", 0.0)
        self.pad = n("hamr:track-padding", 0.0)
        self.NPS = n("hamr:NPS", n("hamr:NFT-to-pole-spacing", 0.0))
        self.fwhm_x, self.fwhm_y = n("hamr:laser-FWHM-x", 20.0), n("hamr:laser-FWHM-y", 20.0)
        self.box_x = min(n("hamr:head-field-x", 20.0), self.Lx)
        self.box_y = min(n("hamr:head-field-y", 20.0), self.Ly)
        self.sigma_x = self.fwhm_x / math.sqrt(8 * math.log(2))
        self.sigma_y = self.fwhm_y / math.sqrt(8 * math.log(2))
        if self.BL <= 0 or self.TW <= 0:
            raise ValueError(f"{run_dir}: a HAMR write needs hamr:bit-size and hamr:track-size in its input")
        # tracks, bits and sequence as src/hamr/initialize.cpp and bit_sequence.cpp decide them (lengths in A there)
        self.bits_per_track = int(math.floor((self.Lx - 0.1) / self.BL))
        max_tracks = int(math.floor((self.Ly - 0.1 - 2 * self.pad) / self.TW))
        n_bits = n("hamr:number-of-bits", 0, "int")
        try_tracks = int(round(n_bits / self.bits_per_track)) if self.bits_per_track else 0
        self.n_tracks = try_tracks if try_tracks < max_tracks else max_tracks
        if inp.str("hamr:bit-sequence-type", "user-defined") == "single-tone-predefined" or not inp.has("hamr:bit-sequence"):
            seq = [-1 if i % 2 == 0 else 1 for i in range(n_bits)]
        else:
            seq = vampire_integers(inp.str("hamr:bit-sequence"))
            seq = seq[:n_bits] if n_bits < len(seq) else seq
        total = self.n_tracks * self.bits_per_track
        self.sequence = (seq + [0] * total)[:total]                    # trimmed, or padded with field-free bits
        vdt = self.v * 1e9 * self.dt                                     # nm per step
        self.step_nm = vdt
        self.nb = int(round(self.BL / vdt))
        self.ne = int(round(self.BL / 2 / vdt)) + int(round(self.NPS / vdt))
        self.nt = self.bits_per_track * self.nb + 2 * self.ne
        self.laser_steps = self.n_tracks * self.nt
        # bits written per track, as the run's log records them. Before VAMPIRE's bit loop was fixed it counted from
        # the end of the lead-in (ne) but stopped at bits_per_track * nb, so a lead-in of a bit or more lost bits at the
        # end of every track (the field stays off over their cells); that model is used when the log has no whole track
        self.written = _bits_written_per_track(segments or [run_dir])
        if self.written is None:
            self.written = max(0, min(self.bits_per_track, -(-(self.bits_per_track * self.nb - self.ne) // self.nb)))
        w, bpt = self.written, self.bits_per_track
        # sequence per bit cell (track * bits_per_track + bit); VAMPIRE takes bit t * written + b of its sequence
        self.sequence = [self.sequence[t * w + b] if b < w else 0 for t in range(self.n_tracks) for b in range(bpt)]

    def params(self):
        """Everything the viewer needs to rebuild the head, as plain numbers."""
        keys = ("dt", "eq", "Tmin", "Tmax", "Hmin", "Hmax", "Hpeak", "field_uz", "rise", "fall", "Lx", "Ly", "Lz", "v", "BL", "TW", "pad", "NPS", "fwhm_x",
                "fwhm_y", "box_x", "box_y", "sigma_x", "sigma_y", "step_nm", "nb", "ne", "bits_per_track", "written", "nt",
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
        return -self.BL / 2 - self.NPS + (j + 1) * self.step_nm, self.track_y(track), on

    def track_y(self, track):
        """Head y (nm) on a track: VAMPIRE takes the track from bits written so far / bits_per_track (integer division),
        so in runs that lost bits per track (see written) later tracks land on earlier ones."""
        return self.TW * (track * self.written // max(self.bits_per_track, 1) + 0.5) + self.pad

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

    def track_of(self, y):
        return np.clip(np.floor((np.asarray(y) - self.pad) / self.TW), 0, self.n_tracks - 1).astype(int)

    def freeze_lag(self, Tc):
        """Distance (nm) behind the laser centre where the track centre line cools through Tc."""
        return self.sigma_x * math.sqrt(2 * math.log((self.Tmax - self.Tmin) / (Tc - self.Tmin)))

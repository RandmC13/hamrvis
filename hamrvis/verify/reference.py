"""Independent reader of VAMPIRE run directories, for checking hamrvis bundles against the raw files.

Written from VAMPIRE's file formats and source (src/vio/read.cpp, src/utility/units.cpp, src/config, src/vio/grain.cpp)
and deliberately sharing no code with hamrvis, so a reading mistake in hamrvis cannot hide itself in its own check.
Simple and slow on purpose: plain loops and dictionaries where hamrvis uses vectorised numpy.
"""
import os
import re

import numpy as np

# the units the test scenarios use, from src/utility/units.cpp: factor to VAMPIRE's internal unit and kind
UNIT = {"A": (1.0, "length"), "nm": (10.0, "length"), "um": (1e4, "length"), "m": (1e10, "length"),
        "s": (1.0, "time"), "ps": (1e-12, "time"), "fs": (1e-15, "time"),
        "T": (1.0, "field"), "kOe": (0.1, "field"), "m/s": (1e10, "velocity")}


def parse_input(path):
    """{key: (value, unit)} as VAMPIRE reads its input: all whitespace removed, '#' starts a comment, then
    key:word=value!unit (src/vio/read.cpp tokenises on ':', '=' and '!'). Later lines replace earlier ones."""
    out = {}
    for raw in open(path):
        line = "".join(raw.split()).split("#")[0]
        if not line:
            continue
        m = re.match(r"([^:=!]+):([^=!]+)(?:=([^!]*))?(?:!(.*))?$", line)
        if m:
            out[f"{m.group(1)}:{m.group(2)}"] = (m.group(3) or "", m.group(4) or "")
    return out


def value(keys, key, kind, default):
    """Internal-unit value of key (A, s, T, A/s); default when absent."""
    if key not in keys:
        return default
    v, u = keys[key]
    if not u:
        return float(v)
    factor, k = UNIT[u]
    assert k == kind, (key, u)
    return float(v) * factor


def input_path(run):
    """The input VAMPIRE opened, from the log line 'Opening main input file "name".'."""
    for line in open(os.path.join(run, "log"), errors="replace"):
        m = re.search(r'Opening main input file "(.*)"', line)
        if m:
            return os.path.join(run, m.group(1))
    return os.path.join(run, "input")


def output_path(run):
    """The output table: 'output' (VAMPIRE's default name), else the file whose second line is '# Output file for
    vampire simulation' (datalog.cpp; a job continued from a checkpoint appends to the table without it)."""
    if os.path.isfile(os.path.join(run, "output")):
        return os.path.join(run, "output")
    for f in sorted(os.listdir(run)):
        p = os.path.join(run, f)
        if os.path.isfile(p) and not f.endswith(".data"):
            with open(p, "rb") as fh:
                lines = fh.read(300).split(b"\n")
            if len(lines) > 1 and lines[1] == b"# Output file for vampire simulation":
                return p
    raise FileNotFoundError(run)


def output_table(run):
    """(column names, float64 rows) of the output table, every number parsed with float()."""
    rows, head = [], None
    for line in open(output_path(run)):
        if line.startswith("#") or not line.strip():
            continue
        if head is None:
            head = line.split()
            continue
        rows.append([float(v) for v in line.split()])
    return head, np.array(rows)


def grain_magnetisation(run, n_grains):
    """float64 [frames, grains, 4] (mx, my, mz, |m|) from grain.txt: the single columns before the magnetisation
    block are skipped by counting them from the input's grain:* lines in order."""
    singles = {"time-steps": 1, "real-time": 1, "temperature": 1, "electron-temperature": 1,
               "phonon-temperature": 1, "applied-field-strength": 1, "applied-field-unit-vector": 3}
    skip = 0
    for k in parse_input(input_path(run)):
        if not k.startswith("grain:"):
            continue
        w = k[6:]
        if w == "magnetisation":
            break
        skip += singles.get(w, 0)
    frames = []
    for line in open(os.path.join(run, "grain.txt")):
        if line.strip() and not line.startswith("#"):
            v = [float(x) for x in line.split()]
            frames.append(np.array(v[skip:skip + 4 * n_grains]).reshape(n_grains, 4))
    return np.array(frames)


def _data_files(meta):
    """Data files named in a config meta file (after 'Number of ... files: N'), in order."""
    text = open(meta).read()
    m = re.search(r"Number of \w+ files:\s*(\d+)", text)
    names = re.findall(r"^(\S+\.data)\s*$", text[m.end():], re.M)[: int(m.group(1))]
    return [os.path.join(os.path.dirname(meta), n) for n in names]


def _binary(path):
    with open(path, "rb") as f:
        return b"\0" in f.read(64)


def coordinates(run):
    """(type, xyz in A) of all atoms, data files concatenated in meta order."""
    types, xyz = [], []
    for p in _data_files(os.path.join(run, "atoms-coords.meta")):
        if _binary(p):
            raw = open(p, "rb").read()
            n = int(np.frombuffer(raw[:8], np.uint64)[0])
            types.append(np.frombuffer(raw[8:8 + 4 * n], np.int32))
            xyz.append(np.frombuffer(raw[8 + 8 * n:8 + 32 * n], np.float64).reshape(n, 3))
        else:
            lines = open(p).read().split("\n")[1:]
            rows = [l.split() for l in lines if l.strip()]
            types.append(np.array([int(r[0]) for r in rows]))
            xyz.append(np.array([[float(c) for c in r[2:5]] for r in rows]))
    return np.concatenate(types), np.concatenate(xyz)


def spins(run, index):
    """float64 [atoms, 3] of complete snapshot index, data files concatenated in meta order."""
    out = []
    for p in _data_files(os.path.join(run, f"spins-{index:08d}.meta")):
        if _binary(p):
            raw = open(p, "rb").read()
            n = int(np.frombuffer(raw[:8], np.uint64)[0])
            out.append(np.frombuffer(raw[8:8 + 24 * n], np.float64).reshape(n, 3))
        else:
            lines = open(p).read().split("\n")[1:]
            out.append(np.array([[float(c) for c in l.split()] for l in lines if l.strip()]))
    return np.concatenate(out)


def column_means(xyz, s, pitch_A, nx, ny, factor=1, sites=None):
    """[3, ny, nx] mean spin of the atoms in each top-down column (NaN: no atom). Atoms sit on sites pitch_A apart on the
    periodic film; a column holds factor x factor sites, so a site's column is (site index mod sites) // factor. sites:
    (x, y) sites on the film, default nx, ny."""
    sx, sy = sites or (nx, ny)
    iy = np.rint(xyz[:, 1] / pitch_A).astype(np.int64) % sy // factor
    ix = np.rint(xyz[:, 0] / pitch_A).astype(np.int64) % sx // factor
    k = iy * nx + ix
    count = np.bincount(k, minlength=nx * ny)
    out = np.full((3, ny, nx), np.nan)
    for c in range(3):
        with np.errstate(invalid="ignore", divide="ignore"):
            m = np.bincount(k, s[:, c], minlength=nx * ny) / count
        out[c] = np.where(count > 0, m, np.nan).reshape(ny, nx)
    return out



def grain_ids(run):
    """int32 grain ID of every atom, data files in meta order (binary files hold them after the coordinates in newer
    VAMPIRE versions, text files in column 6), or None if the files have none."""
    out = []
    for p in _data_files(os.path.join(run, "atoms-coords.meta")):
        if _binary(p):
            raw = open(p, "rb").read()
            n = int(np.frombuffer(raw[:8], np.uint64)[0])
            if len(raw) < 8 + 36 * n:
                return None
            out.append(np.frombuffer(raw[8 + 32 * n:8 + 36 * n], np.int32))
        else:
            rows = [l.split() for l in open(p).read().split("\n")[1:] if l.strip()]
            if len(rows[0]) < 6:
                return None
            out.append(np.array([int(r[5]) for r in rows]))
    return np.concatenate(out)


def magnetic_materials(run):
    """0-based indices of the magnetic materials: every material[n] of the input's material file without
    material[n]:non-magnetic (src/material); None without a material file (then every material counts)."""
    name = parse_input(input_path(run)).get("material:file", ("", ""))[0]
    path = os.path.join(run, name)
    if not name or not os.path.isfile(path):
        return None
    keys = parse_input(path)
    mats = {int(re.match(r"material\[(\d+)\]", k).group(1)) for k in keys if re.match(r"material\[\d+\]", k)}
    off = {int(re.match(r"material\[(\d+)\]", k).group(1)) for k in keys if re.match(r"material\[\d+\]:non-magnetic$", k)}
    return sorted(m - 1 for m in range(1, max(mats) + 1) if m not in off) if mats else None


def snapshot_list(run):
    """[(index, time in s, window)] of the spin snapshots named spins-NNNNNNNN.meta, in index order; window is True for
    a snapshot that holds only the atoms changed since the one before (hamr:active-window, 'Spin data: window')."""
    out = []
    for f in sorted(os.listdir(run)):
        m = re.fullmatch(r"spins-(\d{8})\.meta", f)
        if m:
            text = open(os.path.join(run, f)).read()
            t = float(re.search(r"^Time:\s*(\S+)", text, re.M).group(1))
            w = re.search(r"^Spin data:\s*(\S+)", text, re.M)
            out.append((int(m.group(1)), t, bool(w and w.group(1) == "window")))
    return out


def spin_states(run, wanted=None):
    """Yield (index, float64 [atoms, 3]) for every snapshot in wanted (default all), in order. A window snapshot data
    file is u64 atoms, u64 ranges, ranges x (first position, count) as u64, then the spins of those positions (f64 x 3);
    it changes the spins of the snapshot before it, file by file. The yielded array is updated in place."""
    state = None
    for index, _, window in snapshot_list(run):
        files = _data_files(os.path.join(run, f"spins-{index:08d}.meta"))
        if not window:
            state = []
            for p in files:
                if _binary(p):
                    raw = open(p, "rb").read()
                    n = int(np.frombuffer(raw[:8], np.uint64)[0])
                    state.append(np.frombuffer(raw[8:8 + 24 * n], np.float64).reshape(n, 3).copy())
                else:
                    lines = open(p).read().split("\n")[1:]
                    state.append(np.array([[float(c) for c in l.split()] for l in lines if l.strip()]))
        else:
            assert state is not None, f"window snapshot {index} with no complete one before it"
            for s, p in zip(state, files):
                raw = open(p, "rb").read()
                n, r = (int(v) for v in np.frombuffer(raw[:16], np.uint64))
                assert n == len(s), (p, n, len(s))
                rng = np.frombuffer(raw[16:16 + 16 * r], np.uint64).reshape(r, 2)
                vals = np.frombuffer(raw[16 + 16 * r:], np.float64).reshape(-1, 3)
                k = 0
                for first, count in rng:
                    s[int(first):int(first) + int(count)] = vals[k:k + int(count)]
                    k += int(count)
                assert k == len(vals), p
        if wanted is None or index in wanted:
            yield index, np.concatenate(state) if len(state) > 1 else state[0]


def grain_list(run):
    """float64 [grains, 5] of grain-coordinates.txt: id, x, y (A), atoms, moment."""
    rows = [l.split() for l in open(os.path.join(run, "grain-coordinates.txt")) if l.strip() and not l.startswith("#")]
    return np.array([[float(v) for v in r[:5]] for r in rows])


def dipole_cells(run):
    """float64 [cells, 9] of dipole-field.txt: id, atoms, x, y, z (A), moment, Bx, By, Bz (T)."""
    rows = [l.split() for l in open(os.path.join(run, "dipole-field.txt")) if l.strip() and not l.startswith("#")]
    return np.array([[float(v) for v in r[:9]] for r in rows])


def layer_means(xyz, s, pitch_A, pitch_z_A, z0_A, nz, vx, vy, factor, sites):
    """[3, nz, vy, vx] mean spin per atomic layer on cells of factor x factor sites (NaN: no atom)."""
    sx, sy = sites
    ix = np.rint(xyz[:, 0] / pitch_A).astype(np.int64) % sx // factor
    iy = np.rint(xyz[:, 1] / pitch_A).astype(np.int64) % sy // factor
    iz = np.clip(np.rint((xyz[:, 2] - z0_A) / pitch_z_A).astype(np.int64), 0, nz - 1)
    k = (iz * vy + iy) * vx + ix
    count = np.bincount(k, minlength=nz * vy * vx)
    out = np.full((3, nz, vy, vx), np.nan)
    for c in range(3):
        with np.errstate(invalid="ignore", divide="ignore"):
            m = np.bincount(k, s[:, c], minlength=nz * vy * vx) / count
        out[c] = np.where(count > 0, m, np.nan).reshape(nz, vy, vx)
    return out

class Head:
    """The HAMR head at a step, from the head parameters of a bundle (themselves checked against the input by
    test_head_parameters_in_any_units), written from the definitions in docs/analysis.md and src/hamr/hamr_continuous.cpp:
    the laser starts half a bit plus NPS before the film, moves step_nm per step, writes bits_per_track bits of nb
    steps after a lead-in of ne steps with the field off; scalar on purpose."""

    def __init__(self, p):
        self.p = p
        self.w = p.get("written", p["bits_per_track"])

    def at(self, step):
        """(laser x, laser y, on, track, step in track j)."""
        p = self.p
        k = int(step) - p["eq"]
        on = 1 <= k <= p["n_tracks"] * p["nt"]
        kk = max(k, 1)
        track, j = (kk - 1) // p["nt"], (kk - 1) % p["nt"] + 1
        x = -p["BL"] / 2 - p["NPS"] + (j + 1) * p["step_nm"]
        y = p["TW"] * ((track * self.w) // max(p["bits_per_track"], 1) + 0.5) + p["pad"]
        return x, y, on, track, j

    def field(self, step):
        """Signed head field (T) along the field unit vector."""
        p = self.p
        _, _, on, track, j = self.at(step)
        i = j - 1 - p["ne"]                              # step within the written part of the track
        if not on or i < 0:
            return 0.0
        b, i = divmod(i, p["nb"])
        if b >= p["bits_per_track"]:
            return 0.0
        mag = (p["nb"] - 1 - i) / p["fall"] if i >= p["nb"] - p["fall"] else min((i + 1) / p["rise"], 1.0)
        seq = p["sequence"]
        return mag * p["Hpeak"] * seq[min(track * p["bits_per_track"] + b, len(seq) - 1)]

    def bit(self, step):
        """(track, bit) the field is writing, or None."""
        p = self.p
        _, _, on, track, j = self.at(step)
        i = j - 1 - p["ne"]
        return (track, i // p["nb"]) if on and i >= 0 and i // p["nb"] < self.w else None

    def T(self, x, y, step):
        p = self.p
        hx, hy, on, _, _ = self.at(step)
        if not on:
            return p["Tmin"]
        return p["Tmin"] + (p["Tmax"] - p["Tmin"]) * np.exp(-(x - hx) ** 2 / (2 * p["sigma_x"] ** 2) - (y - hy) ** 2 / (2 * p["sigma_y"] ** 2))

    def in_box(self, x, y, step):
        p = self.p
        hx, hy, on, _, _ = self.at(step)
        return on and abs(x - (hx - p["NPS"])) <= p["box_x"] / 2 and abs(y - hy) <= p["box_y"] / 2

    def peak_T(self, x, y):
        """Hottest laser temperature ever at (x, y): the laser centre comes closest at the nearest point of its path
        (x from the first to the last step of a track) on the nearest track."""
        p = self.p
        x0 = -p["BL"] / 2 - p["NPS"] + 2 * p["step_nm"]
        x1 = -p["BL"] / 2 - p["NPS"] + (p["nt"] + 1) * p["step_nm"]
        dx = min(max(x, x0), x1) - x
        best = max(np.exp(-dx ** 2 / (2 * p["sigma_x"] ** 2) - (y - self.at(p["eq"] + t * p["nt"] + 1)[1]) ** 2 / (2 * p["sigma_y"] ** 2))
                   for t in range(p["n_tracks"]))
        return p["Tmin"] + (p["Tmax"] - p["Tmin"]) * best

    def pass_step(self, x, track):
        """Step at which the laser centre is over x on a track."""
        p = self.p
        return p["eq"] + track * p["nt"] + round((x + p["BL"] / 2 + p["NPS"]) / p["step_nm"]) - 1

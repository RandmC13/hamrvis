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

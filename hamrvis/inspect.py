"""hamrvis inspect: what hamrvis finds in run directories and what each view needs, without reducing anything.

Reads only small files (inputs, logs, meta files, file sizes, the first rows of tables and the first atoms of the
coordinate file), so it is quick on any size of run.
"""
import os

import numpy as np

from . import reduce as red
from . import runs as rn
from . import vampire as vp

# measured on a 1.5M-atom granular write (text snapshots, 2026-11-04): gzip keeps these fractions of the int8 grids;
# reducing ran at ~45 MB/s of spin files per worker
GZ_TOP, GZ_VOL, MB_PER_S = 0.53, 0.68, 45.0


def _lines(path):
    """Data rows of a text table (lines not starting with '#')."""
    n = 0
    with open(path, "rb") as f:
        for line in f:
            n += not line.startswith(b"#") and bool(line.strip())
    return n


def _first_atoms(path, n_max=200000):
    """xyz (A) of the first atoms of a coordinate data file (pitch estimate)."""
    if vp._is_binary(path):
        with open(path, "rb") as f:
            n = int(np.fromfile(f, np.uint64, 1)[0])
            f.seek(8 + 8 * n)                          # type and category, int32 each
            m = min(n, n_max)
            return np.fromfile(f, np.float64, 3 * m).reshape(-1, 3)
    rows = []
    with open(path) as f:
        f.readline()
        for _, line in zip(range(n_max), f):
            rows.append([float(v) for v in line.split()[2:5]])
    return np.array(rows)


def _mb(b):
    return f"{b / 1e6:.1f} MB" if b < 1e9 else f"{b / 1e9:.2f} GB"


def run_report(r, opts):
    d, out = r.dir, []
    say = out.append
    inp = vp.Input(d)
    say(f"{r.id}: {r.type} tab \"{r.title}\"" + (f", {len(r.dirs)} jobs joined: " + ", ".join(
        os.path.basename(x) for x in r.dirs) if len(r.dirs) > 1 else ""))
    say(f"  folder        {d}")
    say(f"  started       {vp.log_start(d) or 'unknown (no log time stamp)'}")
    say(f"  program       {inp.str('sim:program', '(none)')}, integrator {inp.str('sim:integrator', 'llg-heun')}")
    oname = vp.output_name(d)
    names = f"input file '{inp.name}', output table '{oname}'"
    say(f"  files         {names}")
    mat = inp.str("material:file")
    say("  materials     " + (f"{mat}" if mat and os.path.isfile(os.path.join(d, mat)) else
                              f"{mat!r} missing: every material counts as magnetic" if mat else
                              "no material:file: every material counts as magnetic"))
    size, sec = 0, 0.0
    try:
        cols, _ = vp.read_output(d) if oname else (None, None)
        say(f"  output table  {len(cols)} columns ({', '.join(cols[:6])}{', ...' if len(cols) > 6 else ''})")
        size += 8 * len(cols) * _lines(os.path.join(d, oname))
    except (OSError, ValueError, IndexError) as err:
        say(f"  output table  NOT USABLE: {err}")
    # grains
    gfile = os.path.join(d, "grain-coordinates.txt")
    if not os.path.exists(gfile):
        say("  grain views   no: no grain-coordinates.txt (VAMPIRE writes it for granular systems only)")
    else:
        nG = len(vp.read_grain_coordinates(d))
        if not os.path.exists(os.path.join(d, "grain.txt")):
            say(f"  grain views   no: {nG} grains but no grain.txt (add grain:magnetisation and grain:real-time)")
        else:
            try:
                vp.grain_columns(d, nG)
                frames = sum(_lines(os.path.join(x, "grain.txt")) for x in r.dirs)
                fr = -(-frames // opts.get("grain_every", 1))
                gbytes = fr * nG * 4 * (4 if opts.get("grain_precision", 16) == 32 else 2) + fr * 16
                size += gbytes
                say(f"  grain views   yes: {nG} grains x {frames} frames in grain.txt ({_mb(gbytes)} in the bundle)")
            except ValueError as err:
                say(f"  grain views   no: grain.txt not usable: {err}")
    # atoms
    cfiles = vp.coord_files(d)
    snaps = red.run_snapshots(r.dirs, inp.num("sim:time-step", 1e-15, "time")) if cfiles else []
    if not cfiles:
        say("  atom views    no: no atom coordinates (add config:atoms and config:atoms-output-rate)")
    elif not snaps:
        say("  atom views    no: coordinates but no spin snapshots (add config:atoms-output-rate)")
    else:
        meta = os.path.join(d, "atoms-coords.meta")
        n_atoms = int(vp._meta_files(meta)[0].get("Number of atoms", "0").split()[0]) if os.path.exists(meta) else 0
        windowed = sum(1 for *_, w in snaps if w)
        spin_bytes = sum(os.path.getsize(f) for _, _, fs, _ in snaps for f in fs if os.path.exists(f))
        xyz = _first_atoms(cfiles[0])
        g = red.Grid(*[v for v in (inp.num("dimensions:system-size-x"), inp.num("dimensions:system-size-y"))], xyz,
                     None, opts.get("top_max", red.TOP_MAX), opts.get("layer_max", red.LAYER_MAX))
        g.nz = max(1, int(round(inp.num("dimensions:system-size-z") * 10 / g.pz)))   # first atoms may miss layers
        bits = opts.get("precision", 8)
        kept = len(snaps[::opts.get("every", 1)])
        top, vol = 3 * g.nx * g.ny * bits // 8, 3 * g.nz * g.vy * g.vx * bits // 8
        gz = (GZ_TOP, GZ_VOL) if bits == 8 else (1.0, 1.0)
        snap_bytes = kept * (top * gz[0] + (vol * gz[1] if opts.get("volume", True) else 0))
        if opts.get("snapshots", True):
            size += snap_bytes
            sec += spin_bytes / (MB_PER_S * 1e6) / max(1, opts.get("workers") or (os.cpu_count() or 2) - 1)
        say(f"  atom views    yes: {n_atoms or '?'} atoms in {len(cfiles)} coordinate file(s), {len(snaps)} snapshots"
            + (f" ({windowed} window snapshots)" if windowed else "") + f", {_mb(spin_bytes)} of spin files")
        say(f"                grid ~{g.nx} x {g.ny} columns of {g.p * g.tf / 10:.3f} nm, layers ~{g.vx} x {g.vy} x {g.nz}"
            f"; {kept} snapshots kept, ~{_mb(snap_bytes)} in the bundle"
            + ("" if opts.get("snapshots", True) else " (left out: --no-snapshots)"))
        size += 4 * g.nx * g.ny
    # type-specific
    if r.type == "hamr":
        try:
            h = vp.Hamr(d, r.dirs)
            say(f"  head          {h.n_tracks} track(s) x {h.bits_per_track} bits, {h.written} written per track,"
                f" bit {h.BL:g} nm, track {h.TW:g} nm, speed {h.v:g} m/s")
        except (ValueError, ZeroDivisionError) as err:
            say(f"  head          NOT USABLE: {err}")
    last = r.dirs[-1]
    if os.path.exists(os.path.join(last, "dipole-field.txt")):
        say("  field map     yes: dipole-field.txt")
        size += os.path.getsize(os.path.join(last, "dipole-field.txt")) // 6
    elif r.type == "dipole-read":
        say("  field map     no dipole-field.txt (add dipole:output-dipole-field)")
    return out, size, sec


def inspect(paths, order=None, opts=None):
    opts = opts or {}
    found = rn.find(paths, order)
    if not found:
        raise SystemExit(f"no VAMPIRE runs in {paths}")
    print(f"{len(found)} run(s), in tab order ({'as given by --order' if order else 'log start time, then name'}):\n")
    total, secs = 0, 0.0
    for r in found:
        lines, size, sec = run_report(r, opts)
        total, secs = total + size, secs + sec
        print("\n".join(lines) + "\n")
    print(f"estimated bundle ~{_mb(total)}, reduce ~{max(1, round(secs))} s (rough: from a 1.5M-atom run's "
          f"compression and speed)")

"""Make the VAMPIRE test scenarios in tests/scenarios/: python3 tests/make_scenarios.py /path/to/vampire-serial [scenario ...]

file-per-process needs the MPI build: python3 tests/make_scenarios.py /path/to/vampire-parallel file-per-process;
mpi-continue both: VAMPIRE_SERIAL=/path/to/vampire-serial python3 tests/make_scenarios.py /path/to/vampire-parallel mpi-continue.

Every scenario is a folder of runs on the tiny granular film of tests/fixture-run (2,624 atoms), chained by
checkpoints the way HAMR studies chain them. The tests reduce each folder as a whole, so they check that hamrvis
finds whatever runs are there. Runs take seconds each; outputs are committed so the tests need no VAMPIRE.
"""
import os
import re
import shutil
import signal
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "fixture-run")
OUT = os.path.join(HERE, "scenarios")

# the film and output keywords of fixture-run, without its program keywords
FILM = [l for l in open(os.path.join(BASE, "input")) if not l.startswith("#") and not re.match(
    r"(sim:|hamr:|config:atoms-output-rate)", l)]

HAMR = """sim:minimum-temperature = 300
sim:maximum-temperature = 1000
sim:equilibration-temperature = 300
sim:time-step = 1.0 !fs
sim:minimum-applied-field-strength = 0.0
sim:maximum-applied-field-strength = 1.5
sim:applied-field-unit-vector = 0,0,-1
sim:equilibration-time-steps = 200
sim:program = hamr-simulation
sim:integrator = llg-heun
hamr:laser-FWHM-x = 4 !nm
hamr:laser-FWHM-y = 4 !nm
hamr:head-speed = 1500 !m/s
hamr:head-field-x = 3 !nm
hamr:head-field-y = 4 !nm
hamr:NPS = 1 !nm
hamr:bit-size = 3 !nm
hamr:track-size = 5 !nm
hamr:number-of-bits = 3
hamr:bit-sequence-type = single-tone-predefined
hamr:field-rise-time = 0.1e-12 !s
hamr:field-fall-time = 0.1e-12 !s
config:atoms-output-rate = 2000
"""

MC_RELAX = """sim:temperature = 300
sim:equilibration-temperature = 300
sim:time-step = 1.0 !fs
sim:equilibration-time-steps = 0
sim:total-time-steps = 3000
sim:program = time-series
sim:integrator = monte-carlo
config:atoms-output-rate = 1000
"""

LLG_RELAX = MC_RELAX.replace("monte-carlo", "llg-heun")

READ = """sim:temperature = 300
sim:equilibration-temperature = 300
sim:time-step = 1.0 !fs
sim:equilibration-time-steps = 0
sim:total-time-steps = 1
sim:program = time-series
sim:integrator = llg-heun
cells:macro-cell-size = 1.0 !nm
dipole:solver = macrocell
dipole:output-dipole-field
config:atoms-output-rate = 1
"""


def run(vampire, d, text, load=None, save=True, kill_after_snapshot=None, film=None, input_name="input", args=(),
        ranks=1):
    """Write d/<input_name> (film + text), copy the material file and a checkpoint, run VAMPIRE in d (with args, on
    ranks MPI processes if more than 1: vampire is then the parallel build). kill_after_snapshot=k stops VAMPIRE
    once spins-k exists, as a job hitting its wall time would."""
    os.makedirs(d)
    shutil.copy(os.path.join(BASE, "m.mat"), d)
    extra = ""
    if load:
        shutil.copy(os.path.join(load, "vampire0.chk"), d)
        extra += "sim:load-checkpoint = restart\n"
    if save:
        extra += "sim:save-checkpoint = end\n" if kill_after_snapshot is None else \
            "sim:save-checkpoint = continuous\nsim:save-checkpoint-rate = 500\n"
    with open(os.path.join(d, input_name), "w") as f:
        f.write("".join(FILM if film is None else film) + text + extra)
    cmd = (["mpirun", "--oversubscribe", "-np", str(ranks)] if ranks > 1 else []) + [vampire, *args]
    with open(os.path.join(d, "stdout"), "w") as so:
        p = subprocess.Popen(cmd, cwd=d, stdout=so, stderr=subprocess.STDOUT)
        if kill_after_snapshot is not None:
            flag = os.path.join(d, f"spins-{kill_after_snapshot:08d}.meta")
            while p.poll() is None and not os.path.exists(flag):
                time.sleep(0.05)
            time.sleep(0.3)
            p.send_signal(signal.SIGKILL)
        p.wait()
    print(f"  {os.path.relpath(d, HERE)}: exit {p.returncode}")


# a film twice as wide in y: two tracks of three bits, with a user-defined sequence holding field-free (0) bits
WIDE = [l.replace("system-size-y = 6 !nm", "system-size-y = 10 !nm") for l in FILM]
TWO_TRACKS = HAMR.replace("hamr:track-size = 5 !nm", "hamr:track-size = 4.5 !nm").replace(
    "hamr:number-of-bits = 3\nhamr:bit-sequence-type = single-tone-predefined\n",
    "hamr:number-of-bits = 6\nhamr:bit-sequence-type = user-defined\nhamr:bit-sequence = 1,0,-1,0,1,-1\n")

# the same film without grains: a continuous layer (VAMPIRE writes no grain-coordinates.txt for it)
CONTINUOUS = [l for l in FILM if not re.match(r"create:(granular-film|grain-)", l)]

# the write of fixture-run in other units VAMPIRE accepts (same values), and the head speed without a unit, which
# VAMPIRE reads in its internal unit A/s
UNITS = {"dimensions:system-size-x = 12 !nm": "dimensions:system-size-x = 0.012 !um",
         "sim:time-step = 1.0 !fs": "sim:time-step = 0.001 !ps",
         "sim:maximum-applied-field-strength = 1.5": "sim:maximum-applied-field-strength = 15 !kOe",
         "hamr:laser-FWHM-x = 4 !nm": "hamr:laser-FWHM-x = 40 !A",
         "hamr:head-speed = 1500 !m/s": "hamr:head-speed = 1.5e13",
         "hamr:bit-size = 3 !nm": "hamr:bit-size = 3e-9 !m",
         "hamr:field-rise-time = 0.1e-12 !s": "hamr:field-rise-time = 100 !fs"}


def in_units(lines):
    text = "".join(lines)
    for a, b in UNITS.items():
        assert a in text, a
        text = text.replace(a, b)
    return text


def continue_from(vampire, prev, d, text):
    """Second job of a run cut short: copy the checkpoint and resume with load-checkpoint = continue."""
    os.makedirs(d)
    shutil.copy(os.path.join(BASE, "m.mat"), d)
    for f in os.listdir(prev):                    # every process's checkpoint (a job on several processes)
        if f.endswith(".chk"):
            shutil.copy(os.path.join(prev, f), d)
    with open(os.path.join(d, "input"), "w") as f:
        f.write("".join(FILM) + text + "sim:load-checkpoint = continue\nsim:save-checkpoint = end\n")
    with open(os.path.join(d, "stdout"), "w") as so:
        rc = subprocess.call([vampire], cwd=d, stdout=so, stderr=subprocess.STDOUT)
    print(f"  {os.path.relpath(d, HERE)}: exit {rc}")


def scenario(name, vampire):
    root = os.path.join(OUT, name)
    shutil.rmtree(root, ignore_errors=True)
    print(name)
    j = lambda *p: os.path.join(root, *p)
    if name == "chain":                       # MC relax, write, read
        run(vampire, j("relax"), MC_RELAX)
        run(vampire, j("write"), HAMR, load=j("relax"))
        run(vampire, j("read"), READ, load=j("write"), save=False)
    elif name == "write-only":
        run(vampire, j("write"), HAMR, save=False)
    elif name == "split-write":               # one write over two jobs: the first is killed mid-write
        run(vampire, j("job1"), HAMR, kill_after_snapshot=4)
        continue_from(vampire, j("job1"), j("job2"), HAMR)
    elif name == "rewrite":                   # write, read, write a different pattern over it, read; the braces
        # are not VAMPIRE syntax: VAMPIRE reads '{1' as 0, so the first bit of the rewrite is written without field
        seq = "hamr:bit-sequence-type = user-defined\nhamr:bit-sequence = {1,1,-1}\n"
        run(vampire, j("1-write"), HAMR)
        run(vampire, j("2-read"), READ, load=j("1-write"), save=False)
        run(vampire, j("3-rewrite"), HAMR.replace("hamr:bit-sequence-type = single-tone-predefined\n", seq),
            load=j("1-write"))
        run(vampire, j("4-read"), READ, load=j("3-rewrite"), save=False)
    elif name == "two-tracks":                # two tracks, a user-defined sequence with field-free bits
        run(vampire, j("write"), TWO_TRACKS, save=False, film=WIDE)
    elif name == "continuous":                # a write on a film without grains
        run(vampire, j("write"), HAMR, save=False, film=CONTINUOUS)
    elif name == "text-snapshots":            # snapshots in VAMPIRE's text format
        run(vampire, j("write"), HAMR, save=False,
            film=[l.replace("config:output-format = binary", "config:output-format = text") for l in FILM])
    elif name == "file-per-process":          # 2 MPI ranks, one data file per rank (vampire = the parallel build)
        run(vampire, j("write"), HAMR + "config:output-mode = file-per-process\n", save=False, ranks=2)
    elif name == "mpi-continue":              # a write on 2 ranks (one file each) killed mid-write, continued on 1:
        # the second job lists its atoms in another order (vampire = the parallel build, serial = vampire-serial)
        run(vampire, j("job1"), HAMR + "config:output-mode = file-per-process\n", ranks=2, kill_after_snapshot=4)
        continue_from(os.environ.get("VAMPIRE_SERIAL", vampire), j("job1"), j("job2"), HAMR)
    elif name == "custom-names":              # --input-file / --output-file and other units, same write
        run(vampire, j("write"), "", save=False, film=[in_units(FILM + [HAMR])], input_name="write.in",
            args=["--input-file", "write.in", "--output-file", "results.txt"])
    elif name == "llg-relax":                 # an equilibration hamrvis has no specific views for, then a write
        run(vampire, j("anneal"), LLG_RELAX)
        run(vampire, j("write"), HAMR, load=j("anneal"), save=False)
    else:
        raise SystemExit(f"unknown scenario {name}")
    for dp, _, fs in os.walk(root):               # checkpoints are only needed to chain the runs
        for f in fs:
            if f.endswith(".chk"):
                os.remove(os.path.join(dp, f))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    names = sys.argv[2:] or ["chain", "write-only", "split-write", "rewrite", "llg-relax", "two-tracks", "continuous",
                             "text-snapshots", "custom-names"]   # file-per-process: run with the parallel build
    for n in names:
        scenario(n, os.path.abspath(sys.argv[1]))

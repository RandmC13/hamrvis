"""Find the VAMPIRE runs in the given paths, order them, join continued jobs and work out what each run is.

A *run* is one simulation: a run directory, or several when a run was continued from checkpoints over several jobs
(segments). Each run gets one tab in the viewer. Its *type* decides the title and the views specific to that kind
of simulation; everything else (grain views, atom views, field maps) follows from the files the run wrote.
"""
import os
import re

from . import vampire as vp

# Run types, first match wins: (name, tab title, test on the run's Input). Add new kinds of simulation here and
# give them views in viewer/app.js (TYPES there); a run matching none is "generic".
TYPES = [
    ("hamr", "Write", lambda inp: inp.str("sim:program") == "hamr-simulation"),
    ("dipole-read", "Read-back", lambda inp: inp.has("dipole:output-dipole-field")
     or inp.has("dipole:output-atomistic-dipole-field")),
    ("monte-carlo", "Monte Carlo", lambda inp: is_monte_carlo(inp)),
]
TITLES = dict((n, t) for n, t, _ in TYPES)

# input keys that may differ between the jobs of one continued run
_CHECKPOINT_KEYS = re.compile(r"sim:(load|save)-checkpoint")


def is_monte_carlo(inp):
    """True for VAMPIRE's Monte Carlo integrators (monte-carlo, constrained-monte-carlo, hybrid-constrained-...)."""
    return "monte-carlo" in (inp.str("sim:integrator", "") or "")


class Run:
    def __init__(self, dirs):
        self.dirs = dirs                                   # segments, in order
        self.dir = dirs[0]
        inp = vp.Input(self.dir)
        self.type = next((n for n, _, test in TYPES if test(inp)), "generic")
        self.title = TITLES.get(self.type) or (inp.str("sim:program") or "run")
        self.id = None

    def __repr__(self):
        return f"Run({self.id}, {self.type}, {self.dirs})"


def candidates(paths):
    """Run directories among paths and their immediate subfolders. A path given directly must hold a VAMPIRE input
    file; in a scanned folder only subfolders with both an input file and an output table (a run that ran) count,
    so copies of inputs or source trees next to the runs are skipped. File names: see vampire.input_name and
    output_name."""
    runs = []
    for p in paths:
        if os.path.isdir(p) and vp.input_name(p):
            runs.append(p)
            continue
        if not os.path.isdir(p):
            raise SystemExit(f"{p}: not a directory")
        for d in sorted(os.listdir(p)):
            q = os.path.join(p, d)
            if os.path.isdir(q) and vp.input_name(q):
                if vp.output_name(q):
                    runs.append(q)
                else:
                    print(f"skipping {q}: it has an input file but no output table (not a run; if VAMPIRE was run "
                          f"with --output-file, give that name with --output-file)")
    return runs


def _continues(prev, d):
    """True if run directory d continues prev from a checkpoint: same input apart from the checkpoint keys, and it
    resumed at a step prev reached."""
    step = vp.checkpoint_resume_step(d)
    if step is None:
        return False
    strip = lambda k: {a: b for a, b in k.items() if not _CHECKPOINT_KEYS.match(a)}
    if strip(vp.Input(prev).keys) != strip(vp.Input(d).keys):
        return False
    try:
        cols, out = vp.read_output(prev)
    except (OSError, ValueError, IndexError):
        return False
    dt = vp.Input(prev).num("sim:time-step", 1e-15, "time")
    t = out[:, cols.index("Real_time")] if "Real_time" in cols else None
    return t is not None and len(t) and t[-1] / dt >= step - 1


def find(paths, order=None):
    """[Run] in run order. order: run directory names (basenames) in the wanted order; otherwise runs are ordered by
    the time their log was opened (runs chained by checkpoints run one after another), then by name. Consecutive
    directories that continue one another become one run."""
    dirs = [os.path.normpath(d) for d in candidates(paths)]
    if order:
        names = [os.path.basename(d) for d in dirs]
        missing = [o for o in order if o not in names]
        if missing:
            raise SystemExit(f"--order names {missing}, not among the runs found: {names}")
        dirs = [dirs[names.index(o)] for o in order] + [d for d in dirs if os.path.basename(d) not in order]
    else:
        dirs.sort(key=lambda d: (vp.log_start(d) or "~", os.path.basename(d)))
    runs, group = [], []
    for d in dirs:
        if group and _continues(group[-1], d):
            group.append(d)
            continue
        if group:
            runs.append(Run(group))
        group = [d]
    if group:
        runs.append(Run(group))
    used = set()
    for r in runs:                                 # tab id = folder name, made unique and safe for paths and URLs
        base = re.sub(r"[^A-Za-z0-9._-]+", "-", os.path.basename(os.path.abspath(r.dir))) or "run"
        r.id, k = base, 2
        while r.id in used:
            r.id, k = f"{base}-{k}", k + 1
        used.add(r.id)
    return runs

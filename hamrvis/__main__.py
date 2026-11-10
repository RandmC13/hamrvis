"""hamrvis command line: python3 -m hamrvis {inspect,reduce,view,export} ...  (see README.md)"""
import argparse

from . import __version__
from . import reduce as red
from . import vampire as vp


def _run_args(p):
    """Arguments shared by reduce and inspect."""
    p.add_argument("runs", nargs="+", help="run directories, or a folder whose subfolders are run directories")
    p.add_argument("--order", nargs="+", metavar="DIR", help="run folder names in run order (default: the order the "
                                                            "runs started, from their logs, then by name)")
    p.add_argument("--input-file", metavar="NAME", help="name of the runs' input file, if VAMPIRE was run with "
                                                       "--input-file (default: 'input', else the name in the log)")
    p.add_argument("--output-file", metavar="NAME", help="name of the runs' output table, if VAMPIRE was run with "
                                                        "--output-file (default: 'output', else found by its header)")
    p.add_argument("--no-volume", action="store_true", help="skip the per-layer volume (saves ~90%% of the size)")
    p.add_argument("--every", type=int, default=1, help="keep every n-th spin snapshot (the last one is always kept)")
    p.add_argument("--no-snapshots", action="store_true", help="leave the spin snapshots out (grain views only)")
    p.add_argument("--precision", type=int, choices=(8, 16, 32), default=8,
                   help="bits per value of the snapshot grids: 8 (x 1/127, default), 16 (x 1/32767) or 32 (float)")
    p.add_argument("--grain-precision", type=int, choices=(16, 32), default=16,
                   help="bits per value of the grain series: 16 (x 1/32767, default) or 32 (float, as in grain.txt)")
    p.add_argument("--grain-every", type=int, default=1, help="keep every n-th grain.txt frame (the last one is "
                                                              "always kept; the write analysis uses the kept frames)")
    p.add_argument("--top-max", type=float, default=red.TOP_MAX, help="top-down pixels per snapshot before the "
                                                                      "column map is coarsened (default %(default).0f)")
    p.add_argument("--layer-max", type=float, default=red.LAYER_MAX, help="per-layer cells per snapshot before the "
                                                                          "layer grid is coarsened (default %(default).0f)")
    p.add_argument("-j", "--workers", type=int, help="parallel snapshot workers (default: CPUs - 1, as far as --mem allows)")
    p.add_argument("--mem", type=float, metavar="GB", help="memory budget for the snapshot workers (default: what is free, "
                                                           "or left of the job's limit)")


def _opts(a):
    for k in ("every", "grain_every"):
        if getattr(a, k) < 1:
            raise SystemExit(f"--{k.replace('_', '-')} must be 1 or more")
    th = {k: getattr(a, k) for k in red.THRESHOLDS if getattr(a, k, None) is not None}
    return dict(precision=a.precision, grain_precision=a.grain_precision, grain_every=a.grain_every,
                snapshots=not a.no_snapshots, top_max=a.top_max, layer_max=a.layer_max, thresholds=th,
                every=a.every, volume=not a.no_volume, workers=a.workers,
                mem=a.mem * 1e9 if a.mem else None)


def main():
    ap = argparse.ArgumentParser(prog="hamrvis", description="Reduce and view VAMPIRE HAMR simulations")
    ap.add_argument("--version", action="version", version=f"hamrvis {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("inspect", help="list the runs found, what each view needs and whether it is there, and the "
                                       "estimated bundle size (reads small files only)")
    _run_args(i)
    r = sub.add_parser("reduce", help="reduce VAMPIRE runs (run directories or a folder of them) to a viewer bundle, one tab per run")
    _run_args(r)
    r.add_argument("-o", "--out", required=True, help="bundle directory to write")
    r.add_argument("--layer-factor", type=int, help="per-layer cell size in atomic sites (1 = every atom; default 2, "
                                                     "more for large systems)")
    r.add_argument("--tc", type=float, help="Curie temperature (K) for the analysis; default: inferred from the data")
    r.add_argument("--name", help="run name shown in the viewer (default: the runs' common folder name)")
    for k, what in (("decision", "m_z|m| below this in the final direction = not yet decided"),
                    ("heated", "|m| below this at some point = the grain was melted"),
                    ("settled", "|final m_z|m|| below this = still disordered at the end"),
                    ("full_field", "|B| / peak field at the decision above this = full field, below = reversal")):
        r.add_argument(f"--{k.replace('_', '-')}", type=float, help=f"write analysis: {what} (default {red.THRESHOLDS[k]})")
    v = sub.add_parser("view", help="serve the viewer for a bundle and open it in the browser")
    v.add_argument("bundle", help="bundle directory made by 'reduce'")
    v.add_argument("-p", "--port", type=int, default=8737)
    v.add_argument("--no-browser", action="store_true")
    x = sub.add_parser("export", help="write one self-contained HTML page (viewer + thinned bundle) to share")
    x.add_argument("bundle")
    x.add_argument("-o", "--out", required=True, help="HTML file to write")
    x.add_argument("--grain-step", type=int, default=4, help="keep every n-th grain frame (default 4)")
    x.add_argument("--snap-step", type=int, default=1, help="keep every n-th spin snapshot (default 1)")
    x.add_argument("--no-volume", action="store_true", help="leave out the per-layer data")
    x.add_argument("--limit-mb", type=float, default=15.5, help="refuse to write a larger page (default suits sites with a 16 MB limit)")
    x.add_argument("--artifact", action="store_true", help="omit the document skeleton (for hosts that supply their own)")
    a = ap.parse_args()
    if a.cmd in ("reduce", "inspect"):
        vp.NAMES.update(input=a.input_file, output=a.output_file)
    if a.cmd == "inspect":
        from . import inspect
        inspect.inspect(a.runs, a.order, _opts(a))
    elif a.cmd == "reduce":
        red.reduce(a.runs, a.out, volume=not a.no_volume, every=a.every, Tc=a.tc, workers=a.workers, name=a.name,
                   layer_factor=a.layer_factor, order=a.order, opts=_opts(a))
    elif a.cmd == "export":
        from . import export
        export.export(a.bundle, a.out, a.grain_step, a.snap_step, not a.no_volume, a.limit_mb, a.artifact)
    else:
        from . import serve
        serve.serve(a.bundle, a.port, not a.no_browser)


if __name__ == "__main__":
    main()

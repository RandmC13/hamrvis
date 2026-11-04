"""hamrvis command line: python3 -m hamrvis {reduce,view,export} ...  (see README.md)"""
import argparse

from . import reduce as red


def main():
    ap = argparse.ArgumentParser(prog="hamrvis", description="Reduce and view VAMPIRE HAMR simulations")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("reduce", help="reduce VAMPIRE runs (run directories or a folder of them) to a viewer bundle, one tab per run")
    r.add_argument("runs", nargs="+", help="run directories, or a folder whose subfolders are run directories")
    r.add_argument("--order", nargs="+", metavar="DIR", help="run folder names in run order (default: the order the "
                                                            "runs started, from their logs, then by name)")
    r.add_argument("-o", "--out", required=True, help="bundle directory to write")
    r.add_argument("--no-volume", action="store_true", help="skip the per-layer volume (saves ~90%% of the size)")
    r.add_argument("--layer-factor", type=int, help="per-layer cell size in atomic sites (1 = every atom; default 2, "
                                                     "more for large systems)")
    r.add_argument("--every", type=int, default=1, help="keep every n-th spin snapshot (the last one is always kept)")
    r.add_argument("--tc", type=float, help="Curie temperature (K) for the analysis; default: measured from the data")
    r.add_argument("--name", help="run name shown in the viewer (default: the runs' common folder name)")
    r.add_argument("-j", "--workers", type=int, help="parallel snapshot workers (default: CPUs - 1)")
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
    x.add_argument("--limit-mb", type=float, default=15.5, help="refuse to write a larger page (localhost.ai: 16 MB)")
    x.add_argument("--artifact", action="store_true", help="omit the document skeleton (for localhost.ai artifacts)")
    a = ap.parse_args()
    if a.cmd == "reduce":
        red.reduce(a.runs, a.out, volume=not a.no_volume, every=a.every, Tc=a.tc, workers=a.workers, name=a.name,
                   layer_factor=a.layer_factor, order=a.order)
    elif a.cmd == "export":
        from . import export
        export.export(a.bundle, a.out, a.grain_step, a.snap_step, not a.no_volume, a.limit_mb, a.artifact)
    else:
        from . import serve
        serve.serve(a.bundle, a.port, not a.no_browser)


if __name__ == "__main__":
    main()

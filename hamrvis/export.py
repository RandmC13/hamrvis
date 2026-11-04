"""Inline the viewer and a bundle into one self-contained HTML page (to share or publish).

Grain frames are thinned to every n-th frame and spin snapshots to the n-th one (the last always kept) so the page
stays under the size limit; the full bundle is unchanged.
"""
import base64
import copy
import gzip
import json
import os
import re

import numpy as np

from .serve import VIEWER


def _b64gz(data):
    return base64.b64encode(gzip.compress(data, 9)).decode()


def export(bundle, out, grain_step=4, snap_step=1, volume=True, limit_mb=15.5, embedded_page=False):
    man = json.load(open(os.path.join(bundle, "manifest.json")))
    man = copy.deepcopy(man)
    files = {}
    from .reduce import bundle_runs
    for st, e in bundle_runs(man):
        sd = os.path.join(bundle, st)
        for name in ("output.f32", "output.f64", "top_count.u8", "vol_count.u8", "min_m.u8", "grain_map.i16", "dipole.f32",
                     "dipole_atoms.u16"):
            p = os.path.join(sd, name)
            if os.path.exists(p) and (volume or name != "vol_count.u8"):
                files[f"{st}/{name}"] = open(p, "rb").read()
        if "grain_series" in e:
            gs = e["grain_series"]
            nG, nF, nL = len(e["grains"]["id"]), gs["frames"], len(gs["lead"])
            gname, gtype = ("grains.f32", np.float32) if gs.get("dtype") == "f32" else ("grains.i16", np.int16)
            g = np.fromfile(os.path.join(sd, gname), gtype).reshape(nF, nG, 4)
            lead = np.fromfile(os.path.join(sd, "grain_lead.f64")).reshape(nF, nL)
            keep = np.unique(np.r_[np.arange(0, nF, grain_step), nF - 1])
            files[f"{st}/{gname}"] = g[keep].tobytes()
            files[f"{st}/grain_lead.f64"] = lead[keep].tobytes()
            gs["frames"] = len(keep)
        if "snapshots" in e:
            fr = e["snapshots"]["frames"]
            fr = [f for i, f in enumerate(fr) if i % snap_step == 0 or i == len(fr) - 1]
            e["snapshots"]["frames"] = fr
            e["snapshots"]["volume"] = e["snapshots"]["volume"] and volume
            tag = e["snapshots"].get("dtype", "i8")
            for f in fr:
                for kind in ("top", "vol") if e["snapshots"]["volume"] else ("top",):
                    p = os.path.join(sd, kind, f"{f['index']:08d}.{tag}.gz")
                    if os.path.exists(p):
                        files[f"{st}/{kind}/{f['index']:08d}.{tag}.gz"] = gzip.decompress(open(p, "rb").read())
    files["manifest.json"] = json.dumps(man, separators=(",", ":")).encode()
    emb = {k: _b64gz(v) for k, v in files.items()}

    html = open(os.path.join(VIEWER, "index.html")).read()
    style = open(os.path.join(VIEWER, "style.css")).read()
    script = open(os.path.join(VIEWER, "app.js")).read()
    html = html.replace('<link rel="stylesheet" href="style.css">', f"<style>\n{style}</style>")
    html = html.replace('<script src="app.js"></script>',
                        "<script>window.HV_EMBED=" + json.dumps(emb, separators=(",", ":")) + ";</script>\n"
                        f"<script>\n{script}</script>")
    html = re.sub(r"<title>.*?</title>", f"<title>{man['name']} HAMR Run</title>", html)
    if embedded_page:   # localhost.ai artifact: the host supplies the document skeleton
        html = re.sub(r"<!doctype html>|</?html[^>]*>|</?head>|</?body>|<meta charset[^>]*>|<meta name=\"viewport\"[^>]*>", "", html)
    open(out, "w").write(html)
    mb = os.path.getsize(out) / 1e6
    biggest = sorted(((len(v), k) for k, v in emb.items()), reverse=True)[:4]
    print(f"{out}: {mb:.1f} MB (largest: " + ", ".join(f"{k} {n / 1e6:.1f} MB" for n, k in biggest) + ")")
    if mb > limit_mb:
        raise SystemExit(f"{out} is {mb:.1f} MB, over {limit_mb} MB: raise --grain-step or --snap-step, or use --no-volume")

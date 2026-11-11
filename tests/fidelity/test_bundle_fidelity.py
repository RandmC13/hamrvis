"""Bundle layer of the data fidelity checks on the test scenarios: every array and manifest number of their bundles
against the raw files, by hamrvis/verify/check.py (which compares them with hamrvis/verify/reference.py, an
independent reader of VAMPIRE's files; the tolerances are listed there). The same checks run on a bundle of real runs
with: hamrvis verify <bundle>
Run: python3 -m unittest discover tests  (from the repository root)
"""
import ast
import contextlib
import io
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from hamrvis.verify import check as bundle_check         # noqa: E402  (independent of the rest of hamrvis)
import test_hamrvis                                      # noqa: E402  (its window snapshot series; module only, so its tests are not collected here)
from hamrvis import reduce                               # noqa: E402  (only to make the bundles)

SCEN = os.path.join(os.path.dirname(HERE), "scenarios")
CASES = {   # bundle name: (scenario, reduce options)
    "chain": ("chain", {}), "two-tracks": ("two-tracks", {}), "continuous": ("continuous", {}),
    "text-snapshots": ("text-snapshots", {}), "custom-names": ("custom-names", {}),
    "file-per-process": ("file-per-process", {}), "split-write": ("split-write", {}), "rewrite": ("rewrite", {}),
    "p16": ("write-only", dict(precision=16)), "p32": ("write-only", dict(precision=32, grain_precision=32)),
    "coarse": ("chain", dict(top_max=300, layer_max=800)),     # columns of several sites
    "thinned": ("write-only", dict(grain_every=3)),
}


class BundleAgainstRawFiles(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.bundles = {}
        for name, (scen, opts) in CASES.items():
            cls.bundles[name] = os.path.join(cls.tmp, name)
            reduce.reduce([os.path.join(SCEN, scen)], cls.bundles[name], workers=1, opts=opts)
        W = type("W", (test_hamrvis.WindowSnapshots,), {})     # a series of window snapshots, as VAMPIRE writes them
        W.setUpClass()
        cls.window_src = W
        cls.bundles["window"] = os.path.join(cls.tmp, "window")
        reduce.reduce([W.window], cls.bundles["window"], workers=1)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)
        cls.window_src.tearDownClass()

    def test_every_check_passes(self):
        done, full = set(), set()
        for name, b in self.bundles.items():
            man = json.load(open(os.path.join(b, "manifest.json"), encoding="utf-8"))
            for e in man["runs"]:
                res = bundle_check.Results()
                bundle_check.check_run(b, e, res)
                with self.subTest(bundle=name, run=e["id"]):
                    self.assertEqual(res.failed(), [])
                for r in res.rows:
                    done.add((name, r[1].split(" (")[0], r[6]))
                    full.add((name, r[1]))
        what = {w for _, w, _ in done}
        for q in ("output table", "head parameters", "grain list", "grain series", "magnetic materials", "grid",
                  "atoms per column", "atoms per layer cell", "grain map", "mixed columns", "grain centres",
                  "magnetic atoms per grain", "snapshot list", "snapshot columns", "snapshot layers", "lowest column |m|",
                  "dipole field", "atoms per dipole cell"):
            self.assertIn(q, what)
        names = {n for n, _, _ in done}
        self.assertTrue({"p16", "p32", "window", "coarse", "thinned", "split-write", "file-per-process", "text-snapshots"} <= names)
        self.assertIn(("coarse", "grid", "(top map coarsened by an option)"), done)
        self.assertIn(("split-write", "output table (joined jobs)"), full)
        self.assertTrue(any(n == "window" and w.startswith("snapshot columns") for n, w in full))

    def test_a_changed_bundle_fails(self):
        """The checks see a wrong value: one stored column of the chain write changed by one step."""
        import gzip
        src = os.path.join(self.bundles["chain"], "write", "top")
        f = sorted(os.listdir(src))[3]
        raw = bytearray(gzip.open(os.path.join(src, f)).read())
        b = os.path.join(self.tmp, "changed")
        shutil.copytree(self.bundles["chain"], b)
        k = next(i for i in range(len(raw) // 3 * 2, len(raw)) if raw[i] not in (0x80, 0x7f))   # an m_z value, not empty
        raw[k] = (raw[k] + 1) % 256
        with gzip.open(os.path.join(b, "write", "top", f), "wb") as out:
            out.write(bytes(raw))
        e = next(r for r in json.load(open(os.path.join(b, "manifest.json")))["runs"] if r["id"] == "write")
        res = bundle_check.Results()
        bundle_check.check_run(b, e, res)
        self.assertEqual([r[1].split(" (")[0] for r in res.failed()], ["snapshot columns"])


    def verify(self, b, *args):
        """hamrvis verify on bundle b; returns (exit status, printed text)."""
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            st = bundle_check.main([b, *args])
        return st, out.getvalue()

    def test_verify_stores_its_result(self):
        b = os.path.join(self.tmp, "verified")
        shutil.copytree(self.bundles["chain"], b)
        st, txt = self.verify(b)
        v = json.load(open(os.path.join(b, "manifest.json")))["verify"]
        self.assertEqual(st, 0, txt)
        self.assertEqual((v["failed"], v["every"], v["runs"]), (0, 1, None))
        self.assertGreater(v["checks"], 20)
        self.assertIn(f"{v['checks']} checks, 0 failed", txt)
        st, _ = self.verify(b, "--run", "write", "--every", "2")
        v = json.load(open(os.path.join(b, "manifest.json")))["verify"]
        self.assertEqual((st, v["runs"], v["every"]), (0, ["write"], 2))
        self.assertEqual(self.verify(b, "--run", "nosuchrun")[0], 2)
        self.assertEqual(self.verify(os.path.join(self.tmp, "no-such-bundle"))[0], 2)

    def test_verify_fails_on_a_changed_bundle(self):
        b = os.path.join(self.tmp, "changed-output")
        shutil.copytree(self.bundles["chain"], b)
        f = os.path.join(b, "write", "output.f64")
        raw = bytearray(open(f, "rb").read())
        raw[-1] ^= 0x01                                   # last bit of the last value's float64
        open(f, "wb").write(bytes(raw))
        st, txt = self.verify(b, "--run", "write")
        self.assertEqual(st, 1, txt)
        self.assertEqual(json.load(open(os.path.join(b, "manifest.json")))["verify"]["failed"], 1)

    def test_verify_says_when_the_runs_are_missing(self):
        b = os.path.join(self.tmp, "moved")
        shutil.copytree(self.bundles["chain"], b)
        mp = os.path.join(b, "manifest.json")
        man = json.load(open(mp))
        man["runs"][0]["segments"] = ["/nonexistent/run"]
        json.dump(man, open(mp, "w"))
        st, txt = self.verify(b)
        self.assertEqual(st, 2)
        self.assertIn("/nonexistent/run", txt)
        self.assertIn("raw run files", txt)
        self.assertNotIn("verify", json.load(open(mp)))     # nothing stored without a check

    def test_verify_is_independent_of_hamrvis(self):
        """hamrvis/verify imports nothing from the rest of hamrvis but its version."""
        d = os.path.join(os.path.dirname(HERE), "..", "hamrvis", "verify")
        for f in os.listdir(d):
            if f.endswith(".py"):
                for n in ast.walk(ast.parse(open(os.path.join(d, f)).read())):
                    if isinstance(n, ast.ImportFrom) and n.level:
                        ok = n.level == 1 or (n.level == 2 and n.module is None and [a.name for a in n.names] == ["__version__"])
                        self.assertTrue(ok, f"{f}: from {'.' * n.level}{n.module or ''} import {[a.name for a in n.names]}")
                    if isinstance(n, (ast.Import, ast.ImportFrom)):
                        mods = [n.module] if isinstance(n, ast.ImportFrom) else [a.name for a in n.names]
                        self.assertFalse(any(m and m.startswith("hamrvis") for m in mods), f)


if __name__ == "__main__":
    unittest.main()

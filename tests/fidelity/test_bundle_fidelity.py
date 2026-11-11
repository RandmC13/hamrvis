"""Bundle layer of the data fidelity checks on the test scenarios: every array and manifest number of their bundles
against the raw files, by tests/fidelity/bundle_check.py (which compares them with tests/fidelity/reference.py, an
independent reader of VAMPIRE's files; the tolerances are listed there). The same checks run on a bundle of real runs
with: python3 tests/fidelity/bundle_check.py <bundle>
Run: python3 -m unittest discover tests  (from the repository root)
"""
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import bundle_check                                      # noqa: E402  (independent of hamrvis)
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


if __name__ == "__main__":
    unittest.main()

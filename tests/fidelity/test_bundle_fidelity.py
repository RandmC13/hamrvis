"""Bundle layer of the data fidelity checks: every stored array of a bundle against tests/fidelity/reference.py,
an independent reader of the raw VAMPIRE files. The only allowed differences are the stored encodings:
    output table        exact (float64)
    grain series        int16: within half a step (0.5/32767); float32: float32 of grain.txt's value
    snapshot columns    int8: 0.5/127, int16: 0.5/32767, float32: float32 rounding; empty columns marked empty
Run: python3 -m unittest discover tests  (from the repository root)
"""
import gzip
import json
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reference as ref                                  # noqa: E402  (independent of hamrvis)
from hamrvis import reduce                               # noqa: E402  (only to make the bundles)

SCEN = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scenarios")
CASES = {   # bundle name: (scenario, reduce options)
    "chain": ("chain", {}), "two-tracks": ("two-tracks", {}), "continuous": ("continuous", {}),
    "text-snapshots": ("text-snapshots", {}), "custom-names": ("custom-names", {}),
    "file-per-process": ("file-per-process", {}), "split-write": ("split-write", {}),
    "p16": ("write-only", dict(precision=16)), "p32": ("write-only", dict(precision=32, grain_precision=32)),
    "coarse": ("chain", dict(top_max=300, layer_max=800)),     # columns of several sites
}
STEP = {"i8": 127, "i16": 32767, "f32": None}


class BundleAgainstRawFiles(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.man = {}
        for name, (scen, opts) in CASES.items():
            out = os.path.join(cls.tmp, name)
            reduce.reduce([os.path.join(SCEN, scen)], out, workers=1, opts=opts)
            cls.man[name] = json.load(open(os.path.join(out, "manifest.json")))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def runs(self):
        for name, man in self.man.items():
            for e in man["runs"]:
                yield name, e, os.path.join(self.tmp, name, e["id"])

    def test_output_table_exact(self):
        for name, e, d in self.runs():
            cols, rows = None, []
            for seg in e["segments"]:
                c, r = ref.output_table(seg)
                cols = cols or c
                rows.append(r)
            stored = np.fromfile(os.path.join(d, "output.f64")).reshape(-1, len(cols))
            self.assertEqual(e["output"]["columns"], cols, msg=name)
            if len(rows) == 1:
                np.testing.assert_array_equal(stored, rows[0], err_msg=name)
            else:   # joined jobs: every stored row is a row of one of the jobs, in time order
                allrows = {tuple(r) for part in rows for r in part}
                self.assertTrue(all(tuple(r) in allrows for r in stored), msg=name)
                self.assertTrue((np.diff(stored[:, cols.index("Real_time")]) > 0).all(), msg=name)

    def test_grain_series(self):
        checked = 0
        for name, e, d in self.runs():
            if "grain_series" not in e or len(e["segments"]) > 1:
                continue
            nG, nF = len(e["grains"]["id"]), e["grain_series"]["frames"]
            want = ref.grain_magnetisation(e["source"], nG)
            self.assertEqual(len(want), nF, msg=name)
            if e["grain_series"].get("dtype") == "f32":
                got = np.fromfile(os.path.join(d, "grains.f32"), np.float32).reshape(nF, nG, 4)
                np.testing.assert_array_equal(got, want.astype(np.float32), err_msg=name)
            else:
                got = np.fromfile(os.path.join(d, "grains.i16"), np.int16).reshape(nF, nG, 4) / 32767
                np.testing.assert_allclose(got, np.clip(want, -1, 1), atol=0.5 / 32767 + 1e-12, rtol=0, err_msg=name)
            checked += 1
        self.assertGreaterEqual(checked, 7)

    def test_snapshot_columns(self):
        checked, coarse = set(), False
        for name, e, d in self.runs():
            if "snapshots" not in e or len(e["segments"]) > 1:
                continue
            G, tag = e["snapshots"]["grid"], e["snapshots"].get("dtype", "i8")
            typ, xyz = ref.coordinates(e["source"])
            self.assertTrue(set(np.unique(typ)) <= set(e["magnetic_types"]), msg=name)
            xs = np.unique(np.round(xyz[:, 0], 3))
            self.assertAlmostEqual(G["atom_pitch_nm"] * 10, np.diff(xs).min(), places=3, msg=name)
            for fr in e["snapshots"]["frames"]:
                factor = int(round(G["pitch_nm"] / G["atom_pitch_nm"]))
                sites = tuple(int(round(L / G["atom_pitch_nm"])) for L in e["size_nm"][:2])
                want = ref.column_means(xyz, ref.spins(e["source"], fr["index"]), G["atom_pitch_nm"] * 10,
                                        G["nx"], G["ny"], factor, sites)
                coarse = coarse or factor > 1
                raw = gzip.open(os.path.join(d, "top", f"{fr['index']:08d}.{tag}.gz")).read()
                q = np.frombuffer(raw, {"i8": np.int8, "i16": np.int16, "f32": np.float32}[tag]).reshape(want.shape)
                empty = np.isnan(want)
                if STEP[tag]:
                    self.assertTrue((q[empty] == -STEP[tag] - 1).all(), msg=name)
                    np.testing.assert_allclose(q[~empty] / STEP[tag], want[~empty], atol=0.5 / STEP[tag] + 1e-9,
                                               rtol=0, err_msg=f"{name} snapshot {fr['index']}")
                else:
                    self.assertTrue(np.isnan(q[empty]).all(), msg=name)
                    np.testing.assert_allclose(q[~empty], want[~empty], atol=1e-7, rtol=0, err_msg=name)
                checked.add((name, tag))
        self.assertTrue({("p16", "i16"), ("p32", "f32"), ("chain", "i8"), ("file-per-process", "i8"),
                         ("text-snapshots", "i8"), ("continuous", "i8")} <= checked)
        self.assertTrue(coarse, msg="no coarsened grid was checked")

    def test_head_parameters_in_any_units(self):
        for name, e, d in self.runs():
            if "hamr" not in e:
                continue
            keys = ref.parse_input(ref.input_path(e["source"]))
            p = e["hamr"]
            self.assertAlmostEqual(p["v"], ref.value(keys, "hamr:head-speed", "velocity", 30e10) / 1e10, places=9)
            self.assertAlmostEqual(p["BL"], ref.value(keys, "hamr:bit-size", "length", 0) / 10, places=9)
            self.assertAlmostEqual(p["dt"] / ref.value(keys, "sim:time-step", "time", 1e-15), 1, places=12)
            self.assertAlmostEqual(p["Lx"], ref.value(keys, "dimensions:system-size-x", "length", 0) / 10, places=9)
            self.assertAlmostEqual(p["Hmax"], ref.value(keys, "sim:maximum-applied-field-strength", "field", 1), places=12)


if __name__ == "__main__":
    unittest.main()

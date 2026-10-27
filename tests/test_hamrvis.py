"""Tests for hamrvis: python3 -m unittest discover tests  (from the repository root)

fixture-run/ is a tiny granular HAMR write made by VAMPIRE (3 bits, 2624 atoms, 11 binary snapshots, output with
B_applied); regenerate it with tests/make_fixture.sh if the formats change. fixture-leadin/ is the same film with a
lead-in (half a bit + NPS) longer than a bit and two tracks: VAMPIRE reports 6 bits but writes 4 (see its log).
"""
import gzip
import json
import os
import re
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

from hamrvis import export, reduce, vampire as vp

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, "fixture-run")
LEADIN = os.path.join(HERE, "fixture-leadin")
RUN_T0 = os.path.join(HERE, "fixture-t0")    # 0 K time series: snapshot and grain.txt share step 0
APP = os.path.join(HERE, "..", "viewer", "app.js")


class FieldModel(unittest.TestCase):
    def test_matches_vampire_output(self):
        for run in (RUN, LEADIN):
            h = vp.Hamr(run)
            cols, out = vp.read_output(run)
            steps = np.rint(out[:, cols.index("Real_time")] / h.dt).astype(int)
            np.testing.assert_allclose(h.field(steps), out[:, cols.index("B_applied")], atol=1e-12, err_msg=run)

    def test_bits_lost_to_a_long_lead_in(self):
        # the last bit of each track is never written and the second track is written on the first (log: y = 1.25 nm)
        h = vp.Hamr(LEADIN)
        self.assertEqual((h.n_tracks, h.bits_per_track, h.written), (2, 3, 2))
        self.assertEqual(h.sequence, [-1, 1, 0, -1, 1, 0])
        self.assertEqual(h.track_y(1), h.track_y(0))
        self.assertAlmostEqual(h.track_y(0), 1.25)
        self.assertEqual(vp.Hamr(RUN).written, 3)

    def test_default_peak_is_hmax_minus_hmin(self):
        # VAMPIRE's defaults are Hmax = +1 T, Hmin = -1 T, so an input that sets neither writes at 2 T
        with tempfile.TemporaryDirectory() as d:
            text = open(os.path.join(RUN, "input")).read()
            text = re.sub(r"sim:(minimum|maximum)-applied-field-strength.*\n", "", text)
            open(os.path.join(d, "input"), "w").write(text)
            self.assertEqual(vp.Hamr(d).Hpeak, 2.0)


class TrackRules(unittest.TestCase):
    """Tracks, bits and sequence as src/hamr/initialize.cpp and bit_sequence.cpp decide them."""

    def hamr(self, **keys):
        text = open(os.path.join(RUN, "input")).read()
        for k, v in keys.items():
            text = re.sub(rf"hamr:{k}\s*=.*\n", "", text) + (f"hamr:{k} = {v}\n" if v is not None else "")
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "input"), "w").write(text)
            return vp.Hamr(d)

    def test_tracks_fit_the_film(self):
        # 12 x 6 nm film, 3 nm bits: 3 bits per track; 2.5 nm tracks fit twice, so 6 bits give 2 tracks
        h = self.hamr(**{"track-size": "2.5 !nm", "number-of-bits": 6})
        self.assertEqual((h.n_tracks, h.bits_per_track, len(h.sequence)), (2, 3, 6))
        # 9 bits still only fit 2 tracks: the sequence is cut to 6
        h = self.hamr(**{"track-size": "2.5 !nm", "number-of-bits": 9})
        self.assertEqual((h.n_tracks, len(h.sequence)), (2, 6))

    def test_user_sequence_padded_with_zero_bits(self):
        h = self.hamr(**{"bit-sequence-type": "user-defined", "bit-sequence": "1,1", "number-of-bits": 3})
        self.assertEqual(h.sequence, [1, 1, 0])
        self.assertEqual(h.field(np.array([h.eq + h.ne + 2 * h.nb + 5]))[0], 0.0)   # field off in a 0 bit

    def test_grain_columns_follow_the_input(self):
        with tempfile.TemporaryDirectory() as d:
            for f in ("grain.txt", "m.mat"):
                shutil.copy(os.path.join(RUN, f), d)
            text = open(os.path.join(RUN, "input")).read().replace("grain:real-time", "grain:time-steps")
            open(os.path.join(d, "input"), "w").write(text)
            cols, m = vp.read_grain_file(d, 9)
            self.assertEqual(list(cols), ["step"])
            self.assertEqual(m.shape[1:], (9, 4))


class Reduce(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.bundle = os.path.join(cls.tmp, "bundle")
        reduce.reduce([RUN], cls.bundle, workers=2)
        cls.man = json.load(open(os.path.join(cls.bundle, "manifest.json")))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_manifest(self):
        e = self.man["stages"]["write"]
        self.assertEqual(len(e["snapshots"]["frames"]), len(vp.snapshots(RUN)))
        self.assertEqual(len(e["grains"]["id"]), 9)
        self.assertEqual(e["hamr"]["bits_per_track"], 3)
        self.assertNotIn("outline", e["grains"])      # shapes come from the atoms, not grain_shapes.txt

    def test_grain_map_matches_atoms(self):
        # every column of the map holds the grain of its magnetic atoms
        e = self.man["stages"]["write"]
        M = e["grain_map"]
        gm = np.fromfile(os.path.join(self.bundle, "write", "grain_map.i16"), np.int16).reshape(M["ny"], M["nx"])
        typ, xyz, g = vp.read_coords(vp.coord_files(RUN)[0])
        mag = np.isin(typ, e["magnetic_types"])
        ix = np.rint(xyz[mag, 0] / (M["pitch_nm"] * 10)).astype(int) % M["nx"]
        iy = np.rint(xyz[mag, 1] / (M["pitch_nm"] * 10)).astype(int) % M["ny"]
        index = {gid: i for i, gid in enumerate(e["grains"]["id"])}
        np.testing.assert_array_equal(gm[iy, ix], [index[v] for v in g[mag]])
        self.assertEqual(M["mixed_columns"], 0)
        self.assertEqual(int((gm >= 0).sum()), len(set(zip(iy, ix))))

    def test_snapshots_match_grain_txt(self):
        # one step apart (snapshots are written one step after the grain row): small but non-zero differences
        c = self.man["stages"]["write"]["checks"]
        self.assertEqual(c["grain_txt_one_step_earlier"]["n"], len(vp.snapshots(RUN)) - 2)
        self.assertLess(c["column_rounding_max"], 0.5 / 127 + 1e-6)

    def test_snapshot_reading_exact_at_same_step(self):
        with tempfile.TemporaryDirectory() as d:
            reduce.reduce([RUN_T0], os.path.join(d, "b"), workers=1)
            c = json.load(open(os.path.join(d, "b", "manifest.json")))["stages"]["relax"]["checks"]
        self.assertGreaterEqual(c["grain_txt_same_step"]["n"], 1)
        self.assertLess(c["grain_txt_same_step"]["max_diff"], 1e-5)   # grain.txt prints 6 significant digits

    def test_column_average_matches_spins(self):
        e = self.man["stages"]["write"]
        G = e["snapshots"]["grid"]
        idx, _, files = vp.snapshots(RUN)[-1]
        q = np.frombuffer(gzip.open(os.path.join(self.bundle, "write", "top", f"{idx:08d}.i8.gz")).read(), np.int8)
        q = q.reshape(3, G["ny"], G["nx"])
        typ, xyz, _ = vp.read_coords(vp.coord_files(RUN)[0])
        s = vp.read_spins(files[0])
        mag = np.isin(typ, e["magnetic_types"])
        ix = np.rint(xyz[mag, 0] / (G["pitch_nm"] * 10)).astype(int) % G["nx"]
        iy = np.rint(xyz[mag, 1] / (G["pitch_nm"] * 10)).astype(int) % G["ny"]
        n = np.zeros((G["ny"], G["nx"]))
        sz = np.zeros((G["ny"], G["nx"]))
        np.add.at(n, (iy, ix), 1)
        np.add.at(sz, (iy, ix), s[mag, 2])
        filled = n > 0
        np.testing.assert_allclose(q[2][filled] / 127, sz[filled] / n[filled], atol=0.5 / 127 + 1e-9)
        self.assertTrue((q[2][~filled] == -128).all())

    def test_analysis(self):
        a = self.man["stages"]["write"]["analysis"]["grains"]
        self.assertEqual(len(a), 9)
        cats = {g["category"] for g in a}
        self.assertTrue(cats <= {"full-field", "reversal", "zero-field", "unheated", "undecided"})
        for g in a:
            if g["category"] in ("full-field", "reversal", "zero-field"):
                self.assertIn("decision_step", g)

    def test_folder_scan_skips_non_runs(self):
        with tempfile.TemporaryDirectory() as d:
            os.symlink(RUN, os.path.join(d, "write"))
            os.makedirs(os.path.join(d, "src"))
            shutil.copy(os.path.join(RUN, "input"), os.path.join(d, "src"))     # an input copy that never ran
            self.assertEqual([os.path.basename(r) for r in reduce.find_runs([d])], ["write"])

    def test_export(self):
        out = os.path.join(self.tmp, "page.html")
        export.export(self.bundle, out, grain_step=2)
        html = open(out).read()
        self.assertIn("window.HV_EMBED=", html)
        self.assertNotIn('src="app.js"', html)

    @unittest.skipUnless(shutil.which("node"), "node not installed")
    def test_viewer_model_matches_python(self):
        # the viewer's JavaScript port of the head / field / temperature model must agree with vampire.py
        src = open(APP).read()
        model = src[src.index("function Model(p)"):src.index("// ---", src.index("function Model(p)"))]
        for run in (RUN, LEADIN):
            h = vp.Hamr(run)
            steps = list(range(0, h.eq + h.laser_steps + 50, 37))
            js = model + f"""
                const m = Model({json.dumps(h.params())});
                console.log(JSON.stringify({json.dumps(steps)}.map(s => {{ const hd = m.head(s);
                  return [hd.x, hd.y, hd.on ? 1 : 0, m.field(s), m.T(4.0, 2.5, s)]; }})));"""
            res = np.array(json.loads(subprocess.run(["node", "-e", js], capture_output=True, text=True, check=True).stdout))
            x, y, on = h.head(steps)
            np.testing.assert_allclose(res[:, 0][on], x[on], atol=1e-9)
            np.testing.assert_allclose(res[:, 1][on], y[on], atol=1e-9)
            np.testing.assert_array_equal(res[:, 2].astype(bool), on)
            np.testing.assert_allclose(res[:, 3], h.field(steps), atol=1e-12)
            np.testing.assert_allclose(res[:, 4], h.temperature(4.0, 2.5, np.array(steps)), atol=1e-9)


if __name__ == "__main__":
    unittest.main()

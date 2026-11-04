"""Tests for hamrvis: python3 -m unittest discover tests  (from the repository root)

fixture-run/ is a tiny granular HAMR write made by VAMPIRE (3 bits, 2624 atoms, 11 binary snapshots, output with
B_applied); regenerate it with tests/make_fixture.sh if the formats change. fixture-leadin/ is the same film with a
lead-in (half a bit + NPS) longer than a bit and two tracks, made by VAMPIRE before its bit loop was fixed: it reports
6 bits but writes 4 (see its log). fixture-leadin-fixed is the same input run by the fixed VAMPIRE: all 6 bits.
scenarios/ holds folders of chained runs (MC relax, writes, reads, a write split over two jobs, ...) made by
tests/make_scenarios.py; the tests reduce each folder as a whole. They also cover two tracks with field-free bits,
a continuous (grain-less) film, text snapshots, one data file per MPI process, and custom file names with other
units. tests/fidelity/ checks bundles against an independent reader of the raw files.
"""
import gzip
import json
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

from hamrvis import export, reduce, runs, vampire as vp

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, "fixture-run")
LEADIN = os.path.join(HERE, "fixture-leadin")
LEADIN_FIXED = os.path.join(HERE, "fixture-leadin-fixed")
RUN_T0 = os.path.join(HERE, "fixture-t0")    # 0 K time series: snapshot and grain.txt share step 0
APP = os.path.join(HERE, "..", "viewer", "app.js")
SCEN = os.path.join(HERE, "scenarios")


class FieldModel(unittest.TestCase):
    def test_matches_vampire_output(self):
        for run in (RUN, LEADIN, LEADIN_FIXED):
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

    def test_fixed_bit_loop_writes_every_bit(self):
        # the same input run by the fixed VAMPIRE: every bit written, the second track at its own y (log: 3.75 nm)
        h = vp.Hamr(LEADIN_FIXED)
        self.assertEqual((h.n_tracks, h.bits_per_track, h.written), (2, 3, 3))
        self.assertEqual(h.sequence, [-1, 1, -1, 1, -1, 1])
        self.assertAlmostEqual(h.track_y(0), 1.25)
        self.assertAlmostEqual(h.track_y(1), 3.75)

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


class WindowSnapshots(unittest.TestCase):
    """Window snapshots (hamr:active-window): a series built from fixture-run in which each snapshot after the first
    takes only the atoms of a few ranges from the run (the rest keep their spins), written once as full snapshots
    and once as window snapshots holding just those ranges, as VAMPIRE writes them."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.full, cls.window = (os.path.join(cls.tmp, d) for d in ("full", "window"))
        for d in (cls.full, cls.window):
            shutil.copytree(RUN, d, ignore=shutil.ignore_patterns("spins-*"))
        cls.truth, state = [], None
        for k, (idx, _, files, _) in enumerate(vp.snapshots(RUN)):
            new = vp.read_spins(files[0])
            n = len(new)
            ranges = [] if state is None else [(k * 37 % n, min(50, n - k * 37 % n)), (n // 2, n // 3)]
            if state is None:
                state = new.copy()
            for a, c in ranges:
                state[a:a + c] = new[a:a + c]
            cls.truth.append(state.copy())
            meta = open(files[0][:-4] + "meta").read()
            name = f"spins-{idx:08d}"
            for d, window in ((cls.full, False), (cls.window, 0 < k < 10)):
                with open(os.path.join(d, name + ".meta"), "w") as f:
                    f.write(meta + "Spin data: window\n" if window else meta)
                with open(os.path.join(d, name + ".data"), "wb") as f:
                    if window:
                        np.array([n, len(ranges)], np.uint64).tofile(f)
                        np.array(ranges, np.uint64).tofile(f)
                        np.concatenate([state[a:a + c] for a, c in ranges]).tofile(f)
                    else:
                        np.array([n], np.uint64).tofile(f)
                        state.tofile(f)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_series_rebuilds_every_snapshot(self):
        snaps = vp.snapshots(self.window)
        self.assertEqual([w for *_, w in snaps], [False] + [True] * 9 + [False] * (len(snaps) - 10))
        for (idx, _, spins), truth in zip(vp.spin_series(snaps), self.truth):
            np.testing.assert_array_equal(spins[0], truth, err_msg=f"snapshot {idx}")

    def test_reduce_same_as_full_snapshots(self):
        out = {}
        for name, run in (("full", self.full), ("window", self.window)):
            out[name] = os.path.join(self.tmp, "bundle-" + name)
            reduce.reduce([run], out[name], workers=2)
        compared = 0
        full, window = os.path.join(out["full"], "full"), os.path.join(out["window"], "window")   # run folders
        for dp, _, files in os.walk(full):
            for f in files:
                a = os.path.join(dp, f)
                b = os.path.join(window, os.path.relpath(a, full))
                read = (lambda p: gzip.decompress(Path(p).read_bytes())) if f.endswith(".gz") else (lambda p: Path(p).read_bytes())
                self.assertEqual(read(a), read(b), msg=os.path.relpath(a, full))
                compared += 1
        self.assertGreater(compared, 10)


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
        e = self.man["runs"][0]
        self.assertEqual(len(e["snapshots"]["frames"]), len(vp.snapshots(RUN)))
        self.assertEqual(len(e["grains"]["id"]), 9)
        self.assertEqual(e["hamr"]["bits_per_track"], 3)
        self.assertNotIn("outline", e["grains"])      # shapes come from the atoms, not grain_shapes.txt

    def test_grain_map_matches_atoms(self):
        # every column of the map holds the grain of its magnetic atoms
        e = self.man["runs"][0]
        M = e["grain_map"]
        gm = np.fromfile(os.path.join(self.bundle, "fixture-run", "grain_map.i16"), np.int16).reshape(M["ny"], M["nx"])
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
        c = self.man["runs"][0]["checks"]
        self.assertEqual(c["grain_txt_one_step_earlier"]["n"], len(vp.snapshots(RUN)) - 2)
        self.assertLess(c["column_rounding_max"], 0.5 / 127 + 1e-6)

    def test_snapshot_reading_exact_at_same_step(self):
        with tempfile.TemporaryDirectory() as d:
            reduce.reduce([RUN_T0], os.path.join(d, "b"), workers=1)
            c = json.load(open(os.path.join(d, "b", "manifest.json")))["runs"][0]["checks"]
        self.assertGreaterEqual(c["grain_txt_same_step"]["n"], 1)
        self.assertLess(c["grain_txt_same_step"]["max_diff"], 1e-5)   # grain.txt prints 6 significant digits

    def test_column_average_matches_spins(self):
        e = self.man["runs"][0]
        G = e["snapshots"]["grid"]
        idx, _, files, _ = vp.snapshots(RUN)[-1]
        q = np.frombuffer(gzip.open(os.path.join(self.bundle, "fixture-run", "top", f"{idx:08d}.i8.gz")).read(), np.int8)
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
        a = self.man["runs"][0]["analysis"]["grains"]
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
            self.assertEqual([os.path.basename(r) for r in runs.candidates([d])], ["write"])

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


class RunDiscovery(unittest.TestCase):
    """Runs are found, ordered and typed from their own files, whatever mix a folder holds."""

    def found(self, name):
        return [(r.id, r.type, [os.path.basename(d) for d in r.dirs]) for r in runs.find([os.path.join(SCEN, name)])]

    def test_types_and_order(self):
        self.assertEqual(self.found("chain"), [("relax", "monte-carlo", ["relax"]), ("write", "hamr", ["write"]),
                                               ("read", "dipole-read", ["read"])])
        self.assertEqual(self.found("write-only"), [("write", "hamr", ["write"])])
        self.assertEqual([t for _, t, _ in self.found("rewrite")], ["hamr", "dipole-read", "hamr", "dipole-read"])
        self.assertEqual(self.found("llg-relax"), [("anneal", "generic", ["anneal"]), ("write", "hamr", ["write"])])

    def test_continued_jobs_are_one_run(self):
        self.assertEqual(self.found("split-write"), [("job1", "hamr", ["job1", "job2"])])

    def test_order_follows_the_logs_not_the_names(self):
        # names sort as read, relax, write; the logs say relax ran first
        self.assertEqual([i for i, _, _ in self.found("chain")], ["relax", "write", "read"])
        order = [r.id for r in runs.find([os.path.join(SCEN, "chain")], order=["read", "write"])]
        self.assertEqual(order, ["read", "write", "relax"])

    def test_resumed_log_counts_bits_once(self):
        # the second job replays its loop and logs bits 1-3 again
        d = [os.path.join(SCEN, "split-write", j) for j in ("job1", "job2")]
        self.assertEqual(vp.Hamr(d[0], d).written, 3)


class Bundles(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.man = {}
        for name in ("chain", "write-only", "split-write", "rewrite", "llg-relax"):
            reduce.reduce([os.path.join(SCEN, name)], os.path.join(cls.tmp, name), workers=1)
            cls.man[name] = json.load(open(os.path.join(cls.tmp, name, "manifest.json")))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_joined_jobs_equal_one_job(self):
        # the same write run uninterrupted and as two jobs (the first killed mid-write, the second continued from
        # its checkpoint, overlapping it) gives the same bundle, byte for byte
        a, b = os.path.join(self.tmp, "write-only", "write"), os.path.join(self.tmp, "split-write", "job1")
        compared = 0
        for dp, _, files in os.walk(a):
            for f in files:
                p = os.path.join(dp, f)
                q = os.path.join(b, os.path.relpath(p, a))
                read = (lambda x: gzip.decompress(Path(x).read_bytes())) if f.endswith(".gz") else (lambda x: Path(x).read_bytes())
                self.assertEqual(read(p), read(q), msg=os.path.relpath(p, a))
                compared += 1
        self.assertGreater(compared, 20)
        ea, eb = self.man["write-only"]["runs"][0], self.man["split-write"]["runs"][0]
        for k in ("hamr", "analysis", "Tc", "output", "grain_series", "snapshots"):
            self.assertEqual(ea[k], eb[k], msg=k)
        self.assertEqual(len(eb["segments"]), 2)

    def test_reads_compare_with_the_write_before_them(self):
        r = self.man["rewrite"]["runs"]
        self.assertEqual([e.get("ref") for e in r], [None, "1-write", None, "3-rewrite"])
        self.assertNotIn("ref", self.man["chain"]["runs"][0])     # the relax comes before any write
        self.assertTrue(all("dipole" in e for e in r if e["type"] == "dipole-read"))

    def test_time_axis(self):
        self.assertEqual([e["monte_carlo"] for e in self.man["chain"]["runs"]], [True, False, False])

    def test_freezing_temperature_above_laser_minimum(self):
        # decisions far from the laser (thermal flips of small grains) used to pull the estimate to Tmin
        e = self.man["chain"]["runs"][1]
        self.assertGreater(e["Tc"], e["hamr"]["Tmin"] + 50)

    def test_reduce_again_replaces_the_bundle(self):
        out = os.path.join(self.tmp, "again")
        reduce.reduce([os.path.join(SCEN, "rewrite")], out, workers=1)
        reduce.reduce([os.path.join(SCEN, "write-only")], out, workers=1)
        self.assertEqual(sorted(os.listdir(out)), ["manifest.json", "write"])
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "notes.txt"), "w").write("not a bundle")
            with self.assertRaises(SystemExit):
                reduce.reduce([os.path.join(SCEN, "write-only")], d, workers=1)
            self.assertTrue(os.path.exists(os.path.join(d, "notes.txt")))

    def test_export_takes_every_run(self):
        out = os.path.join(self.tmp, "rewrite.html")
        export.export(os.path.join(self.tmp, "rewrite"), out)
        html = open(out).read()
        for rid in ("1-write", "2-read", "3-rewrite", "4-read"):
            self.assertIn(f'"{rid}/output.f64"', html)

    @unittest.skipUnless(shutil.which("node"), "node not installed")
    def test_viewer_reads_version_1_manifests(self):
        # bundles reduced before the run list keep their three tabs, titles and links (#relax, #write, #read)
        src = open(APP).read()
        fn = src[src.index("function normalise(M)"):src.index("\n}\n", src.index("function normalise(M)")) + 3]
        v1 = {"version": 1, "stages": {
            "read": {"input": {"sim:integrator": "llg-heun"}, "dipole": {}},
            "relax": {"input": {"sim:integrator": "monte-carlo"}},
            "write": {"input": {"sim:integrator": "llg-heun"}, "hamr": {}}}}
        js = fn + f"""const M = normalise({json.dumps(v1)});
            console.log(JSON.stringify(M.runs.map(e => [e.id, e.title, e.type, e.ref || null, e.monte_carlo])));"""
        res = json.loads(subprocess.run(["node", "-e", js], capture_output=True, text=True, check=True).stdout)
        self.assertEqual(res, [["relax", "Relaxation", "monte-carlo", None, True], ["write", "Write", "hamr", None, False],
                               ["read", "Read-back", "dipole-read", "write", False]])


def copy_run(src, d, edit=None):
    """Copy run directory src to d (files only); edit(input text) -> new input text."""
    shutil.copytree(src, d)
    if edit:
        name = vp.input_name(d)
        text = open(os.path.join(d, name)).read()
        open(os.path.join(d, name), "w").write(edit(text))
    return d


class Parsing(unittest.TestCase):
    """VAMPIRE's file names, units and missing or damaged files."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)
        vp.NAMES.update(input=None, output=None)

    def write(self, text):
        d = os.path.join(self.tmp, "r")
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, "input"), "w").write(text)
        return vp.Input(d)

    def test_units_as_vampire_converts_them(self):
        inp = self.write("a:x = 0.012 !um\nb:t = 0.001 !ps\nc:h = 15 !kOe\nd:v = 1.5e13\ne:v = 1500 !m/s\n"
                         "f:x = 30\ng:n = 7 !fs\nh:t = 300\n")
        self.assertAlmostEqual(inp.num("a:x"), 12.0)                       # nm
        self.assertAlmostEqual(inp.num("b:t", kind="time"), 1e-15)
        self.assertAlmostEqual(inp.num("c:h", kind="field"), 1.5)
        self.assertAlmostEqual(inp.num("d:v", kind="velocity"), 1500.0)    # no unit: VAMPIRE's A/s
        self.assertAlmostEqual(inp.num("e:v", kind="velocity"), 1500.0)
        self.assertAlmostEqual(inp.num("f:x"), 3.0)                        # no unit: A
        self.assertEqual(inp.num("g:n", kind="int"), 7)                    # VAMPIRE ignores units on integers
        self.assertEqual(inp.num("h:t", kind="none"), 300.0)

    def test_bad_units_say_what_is_wrong(self):
        inp = self.write("a:x = 3 !furlong\nb:x = 3 !ns\nc:t = 300 !K\n")
        with self.assertRaisesRegex(ValueError, "not one of VAMPIRE's units"):
            inp.num("a:x")
        with self.assertRaisesRegex(ValueError, "reads b:x as a length, but 'ns' is a time unit"):
            inp.num("b:x")
        with self.assertRaisesRegex(ValueError, "not one of VAMPIRE's units"):   # VAMPIRE rejects !K too
            inp.num("c:t", kind="none")

    def test_keys_as_vampire_reads_them(self):
        inp = self.write("sim : time-step\t= 1 !fs # comment\nhamr:bit-size = 1 !nm\nhamr:bit-size = 2 !nm\n")
        self.assertEqual(inp.keys["sim:time-step"], ("1", "fs"))
        self.assertEqual(inp.num("hamr:bit-size"), 2.0)
        self.assertEqual(list(inp.keys)[-1], "hamr:bit-size")

    def test_custom_file_names(self):
        d = os.path.join(SCEN, "custom-names", "write")
        self.assertEqual((vp.input_name(d), vp.output_name(d)), ("write.in", "results.txt"))   # log, header
        self.assertEqual([os.path.basename(r) for r in runs.candidates([os.path.dirname(d)])], ["write"])
        vp.NAMES.update(input="nothing.in")
        self.assertEqual(vp.input_name(d), "write.in")         # a name that is not there: found as before

    def test_continued_jobs_with_a_custom_output_name(self):
        # VAMPIRE appends a continued job's rows to its table without the header, so only --output-file finds it
        root = os.path.join(self.tmp, "split")
        shutil.copytree(os.path.join(SCEN, "split-write"), root)
        for j in ("job1", "job2"):
            os.rename(os.path.join(root, j, "output"), os.path.join(root, j, "results.txt"))
        self.assertEqual([r.dirs for r in runs.find([root])], [[os.path.join(root, "job1")]])
        vp.NAMES.update(output="results.txt")
        found = runs.find([root])
        self.assertEqual([[os.path.basename(d) for d in r.dirs] for r in found], [["job1", "job2"]])
        a, b = os.path.join(self.tmp, "a"), os.path.join(self.tmp, "b")
        reduce.reduce([root], a, workers=1)
        vp.NAMES.update(output=None)
        reduce.reduce([os.path.join(SCEN, "split-write")], b, workers=1)
        for f in ("output.f64", "grains.i16", "grain_lead.f64"):
            self.assertEqual(Path(a, "job1", f).read_bytes(), Path(b, "job1", f).read_bytes(), msg=f)

    def test_custom_names_and_units_give_the_same_bundle(self):
        # the same write as write-only, run with --input-file / --output-file and the lengths, times, fields and
        # speed in other units (VAMPIRE's own output is identical)
        a, b = os.path.join(self.tmp, "a"), os.path.join(self.tmp, "b")
        reduce.reduce([os.path.join(SCEN, "write-only")], a, workers=1)
        reduce.reduce([os.path.join(SCEN, "custom-names")], b, workers=1)
        ea, eb = (json.load(open(os.path.join(x, "manifest.json")))["runs"][0] for x in (a, b))
        for k in ("hamr", "analysis", "Tc", "dt", "size_nm", "snapshots", "grain_series", "output", "checks"):
            self.assertEqual(ea[k], eb[k], msg=k)
        for f in os.listdir(os.path.join(a, "write")):
            if os.path.isfile(os.path.join(a, "write", f)):
                self.assertEqual(Path(a, "write", f).read_bytes(), Path(b, "write", f).read_bytes(), msg=f)

    def test_missing_material_file(self):
        d = copy_run(os.path.join(SCEN, "write-only", "write"), os.path.join(self.tmp, "w"))
        os.remove(os.path.join(d, "m.mat"))
        self.assertIsNone(reduce.magnetic_types(d))
        reduce.reduce([d], os.path.join(self.tmp, "b"), workers=1)
        e = json.load(open(os.path.join(self.tmp, "b", "manifest.json")))["runs"][0]
        self.assertIn("snapshots", e)
        self.assertEqual(len(e["grains"]["magnetic_atoms"]), 9)

    def test_continuous_film_has_atom_views_only(self):
        reduce.reduce([os.path.join(SCEN, "continuous")], os.path.join(self.tmp, "b"), workers=1)
        e = json.load(open(os.path.join(self.tmp, "b", "manifest.json")))["runs"][0]
        self.assertEqual(e["grains"]["id"], [])
        self.assertNotIn("grain_series", e)
        self.assertNotIn("grain_map", e)
        self.assertEqual(len(e["snapshots"]["frames"]), 11)
        self.assertTrue(e["notes"])

    def test_unusable_grain_file_skips_grain_views(self):
        edit = lambda t: t.replace("grain:magnetisation", "grain:material-magnetisation\ngrain:magnetisation")
        d = copy_run(os.path.join(SCEN, "write-only", "write"), os.path.join(self.tmp, "w"), edit)
        reduce.reduce([d], os.path.join(self.tmp, "b"), workers=1)
        e = json.load(open(os.path.join(self.tmp, "b", "manifest.json")))["runs"][0]
        self.assertNotIn("grain_series", e)
        self.assertIn("snapshots", e)
        self.assertIn("material-magnetisation", " ".join(e["notes"]))

    def test_truncated_text_snapshot(self):
        d = copy_run(os.path.join(SCEN, "text-snapshots", "write"), os.path.join(self.tmp, "w"))
        p = os.path.join(d, "spins-00000005.data")
        data = open(p).read()
        open(p, "w").write(data[: len(data) // 2])
        with self.assertRaisesRegex(IOError, "truncated or damaged text file"):
            vp.read_spins(p)
        reduce.reduce([d], os.path.join(self.tmp, "b"), workers=1)
        e = json.load(open(os.path.join(self.tmp, "b", "manifest.json")))["runs"][0]
        self.assertEqual(e["snapshots"]["skipped"], [5])

    def test_half_written_last_row(self):
        d = copy_run(os.path.join(SCEN, "write-only", "write"), os.path.join(self.tmp, "w"))
        for f in ("output", "grain.txt"):
            text = open(os.path.join(d, f)).read()
            open(os.path.join(d, f), "w").write(text[: -len(text.splitlines()[-1]) // 2 - 1])
        cols, out = vp.read_output(d)
        full = vp.read_output(os.path.join(SCEN, "write-only", "write"))[1]
        np.testing.assert_array_equal(out, full[:-1])
        lead, m = vp.read_grain_file(d, 9)
        self.assertEqual(len(m), len(vp.read_grain_file(os.path.join(SCEN, "write-only", "write"), 9)[1]) - 1)

    def test_dipole_cells_default_to_vampire_size(self):
        # cells:macro-cell-size defaults to 10 A in VAMPIRE; the read was run with 1 nm cells written out
        edit = lambda t: re.sub(r"cells:macro-cell-size.*\n", "", t)
        d = copy_run(os.path.join(SCEN, "chain", "read"), os.path.join(self.tmp, "r"), edit)
        e = dict(size_nm=[12.0, 6.0, 1.2])
        os.makedirs(os.path.join(self.tmp, "out"))
        reduce.reduce_dipole(d, os.path.join(self.tmp, "out"), e)
        self.assertEqual(e["dipole"], dict(nx=13, ny=7, nz=2, cell_nm=1.0))


class ReduceOptions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def make(self, name, **opts):
        out = os.path.join(self.tmp, name)
        reduce.reduce([os.path.join(SCEN, "write-only")], out, workers=1, opts=opts)
        return out, json.load(open(os.path.join(out, "manifest.json")))["runs"][0]

    def test_defaults_add_nothing(self):
        _, e = self.make("default")
        self.assertNotIn("dtype", e["snapshots"])
        self.assertNotIn("dtype", e["grain_series"])
        self.assertNotIn("every", e["grain_series"])

    def test_precision(self):
        for bits, tag in ((16, "i16"), (32, "f32")):
            out, e = self.make(f"p{bits}", precision=bits, grain_precision=32 if bits == 32 else 16)
            self.assertEqual(e["snapshots"]["dtype"], tag)
            self.assertTrue(os.path.exists(os.path.join(out, "write", "top", f"00000000.{tag}.gz")))
            self.assertLess(e["checks"]["column_rounding_max"], 0.5 / 32767 + 1e-7)
        self.assertEqual(e["grain_series"]["dtype"], "f32")
        self.assertTrue(os.path.exists(os.path.join(out, "write", "grains.f32")))
        export.export(out, os.path.join(self.tmp, "p32.html"))
        html = open(os.path.join(self.tmp, "p32.html")).read()
        self.assertIn('"write/grains.f32"', html)
        self.assertIn('"write/top/00000000.f32.gz"', html)

    def test_grain_every_and_no_snapshots(self):
        _, full = self.make("full")
        _, e = self.make("thin", grain_every=4, snapshots=False)
        n = full["grain_series"]["frames"]
        self.assertEqual(e["grain_series"]["frames"], len(range(0, n - 1, 4)) + 1)
        self.assertNotIn("snapshots", e)
        self.assertIn("grain_map", e)                       # grain shapes still come from the atoms

    def test_thresholds_and_grid_limits(self):
        _, e = self.make("th", thresholds=dict(heated=0.0), top_max=200)
        self.assertEqual(e["analysis"]["thresholds"]["heated"], 0.0)
        self.assertEqual({g["category"] for g in e["analysis"]["grains"]}, {"unheated"})
        self.assertEqual(e["snapshots"]["grid"]["nx"], 20)  # 40 x 20 columns > 200 pixels: 2 x 2 columns per pixel

    def test_inspect_every_scenario(self):
        import contextlib, io
        from hamrvis import inspect
        for name in sorted(os.listdir(SCEN)):
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                inspect.inspect([os.path.join(SCEN, name)])
            text = buf.getvalue()
            self.assertIn("estimated bundle", text, msg=name)
            if name == "continuous":
                self.assertIn("grain views   no", text)
            if name == "custom-names":
                self.assertIn("'write.in', output table 'results.txt'", text)

    def test_serve_refuses_a_busy_port(self):
        import socket
        from hamrvis import serve
        out, _ = self.make("serve")
        with socket.socket() as s:
            s.bind(("127.0.0.1", 0))
            s.listen()
            with self.assertRaisesRegex(SystemExit, "another viewer or program is using it"):
                serve.serve(out, s.getsockname()[1], browser=False)

    @unittest.skipUnless(shutil.which("node"), "node not installed")
    def test_viewer_decodes_every_precision(self):
        src = open(APP).read()
        code = src[src.index("const SNAP = "):src.index("// one snapshot array")]
        js = code + """const raw = { i8: [0, 127, -127, -128], i16: [0, 32767, -32767, -32768], f32: [0, 0.25, -1, NaN] };
            const out = Object.entries(SNAP).map(([k, [T, sc, E]]) => { const a = new T(raw[k]); a.S = sc; a.E = E;
              return [0, 1, 2, 3].map(i => snapVal(a, i)); });
            console.log(JSON.stringify(out.map(r => r.map(v => v !== v ? null : v))));"""
        res = json.loads(subprocess.run(["node", "-e", js], capture_output=True, text=True, check=True).stdout)
        self.assertEqual(res, [[0, 1, -1, None], [0, 1, -1, None], [0, 0.25, -1, None]])


class NewScenarios(unittest.TestCase):
    def test_head_field_matches_vampire_in_every_write(self):
        h = vp.Hamr(os.path.join(SCEN, "two-tracks", "write"))
        self.assertEqual((h.n_tracks, h.bits_per_track, h.written, h.sequence), (2, 3, 3, [1, 0, -1, 0, 1, -1]))
        # '{1,1,-1}' is not VAMPIRE syntax: VAMPIRE reads the first bit as 0 (its stdout: "0 1 -1")
        self.assertEqual(vp.Hamr(os.path.join(SCEN, "rewrite", "3-rewrite")).sequence, [0, 1, -1])
        n = 0
        for r in (r for name in os.listdir(SCEN) for r in runs.find([os.path.join(SCEN, name)]) if r.type == "hamr"):
            h = vp.Hamr(r.dir, r.dirs)
            for d in r.dirs:
                cols, out = vp.read_output(d)
                steps = np.rint(out[:, cols.index("Real_time")] / h.dt).astype(int)
                np.testing.assert_allclose(h.field(steps), out[:, cols.index("B_applied")], atol=1e-12, err_msg=d)
            n += 1
        self.assertGreaterEqual(n, 10)

    def test_file_per_process(self):
        d = os.path.join(SCEN, "file-per-process", "write")
        self.assertEqual(len(vp.coord_files(d)), 2)
        self.assertEqual(len(vp.snapshots(d)[0][2]), 2)


if __name__ == "__main__":
    unittest.main()

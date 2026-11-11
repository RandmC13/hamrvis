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
import contextlib
import functools
import gzip
import http.server
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np

from hamrvis import export, reduce, runs, serve, vampire as vp

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, "fixture-run")
LEADIN = os.path.join(HERE, "fixture-leadin")
LEADIN_FIXED = os.path.join(HERE, "fixture-leadin-fixed")
RUN_T0 = os.path.join(HERE, "fixture-t0")    # 0 K time series: snapshot and grain.txt share step 0
APP = os.path.join(HERE, "..", "hamrvis", "viewer", "app.js")
SCEN = os.path.join(HERE, "scenarios")


def read_coords(path):
    """(type, xyz, grain) of a whole coordinate file."""
    parts = list(vp.coord_chunks(path))
    return tuple(np.concatenate([p[i] for p in parts]) for i in range(3))


def read_grain_file(run_dir, n_grains):
    """({name: column}, float64 [frames, grains, 4]) of a whole grain.txt."""
    names, rows = vp.grain_rows(run_dir, n_grains)
    d = np.array(list(rows))
    return {n: d[:, i] for i, n in enumerate(names)}, d[:, len(names):].reshape(len(d), n_grains, 4)


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
            cols, m = read_grain_file(d, 9)
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
        typ, xyz, g = read_coords(vp.coord_files(RUN)[0])
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
        typ, xyz, _ = read_coords(vp.coord_files(RUN)[0])
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
        self.assertTrue(cats <= {"full-field", "reversal", "zero-field", "unheated", "undecided", "empty"})
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
        html = open(out, encoding="utf-8").read()
        self.assertIn("window.HV_EMBED=", html)
        self.assertNotIn('src="app.js"', html)

    @unittest.skipUnless(shutil.which("node"), "node not installed")
    def test_viewer_model_matches_python(self):
        # the viewer's JavaScript port of the head / field / temperature model must agree with vampire.py
        src = open(APP, encoding="utf-8").read()
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

    def test_head_checked_against_vampire_at_reduce_time(self):
        n = 0
        for man in self.man.values():
            for e in (e for e in man["runs"] if e["type"] == "hamr"):
                hc = e["head_check"]
                self.assertEqual((hc["B_applied"]["bad"], hc["positions"]["bad"]), (0, 0), e["source"])
                self.assertGreater(hc["positions"]["lines"], 0)
                self.assertEqual(e["hamr"]["written_source"], "log")
                n += 1
        self.assertGreaterEqual(n, 5)

    @unittest.skipUnless(shutil.which("node"), "node not installed")
    def test_viewer_scores_grains_as_the_reducer(self):
        # the viewer rescores the write at the lag in use; at the bundle's own lag it must give the reducer's numbers
        manifests = list(self.man.values())
        for name in ("h3", "h3-v1", "2x11", "h3-scale-window", "h3-scale-full"):   # real-run bundles in bundles/ (git-ignored), when present
            mp = os.path.join(HERE, "..", "bundles", name, "manifest.json")
            if os.path.exists(mp):
                manifests.append(json.load(open(mp, encoding="utf-8")))
        n = 0
        for man in manifests:
            for e in man.get("runs") or man["stages"].values():
                if e.get("analysis"):
                    self.assertEqual(viewer_scores(e, e["freeze_lag_nm"]), [[a["bit"], a["intended"], a["correct"]] for a in e["analysis"]["grains"]])
                    n += 1
        self.assertGreaterEqual(n, 5)

    @unittest.skipUnless(shutil.which("node"), "node not installed")
    def test_viewer_rescores_at_another_lag(self):
        # at another lag the viewer's scoring must equal write_analysis run again with that lag (here via --tc)
        h = vp.Hamr(os.path.join(SCEN, "write-only", "write"))
        e = next(e for e in self.man["write-only"]["runs"] if e["type"] == "hamr")
        Tc = h.Tmin + (h.Tmax - h.Tmin) * 0.9
        out = os.path.join(self.tmp, "tc")
        reduce.reduce([os.path.join(SCEN, "write-only")], out, workers=1, Tc=Tc, opts=dict(snapshots=False))
        f = next(f for f in json.load(open(os.path.join(out, "manifest.json")))["runs"] if f["type"] == "hamr")
        self.assertNotEqual(f["freeze_lag_nm"], e["freeze_lag_nm"])
        self.assertEqual(f["Tc_kind"], "given")
        self.assertEqual(viewer_scores(e, f["freeze_lag_nm"]), [[a["bit"], a["intended"], a["correct"]] for a in f["analysis"]["grains"]])

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
        html = open(out, encoding="utf-8").read()
        for rid in ("1-write", "2-read", "3-rewrite", "4-read"):
            self.assertIn(f'"{rid}/output.f64"', html)

    @unittest.skipUnless(shutil.which("node"), "node not installed")
    def _run_with_input(self, tmp, edit):
        """A copy of a write's input with edit(lines) applied, as a run directory."""
        os.makedirs(tmp, exist_ok=True)
        lines = open(os.path.join(SCEN, "write-only", "write", "input")).read().split("\n")
        with open(os.path.join(tmp, "input"), "w") as f:
            f.write("\n".join(edit(lines)))
        return tmp

    def test_defaults_are_vampires(self):
        # sim::Tmax starts at 300 K (src/simulate/sim.cpp): a write that does not set it is heated to 300 K
        with tempfile.TemporaryDirectory() as tmp:
            h = vp.Hamr(self._run_with_input(tmp, lambda ls: [l for l in ls if "maximum-temperature" not in l]))
            self.assertEqual(h.Tmax, 300.0)
            self.assertEqual(h.Tmin, vp.Hamr(os.path.join(SCEN, "write-only", "write")).Tmin)

    def test_system_size_sets_every_axis(self):
        # dimensions:system-size sets all three, -x/-y/-z one, the last line wins (src/create/interface.cpp)
        no_size = lambda ls: [l for l in ls if "system-size" not in l]
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(vp.Input(self._run_with_input(tmp, no_size)).system_size(), [10.0, 10.0, 10.0])
            h = vp.Hamr(self._run_with_input(tmp, lambda ls: no_size(ls) + ["dimensions:system-size = 12 !nm"]))
            self.assertEqual((h.Lx, h.Ly, h.Lz), (12.0, 12.0, 12.0))
            h = vp.Hamr(self._run_with_input(tmp, lambda ls: no_size(ls) + ["dimensions:system-size = 12 !nm",
                                                                            "dimensions:system-size-z = 3 !nm"]))
            self.assertEqual((h.Lx, h.Ly, h.Lz), (12.0, 12.0, 3.0))
            h = vp.Hamr(self._run_with_input(tmp, lambda ls: no_size(ls) + ["dimensions:system-size-z = 3 !nm",
                                                                            "dimensions:system-size = 12 !nm"]))
            self.assertEqual((h.Lx, h.Ly, h.Lz), (12.0, 12.0, 12.0))

    def test_write_without_bit_or_track_size_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self._run_with_input(tmp, lambda ls: [l for l in ls if "hamr:bit-size" not in l])
            with self.assertRaisesRegex(ValueError, "hamr:bit-size"):
                vp.Hamr(run)

    def test_empty_output_table_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self._run_with_input(tmp, lambda ls: ls)
            with open(os.path.join(run, "output"), "w") as f:
                f.write("#----\n# Output file for vampire simulation\n")
            with self.assertRaisesRegex(ValueError, "no rows"):
                vp.read_output(run)

    def test_server_stays_inside_the_bundle(self):
        h = serve.Handler.__new__(serve.Handler)
        with tempfile.TemporaryDirectory() as tmp:
            h.bundle = os.path.abspath(tmp)
            self.assertEqual(h.translate_path("/bundle/write/output.f64"), os.path.join(h.bundle, "write", "output.f64"))
            for bad in ("/bundle/../secret", "/bundle/a/../../secret", "/bundle/" + os.path.abspath(os.sep) + "secret",
                        "/bundle/C:/secret"):
                path = h.translate_path(bad)
                self.assertTrue(path.startswith(h.bundle + os.sep), path)
                self.assertFalse(os.path.exists(path))

    def test_test_data_has_no_machine_details(self):
        for root in [SCEN] + [os.path.join(HERE, d) for d in os.listdir(HERE) if d.startswith("fixture-")]:
            for d, _, files in os.walk(root):
                for f in files:
                    if f in ("log", "stdout", "output", "results.txt"):
                        data = open(os.path.join(d, f), "rb").read()
                        self.assertNotIn(b"/home/", data, os.path.join(d, f))

    def test_version_is_the_same_everywhere(self):
        import hamrvis
        self.assertIn(f'HAMRVIS_VERSION = "{hamrvis.__version__}"', open(APP, encoding="utf-8").read())
        self.assertIn(f"version: {hamrvis.__version__}", open(os.path.join(HERE, "..", "CITATION.cff")).read())
        with tempfile.TemporaryDirectory() as tmp:
            reduce.reduce([os.path.join(SCEN, "write-only")], tmp + "/b", workers=1)
            self.assertEqual(json.load(open(tmp + "/b/manifest.json"))["hamrvis_version"], hamrvis.__version__)

    def test_viewer_reads_version_1_manifests(self):
        # bundles reduced before the run list keep their three tabs, titles and links (#relax, #write, #read)
        src = open(APP, encoding="utf-8").read()
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


def viewer_scores(e, lag):
    """[bit, intended, correct] per grain from the viewer's own scoring code (app.js scoreGrains, run in node)."""
    src = open(APP, encoding="utf-8").read()
    a = src.index("function scoreGrains(e, lag)")
    js = src[a:src.index("// ---", a)] + f"""
        const e = {json.dumps(dict(hamr=e["hamr"], grains=dict(x=e["grains"]["x"]), analysis=e["analysis"]))};
        console.log(JSON.stringify(scoreGrains(e, {json.dumps(lag)}).map(a => [a.bit, a.intended, a.correct])));"""
    return json.loads(subprocess.run(["node", "-e", js], capture_output=True, text=True, check=True).stdout)


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
        self.assertIn("every material counted as magnetic", " ".join(e["notes"]))   # said on the page, not only in the log

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
        lead, m = read_grain_file(d, 9)
        self.assertEqual(len(m), len(read_grain_file(os.path.join(SCEN, "write-only", "write"), 9)[1]) - 1)

    def test_dipole_cells_default_to_vampire_size(self):
        # cells:macro-cell-size defaults to 10 A in VAMPIRE; the read was run with 1 nm cells written out
        edit = lambda t: re.sub(r"cells:macro-cell-size.*\n", "", t)
        d = copy_run(os.path.join(SCEN, "chain", "read"), os.path.join(self.tmp, "r"), edit)
        e = dict(size_nm=[12.0, 6.0, 1.2])
        os.makedirs(os.path.join(self.tmp, "out"))
        reduce.reduce_dipole(d, os.path.join(self.tmp, "out"), e)
        self.assertEqual(e["dipole"], dict(nx=13, ny=7, nz=2, cell_nm=1.0))


class EmptyGrains(unittest.TestCase):
    """A grain with no magnetic atoms: VAMPIRE writes 0 for it in grain.txt, which must not read as a melted grain."""

    def reduce_with_empty_grain(self, tmp, g, keep_coordinates):
        d = copy_run(os.path.join(SCEN, "write-only", "write"), os.path.join(tmp, "write"))
        lines = open(os.path.join(d, "grain.txt")).read().splitlines()
        rows = [l.split() if l.strip() and not l.startswith("#") else None for l in lines]
        for r in rows:
            if r:
                r[1 + 4 * g:5 + 4 * g] = ["0"] * 4          # one lead column (time), then mx my mz |m| per grain
        open(os.path.join(d, "grain.txt"), "w").write("\n".join(" ".join(r) if r else l for r, l in zip(rows, lines)) + "\n")
        if keep_coordinates:                                # the grain's atoms get an ID of no listed grain
            path = os.path.join(d, "atoms-coords.data")
            n = int(np.fromfile(path, np.uint64, 1)[0])
            ids = np.memmap(path, np.int32, "r+", offset=8 + 32 * n, shape=(n,))
            ids[ids == g] = 99
            ids.flush()
        else:                                               # no coordinate files: |m| = 0 throughout marks it
            for f in os.listdir(d):
                if f.startswith(("atoms-coords", "spins-")):
                    os.remove(os.path.join(d, f))
        reduce.reduce([d], os.path.join(tmp, "b"), workers=1)
        return json.load(open(os.path.join(tmp, "b", "manifest.json")))["runs"][0]

    def test_empty_grain_is_left_out(self):
        with tempfile.TemporaryDirectory() as tmp:
            reduce.reduce([os.path.join(SCEN, "write-only")], os.path.join(tmp, "ref"), workers=1)
            ref = json.load(open(os.path.join(tmp, "ref", "manifest.json")))["runs"][0]
        a = ref["analysis"]["grains"]
        g = 1
        self.assertEqual((a[g]["category"], a[g]["correct"]), ("full-field", True))   # melted and scored before
        for keep in (True, False):
            with self.subTest(coordinates=keep), tempfile.TemporaryDirectory() as tmp:
                e = self.reduce_with_empty_grain(tmp, g, keep)
                if keep:
                    self.assertEqual(e["grains"]["magnetic_atoms"][g], 0)
                b = e["analysis"]["grains"]
                self.assertEqual((b[g]["category"], b[g]["correct"]), ("empty", False))
                melted = lambda an: [r for r in an if r["category"] not in ("unheated", "empty") and r["intended"]]
                self.assertEqual(len(melted(b)), len(melted(a)) - 1)
                self.assertEqual([r["correct"] for k, r in enumerate(b) if k != g], [r["correct"] for k, r in enumerate(a) if k != g])
                if shutil.which("node"):                    # the viewer's scoring agrees
                    self.assertEqual(viewer_scores(e, e["freeze_lag_nm"]), [[r["bit"], r["intended"], r["correct"]] for r in b])


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
        html = open(os.path.join(self.tmp, "p32.html"), encoding="utf-8").read()
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
        src = open(APP, encoding="utf-8").read()
        code = src[src.index("const SNAP = "):src.index("// one snapshot array")]
        js = code + """const raw = { i8: [0, 127, -127, -128], i16: [0, 32767, -32767, -32768], f32: [0, 0.25, -1, NaN] };
            const out = Object.entries(SNAP).map(([k, [T, sc, E]]) => { const a = new T(raw[k]); a.S = sc; a.E = E;
              return [0, 1, 2, 3].map(i => snapVal(a, i)); });
            console.log(JSON.stringify(out.map(r => r.map(v => v !== v ? null : v))));"""
        res = json.loads(subprocess.run(["node", "-e", js], capture_output=True, text=True, check=True).stdout)
        self.assertEqual(res, [[0, 1, -1, None], [0, 1, -1, None], [0, 0.25, -1, None]])


class HeadCheck(unittest.TestCase):
    def test_a_wrong_model_is_caught(self):
        run = os.path.join(SCEN, "write-only", "write")
        cols, out = vp.read_output(run)
        h = vp.Hamr(run)
        self.assertEqual({k: v["bad"] for k, v in h.check([run], cols, out).items()}, {"B_applied": 0, "positions": 0})
        h.Hpeak *= 1.001
        h.step_nm *= 1.001
        bad = h.check([run], cols, out)
        self.assertGreater(bad["B_applied"]["bad"], 0)
        self.assertGreater(bad["positions"]["bad"], 0)

    def test_nothing_to_check_against(self):
        h = vp.Hamr(os.path.join(SCEN, "write-only", "write"))
        self.assertEqual(h.check([], ["Real_time"], np.zeros((0, 1))), {})


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


class Scale(unittest.TestCase):
    """Generated run of known content (tests/make_large.py) and what makes the reducer work at any size."""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, HERE)
        import make_large
        cls.ml = make_large
        cls.tmp = tempfile.mkdtemp()
        cls.rundir = os.path.join(cls.tmp, "run")
        with contextlib.redirect_stdout(io.StringIO()):
            make_large.make(cls.rundir, 100, 30, 2e5, 60, 400, 3)
        cls.h = vp.Hamr(cls.rundir)
        cls.nF, cls.film = 0, make_large.Film(100, 30, 2e5, 60)
        cls.x, cls.y, cls.mags = cls.film.centres()
        cls.state = make_large.State(cls.h, cls.x, cls.y)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def bundle(self, name, **kw):
        out = os.path.join(self.tmp, name)
        reduce.reduce([self.rundir], out, workers=1, **kw)
        return out, json.load(open(os.path.join(out, "manifest.json")))["runs"][0]

    def test_grain_series_and_overview_hold_the_known_answer(self):
        old = reduce.GRAIN_WHOLE_MAX
        reduce.GRAIN_WHOLE_MAX = 0                        # as a long series
        try:
            out, e = self.bundle("long")
        finally:
            reduce.GRAIN_WHOLE_MAX = old
        gs = e["grain_series"]
        nG, nF, bins = len(e["grains"]["id"]), gs["frames"], gs["overview"]["bins"]
        steps = np.rint(np.fromfile(os.path.join(out, "run", "grain_lead.f64")).reshape(nF, 2)[:, 0] / self.h.dt).astype(int)
        g = np.fromfile(os.path.join(out, "run", "grains.i16"), np.int16).reshape(nF, nG, 4) / 32767
        for f in range(0, nF, 37):
            d, m = self.state.at(steps[f])
            np.testing.assert_allclose(g[f, :, 2], d, atol=1e-4)
            np.testing.assert_allclose(g[f, :, 3], m, atol=1e-4)
        ov = np.fromfile(os.path.join(out, "run", "grains_overview.i16"), np.int16).reshape(bins, nG, 4) / 32767
        mz = g[:, :, 2] * g[:, :, 3]
        for b in range(bins):
            part, ml = mz[b * nF // bins:(b + 1) * nF // bins], g[b * nF // bins:(b + 1) * nF // bins, :, 3]
            np.testing.assert_allclose(ov[b, :, 0], part.min(0), atol=1 / 32767)
            np.testing.assert_allclose(ov[b, :, 1], part.max(0), atol=1 / 32767)
            np.testing.assert_allclose(ov[b, :, 2], ml.min(0), atol=1 / 32767)
            np.testing.assert_allclose(ov[b, :, 3], ml.max(0), atol=1 / 32767)
        # export holds the thinned series whole, so without the overview
        html = os.path.join(self.tmp, "long.html")
        with contextlib.redirect_stdout(io.StringIO()):
            export.export(out, html, grain_step=20, volume=False, limit_mb=1000)
        self.assertNotIn('"overview"', re.search(r'manifest.json":"([^"]*)"', open(html, encoding="utf-8").read()).group(0))

    def test_snapshots_hold_the_known_answer(self):
        out, e = self.bundle("snaps")
        G = e["snapshots"]["grid"]
        self.assertEqual((G["nx"], G["ny"]), (self.film.nx, self.film.ny))
        fr = e["snapshots"]["frames"][-1]
        d, m = self.state.at(round(fr["time"] / self.h.dt))
        q = np.frombuffer(gzip.open(os.path.join(out, "run", "top", f"{fr['index']:08d}.i8.gz")).read(), np.int8).reshape(3, G["ny"], G["nx"])
        gm = np.fromfile(os.path.join(out, "run", "grain_map.i16"), np.int16).reshape(G["ny"], G["nx"])
        ok = gm >= 0
        self.assertEqual(int(ok.sum()), int((self.mags / self.film.nz).sum()))
        np.testing.assert_allclose(q[2][ok] / 127, (d * m)[gm[ok]], atol=0.5 / 127 + 1e-9)
        self.assertTrue((q[2][~ok] == -128).all())
        self.assertLess(e["checks"]["grain_txt_one_step_earlier"]["max_diff"], 2e-3)   # spins against grain.txt
        self.assertGreaterEqual(e["checks"]["grain_txt_one_step_earlier"]["n"], 2)

    def test_pieces_change_nothing(self):
        out_a, a = self.bundle("whole")
        old = vp.coord_chunks.__defaults__
        vp.coord_chunks.__defaults__ = (777,)
        try:
            out_b, b = self.bundle("pieces")
        finally:
            vp.coord_chunks.__defaults__ = old
        self.assertEqual(open(os.path.join(out_a, "run", "grain_map.i16"), "rb").read(), open(os.path.join(out_b, "run", "grain_map.i16"), "rb").read())
        self.assertEqual(open(os.path.join(out_a, "run", "top_count.u8"), "rb").read(), open(os.path.join(out_b, "run", "top_count.u8"), "rb").read())
        for k in ("x", "y", "magnetic_atoms"):
            self.assertEqual(a["grains"][k], b["grains"][k])
        f = vp.snapshots(self.rundir)[1][2][0]
        np.testing.assert_array_equal(np.concatenate(list(vp.spin_chunks(f, 777))), vp.read_spins(f))

    def test_interrupted_reduce_carries_on(self):
        def tree(out):
            return {os.path.relpath(os.path.join(d, f), out): gzip.decompress(open(os.path.join(d, f), "rb").read()) if f.endswith(".gz") else open(os.path.join(d, f), "rb").read()
                    for d, _, fs in os.walk(out) for f in fs if f != "manifest.json"}
        whole, _ = self.bundle("uninterrupted")
        # during the snapshots: the first save of progress is the last thing the job does
        saver, every = reduce._save_progress, reduce.SAVE_EVERY

        def killed(*a):
            saver(*a)
            raise KeyboardInterrupt
        reduce._save_progress, reduce.SAVE_EVERY = killed, -1
        out = os.path.join(self.tmp, "killed")
        try:
            with self.assertRaises(KeyboardInterrupt):
                reduce.reduce([self.rundir], out, workers=1)
        finally:
            reduce._save_progress, reduce.SAVE_EVERY = saver, every
        self.assertTrue(os.path.exists(os.path.join(out, reduce.WORK, "run", "progress.npz")))
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            reduce.reduce([self.rundir], out, workers=1)
        self.assertIn("carrying on", buf.getvalue())
        self.assertIn("already parsed", buf.getvalue())
        self.assertFalse(os.path.exists(os.path.join(out, reduce.WORK)))
        self.assertEqual(tree(whole), tree(out))
        ma, mb = (json.load(open(os.path.join(p, "manifest.json")))["runs"] for p in (whole, out))
        self.assertEqual(ma, mb)

    def test_interrupted_between_runs_and_other_options_start_again(self):
        out = os.path.join(self.tmp, "between")
        scen = os.path.join(SCEN, "chain")
        atoms = reduce.reduce_atoms

        def stop_at_write(dirs, *a, **k):
            if dirs[0].endswith("write"):
                raise KeyboardInterrupt
            return atoms(dirs, *a, **k)
        reduce.reduce_atoms = stop_at_write
        try:
            with self.assertRaises(KeyboardInterrupt):
                reduce.reduce([scen], out, workers=1)
        finally:
            reduce.reduce_atoms = atoms
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            reduce.reduce([scen], out, workers=1)
        self.assertIn("1 of 3 runs were complete", buf.getvalue())
        self.assertEqual([e["id"] for e in json.load(open(os.path.join(out, "manifest.json")))["runs"]], ["relax", "write", "read"])
        self.assertEqual(json.load(open(os.path.join(out, "manifest.json")))["runs"][2]["ref"], "write")
        # nothing to carry on after a complete reduce, and a different request starts again
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            reduce.reduce([scen], out, every=2, workers=1)
        self.assertNotIn("carrying on", buf.getvalue())

    def test_changed_options_leave_nothing_of_the_unfinished_reduce(self):
        out = os.path.join(self.tmp, "changed")
        atoms = reduce.reduce_atoms

        def stop(*a, **k):
            raise KeyboardInterrupt
        reduce.reduce_atoms = stop
        try:
            with self.assertRaises(KeyboardInterrupt):
                reduce.reduce([self.rundir], out, workers=1, opts=dict(grain_precision=32))
        finally:
            reduce.reduce_atoms = atoms
        self.assertTrue(os.path.exists(os.path.join(out, "run", "grains.f32")))
        reduce.reduce([self.rundir], out, workers=1)      # different options: starts again
        self.assertFalse(os.path.exists(os.path.join(out, "run", "grains.f32")))
        self.assertFalse(os.path.exists(os.path.join(out, reduce.WORK)))

    def test_accumulator_matches_the_whole_series(self):
        rng = np.random.default_rng(7)
        m = rng.uniform(-1, 1, (200, 40, 4)).astype(np.float32)
        m[:, :, 3] = np.abs(m[:, :, 3])
        m[-1, 3, 2] = 0                                   # final m_z exactly 0
        m[150:, 5, 2] = np.sign(m[150, 5, 2]) or 1        # settled long before the end
        m[:, 6, 2], m[:, 6, 3] = 1, .9                    # never below the threshold
        acc = reduce.GrainAccumulator(40, reduce.DECISION)
        for f in m:
            acc.add(f)
        final, minM, anyb, dec, m_dec = acc.outcome()
        # the formulas on the whole array
        mz = m[:, :, 2] * m[:, :, 3]
        want_final = mz[-1]
        below = mz * np.sign(want_final) < reduce.DECISION
        want_dec = np.minimum(np.where(below.any(0), len(m) - np.argmax(below[::-1], axis=0), 0), len(m) - 1)
        np.testing.assert_array_equal(final, want_final)
        np.testing.assert_array_equal(minM, m[:, :, 3].min(0))
        np.testing.assert_array_equal(anyb, below.any(0))
        np.testing.assert_array_equal(dec, want_dec)
        got = below.any(0)
        np.testing.assert_array_equal(m_dec[got], m[want_dec, np.arange(40), 3][got])

    def test_outlines_match_the_plain_search(self):
        coords = np.array(vp.read_grain_coordinates(RUN))
        shapes = vp.read_grain_shapes(RUN)
        res = vp.grain_outlines(coords, shapes)
        cen = np.array([s.mean(axis=0) for s in shapes])
        for (x, y), o in zip(coords[:, 1:3], res):
            inside = [s for s in shapes if vp._contains(s, x, y)]
            want = inside[0] if len(inside) == 1 else shapes[int(np.argmin(np.hypot(*(cen - (x, y)).T)))]
            np.testing.assert_array_equal(o, want)

    def test_worker_count_follows_the_memory_budget(self):
        self.assertEqual(reduce.worker_count(100, 1e9, 3, 1e9)[0], 3)              # -j wins
        self.assertEqual(reduce.worker_count(100, 1e9, None, 2.5e9)[0], min(2, reduce._cpus() - 1 or 1))
        self.assertEqual(reduce.worker_count(100, 1e12, None, 1e9)[0], 1)
        self.assertEqual(reduce.worker_count(1, 1e6, None, 1e12)[0], 1)            # never more than jobs

    def test_server_answers_byte_ranges(self):
        out, _ = self.bundle("ranges")
        srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(serve.Handler, bundle=out))
        threading.Thread(target=srv.serve_forever, daemon=True).start()
        try:
            url = f"http://127.0.0.1:{srv.server_address[1]}/bundle/run/grains.i16"
            raw = open(os.path.join(out, "run", "grains.i16"), "rb").read()
            r = urllib.request.urlopen(urllib.request.Request(url, headers={"Range": "bytes=800-815"}))
            self.assertEqual((r.status, r.read(), r.headers["Content-Range"]), (206, raw[800:816], f"bytes 800-815/{len(raw)}"))
            self.assertEqual(urllib.request.urlopen(url).read(), raw)                 # no Range: whole file
            r = urllib.request.urlopen(urllib.request.Request(url, headers={"Range": f"bytes={len(raw) - 2}-{len(raw) + 9}"}))
            self.assertEqual(r.read(), raw[-2:])                                     # past the end: what is there
            with self.assertRaises(urllib.error.HTTPError) as ctx:
                urllib.request.urlopen(urllib.request.Request(url, headers={"Range": f"bytes={len(raw) + 5}-{len(raw) + 9}"}))
            self.assertEqual(ctx.exception.code, 416)
            with self.assertRaises(urllib.error.HTTPError):                          # never outside the bundle
                urllib.request.urlopen(urllib.request.Request(url.replace("run/grains.i16", "../../input"), headers={"Range": "bytes=0-5"}))
        finally:
            srv.shutdown()


if __name__ == "__main__":
    unittest.main()

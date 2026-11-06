"""Every scenario's bundle in a real browser: no page errors, every lane drawn. Skipped without node and Playwright."""
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import unittest
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCEN = os.path.join(HERE, "scenarios")
sys.path.insert(0, ROOT)
from hamrvis import reduce  # noqa: E402


def have_playwright():
    return shutil.which("node") and subprocess.run(["node", "-e", "require('playwright')"], cwd=HERE,
                                                   capture_output=True).returncode == 0


@unittest.skipUnless(have_playwright() or os.environ.get("HAMRVIS_REQUIRE_VIEWER_TEST"), "needs node and Playwright")
class ViewerSmoke(unittest.TestCase):
    def test_every_scenario_draws(self):
        for name in sorted(os.listdir(SCEN)):
            with self.subTest(name), tempfile.TemporaryDirectory() as tmp:
                out = os.path.join(tmp, "bundle")
                reduce.reduce([os.path.join(SCEN, name)], out, workers=1)
                runs = [r["id"] for r in json.load(open(os.path.join(out, "manifest.json")))["runs"]]
                with socket.socket() as s:
                    s.bind(("127.0.0.1", 0))
                    port = s.getsockname()[1]
                view = subprocess.Popen([sys.executable, "-m", "hamrvis", "view", out, "-p", str(port), "--no-browser"],
                                        cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                try:
                    url = f"http://127.0.0.1:{port}/"
                    for _ in range(100):
                        try:
                            urllib.request.urlopen(url)
                            break
                        except OSError:
                            time.sleep(0.1)
                    r = subprocess.run(["node", "viewer_smoke.js", url] + runs, cwd=HERE, capture_output=True, text=True)
                    self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
                finally:
                    view.kill()
                    view.wait()


if __name__ == "__main__":
    unittest.main()

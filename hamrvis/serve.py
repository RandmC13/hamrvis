"""Serve the viewer and one bundle on localhost: / is hamrvis/viewer/, /bundle/ is the bundle directory."""
import functools
import http.server
import io
import os
import re
import threading
import webbrowser

VIEWER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "viewer")


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, bundle, **k):
        self.bundle = bundle
        super().__init__(*a, directory=VIEWER, **k)

    def translate_path(self, path):
        path = path.split("?", 1)[0].split("#", 1)[0]
        if path.startswith("/bundle/"):
            rel = os.path.normpath(path[len("/bundle/"):]).lstrip(os.sep)
            full = os.path.abspath(os.path.join(self.bundle, rel))    # abspath: "C:/..." would replace the bundle on Windows
            if rel.startswith("..") or not full.startswith(self.bundle + os.sep):
                return os.path.join(self.bundle, "__outside__")       # 404, never leave the bundle
            return full
        return super().translate_path(path)

    def send_head(self):
        """A plain file, or the byte range the viewer asked for: long grain series are read in blocks."""
        m = re.fullmatch(r"bytes=(\d+)-(\d+)", self.headers.get("Range", ""))
        path = self.translate_path(self.path)
        if not m or not os.path.isfile(path):
            return super().send_head()
        size = os.path.getsize(path)
        a, b = int(m.group(1)), min(int(m.group(2)), size - 1)
        if a > b:
            self.send_error(416)
            return None
        with open(path, "rb") as f:
            f.seek(a)
            data = f.read(b - a + 1)
        self.send_response(206)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Range", f"bytes {a}-{b}/{size}")
        self.send_header("Content-Length", str(b - a + 1))
        self.end_headers()
        return io.BytesIO(data)

    def end_headers(self):
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Cache-Control", "no-store")                 # a rebuilt bundle shows on reload
        super().end_headers()

    def log_message(self, fmt, *args):
        if args and str(args[1])[0] in "45":                           # only report failed requests
            super().log_message(fmt, *args)


def serve(bundle, port=8737, browser=True):
    bundle = os.path.abspath(bundle)
    if not os.path.exists(os.path.join(bundle, "manifest.json")):
        raise SystemExit(f"{bundle} has no manifest.json: make it with 'python3 -m hamrvis reduce'")
    try:
        srv = http.server.ThreadingHTTPServer(("127.0.0.1", port), functools.partial(Handler, bundle=bundle))
    except OSError as err:
        raise SystemExit(f"cannot serve on port {port} ({err.strerror}): another viewer or program is using it; "
                         f"stop it or choose another port with -p") from None
    url = f"http://127.0.0.1:{port}/"
    print(f"Viewing {bundle}\nOpen {url}  (Ctrl+C to stop)", flush=True)
    if browser:
        threading.Timer(0.5, webbrowser.open, (url,)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass

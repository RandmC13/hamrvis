"""Serve the viewer and one bundle on localhost: / is viewer/, /bundle/ is the bundle directory."""
import functools
import http.server
import os
import threading
import webbrowser

VIEWER = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "viewer")


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, bundle, **k):
        self.bundle = bundle
        super().__init__(*a, directory=VIEWER, **k)

    def translate_path(self, path):
        path = path.split("?", 1)[0].split("#", 1)[0]
        if path.startswith("/bundle/"):
            rel = os.path.normpath(path[len("/bundle/"):]).lstrip(os.sep)
            if rel.startswith(".."):
                return os.path.join(self.bundle, "__outside__")       # 404, never leave the bundle
            return os.path.join(self.bundle, rel)
        return super().translate_path(path)

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")                 # a rebuilt bundle shows on reload
        super().end_headers()

    def log_message(self, fmt, *args):
        if args and str(args[1])[0] in "45":                           # only report failed requests
            super().log_message(fmt, *args)


def serve(bundle, port=8737, browser=True):
    bundle = os.path.abspath(bundle)
    if not os.path.exists(os.path.join(bundle, "manifest.json")):
        raise SystemExit(f"{bundle} has no manifest.json: make it with 'python3 -m hamrvis reduce'")
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", port), functools.partial(Handler, bundle=bundle))
    url = f"http://127.0.0.1:{port}/"
    print(f"Viewing {bundle}\nOpen {url}  (Ctrl+C to stop)", flush=True)
    if browser:
        threading.Timer(0.5, webbrowser.open, (url,)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass

"""Minimal JSON HTTP service (standard library only), intended as a backend for a
browser interface.

POST /calculate  body: {"canrisk": "<CanRisk file text>"} or {"pedigree": {...}},
                 optional "settings": {"model": "BC"|"OC", "population": "UK",
                 "mutation_frequencies": "UK", "sensitivities": {...}}
                 -> list of results (one per pedigree)
GET  /populations -> list of available incidence populations
"""
from dataclasses import fields
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json

from .api import Settings, calculate
from .incidence import available_populations
from .pedigree import Pedigree, PedigreeError, read_canrisk

MAX_BODY = 2_000_000


def handle_request(payload):
    """Process a decoded JSON request; returns a JSON-serialisable object."""
    allowed = {f.name for f in fields(Settings)}
    opts = {k: v for k, v in (payload.get("settings") or {}).items() if k in allowed}
    settings = Settings(**opts)
    if "canrisk" in payload:
        peds = read_canrisk(payload["canrisk"])
    elif "pedigree" in payload:
        peds = [Pedigree.from_dict(payload["pedigree"])]
    else:
        raise PedigreeError("request must contain 'canrisk' or 'pedigree'")
    return [calculate(p, settings).to_dict() for p in peds]


class _Handler(BaseHTTPRequestHandler):
    cors = False

    def _send(self, code, obj):
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        if self.cors:
            self.send_header("Access-Control-Allow-Origin", "*")
            self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()
        self.wfile.write(body)

    def do_OPTIONS(self):
        self._send(204, {})

    def do_GET(self):
        if self.path == "/populations":
            self._send(200, available_populations())
        else:
            self._send(404, {"error": "not found"})

    def do_POST(self):
        if self.path != "/calculate":
            self._send(404, {"error": "not found"})
            return
        n = int(self.headers.get("Content-Length", 0))
        if n > MAX_BODY:
            self._send(413, {"error": "request too large"})
            return
        try:
            payload = json.loads(self.rfile.read(n) or b"{}")
            self._send(200, handle_request(payload))
        except (PedigreeError, ValueError, KeyError, TypeError) as e:
            self._send(400, {"error": str(e)})


def serve(host="127.0.0.1", port=8080, cors=False):
    _Handler.cors = cors
    httpd = ThreadingHTTPServer((host, port), _Handler)
    print(f"openboadicea service on http://{host}:{port}")
    httpd.serve_forever()

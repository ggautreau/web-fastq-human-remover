#!/usr/bin/env python3
"""Plain static server — no special headers needed, same as GitHub Pages."""
import http.server, socketserver, os
DIR = os.path.dirname(os.path.abspath(__file__))
class H(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw): super().__init__(*a, directory=DIR, **kw)
    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()
    def log_message(self, *a): pass
socketserver.TCPServer.allow_reuse_address = True
with socketserver.TCPServer(("127.0.0.1", 8767), H) as httpd:
    print("http://127.0.0.1:8767/")
    httpd.serve_forever()

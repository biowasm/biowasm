#!/usr/bin/python3
"""
Simple CORS enabled, but open HTTP server for test purposes.
Maps `.wasm` extension files to MIME type `application/wasm`.
Opens on TCP port 80.
"""
from http.server import SimpleHTTPRequestHandler
import socketserver

class CORSRequestHandler (SimpleHTTPRequestHandler):
    def end_headers (self):
        self.send_header('Access-Control-Allow-Origin', '*')
        SimpleHTTPRequestHandler.end_headers(self)

handler = CORSRequestHandler
handler.extensions_map['.wasm'] = 'application/wasm'

httpd = socketserver.TCPServer(('', 80), handler)
httpd.serve_forever()

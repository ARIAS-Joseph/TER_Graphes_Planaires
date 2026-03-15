#!/usr/bin/env python3
import json, os, platform, subprocess, threading, webbrowser
from http.server import SimpleHTTPRequestHandler, HTTPServer

PORT = 8080
GRAPH_FILE = "graph.txt"

system = platform.system()
if system == "Windows":
    EXE = os.path.join(os.path.dirname(__file__), "graph_tool.exe")
elif system == "Darwin":
    EXE = os.path.join(os.path.dirname(__file__), "graph_tool_mac")
else:
    EXE = os.path.join(os.path.dirname(__file__), "graph_tool_linux")


class Handler(SimpleHTTPRequestHandler):

    def do_OPTIONS(self):
        self.send_response(200)
        self._cors()
        self.end_headers()

    def do_POST(self):
        length = int(self.headers.get("Content-Length", 0))
        body = json.loads(self.rfile.read(length))

        action = body.get("action", "")
        graph  = body.get("graph", None)   # None pour generate_graph
        args   = body.get("args", [])

        # Écrire le graphe dans le fichier SEULEMENT si fourni (pas pour generate_graph)
        if graph is not None:
            with open(GRAPH_FILE, "w", encoding="utf-8") as f:
                f.write(graph)

        # Construire la commande
        cmd = [EXE, action, GRAPH_FILE] + [str(a) for a in args]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            stderr_out = result.stderr.strip()
            if result.stdout.strip():
                print(f"[C stdout]\n{result.stdout}")
            if stderr_out:
                print(f"[C stderr] {stderr_out}")
            if result.returncode != 0:
                self._json_error(f"Programme C a retourne code {result.returncode}: {stderr_out}")
                return
        except FileNotFoundError:
            self._json_error(f"Executable introuvable : {EXE}")
            return
        except subprocess.TimeoutExpired:
            self._json_error("Timeout : le programme C a pris trop de temps")
            return

        # Lire le graphe mis à jour
        try:
            with open(GRAPH_FILE, encoding="utf-8") as f:
                updated_graph = f.read()
        except FileNotFoundError:
            self._json_error("Le programme C n'a pas produit de fichier de sortie")
            return

        if not updated_graph.strip():
            self._json_error("Le programme C a produit un fichier vide")
            return

        self.send_response(200)
        self._cors()
        self.send_header("Content-Type", "application/json")
        self.end_headers()
        self.wfile.write(json.dumps({
            "graph": updated_graph,
            "stdout": result.stdout,
            "ok": True
        }).encode())

    def _cors(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Allow-Methods", "POST, GET, OPTIONS")

    def _json_error(self, msg):
        self.send_response(500)
        self._cors()
        self.send_header("Content-Type", "application/json")
        self.end_headers()
        self.wfile.write(json.dumps({"error": msg}).encode())
        print(f"[ERREUR] {msg}")

    def log_message(self, fmt, *args):
        if "POST" in (fmt % args):
            print(f"[API] {fmt % args}")


if __name__ == "__main__":
    if not os.path.exists(EXE):
        print(f"⚠  Executable '{EXE}' introuvable.")
        print("   Compilez graph_tool et placez l'exe ici.\n")

    server = HTTPServer(("localhost", PORT), Handler)
    url = f"http://localhost:{PORT}/visualizer2.html"
    print(f"✓  Serveur demarre → {url}")
    print("   Ctrl+C pour arreter\n")

    threading.Timer(1.0, lambda: webbrowser.open(url)).start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nServeur arrete.")
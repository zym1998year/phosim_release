#!/usr/bin/env python
##
## @package phosim
## @file phosim_gui2.py
## @brief compact, reliable Tk GUI front-end for phosim.py
##
## A small alternative to the legacy PAGE-generated phosim_gui, focused on
## running reliably on Windows/Linux/macOS with the Python stdlib only:
##   - resizable layout, file/dir browse dialogs, examples/ dropdowns
##   - builds an explicit argv (never shell=True) and shows it before running
##   - streams phosim output live into the window and reports the exit code
##   - Stop button terminates the whole process tree
##   - remembers the last-used settings in ~/.phosim_gui2.json
##
## Usage:  python tools/phosim_gui2.py        (phosim_gui2.bat on Windows)
##
## @warning This code is not fully validated
## and not ready for full release.  Please
## treat results with caution.
##

import glob
import json
import os
import queue
import shlex
import signal
import subprocess
import sys
import threading
import time
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinter.scrolledtext import ScrolledText

PHOSIM_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".phosim_gui2.json")


# ---------------------------------------------------------------------------
# repo introspection helpers
# ---------------------------------------------------------------------------

def list_instruments():
    """Instrument = subdirectory of data/ containing focalplanelayout.txt."""
    out = []
    dataDir = os.path.join(PHOSIM_ROOT, "data")
    try:
        for name in sorted(os.listdir(dataDir)):
            if os.path.isfile(os.path.join(dataDir, name, "focalplanelayout.txt")):
                out.append(name)
    except OSError:
        pass
    return out or ["generic"]


def list_examples():
    """Plain files directly under examples/ (instance catalogs and command files)."""
    out = []
    exDir = os.path.join(PHOSIM_ROOT, "examples")
    try:
        for name in sorted(os.listdir(exDir)):
            if os.path.isfile(os.path.join(exDir, name)):
                out.append("examples/" + name)
    except OSError:
        pass
    return out


def split_extra_args(text):
    """Split the free-form extra-args string into argv tokens.

    A POSIX shlex with the escape character disabled gives the same behavior on
    every platform: double/single quotes group and are stripped, while
    backslashes (Windows paths) are kept literal.
    """
    lex = shlex.shlex(text, posix=True)
    lex.whitespace_split = True
    lex.escape = ""
    return list(lex)


def child_environment():
    """Environment for the phosim child process.

    On Windows the raytrace/e2adc executables need the cfitsio/fftw/zlib DLLs;
    when running inside a conda env they live in <prefix>\\Library\\bin, which is
    not always on PATH. Prepend it so a plain double-click launch still works.
    """
    env = os.environ.copy()
    if os.name == "nt":
        dllDir = os.path.join(sys.prefix, "Library", "bin")
        if os.path.isdir(dllDir):
            env["PATH"] = dllDir + os.pathsep + env.get("PATH", "")
    return env


# ---------------------------------------------------------------------------
# subprocess runner (worker thread; all UI updates flow through a queue)
# ---------------------------------------------------------------------------

class PhosimRunner:
    """Runs one phosim.py invocation and streams its output line by line.

    Thread-safety contract: the worker thread only ever puts items on
    self.events; the Tk main loop drains the queue via root.after(). Events:
    ("line", text) | ("done", returncode).
    """

    def __init__(self):
        self.events = queue.Queue()
        self.proc = None
        self._thread = None

    @property
    def running(self):
        return self.proc is not None and self.proc.poll() is None

    def start(self, argv, cwd):
        # Windows: CREATE_NO_WINDOW avoids a console popup and taskkill /T can
        # still reach the tree. POSIX: own session so stop() can killpg().
        creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        self.proc = subprocess.Popen(
            argv,
            cwd=cwd,
            env=child_environment(),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
            creationflags=creationflags,
            start_new_session=(os.name != "nt"),
        )
        self._thread = threading.Thread(target=self._pump, daemon=True)
        self._thread.start()

    def _pump(self):
        try:
            for line in self.proc.stdout:
                self.events.put(("line", line.rstrip("\n")))
        except (ValueError, OSError):
            pass  # pipe closed by stop() while we were blocked reading
        except Exception as e:  # reader must never take the app down
            self.events.put(("line", "[gui] output reader error: %r" % (e,)))
        rc = self.proc.wait()
        self.events.put(("done", rc))

    def stop(self):
        """Terminate phosim.py and its raytrace/e2adc children."""
        if not self.running:
            return
        try:
            if os.name == "nt":
                subprocess.call(
                    ["taskkill", "/T", "/F", "/PID", str(self.proc.pid)],
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                    creationflags=subprocess.CREATE_NO_WINDOW,
                )
            else:
                os.killpg(os.getpgid(self.proc.pid), signal.SIGTERM)
        except Exception as e:
            self.events.put(("line", "[gui] stop failed: %r" % (e,)))
        # The tree kill can race with a child being spawned; a survivor holding
        # the pipe's write end would leave _pump blocked forever. Closing our
        # read end unblocks it so the "done" event is always delivered.
        try:
            if self.proc.stdout:
                self.proc.stdout.close()
        except OSError:
            pass


# ---------------------------------------------------------------------------
# main application
# ---------------------------------------------------------------------------

class App:
    POLL_MS = 100

    def __init__(self, root):
        self.root = root
        self.runner = PhosimRunner()
        self._t0 = None
        self._stop_requested = False
        self._done_handled = True   # no run yet
        self._dead_since = None
        root.title("PhoSim launcher (phosim_gui2)")
        root.minsize(720, 560)
        root.protocol("WM_DELETE_WINDOW", self.on_close)

        self._build_vars()
        self._build_ui()
        self._load_settings()
        self._refresh_preview()
        self._poll()

    # ---- state ----------------------------------------------------------

    def _build_vars(self):
        v = tk.StringVar
        self.catalog = v(value="examples/star")
        self.command = v(value="examples/nobackground")
        self.instrument = v(value="generic")
        self.sensor = v(value="")
        self.processors = tk.IntVar(value=1)
        self.threads = tk.IntVar(value=1)
        self.run_e2adc = tk.BooleanVar(value=True)
        self.outputDir = v(value="output")
        self.workDir = v(value="work")
        self.extraArgs = v(value="")
        self.status = v(value="Ready.")
        for var in (self.catalog, self.command, self.instrument, self.sensor,
                    self.outputDir, self.workDir, self.extraArgs):
            var.trace_add("write", lambda *_: self._refresh_preview())
        for var in (self.processors, self.threads, self.run_e2adc):
            var.trace_add("write", lambda *_: self._refresh_preview())

    # ---- layout ---------------------------------------------------------

    def _build_ui(self):
        pad = dict(padx=6, pady=3)
        outer = ttk.Frame(self.root, padding=8)
        outer.grid(sticky="nsew")
        self.root.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)
        outer.columnconfigure(0, weight=1)
        # log pane row stretches
        outer.rowconfigure(4, weight=1)

        # -- inputs
        box = ttk.LabelFrame(outer, text="Simulation input", padding=6)
        box.grid(row=0, column=0, sticky="ew", **pad)
        box.columnconfigure(1, weight=1)
        examples = list_examples()

        ttk.Label(box, text="Instance catalog").grid(row=0, column=0, sticky="w")
        ttk.Combobox(box, textvariable=self.catalog, values=examples)\
            .grid(row=0, column=1, sticky="ew", **pad)
        ttk.Button(box, text="Browse...", command=lambda: self._pick_file(self.catalog))\
            .grid(row=0, column=2, **pad)

        ttk.Label(box, text="Command file (optional)").grid(row=1, column=0, sticky="w")
        ttk.Combobox(box, textvariable=self.command, values=[""] + examples)\
            .grid(row=1, column=1, sticky="ew", **pad)
        ttk.Button(box, text="Browse...", command=lambda: self._pick_file(self.command))\
            .grid(row=1, column=2, **pad)

        # -- options
        box = ttk.LabelFrame(outer, text="Options", padding=6)
        box.grid(row=1, column=0, sticky="ew", **pad)
        for c in (1, 3, 5):
            box.columnconfigure(c, weight=1)

        ttk.Label(box, text="Instrument (-i)").grid(row=0, column=0, sticky="w")
        ttk.Combobox(box, textvariable=self.instrument, values=list_instruments())\
            .grid(row=0, column=1, sticky="ew", **pad)
        ttk.Label(box, text="Sensor (-s)").grid(row=0, column=2, sticky="w")
        ttk.Entry(box, textvariable=self.sensor).grid(row=0, column=3, sticky="ew", **pad)
        ttk.Checkbutton(box, text="Run e2adc (-e)", variable=self.run_e2adc)\
            .grid(row=0, column=4, columnspan=2, sticky="w", **pad)

        ttk.Label(box, text="Processes (-p)").grid(row=1, column=0, sticky="w")
        ttk.Spinbox(box, from_=1, to=max(1, os.cpu_count() or 1),
                    textvariable=self.processors, width=6)\
            .grid(row=1, column=1, sticky="w", **pad)
        ttk.Label(box, text="Threads (-t)").grid(row=1, column=2, sticky="w")
        ttk.Spinbox(box, from_=1, to=max(1, os.cpu_count() or 1),
                    textvariable=self.threads, width=6)\
            .grid(row=1, column=3, sticky="w", **pad)
        ttk.Label(box, text="Extra args").grid(row=1, column=4, sticky="w")
        ttk.Entry(box, textvariable=self.extraArgs).grid(row=1, column=5, sticky="ew", **pad)

        # -- directories
        box = ttk.LabelFrame(outer, text="Directories (relative to the phosim root)", padding=6)
        box.grid(row=2, column=0, sticky="ew", **pad)
        box.columnconfigure(1, weight=1)
        box.columnconfigure(4, weight=1)
        ttk.Label(box, text="Output (-o)").grid(row=0, column=0, sticky="w")
        ttk.Entry(box, textvariable=self.outputDir).grid(row=0, column=1, sticky="ew", **pad)
        ttk.Button(box, text="Browse...", command=lambda: self._pick_dir(self.outputDir))\
            .grid(row=0, column=2, **pad)
        ttk.Label(box, text="Work (-w)").grid(row=0, column=3, sticky="w")
        ttk.Entry(box, textvariable=self.workDir).grid(row=0, column=4, sticky="ew", **pad)
        ttk.Button(box, text="Browse...", command=lambda: self._pick_dir(self.workDir))\
            .grid(row=0, column=5, **pad)

        # -- command preview + actions
        box = ttk.Frame(outer)
        box.grid(row=3, column=0, sticky="ew", **pad)
        box.columnconfigure(0, weight=1)
        self.preview = ttk.Entry(box, state="readonly")
        self.preview.grid(row=0, column=0, sticky="ew", **pad)
        self.runBtn = ttk.Button(box, text="Run PhoSim", command=self.on_run)
        self.runBtn.grid(row=0, column=1, **pad)
        self.stopBtn = ttk.Button(box, text="Stop", command=self.on_stop, state="disabled")
        self.stopBtn.grid(row=0, column=2, **pad)
        ttk.Button(box, text="Open output", command=self.on_open_output).grid(row=0, column=3, **pad)

        # -- log
        self.log = ScrolledText(outer, height=16, state="disabled", wrap="none",
                                font=("Consolas" if os.name == "nt" else "Courier", 9))
        self.log.grid(row=4, column=0, sticky="nsew", **pad)

        # -- status bar
        ttk.Label(outer, textvariable=self.status, anchor="w", relief="sunken")\
            .grid(row=5, column=0, sticky="ew")

    # ---- small helpers --------------------------------------------------

    def _pick_file(self, var):
        path = filedialog.askopenfilename(initialdir=PHOSIM_ROOT)
        if path:
            var.set(self._relativize(path))

    def _pick_dir(self, var):
        path = filedialog.askdirectory(initialdir=PHOSIM_ROOT)
        if path:
            var.set(self._relativize(path))

    @staticmethod
    def _relativize(path):
        """Prefer repo-relative paths so saved settings stay portable."""
        try:
            rel = os.path.relpath(path, PHOSIM_ROOT)
            if not rel.startswith(".."):
                return rel.replace("\\", "/")
        except ValueError:  # different drive on Windows
            pass
        return path

    MAX_LOG_LINES = 5000

    def _append_log(self, text):
        self.log.configure(state="normal")
        self.log.insert("end", text + "\n")
        # bound memory/redraw cost on very long runs
        lines = int(self.log.index("end-1c").split(".")[0])
        if lines > self.MAX_LOG_LINES:
            self.log.delete("1.0", "%d.0" % (lines - self.MAX_LOG_LINES + 1))
        self.log.see("end")
        self.log.configure(state="disabled")

    # ---- command construction -------------------------------------------

    def build_argv(self):
        argv = [sys.executable, os.path.join(PHOSIM_ROOT, "phosim.py"),
                self.catalog.get().strip()]
        cmd = self.command.get().strip()
        if cmd:
            argv += ["-c", cmd]
        inst = self.instrument.get().strip()
        if inst and inst != "generic":
            argv += ["-i", inst]
        sensor = self.sensor.get().strip()
        if sensor:
            argv += ["-s", sensor]
        argv += ["-p", str(max(1, self.processors.get())),
                 "-t", str(max(1, self.threads.get()))]
        if not self.run_e2adc.get():
            argv += ["-e", "0"]
        out = self.outputDir.get().strip()
        if out and out != "output":
            argv += ["-o", out]
        work = self.workDir.get().strip()
        if work and work != "work":
            argv += ["-w", work]
        extra = self.extraArgs.get().strip()
        if extra:
            argv += split_extra_args(extra)
        return argv

    def _refresh_preview(self):
        try:
            argv = self.build_argv()
        except tk.TclError:  # spinbox mid-edit
            return
        shown = argv[:]
        shown[0] = os.path.basename(shown[0])
        shown[1] = "phosim.py"
        text = " ".join(a if " " not in a else '"%s"' % a for a in shown)
        self.preview.configure(state="normal")
        self.preview.delete(0, "end")
        self.preview.insert(0, text)
        self.preview.configure(state="readonly")

    def _resolve(self, path):
        return path if os.path.isabs(path) else os.path.join(PHOSIM_ROOT, path)

    def _validate(self):
        try:
            self.processors.get()
            self.threads.get()
        except tk.TclError:
            return "Processes (-p) and threads (-t) must be integers."
        catalog = self.catalog.get().strip()
        if not catalog:
            return "Please choose an instance catalog."
        if not os.path.isfile(self._resolve(catalog)):
            return "Instance catalog not found: %s" % catalog
        cmd = self.command.get().strip()
        if cmd and not os.path.isfile(self._resolve(cmd)):
            return "Command file not found: %s" % cmd
        for label, var in (("output", self.outputDir), ("work", self.workDir)):
            d = var.get().strip()
            if d:
                try:
                    os.makedirs(self._resolve(d), exist_ok=True)
                except OSError as e:
                    return "Cannot create %s directory %s: %s" % (label, d, e)
        return None

    # ---- actions ---------------------------------------------------------

    def on_run(self):
        if self.runner.running:
            return
        err = self._validate()
        if err:
            messagebox.showerror("PhoSim", err, parent=self.root)
            return
        self._save_settings()
        argv = self.build_argv()
        self._append_log("$ " + " ".join(argv))
        self._stop_requested = False
        self._done_handled = False
        self._dead_since = None
        self._t0 = time.time()
        try:
            self.runner.start(argv, cwd=PHOSIM_ROOT)
        except OSError as e:
            self._append_log("[gui] failed to start: %s" % e)
            self.status.set("Failed to start.")
            return
        self.runBtn.configure(state="disabled")
        self.stopBtn.configure(state="normal")
        self.status.set("Running (pid %d)..." % self.runner.proc.pid)

    def on_stop(self):
        if not self.runner.running:
            return  # already finished; never overwrite the final status
        self._stop_requested = True
        self.runner.stop()
        self.status.set("Stopping...")

    def on_open_output(self):
        path = self._resolve(self.outputDir.get().strip() or "output")
        if not os.path.isdir(path):
            messagebox.showinfo("PhoSim", "Output directory does not exist yet.",
                                parent=self.root)
            return
        if os.name == "nt":
            os.startfile(path)
        elif sys.platform == "darwin":
            subprocess.Popen(["open", path])
        else:
            subprocess.Popen(["xdg-open", path])

    def on_close(self):
        if self.runner.running:
            if not messagebox.askyesno("PhoSim", "A run is in progress. Stop it and exit?",
                                       parent=self.root):
                return
            self.runner.stop()
        self._save_settings()
        self.root.destroy()

    # ---- event pump -------------------------------------------------------

    def _poll(self):
        try:
            while True:
                kind, payload = self.runner.events.get_nowait()
                if kind == "line":
                    self._append_log(payload)
                elif kind == "done" and not self._done_handled:
                    self._done_handled = True
                    self._on_done(payload)
        except queue.Empty:
            pass
        if self.runner.running and self._t0 is not None:
            self.status.set("Running (pid %d, %ds)..."
                            % (self.runner.proc.pid, int(time.time() - self._t0)))
        elif not self._done_handled and self.runner.proc is not None:
            # Watchdog: the process is dead but no "done" arrived. This happens
            # when the tree kill raced a child spawn and a survivor still holds
            # the stdout write end, keeping the reader thread blocked. Give the
            # normal path a grace period, then synthesize completion.
            if self._dead_since is None:
                self._dead_since = time.time()
            elif time.time() - self._dead_since > 2.0:
                self._done_handled = True
                self._on_done(self.runner.proc.returncode)
        self.root.after(self.POLL_MS, self._poll)

    def _on_done(self, rc):
        elapsed = int(time.time() - self._t0) if self._t0 else 0
        self.runBtn.configure(state="normal")
        self.stopBtn.configure(state="disabled")
        if getattr(self, "_stop_requested", False):
            self._append_log("[gui] run stopped by user after %ds" % elapsed)
            self.status.set("Stopped by user (%ds)." % elapsed)
            return
        if rc == 0:
            produced = sorted(glob.glob(
                os.path.join(self._resolve(self.outputDir.get().strip() or "output"),
                             "*.fits*")))
            self._append_log("[gui] finished OK in %ds; %d FITS file(s) in output:"
                             % (elapsed, len(produced)))
            for p in produced:
                self._append_log("[gui]   " + os.path.basename(p))
            self.status.set("Finished OK (%ds, %d FITS)." % (elapsed, len(produced)))
        else:
            self._append_log("[gui] phosim exited with code %s after %ds" % (rc, elapsed))
            self.status.set("Failed (exit %s)." % rc)

    # ---- settings ---------------------------------------------------------

    _PERSISTED = ("catalog", "command", "instrument", "sensor", "processors",
                  "threads", "run_e2adc", "outputDir", "workDir", "extraArgs")

    def _save_settings(self):
        data = {}
        for name in self._PERSISTED:
            try:
                data[name] = getattr(self, name).get()
            except tk.TclError:
                pass
        try:
            with open(SETTINGS_FILE, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=1)
        except OSError:
            pass

    def _load_settings(self):
        try:
            with open(SETTINGS_FILE, encoding="utf-8") as f:
                data = json.load(f)
        except (OSError, ValueError):
            return
        for name in self._PERSISTED:
            if name in data:
                try:
                    getattr(self, name).set(data[name])
                except tk.TclError:
                    pass


def main():
    root = tk.Tk()
    try:
        ttk.Style().theme_use("vista" if os.name == "nt" else "clam")
    except tk.TclError:
        pass
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()

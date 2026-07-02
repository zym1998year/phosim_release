#!/usr/bin/env python
##
## @package phosim
## @file phosim_gui2.py
## @brief compact, reliable Tk GUI front-end for phosim.py
##
## A stdlib-only alternative to the legacy PAGE-generated phosim_gui.
## Parameter coverage follows the official documentation split
## (www.phosim.org/documentation/physics-commands): the "Physics Module
## Switches" group is exposed as controls, the ~200 "Advanced Commands"
## stay reachable through the free-form custom-commands box (appended last,
## so they override everything else - PhoSim applies later lines over
## earlier ones). Observation parameters are written to a generated
## instance catalog only when the user sets them; empty means PhoSim's
## own default/random behavior (the official docs state no key is required).
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
import re
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
# palette / fonts
# ---------------------------------------------------------------------------

BG = "#eceff4"        # window background
CARD = "#ffffff"      # section background
ACCENT = "#2563eb"    # primary action
ACCENT_DARK = "#1d4ed8"
TXT = "#111827"
SUB = "#6b7280"
OK_GREEN = "#15803d"
ERR_RED = "#b91c1c"
RUN_ORANGE = "#b45309"
LOG_BG = "#16181d"
LOG_FG = "#d7dae0"

FONT = ("Segoe UI", 10) if os.name == "nt" else ("Helvetica", 11)
FONT_BOLD = (FONT[0], FONT[1], "bold")
FONT_SMALL = (FONT[0], FONT[1] - 1)
FONT_MONO = ("Consolas", 9) if os.name == "nt" else ("Courier", 10)


def apply_style(root):
    style = ttk.Style(root)
    try:
        style.theme_use("clam")
    except tk.TclError:
        pass
    root.configure(background=BG)
    style.configure(".", font=FONT, background=BG, foreground=TXT)
    style.configure("TFrame", background=BG)
    style.configure("Card.TFrame", background=CARD)
    style.configure("TLabel", background=BG, foreground=TXT)
    style.configure("Card.TLabel", background=CARD, foreground=TXT)
    style.configure("Sub.TLabel", background=CARD, foreground=SUB, font=FONT_SMALL)
    style.configure("Head.TLabel", background=CARD, foreground=TXT, font=FONT_BOLD)
    style.configure("Status.TLabel", background=BG, foreground=SUB, font=FONT_SMALL)

    style.configure("Card.TLabelframe", background=CARD, borderwidth=1,
                    relief="solid", bordercolor="#d8dce3")
    style.configure("Card.TLabelframe.Label", background=CARD,
                    foreground=TXT, font=FONT_BOLD)

    style.configure("TNotebook", background=BG, borderwidth=0, tabmargins=(8, 6, 8, 0))
    style.configure("TNotebook.Tab", font=FONT, padding=(16, 7))
    style.map("TNotebook.Tab",
              background=[("selected", CARD), ("!selected", BG)],
              foreground=[("selected", ACCENT), ("!selected", SUB)])

    style.configure("TButton", padding=(10, 5))
    style.configure("Accent.TButton", background=ACCENT, foreground="#ffffff",
                    padding=(18, 6), font=FONT_BOLD, borderwidth=0)
    style.map("Accent.TButton",
              background=[("active", ACCENT_DARK), ("disabled", "#9db4e8")])

    for w in ("TCheckbutton", "TRadiobutton"):
        style.configure("Card.%s" % w, background=CARD, foreground=TXT)
        style.map("Card.%s" % w, background=[("active", CARD)])

    style.configure("TEntry", fieldbackground="#ffffff")
    style.configure("TCombobox", fieldbackground="#ffffff")
    style.configure("TSpinbox", fieldbackground="#ffffff")


# ---------------------------------------------------------------------------
# repo introspection helpers
# ---------------------------------------------------------------------------

def list_instruments():
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
    """POSIX shlex with escapes disabled: quotes strip on every platform,
    Windows backslash paths stay literal."""
    lex = shlex.shlex(text, posix=True)
    lex.whitespace_split = True
    lex.escape = ""
    return list(lex)


def child_environment():
    """conda on Windows keeps the cfitsio/fftw/zlib DLLs in <prefix>\\Library\\bin;
    prepend it so raytrace/e2adc resolve them even from a double-click launch."""
    env = os.environ.copy()
    if os.name == "nt":
        dllDir = os.path.join(sys.prefix, "Library", "bin")
        if os.path.isdir(dllDir):
            env["PATH"] = dllDir + os.pathsep + env.get("PATH", "")
    return env


# ---------------------------------------------------------------------------
# subprocess runner (worker thread; UI updated only from the Tk main loop)
# ---------------------------------------------------------------------------

class PhosimRunner:
    def __init__(self):
        self.events = queue.Queue()
        self.proc = None
        self._thread = None

    @property
    def running(self):
        return self.proc is not None and self.proc.poll() is None

    def start(self, argv, cwd):
        creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        self.proc = subprocess.Popen(
            argv, cwd=cwd, env=child_environment(),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace", bufsize=1,
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
            pass  # pipe closed by stop()
        except Exception as e:
            self.events.put(("line", "[gui] output reader error: %r" % (e,)))
        rc = self.proc.wait()
        self.events.put(("done", rc))

    def stop(self):
        if not self.running:
            return
        try:
            if os.name == "nt":
                subprocess.call(["taskkill", "/T", "/F", "/PID", str(self.proc.pid)],
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                creationflags=subprocess.CREATE_NO_WINDOW)
            else:
                os.killpg(os.getpgid(self.proc.pid), signal.SIGTERM)
        except Exception as e:
            self.events.put(("line", "[gui] stop failed: %r" % (e,)))
        try:
            if self.proc.stdout:
                self.proc.stdout.close()   # unblock the reader thread
        except OSError:
            pass


# ---------------------------------------------------------------------------
# parameter tables (ranges follow the legacy GUI's validation, [min, max))
# ---------------------------------------------------------------------------

# key, label, unit/hint, (lo, hi) or None, is_float
OBS_FIELDS = [
    ("rightascension", "Right ascension", "deg", (0.0, 360.0), True),
    ("declination", "Declination", "deg", (-90.0, 90.0), True),
    ("azimuth", "Azimuth", "deg", (0.0, 360.0), True),
    ("altitude", "Altitude", "deg", (0.0, 90.0), True),
    ("rottelpos", "Telescope rotation", "deg", (-360.0, 360.0), True),
    ("mjd", "MJD / 'now' / 'tonight'", "", None, False),
    ("vistime", "Visit exposure time", "s", (0.0001, 10000.0), True),
    ("nsnap", "Number of snaps", "", (1, 1000000), False),
    ("filter", "Filter index", "0-5", (0, 8), False),
    ("seed", "Random seed (obsseed)", "", (1, 100000000), False),
    ("seeing", "Seeing @500nm", "arcsec", (0.0, 10.0), True),
    ("pressure", "Air pressure", "mmHg", (0.0, 2000.0), True),
    ("temperature", "Air temperature", "degC", (-274.0, 400.0), True),
]

# Physics Module Switches (official docs group). kind:
#   sw1  default-on switch, write "key 0" when unchecked
#   sw0  default-off switch, write "key 1" when checked
#   bare default-off clear command, write bare "key" when checked
DET_SWITCHES = [
    ("detectormode", "Charge diffusion", "sw1"),
    ("chargesharing", "Charge sharing", "sw1"),
    ("pixelerror", "Pixel boundary errors", "sw1"),
    ("fringing", "Fringing", "sw1"),
    ("fieldanisotropy", "Field anisotropy", "sw1"),
    ("impurityvariation", "Impurity variation", "sw1"),
    ("saturation", "Saturation", "sw1"),
    ("blooming", "Blooming", "sw1"),
    ("cleardefects", "Clear detector defects", "bare"),
]
ATM_SWITCHES = [
    ("atmosphericdispersion", "Atmospheric dispersion", "sw1"),
    ("clearturbulence", "Clear turbulence", "bare"),
    ("clearclouds", "Clear clouds", "bare"),
    ("clearopacity", "Clear opacity", "bare"),
]
TEL_SWITCHES = [
    ("trackingmode", "Tracking errors", "sw1"),
    ("windshake", "Wind shake", "sw1"),
    ("coatingmode", "Optics coatings", "sw1"),
    ("contaminationmode", "Optics contamination", "sw1"),
    ("clearperturbations", "Clear perturbations", "bare"),
    ("cleartracking", "Clear tracking", "bare"),
]


class App:
    POLL_MS = 100
    MAX_LOG_LINES = 5000

    def __init__(self, root):
        self.root = root
        self.runner = PhosimRunner()
        self._t0 = None
        self._stop_requested = False
        self._done_handled = True
        self._dead_since = None
        root.title("PhoSim launcher")
        root.minsize(940, 700)
        apply_style(root)
        root.protocol("WM_DELETE_WINDOW", self.on_close)

        self._build_vars()
        self._build_ui()
        self._load_settings()
        self._sync_modes()
        self._refresh_preview()
        self._poll()

    # ---- state ----------------------------------------------------------

    def _build_vars(self):
        v = tk.StringVar
        # execution
        self.instrument = v(value="generic")
        self.sensor = v(value="")
        self.processors = tk.IntVar(value=1)
        self.threads = tk.IntVar(value=1)
        self.run_e2adc = tk.BooleanVar(value=True)
        self.run_ds9 = tk.BooleanVar(value=False)
        self.outputDir = v(value="output")
        self.workDir = v(value="work")
        self.sedDir = v(value="")
        self.extraArgs = v(value="")
        # input mode
        self.inputMode = v(value="file")          # file | build
        self.catalog = v(value="examples/star")
        # observation (build mode); empty string = PhoSim default
        self.obshistid = v(value="9999")
        self.obs = {key: v(value="") for key, *_ in OBS_FIELDS}
        # source
        self.sourceMode = v(value="star")         # star|grid|stars|galaxies|include|custom
        self.srcStarMag = v(value="20")
        self.srcGrid = [v(value="1"), v(value="1"), v(value="20")]   # spacing width mag
        self.srcStars = [v(value="15"), v(value="20"), v(value="1")]  # minM maxM fov
        self.srcGalaxies = [v(value="15"), v(value="20"), v(value="1")]
        self.srcInclude = v(value="")
        # physics
        self.cmdMode = v(value="none")            # none | file
        self.cmdFile = v(value="examples/nobackground")
        self.background = v(value="normal")       # normal|quick|single|off
        self.atmosphereMode = v(value="2")        # 2|1|0
        self.perfectTelescope = tk.BooleanVar(value=False)
        self.opticsOnly = tk.BooleanVar(value=False)
        self.clearEverything = tk.BooleanVar(value=False)
        self.diffraction = v(value="1")           # 1 MC | 2 FFT | 0 off
        self.pertControl = tk.BooleanVar(value=True)
        self.pertEnv = tk.BooleanVar(value=True)
        self.switches = {}
        for key, _label, kind in DET_SWITCHES + ATM_SWITCHES + TEL_SWITCHES:
            self.switches[key] = tk.BooleanVar(value=(kind == "sw1"))
        self.status = v(value="Ready.")

        watched = ([self.instrument, self.sensor, self.extraArgs, self.outputDir,
                    self.workDir, self.sedDir, self.inputMode, self.catalog,
                    self.obshistid, self.sourceMode, self.cmdMode, self.cmdFile,
                    self.background, self.atmosphereMode, self.diffraction,
                    self.processors, self.threads, self.run_e2adc, self.run_ds9,
                    self.perfectTelescope, self.opticsOnly, self.clearEverything,
                    self.pertControl, self.pertEnv]
                   + list(self.obs.values()) + list(self.switches.values()))
        for var in watched:
            var.trace_add("write", lambda *_: self._on_state_change())

    def _on_state_change(self):
        self._sync_modes()
        self._refresh_preview()

    # ---- layout ----------------------------------------------------------

    def _card(self, parent, title):
        f = ttk.LabelFrame(parent, text=" %s " % title, style="Card.TLabelframe",
                           padding=(10, 6, 10, 8))
        return f

    def _build_ui(self):
        outer = ttk.Frame(self.root, padding=(10, 8, 10, 6))
        outer.grid(sticky="nsew")
        self.root.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)
        outer.columnconfigure(0, weight=1)
        outer.rowconfigure(1, weight=3)   # notebook
        outer.rowconfigure(3, weight=2)   # log

        # header
        head = ttk.Frame(outer)
        head.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        head.columnconfigure(1, weight=1)
        ttk.Label(head, text="PhoSim", font=(FONT[0], 15, "bold")).grid(row=0, column=0, sticky="w")
        ttk.Label(head, text="  photon Monte Carlo image simulator",
                  foreground=SUB).grid(row=0, column=1, sticky="w")
        self.statusLbl = ttk.Label(head, textvariable=self.status, style="Status.TLabel")
        self.statusLbl.grid(row=0, column=2, sticky="e")

        nb = ttk.Notebook(outer)
        nb.grid(row=1, column=0, sticky="nsew")
        self.nb = nb
        tabObs = ttk.Frame(nb, style="Card.TFrame", padding=10)
        tabPhy = ttk.Frame(nb, style="Card.TFrame", padding=10)
        tabExe = ttk.Frame(nb, style="Card.TFrame", padding=10)
        nb.add(tabObs, text="Observation")
        nb.add(tabPhy, text="Physics")
        nb.add(tabExe, text="Execution")
        self._build_tab_observation(tabObs)
        self._build_tab_physics(tabPhy)
        self._build_tab_execution(tabExe)

        # action row
        act = ttk.Frame(outer)
        act.grid(row=2, column=0, sticky="ew", pady=6)
        act.columnconfigure(0, weight=1)
        self.preview = ttk.Entry(act, state="readonly", font=FONT_MONO)
        self.preview.grid(row=0, column=0, sticky="ew", padx=(0, 8), ipady=3)
        self.runBtn = ttk.Button(act, text="▶  Run PhoSim", style="Accent.TButton",
                                 command=self.on_run)
        self.runBtn.grid(row=0, column=1, padx=(0, 6))
        self.stopBtn = ttk.Button(act, text="Stop", command=self.on_stop, state="disabled")
        self.stopBtn.grid(row=0, column=2, padx=(0, 6))
        ttk.Button(act, text="Open output", command=self.on_open_output).grid(row=0, column=3)

        # log console
        self.log = ScrolledText(outer, height=11, state="disabled", wrap="none",
                                font=FONT_MONO, background=LOG_BG, foreground=LOG_FG,
                                insertbackground=LOG_FG, relief="flat",
                                padx=8, pady=6)
        self.log.grid(row=3, column=0, sticky="nsew")
        self.log.tag_configure("gui", foreground="#7aa2f7")
        self.log.tag_configure("err", foreground="#f7768e")

    # -- Observation tab

    def _build_tab_observation(self, tab):
        tab.columnconfigure(0, weight=1, uniform="col")
        tab.columnconfigure(1, weight=1, uniform="col")
        tab.rowconfigure(1, weight=1)

        # input mode
        mode = self._card(tab, "Instance catalog")
        mode.grid(row=0, column=0, columnspan=2, sticky="ew", pady=(0, 8))
        mode.columnconfigure(2, weight=1)
        ttk.Radiobutton(mode, text="Use an existing catalog file", value="file",
                        variable=self.inputMode, style="Card.TRadiobutton")\
            .grid(row=0, column=0, sticky="w", padx=(0, 10))
        ttk.Radiobutton(mode, text="Build the observation here", value="build",
                        variable=self.inputMode, style="Card.TRadiobutton")\
            .grid(row=0, column=1, sticky="w")
        self.catalogCombo = ttk.Combobox(mode, textvariable=self.catalog, values=list_examples())
        self.catalogCombo.grid(row=1, column=0, columnspan=3, sticky="ew", pady=(6, 0))
        self.catalogBrowse = ttk.Button(mode, text="Browse...",
                                        command=lambda: self._pick_file(self.catalog))
        self.catalogBrowse.grid(row=1, column=3, padx=(8, 0), pady=(6, 0))

        # observation parameters
        obs = self._card(tab, "Observation parameters (empty = PhoSim default / random)")
        obs.grid(row=1, column=0, sticky="nsew", padx=(0, 8))
        obs.columnconfigure(1, weight=1)
        self.obsWidgets = []
        ttk.Label(obs, text="Observation ID", style="Card.TLabel").grid(row=0, column=0, sticky="w", pady=2)
        w = ttk.Entry(obs, textvariable=self.obshistid, width=14)
        w.grid(row=0, column=1, sticky="w", padx=8, pady=2)
        self.obsWidgets.append(w)
        for i, (key, label, unit, rng, _isf) in enumerate(OBS_FIELDS, start=1):
            ttk.Label(obs, text=label, style="Card.TLabel").grid(row=i, column=0, sticky="w", pady=2)
            w = ttk.Entry(obs, textvariable=self.obs[key], width=14)
            w.grid(row=i, column=1, sticky="w", padx=8, pady=2)
            self.obsWidgets.append(w)
            hint = unit
            if rng:
                hint = ("%s  [%g, %g)" % (unit, rng[0], rng[1])).strip()
            ttk.Label(obs, text=hint, style="Sub.TLabel").grid(row=i, column=2, sticky="w")

        # sources
        src = self._card(tab, "Sources (at least one is needed to make photons)")
        src.grid(row=1, column=1, sticky="nsew")
        src.columnconfigure(1, weight=1)
        self.srcWidgets = []

        def srow(r, value, text, entries, hints):
            rb = ttk.Radiobutton(src, text=text, value=value, variable=self.sourceMode,
                                 style="Card.TRadiobutton")
            rb.grid(row=r, column=0, sticky="w", pady=3)
            self.srcWidgets.append(rb)
            fr = ttk.Frame(src, style="Card.TFrame")
            fr.grid(row=r, column=1, sticky="w")
            for j, (var, hint) in enumerate(zip(entries, hints)):
                e = ttk.Entry(fr, textvariable=var, width=7)
                e.grid(row=0, column=2 * j, padx=(6, 2))
                self.srcWidgets.append(e)
                ttk.Label(fr, text=hint, style="Sub.TLabel").grid(row=0, column=2 * j + 1)

        srow(0, "star", "Single star", [self.srcStarMag], ["mag"])
        srow(1, "grid", "Star grid", self.srcGrid, ["deg", "deg", "mag"])
        srow(2, "stars", "Random stars", self.srcStars, ["minM", "maxM", "fov°"])
        srow(3, "galaxies", "Random galaxies", self.srcGalaxies, ["minM", "maxM", "fov°"])
        rb = ttk.Radiobutton(src, text="Include catalog file", value="include",
                             variable=self.sourceMode, style="Card.TRadiobutton")
        rb.grid(row=4, column=0, sticky="w", pady=3)
        self.srcWidgets.append(rb)
        fr = ttk.Frame(src, style="Card.TFrame")
        fr.grid(row=4, column=1, sticky="ew")
        fr.columnconfigure(0, weight=1)
        e = ttk.Entry(fr, textvariable=self.srcInclude)
        e.grid(row=0, column=0, sticky="ew", padx=(6, 4))
        self.srcWidgets.append(e)
        b = ttk.Button(fr, text="...", width=3, command=lambda: self._pick_file(self.srcInclude))
        b.grid(row=0, column=1)
        self.srcWidgets.append(b)
        rb = ttk.Radiobutton(src, text="Custom object lines (below)", value="custom",
                             variable=self.sourceMode, style="Card.TRadiobutton")
        rb.grid(row=5, column=0, sticky="w", pady=3)
        self.srcWidgets.append(rb)
        src.rowconfigure(6, weight=1)
        self.objText = tk.Text(src, height=5, font=FONT_MONO, relief="solid",
                               borderwidth=1, highlightthickness=0)
        self.objText.grid(row=6, column=0, columnspan=2, sticky="nsew", pady=(4, 0))
        self.srcWidgets.append(self.objText)
        ttk.Label(src, text="e.g.  object 0 187.5 -0.6 20 ../sky/sed_flat.txt "
                            "0 0 0 0 0 0 point none none",
                  style="Sub.TLabel").grid(row=7, column=0, columnspan=2, sticky="w")

    # -- Physics tab

    def _build_tab_physics(self, tab):
        for c in range(3):
            tab.columnconfigure(c, weight=1, uniform="p")
        tab.rowconfigure(2, weight=1)

        # command file + presets row
        top = self._card(tab, "Physics command file (lines later in the file override earlier ones)")
        top.grid(row=0, column=0, columnspan=3, sticky="ew", pady=(0, 8))
        top.columnconfigure(2, weight=1)
        ttk.Radiobutton(top, text="None (switches below only)", value="none",
                        variable=self.cmdMode, style="Card.TRadiobutton")\
            .grid(row=0, column=0, sticky="w", padx=(0, 12))
        ttk.Radiobutton(top, text="Start from file:", value="file",
                        variable=self.cmdMode, style="Card.TRadiobutton")\
            .grid(row=0, column=1, sticky="w")
        self.cmdCombo = ttk.Combobox(top, textvariable=self.cmdFile, values=list_examples())
        self.cmdCombo.grid(row=0, column=2, sticky="ew", padx=8)
        self.cmdBrowse = ttk.Button(top, text="Browse...",
                                    command=lambda: self._pick_file(self.cmdFile))
        self.cmdBrowse.grid(row=0, column=3)

        # master modes
        master = self._card(tab, "Master modes")
        master.grid(row=1, column=0, columnspan=3, sticky="ew", pady=(0, 8))
        ttk.Label(master, text="Sky background:", style="Card.TLabel").grid(row=0, column=0, sticky="w")
        for j, (val, txt) in enumerate([("normal", "Normal"), ("quick", "Quick"),
                                        ("single", "Single-photon"), ("off", "Off")]):
            ttk.Radiobutton(master, text=txt, value=val, variable=self.background,
                            style="Card.TRadiobutton").grid(row=0, column=1 + j, sticky="w", padx=6)
        ttk.Label(master, text="Atmosphere:", style="Card.TLabel").grid(row=1, column=0, sticky="w", pady=(4, 0))
        for j, (val, txt) in enumerate([("2", "Full"), ("1", "Diffraction only"), ("0", "Off")]):
            ttk.Radiobutton(master, text=txt, value=val, variable=self.atmosphereMode,
                            style="Card.TRadiobutton").grid(row=1, column=1 + j, sticky="w", padx=6, pady=(4, 0))
        ttk.Label(master, text="Diffraction:", style="Card.TLabel").grid(row=2, column=0, sticky="w", pady=(4, 0))
        for j, (val, txt) in enumerate([("1", "Monte Carlo"), ("2", "FFT"), ("0", "Off")]):
            ttk.Radiobutton(master, text=txt, value=val, variable=self.diffraction,
                            style="Card.TRadiobutton").grid(row=2, column=1 + j, sticky="w", padx=6, pady=(4, 0))
        cb1 = ttk.Checkbutton(master, text="Perfect telescope", variable=self.perfectTelescope,
                              style="Card.TCheckbutton")
        cb1.grid(row=0, column=6, sticky="w", padx=(20, 0))
        cb2 = ttk.Checkbutton(master, text="Optics only", variable=self.opticsOnly,
                              style="Card.TCheckbutton")
        cb2.grid(row=1, column=6, sticky="w", padx=(20, 0), pady=(4, 0))
        cb3 = ttk.Checkbutton(master, text="Clear everything", variable=self.clearEverything,
                              style="Card.TCheckbutton")
        cb3.grid(row=2, column=6, sticky="w", padx=(20, 0), pady=(4, 0))

        # switch groups
        def group(col, title, items, extra=None):
            g = self._card(tab, title)
            g.grid(row=2, column=col, sticky="nsew", padx=(0 if col == 0 else 8, 0))
            for r, (key, label, kind) in enumerate(items):
                ttk.Checkbutton(g, text=label, variable=self.switches[key],
                                style="Card.TCheckbutton").grid(row=r, column=0, sticky="w", pady=2)
            if extra:
                extra(g, len(items))
            return g

        group(0, "Detector", DET_SWITCHES)
        group(1, "Atmosphere", ATM_SWITCHES)

        def tel_extra(g, r0):
            ttk.Separator(g).grid(row=r0, column=0, sticky="ew", pady=4)
            ttk.Label(g, text="Perturbations:", style="Card.TLabel").grid(row=r0 + 1, column=0, sticky="w")
            ttk.Checkbutton(g, text="Control-system errors", variable=self.pertControl,
                            style="Card.TCheckbutton").grid(row=r0 + 2, column=0, sticky="w", pady=2)
            ttk.Checkbutton(g, text="Environment & fabrication", variable=self.pertEnv,
                            style="Card.TCheckbutton").grid(row=r0 + 3, column=0, sticky="w", pady=2)
        group(2, "Telescope & optics", TEL_SWITCHES, tel_extra)

        # custom commands
        cust = self._card(tab, "Custom physics commands (appended last = highest priority; "
                               "any of the ~200 advanced keys)")
        cust.grid(row=3, column=0, columnspan=3, sticky="ew", pady=(8, 0))
        cust.columnconfigure(0, weight=1)
        self.customCmds = tk.Text(cust, height=3, font=FONT_MONO, relief="solid",
                                  borderwidth=1, highlightthickness=0)
        self.customCmds.grid(row=0, column=0, sticky="ew")

    # -- Execution tab

    def _build_tab_execution(self, tab):
        tab.columnconfigure(0, weight=1, uniform="e")
        tab.columnconfigure(1, weight=1, uniform="e")

        run = self._card(tab, "Run")
        run.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
        run.columnconfigure(1, weight=1)
        ttk.Label(run, text="Instrument (-i)", style="Card.TLabel").grid(row=0, column=0, sticky="w", pady=3)
        ttk.Combobox(run, textvariable=self.instrument, values=list_instruments())\
            .grid(row=0, column=1, sticky="ew", padx=8)
        ttk.Label(run, text="Sensor (-s)", style="Card.TLabel").grid(row=1, column=0, sticky="w", pady=3)
        ttk.Entry(run, textvariable=self.sensor).grid(row=1, column=1, sticky="ew", padx=8)
        ttk.Label(run, text="chip1|chip2 (empty = all)", style="Sub.TLabel").grid(row=2, column=1, sticky="w", padx=8)
        ttk.Label(run, text="Processes (-p)", style="Card.TLabel").grid(row=3, column=0, sticky="w", pady=3)
        ttk.Spinbox(run, from_=1, to=max(1, os.cpu_count() or 1), textvariable=self.processors,
                    width=6).grid(row=3, column=1, sticky="w", padx=8)
        ttk.Label(run, text="Threads (-t)", style="Card.TLabel").grid(row=4, column=0, sticky="w", pady=3)
        ttk.Spinbox(run, from_=1, to=max(1, os.cpu_count() or 1), textvariable=self.threads,
                    width=6).grid(row=4, column=1, sticky="w", padx=8)
        ttk.Checkbutton(run, text="Run e2adc amplifier images (-e)", variable=self.run_e2adc,
                        style="Card.TCheckbutton").grid(row=5, column=0, columnspan=2, sticky="w", pady=3)
        ttk.Checkbutton(run, text="Open result in ds9 (--ds9)", variable=self.run_ds9,
                        style="Card.TCheckbutton").grid(row=6, column=0, columnspan=2, sticky="w", pady=3)
        ttk.Label(run, text="Extra phosim.py args", style="Card.TLabel").grid(row=7, column=0, sticky="w", pady=3)
        ttk.Entry(run, textvariable=self.extraArgs).grid(row=7, column=1, sticky="ew", padx=8)
        ttk.Label(run, text="e.g. --keepscreens=1 --sensor ...", style="Sub.TLabel")\
            .grid(row=8, column=1, sticky="w", padx=8)

        dirs = self._card(tab, "Directories (relative to the phosim root)")
        dirs.grid(row=0, column=1, sticky="nsew")
        dirs.columnconfigure(1, weight=1)
        for r, (label, var) in enumerate([("Output (-o)", self.outputDir),
                                          ("Work (-w)", self.workDir),
                                          ("SED (--sed)", self.sedDir)]):
            ttk.Label(dirs, text=label, style="Card.TLabel").grid(row=r, column=0, sticky="w", pady=3)
            ttk.Entry(dirs, textvariable=var).grid(row=r, column=1, sticky="ew", padx=8)
            ttk.Button(dirs, text="Browse...",
                       command=lambda v=var: self._pick_dir(v)).grid(row=r, column=2)
        ttk.Label(dirs, text="SED empty = data/SEDs default", style="Sub.TLabel")\
            .grid(row=3, column=1, sticky="w", padx=8)

    # ---- small helpers ---------------------------------------------------

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
        try:
            rel = os.path.relpath(path, PHOSIM_ROOT)
            if not rel.startswith(".."):
                return rel.replace("\\", "/")
        except ValueError:
            pass
        return path

    def _resolve(self, path):
        return path if os.path.isabs(path) else os.path.join(PHOSIM_ROOT, path)

    def _append_log(self, text, tag=None):
        self.log.configure(state="normal")
        self.log.insert("end", text + "\n", tag or ())
        lines = int(self.log.index("end-1c").split(".")[0])
        if lines > self.MAX_LOG_LINES:
            self.log.delete("1.0", "%d.0" % (lines - self.MAX_LOG_LINES + 1))
        self.log.see("end")
        self.log.configure(state="disabled")

    def _sync_modes(self):
        """Enable/disable widget groups according to the selected modes."""
        build = self.inputMode.get() == "build"
        try:
            self.catalogCombo.configure(state="normal" if not build else "disabled")
            self.catalogBrowse.configure(state="normal" if not build else "disabled")
            for w in self.obsWidgets + self.srcWidgets:
                if isinstance(w, tk.Text):
                    w.configure(state="normal" if build else "disabled",
                                background="#ffffff" if build else "#f3f4f6")
                else:
                    w.configure(state="normal" if build else "disabled")
            useFile = self.cmdMode.get() == "file"
            self.cmdCombo.configure(state="normal" if useFile else "disabled")
            self.cmdBrowse.configure(state="normal" if useFile else "disabled")
        except AttributeError:
            pass  # during construction

    # ---- validation ------------------------------------------------------

    def _validate_obs_ranges(self):
        for key, label, _unit, rng, isf in OBS_FIELDS:
            raw = self.obs[key].get().strip()
            if not raw or rng is None:
                continue
            if key == "mjd" and raw in ("now", "tonight"):
                continue
            try:
                val = float(raw) if isf else int(raw)
            except ValueError:
                return "%s: '%s' is not a number." % (label, raw)
            lo, hi = rng
            if not (lo <= val < hi):
                return "%s must be in [%g, %g)." % (label, lo, hi)
        ra = self.obs["rightascension"].get().strip()
        dec = self.obs["declination"].get().strip()
        az = self.obs["azimuth"].get().strip()
        alt = self.obs["altitude"].get().strip()
        if bool(ra) != bool(dec):
            return "Right ascension and declination must be set together."
        if bool(az) != bool(alt):
            return "Azimuth and altitude must be set together."
        return None

    def _validate(self):
        try:
            self.processors.get()
            self.threads.get()
        except tk.TclError:
            return "Processes (-p) and threads (-t) must be integers."
        obsid = self.obshistid.get().strip()
        if obsid and not re.fullmatch(r"[A-Za-z0-9_-]+", obsid):
            return "Observation ID may only contain letters, digits, '_' and '-' " \
                   "(it becomes part of the output file names)."
        if self.inputMode.get() == "file":
            catalog = self.catalog.get().strip()
            if not catalog:
                return "Please choose an instance catalog file."
            if not os.path.isfile(self._resolve(catalog)):
                return "Instance catalog not found: %s" % catalog
        else:
            err = self._validate_obs_ranges()
            if err:
                return err
            if self.sourceMode.get() == "include":
                inc = self.srcInclude.get().strip()
                if not inc or not os.path.isfile(self._resolve(inc)):
                    return "Include catalog file not found: %s" % (inc or "(empty)")
            if self.sourceMode.get() == "custom" and not self.objText.get("1.0", "end").strip():
                return "Custom source mode selected but no object lines were given."
        if self.cmdMode.get() == "file":
            cf = self.cmdFile.get().strip()
            if not cf or not os.path.isfile(self._resolve(cf)):
                return "Physics command file not found: %s" % (cf or "(empty)")
        for label, var in (("output", self.outputDir), ("work", self.workDir)):
            d = var.get().strip()
            if d:
                try:
                    os.makedirs(self._resolve(d), exist_ok=True)
                except OSError as e:
                    return "Cannot create %s directory %s: %s" % (label, d, e)
        sed = self.sedDir.get().strip()
        if sed and not os.path.isdir(self._resolve(sed)):
            return "SED directory not found: %s" % sed
        return None

    # ---- file generation ---------------------------------------------------

    def _generated_dir(self):
        d = os.path.join(self._resolve(self.workDir.get().strip() or "work"))
        os.makedirs(d, exist_ok=True)
        return d

    def generate_instance_lines(self):
        """Instance catalog content for build mode (white-list: only set keys)."""
        lines = ["# generated by phosim_gui2"]
        obsid = self.obshistid.get().strip() or "9999"
        lines.append("obshistid %s" % obsid)
        for key, *_ in OBS_FIELDS:
            val = self.obs[key].get().strip()
            if val:
                lines.append("%s %s" % (key, val))
        mode = self.sourceMode.get()
        if mode == "star":
            mag = self.srcStarMag.get().strip() or "20"
            lines.append("stargrid 1 1 %s" % mag)
        elif mode == "grid":
            lines.append("stargrid %s" % " ".join(v.get().strip() or "1" for v in self.srcGrid))
        elif mode == "stars":
            lines.append("stars %s" % " ".join(v.get().strip() or "15" for v in self.srcStars))
        elif mode == "galaxies":
            lines.append("galaxies %s" % " ".join(v.get().strip() or "15" for v in self.srcGalaxies))
        elif mode == "include":
            lines.append("includeobj %s" % self.srcInclude.get().strip())
        elif mode == "custom":
            lines.extend(l for l in self.objText.get("1.0", "end").splitlines() if l.strip())
        return lines

    def _base_command_lines(self):
        """Content of the base physics file, cached by (path, mtime) so the
        per-keystroke preview refresh does not re-read large files."""
        if self.cmdMode.get() != "file":
            return []
        rel = self.cmdFile.get().strip()
        path = self._resolve(rel)
        try:
            mtime = os.path.getmtime(path)
            cached = getattr(self, "_base_cache", None)
            if cached and cached[0] == path and cached[1] == mtime:
                return list(cached[2])
            with open(path, encoding="utf-8", errors="replace") as f:
                lines = ["# --- from %s" % rel] + [l.rstrip("\n") for l in f]
            self._base_cache = (path, mtime, list(lines))
            return lines
        except OSError as e:
            return ["# could not read base file: %s" % e]

    def generate_override_lines(self):
        """GUI-driven physics lines only (without the base file). Order matters
        (later overrides earlier): background preset -> master modes ->
        switches -> custom."""
        lines = []
        bg = self.background.get()
        if bg == "quick":
            lines.append("quickbackground")
        elif bg == "single":
            lines.append("singlephotonbackground")
        elif bg == "off":
            lines.append("backgroundmode 0")
        if self.clearEverything.get():
            lines.append("cleareverything")
        if self.atmosphereMode.get() != "2":
            lines.append("atmospheremode %s" % self.atmosphereMode.get())
        if self.perfectTelescope.get():
            lines.append("telescopemode 0")
        if self.opticsOnly.get():
            lines.append("opticsonlymode 1")
        if self.diffraction.get() != "1":
            lines.append("diffractionmode %s" % self.diffraction.get())
        pert = (1 if self.pertControl.get() else 0) + (2 if self.pertEnv.get() else 0)
        if pert != 3:
            lines.append("perturbationmode %d" % pert)
        for key, _label, kind in DET_SWITCHES + ATM_SWITCHES + TEL_SWITCHES:
            on = self.switches[key].get()
            if kind == "sw1" and not on:
                lines.append("%s 0" % key)
            elif kind == "bare" and on:
                lines.append(key)
        custom = self.customCmds.get("1.0", "end").strip()
        if custom:
            lines.append("# --- custom")
            lines.extend(l for l in custom.splitlines() if l.strip())
        return lines

    def generate_command_lines(self):
        return self._base_command_lines() + self.generate_override_lines()

    def build_argv(self, write_files=False):
        """argv for phosim.py; generated files are only written when running."""
        obsid = self.obshistid.get().strip() or "9999"
        # the id goes into generated file names; keep it filesystem-safe
        safeid = re.sub(r"[^A-Za-z0-9_-]", "_", obsid) or "9999"
        if self.inputMode.get() == "file":
            catalog = self.catalog.get().strip()
        else:
            catalog = os.path.join(self.workDir.get().strip() or "work",
                                   "gui_%s.inst" % safeid).replace("\\", "/")
            if write_files:
                with open(self._resolve(catalog), "w", encoding="utf-8", newline="\n") as f:
                    f.write("\n".join(self.generate_instance_lines()) + "\n")
        overrides = self.generate_override_lines()
        cmdArg = None
        if self.cmdMode.get() == "file" and not overrides:
            # base file only, no GUI overrides: pass it through untouched
            cmdArg = self.cmdFile.get().strip()
        elif overrides or self.cmdMode.get() == "file":
            cmdArg = os.path.join(self.workDir.get().strip() or "work",
                                  "gui_%s.cmd" % safeid).replace("\\", "/")
            if write_files:
                with open(self._resolve(cmdArg), "w", encoding="utf-8", newline="\n") as f:
                    f.write("\n".join(self._base_command_lines() + overrides) + "\n")

        argv = [sys.executable, os.path.join(PHOSIM_ROOT, "phosim.py"), catalog]
        if cmdArg:
            argv += ["-c", cmdArg]
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
        if self.run_ds9.get():
            argv += ["--ds9"]
        out = self.outputDir.get().strip()
        if out and out != "output":
            argv += ["-o", out]
        work = self.workDir.get().strip()
        if work and work != "work":
            argv += ["-w", work]
        sed = self.sedDir.get().strip()
        if sed:
            argv += ["--sed", sed]
        extra = self.extraArgs.get().strip()
        if extra:
            argv += split_extra_args(extra)
        return argv

    def _refresh_preview(self):
        try:
            argv = self.build_argv(write_files=False)
        except tk.TclError:
            return
        shown = argv[:]
        shown[0] = os.path.basename(shown[0])
        shown[1] = "phosim.py"
        text = " ".join(a if " " not in a else '"%s"' % a for a in shown)
        self.preview.configure(state="normal")
        self.preview.delete(0, "end")
        self.preview.insert(0, text)
        self.preview.configure(state="readonly")

    # ---- actions -----------------------------------------------------------

    def on_run(self):
        if self.runner.running:
            return
        err = self._validate()
        if err:
            messagebox.showerror("PhoSim", err, parent=self.root)
            return
        self._save_settings()
        try:
            argv = self.build_argv(write_files=True)
        except OSError as e:
            messagebox.showerror("PhoSim", "Could not write generated files: %s" % e,
                                 parent=self.root)
            return
        self._append_log("$ " + " ".join(argv), "gui")
        if self.inputMode.get() == "build":
            self._append_log("[gui] instance catalog: %s" % argv[2], "gui")
        self._stop_requested = False
        self._done_handled = False
        self._dead_since = None
        self._t0 = time.time()
        try:
            self.runner.start(argv, cwd=PHOSIM_ROOT)
        except OSError as e:
            self._append_log("[gui] failed to start: %s" % e, "err")
            self._set_status("Failed to start.", ERR_RED)
            return
        self.runBtn.configure(state="disabled")
        self.stopBtn.configure(state="normal")
        self._set_status("Running (pid %d)..." % self.runner.proc.pid, RUN_ORANGE)

    def on_stop(self):
        if not self.runner.running:
            return
        self._stop_requested = True
        self.runner.stop()
        self._set_status("Stopping...", RUN_ORANGE)

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

    def _set_status(self, text, color=SUB):
        self.status.set(text)
        self.statusLbl.configure(foreground=color)

    # ---- event pump ----------------------------------------------------------

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
            self._set_status("Running (pid %d, %ds)..."
                             % (self.runner.proc.pid, int(time.time() - self._t0)),
                             RUN_ORANGE)
        elif not self._done_handled and self.runner.proc is not None:
            # watchdog: tree-kill raced a child spawn and a survivor holds the
            # stdout pipe; synthesize completion after a grace period
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
        if self._stop_requested:
            self._append_log("[gui] run stopped by user after %ds" % elapsed, "gui")
            self._set_status("Stopped by user (%ds)." % elapsed, RUN_ORANGE)
            return
        if rc == 0:
            produced = sorted(glob.glob(
                os.path.join(self._resolve(self.outputDir.get().strip() or "output"),
                             "*.fits*")))
            self._append_log("[gui] finished OK in %ds; %d FITS file(s) in output:"
                             % (elapsed, len(produced)), "gui")
            for p in produced:
                self._append_log("[gui]   " + os.path.basename(p), "gui")
            self._set_status("Finished OK (%ds, %d FITS)." % (elapsed, len(produced)),
                             OK_GREEN)
        else:
            self._append_log("[gui] phosim exited with code %s after %ds" % (rc, elapsed),
                             "err")
            self._set_status("Failed (exit %s)." % rc, ERR_RED)

    # ---- settings --------------------------------------------------------------

    _PERSISTED = ("instrument", "sensor", "processors", "threads", "run_e2adc",
                  "run_ds9", "outputDir", "workDir", "sedDir", "extraArgs",
                  "inputMode", "catalog", "obshistid", "sourceMode", "srcStarMag",
                  "srcInclude", "cmdMode", "cmdFile", "background",
                  "atmosphereMode", "diffraction", "perfectTelescope",
                  "opticsOnly", "clearEverything", "pertControl", "pertEnv")

    def _save_settings(self):
        data = {}
        for name in self._PERSISTED:
            try:
                data[name] = getattr(self, name).get()
            except tk.TclError:
                pass
        data["obs"] = {k: v.get() for k, v in self.obs.items()}
        data["switches"] = {k: v.get() for k, v in self.switches.items()}
        data["srcGrid"] = [v.get() for v in self.srcGrid]
        data["srcStars"] = [v.get() for v in self.srcStars]
        data["srcGalaxies"] = [v.get() for v in self.srcGalaxies]
        data["customCmds"] = self.customCmds.get("1.0", "end").strip()
        data["objText"] = self.objText.get("1.0", "end").strip()
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
        for k, val in (data.get("obs") or {}).items():
            if k in self.obs:
                self.obs[k].set(val)
        for k, val in (data.get("switches") or {}).items():
            if k in self.switches:
                self.switches[k].set(bool(val))
        for attr in ("srcGrid", "srcStars", "srcGalaxies"):
            vals = data.get(attr)
            if isinstance(vals, list):
                for var, val in zip(getattr(self, attr), vals):
                    var.set(val)
        for attr, widget in (("customCmds", getattr(self, "customCmds", None)),
                             ("objText", getattr(self, "objText", None))):
            if widget is not None and data.get(attr):
                # the inputMode trace may have already disabled the widget, and
                # Text.insert on a disabled widget is a silent no-op; force it
                # writable for the restore (__init__ re-syncs states afterwards)
                widget.configure(state="normal")
                widget.delete("1.0", "end")
                widget.insert("1.0", data[attr])


def main():
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()

"""A window for looking at XL-MS crosslink data as a matrix.
See docs/usage.md for how to use it.

The data comes from the same code the command line tool uses
(scripts/print_matrix.py). The matrix is drawn either through an embedded
pygame window, which is faster, or without pygame as a plain numpy/PIL image.
When you zoom out far enough, neighbouring rows and columns are merged into
blocks and their scores are combined."""
from __future__ import annotations

import json
import math
import os
import sys
import threading
import traceback
import tkinter as tk
import tkinter.filedialog as filedialog
import tkinter.font as tkfont
import tkinter.messagebox as messagebox
import tkinter.ttk as ttk
import urllib.request

import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageTk

try:
    import pygame as _pygame
except Exception:
    _pygame = None  # type: ignore[assignment]

try:
    import matplotlib
    _MPL = True
except Exception:
    _MPL = False

_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_DIR, "..", "scripts"))
sys.path.insert(0, os.path.join(_DIR, "..", "src"))

from print_matrix import (  # type: ignore[import]
    COL_W,
    ORDER_MODES,
    ROW_W,
    add_unlinked_residues,
    build_matrix,
    build_residue_matrix,
    protein_of,
    short_accession,
    sort_proteins,
    sort_residues,
)
from xlms.io import read_fasta  # type: ignore[import]


# -----------------------------------------------------------------------------
# Settings
#
# Cell sizes, how zooming out into blocks behaves, the mouse, fonts and
# colours.
# -----------------------------------------------------------------------------

CELL_PX_MIN = 4
CELL_PX_MAX = 32
CELL_PX_DEFAULT = 10
SCORE_THRESHOLD_PX = 18   # from this cell size on, the score is written into each cell

BLOCK_ZOOM_STEPS = 20     # at least this many slider steps below CELL_PX_MIN are used for blocks
BLOCK_ZOOM_GROWTH = 1.3   # each step makes blocks this much bigger (feels like a log zoom)
ZOOM_OUT_FIT_PX = 200     # big matrices get more block steps, until the whole matrix fits into this many pixels
OVERVIEW_PAD_PX = 12      # zoomed out into blocks: empty space around the matrix
_AGG_METHODS = {"Geometric Mean": "geomean", "Max Score": "max",
                "Root Sum of Squares": "rss", "Count": "count"}                    # dropdown -> key
_QUADRANTS = ["TT", "TD", "DT", "DD"]   # row target/decoy + column target/decoy; index = 2*row_decoy + col_decoy
_AGG_NAMES = {"geomean": "Geometric mean", "max": "Max", "rss": "Root sum of squares"}  # key -> info text

DRAG_THRESHOLD_PX = 4     # moving the mouse further than this with the button down is a drag, not a click
WHEEL_ZOOM_STEP = 0.15    # one notch of the mouse wheel changes the cell size by about this share

LABEL_FONT_SIZE = 9
LABEL_ANGLE = 45          # the column names go up to the right at this angle (degrees)
LABEL_GAP_PX = 3          # space between the column names and the line under them
SIDEBAR_W = 270          # width of the control panel on the left, in pixels
BG = "white"
_BG = (255, 255, 255)
_FG = (0, 0, 0)
_SEP = (153, 153, 153)
_HL = (0, 120, 215)  # outline of the selected cell

# The first of these fonts that exists is used (Windows, then Linux, then
# macOS). If none is found, Pillow's own font is used.
_MONO_FONTS = [
    "C:/Windows/Fonts/consola.ttf",
    "C:/Windows/Fonts/cour.ttf",
    "C:/Windows/Fonts/lucon.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",           # Debian, Ubuntu
    "/usr/share/fonts/dejavu-sans-mono-fonts/DejaVuSansMono.ttf",    # Fedora
    "/usr/share/fonts/TTF/DejaVuSansMono.ttf",                       # Arch
    "/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf",
    "/System/Library/Fonts/Menlo.ttc",
]


# -----------------------------------------------------------------------------
# Small drawing helpers
#
# Fonts, colormaps and the maths for turning a score into a colour, plus the
# different ways to combine the scores inside a block.
# -----------------------------------------------------------------------------

def _find_pil_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    for path in _MONO_FONTS:
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def _make_lut(cmap_name: str) -> np.ndarray:
    """A table of 256 RGB colours taken from a matplotlib colormap."""
    if not _MPL:
        return np.zeros((256, 3), dtype=np.uint8)
    cmap = matplotlib.colormaps[cmap_name]
    return (cmap(np.linspace(0, 1, 256))[:, :3] * 255).astype(np.uint8)


def _score_to_bg(score: float, lut: np.ndarray, score_sorted: np.ndarray) -> tuple[int, int, int]:
    """The colormap colour for a score. The colour follows the score's rank
    among all scores (quantile): the lowest score gets the first colour of
    the colormap, the median the middle one, the highest the last."""
    if len(score_sorted) > 0:
        idx = int(np.searchsorted(score_sorted, score) / len(score_sorted) * 255)
    else:
        idx = 128
    r, g, b = lut[max(0, min(255, idx))]
    return int(r), int(g), int(b)


def _auto_fg(bg_rgb: tuple[int, int, int]) -> tuple[int, int, int]:
    """Black or white text, whichever is easier to read on this background."""
    lum = 0.299 * bg_rgb[0] + 0.587 * bg_rgb[1] + 0.114 * bg_rgb[2]
    return (0, 0, 0) if lum > 128 else (255, 255, 255)


def _aggregate_groups(scores: np.ndarray, starts: list[int], method: str) -> np.ndarray:
    """Combine the scores of many groups at once, one value per group.
    Group k is scores[starts[k]:starts[k+1]]. The method is "geomean"
    (geometric mean), "max" (highest score), "rss" (square root of the sum
    of the squared scores) or "count" (number of scores)."""
    if len(starts) < 2:
        return np.zeros(0)
    first, sizes = starts[:-1], np.diff(starts)
    if method == "count":
        return sizes.astype(float)
    if method == "max":
        return np.maximum.reduceat(scores, first)
    if method == "rss":
        return np.sqrt(np.add.reduceat(scores * scores, first))
    # Geometric mean = exp(mean of the logarithms). Logarithms only exist for
    # positive numbers, so scores of 0 or below are left out. A group without
    # any positive score gets the normal mean instead, because dropping it
    # would make real crosslinks look like an empty cell.
    positive = scores > 0
    logs = np.log(np.where(positive, scores, 1.0))   # the 1.0 is never used, it only avoids log(0)
    n_positive = np.add.reduceat(positive.astype(int), first)
    log_sum = np.add.reduceat(np.where(positive, logs, 0.0), first)
    mean = np.add.reduceat(scores, first) / sizes
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(n_positive > 0, np.exp(log_sum / n_positive), mean)


# -----------------------------------------------------------------------------
# A slider with two handles
#
# Tkinter doesn't have one, so this draws its own on a canvas. The left
# handle sets the lower limit and the right handle the upper limit; the two
# can't pass each other.
# -----------------------------------------------------------------------------

class RangeSlider(tk.Canvas):
    """Pick a range [lo, hi] inside [min, max]. command(lo, hi) is called
    while you drag, release_command(lo, hi) when you let go."""
    PAD = 9       # space between the canvas edge and the end of the track
    HANDLE = 5    # half the width of a handle

    def __init__(self, master, command, release_command, **kw) -> None:
        kw.setdefault("height", 26)
        kw.setdefault("highlightthickness", 0)
        super().__init__(master, **kw)
        self._command = command
        self._release_command = release_command
        self._min, self._max = 0.0, 1.0
        self._lo, self._hi = 0.0, 1.0
        self._drag: str | None = None   # which handle is being dragged: "lo", "hi" or None
        self.bind("<Configure>", lambda _e: self._redraw())
        self.bind("<ButtonPress-1>", self._on_press)
        self.bind("<B1-Motion>", self._on_drag)
        self.bind("<ButtonRelease-1>", self._on_release)

    def set_bounds(self, lo: float, hi: float) -> None:
        """Set the full range and move both handles to its ends."""
        self._min, self._max = lo, hi
        self._lo, self._hi = lo, hi
        self._redraw()

    def values(self) -> tuple[float, float]:
        return self._lo, self._hi

    def _x(self, value: float) -> float:
        span = self._max - self._min
        frac = (value - self._min) / span if span > 0 else 0.0
        return self.PAD + frac * (self.winfo_width() - 2 * self.PAD)

    def _value(self, x: float) -> float:
        frac = (x - self.PAD) / max(1, self.winfo_width() - 2 * self.PAD)
        if frac <= 0:
            return self._min   # return the exact ends, so the lowest/highest score stays included
        if frac >= 1:
            return self._max
        return self._min + frac * (self._max - self._min)

    def _redraw(self) -> None:
        self.delete("all")
        y = self.winfo_height() // 2
        xl, xh = self._x(self._lo), self._x(self._hi)
        self.create_line(self.PAD, y, self.winfo_width() - self.PAD, y, fill="#bbbbbb", width=4)
        self.create_line(xl, y, xh, y, fill="#1a6fb5", width=4)
        for x in (xl, xh):
            self.create_rectangle(x - self.HANDLE, y - 8, x + self.HANDLE, y + 8,
                                  fill="white", outline="#1a6fb5", width=2)

    def _on_press(self, event: tk.Event) -> None:
        if self._max <= self._min:
            return
        d_lo, d_hi = abs(event.x - self._x(self._lo)), abs(event.x - self._x(self._hi))
        if d_lo != d_hi:
            self._drag = "lo" if d_lo < d_hi else "hi"
        elif self._hi >= self._max:      # both handles sit on top of each other at the right end
            self._drag = "lo"
        elif self._lo <= self._min:      # ... or at the left end
            self._drag = "hi"
        else:
            self._drag = "lo" if event.x < self._x(self._lo) else "hi"
        self._on_drag(event)

    def _on_drag(self, event: tk.Event) -> None:
        if self._drag is None:
            return
        value = self._value(event.x)
        if self._drag == "lo":
            self._lo = min(value, self._hi)
        else:
            self._hi = max(value, self._lo)
        self._redraw()
        self._command(self._lo, self._hi)

    def _on_release(self, _event: tk.Event) -> None:
        if self._drag is not None:
            self._drag = None
            self._release_command(self._lo, self._hi)


# -----------------------------------------------------------------------------
# The application window
# -----------------------------------------------------------------------------

class MatrixApp(tk.Tk):

    # -------------------------------------------------------------------------
    # Setting up the window
    #
    # All the state the app keeps, then the control panel, the matrix area
    # and (if available) the pygame renderer.
    # -------------------------------------------------------------------------

    def __init__(self, initial_csv: str | None = None) -> None:
        super().__init__()
        self.title("XL-MS Matrix Viewer")
        self.geometry("1200x800")
        self.minsize(700, 450)

        # The loaded data. Axis items are protein names, or residues in residue mode.
        self._level: str = "protein"
        self._proteins: list = []
        self._sparse: dict[tuple, float] = {}
        self._decoy: dict = {}
        self._section: dict | None = None
        self._seqs: dict[str, str] = {}
        self._split: int = 0              # index of the first decoy
        self._col_labels: list[str] = []
        self._row_labels: list[str] = []
        self._item_index: dict = {}
        # Every crosslink once more as arrays (row item, column item, score),
        # in the same order as self._sparse. numpy is much faster on these.
        self._link_rows: np.ndarray = np.zeros(0, dtype=np.int64)
        self._link_cols: np.ndarray = np.zeros(0, dtype=np.int64)
        self._link_scores: np.ndarray = np.zeros(0)

        # Blocks. All drawing code works on blocks instead of single items.
        # Without zooming out every block holds exactly one item.
        self._cell_px: int = CELL_PX_DEFAULT
        self._block_size: int = 1
        self._blocks: list[tuple[int, int]] = []    # (start, end) of each block
        # (start, end) of the sections in axis order, also cut at the
        # target/decoy line; blocks are built from these
        self._runs: list[tuple[int, int]] = []
        # Zoomed out into blocks the labels can't be read, so the headers
        # and row labels are hidden and the matrix gets the whole area
        self._overview: bool = False
        self._item_to_block: list[int] = []
        self._split_block: int = 0
        self._agg_method: str = "geomean"
        # Folded, the crosslinks A→B and B→A share one cell below the
        # diagonal and the matrix becomes a triangle
        self._folded: bool = False
        # The crosslinks sorted into cells: _group_scores holds the scores cell
        # by cell, and cell k owns the part _group_starts[k]:_group_starts[k+1].
        self._group_scores: np.ndarray = np.zeros(0)
        self._group_starts: list[int] = [0]

        # Every cell that has a value, row by row, ready for drawing
        self._cell_rows: np.ndarray = np.zeros(0, dtype=np.int64)    # row block
        self._cell_cols: np.ndarray = np.zeros(0, dtype=np.int64)    # column block
        self._cell_scores: np.ndarray = np.zeros(0)
        self._cell_rgb: np.ndarray = np.zeros((0, 3), dtype=np.uint8)

        # Colours
        self._cmap_name: str = "coolwarm"
        self._cmap_lut: np.ndarray | None = None   # None means greyscale
        self._min_score: float = 0.0
        self._max_score: float = 1.0
        self._score_sorted: np.ndarray = np.array([])
        # What the cell colours are measured against: all scores, or with
        # "Count" all the counts on screen. Sorted, for the quantile.
        self._colour_ref: np.ndarray = np.array([])

        # Score cutoff: cells with a score outside [lo, hi] are hidden
        self._cut_lo: float = float("-inf")
        self._cut_hi: float = float("inf")

        # Redrawing and the current selection
        self._pending: bool = False
        self._selected: tuple[int, int] | None = None   # (row block, column block)
        self._selected_items: tuple[int, int] | None = None   # the clicked row and column item, to find the cell again after zooming
        self._sel_token: int = 0                        # so an old UniProt answer can't overwrite a newer click
        self._uniprot_cache: dict[str, tuple[str, int]] = {}

        # Mouse: where a press started (mouse x, y and scroll offset x, y) and
        # whether it has turned into a drag
        self._drag_start: tuple[int, int, float, float] | None = None
        self._dragging: bool = False
        self._zoom_value: int = CELL_PX_DEFAULT    # the zoom slider value that is shown right now
        self._zoom_min: int = CELL_PX_MIN - BLOCK_ZOOM_STEPS   # lowest slider value; lower for big matrices

        # Fonts: the Tk font is only used for measuring, the PIL fonts for drawing
        self._lf: tkfont.Font | None = None
        self._cw: int = 8
        self._ch: int = 13
        self._pil_lf: ImageFont.FreeTypeFont | ImageFont.ImageFont | None = None
        self._pil_cf: ImageFont.FreeTypeFont | ImageFont.ImageFont | None = None
        # How far the slanted column names reach up and to the right of the
        # point where they start (the longest name counts)
        self._slant_h: int = 0
        self._slant_w: int = 0

        # Tk forgets images that nobody holds on to, so we keep them here
        self._mx_photo: ImageTk.PhotoImage | None = None
        self._ch_photo: ImageTk.PhotoImage | None = None
        self._rl_photo: ImageTk.PhotoImage | None = None
        # Writing text is slow, so every name and every score is only written once
        self._label_cache: dict[str, Image.Image] = {}
        self._slant_cache: dict[str, tuple[Image.Image, int, int]] = {}
        self._text_cache: dict[str, tuple[Image.Image, int, int]] = {}

        # pygame renderer (only used if pygame is installed)
        self._pg_host: tk.Frame | None = None       # the window pygame draws into
        self._pg_screen = None
        self._pg_size: tuple[int, int] = (0, 0)
        self._pg_font = None
        self._pg_cache: dict[tuple, object] = {}   # (colour, text) -> ready-made cell
        self._pg_frame = None                        # finished picture, shown by _present()
        self._last_render_error: str | None = None

        self._build_sidebar()
        self._build_matrix_area()

        self.update_idletasks()
        self._init_fonts()
        if _pygame is not None:
            self._init_pygame()

        if initial_csv:
            self._csv_var.set(initial_csv)
            self._load()

    # -------------------------------------------------------------------------
    # The control panel on the left
    #
    # Four groups of controls (data, display, score cutoff, selection) in a
    # panel that scrolls when the window is too short to show all of them.
    # -------------------------------------------------------------------------

    def _build_sidebar(self) -> None:
        outer = tk.Frame(self, width=SIDEBAR_W)
        outer.pack(side=tk.LEFT, fill=tk.Y, padx=(4, 0), pady=4)
        outer.pack_propagate(False)
        footer = tk.Frame(outer)
        footer.pack(side=tk.BOTTOM, fill=tk.X)
        side = self._scrollable(outer)
        wrap = SIDEBAR_W - 50   # line width for labels that can get long

        def group(title: str) -> tk.LabelFrame:
            frame = tk.LabelFrame(side, text=title, padx=6, pady=4)
            frame.pack(fill=tk.X, pady=(0, 6))
            frame.columnconfigure(1, weight=1)
            return frame

        def field(frame, row: int, label, widget, button=None) -> None:
            """Put a label, a widget and optionally a button on one row.
            The label can also be a StringVar, so it can change later."""
            if isinstance(label, tk.Variable):
                tk.Label(frame, textvariable=label).grid(row=row, column=0, sticky="w")
            else:
                tk.Label(frame, text=label).grid(row=row, column=0, sticky="w")
            widget.grid(row=row, column=1, columnspan=1 if button else 2, sticky="ew",
                        padx=(4, 0), pady=1)
            if button:
                button.grid(row=row, column=2, padx=(2, 0))

        # Data: which files to load and how
        data = group("Data")
        self._csv_var = tk.StringVar()
        field(data, 0, "CSV", tk.Entry(data, textvariable=self._csv_var, width=10),
              tk.Button(data, text="…", width=2, command=self._browse_csv))
        self._fasta_var = tk.StringVar()
        field(data, 1, "FASTA", tk.Entry(data, textvariable=self._fasta_var, width=10),
              tk.Button(data, text="…", width=2, command=self._browse_fasta))
        self._level_var = tk.StringVar(value="protein")
        field(data, 2, "Level", ttk.Combobox(data, textvariable=self._level_var, width=10,
                                             state="readonly", values=["protein", "residue"]))
        self._order_var = tk.StringVar(value="confidence")
        field(data, 3, "Order", ttk.Combobox(data, textvariable=self._order_var, width=10,
                                             state="readonly", values=ORDER_MODES))
        self._species_var = tk.StringVar(value="9606")
        field(data, 4, "Species", tk.Entry(data, textvariable=self._species_var, width=10))

        buttons = tk.Frame(data)
        buttons.grid(row=5, column=0, columnspan=3, sticky="ew", pady=(4, 0))
        self._load_btn = tk.Button(buttons, text="Load", command=self._load, width=8)
        self._load_btn.pack(side=tk.LEFT)
        tk.Button(buttons, text="Export PNG…", command=self._export).pack(side=tk.LEFT, padx=4)
        self._status_var = tk.StringVar(value="No data loaded.")
        tk.Label(data, textvariable=self._status_var, fg="gray", wraplength=wrap,
                 justify=tk.LEFT).grid(row=6, column=0, columnspan=3, sticky="w")

        # Display: colours, zoom and how blocks are combined
        display = group("Display")
        self._cmap_var = tk.StringVar(value="(none)")
        field(display, 0, "Colormap", ttk.Combobox(
            display, textvariable=self._cmap_var, width=10, state="readonly",
            values=["(none)", "coolwarm"] if _MPL else ["(none)"],
        ))
        self._cmap_var.trace_add("write", lambda *_: self._on_colormap_change())
        if not _MPL:
            tk.Label(display, text="(matplotlib not available)", fg="red").grid(
                row=1, column=0, columnspan=3, sticky="w")

        self._zoom_label_var = tk.StringVar(value=f"{CELL_PX_DEFAULT}px")
        self._zoom = tk.Scale(
            display, from_=CELL_PX_MIN - BLOCK_ZOOM_STEPS, to=CELL_PX_MAX,
            orient=tk.HORIZONTAL, showvalue=False, command=self._on_zoom,
        )
        self._zoom.set(CELL_PX_DEFAULT)
        field(display, 2, "Zoom", self._zoom)
        tk.Label(display, textvariable=self._zoom_label_var, fg="gray").grid(
            row=3, column=1, columnspan=2, sticky="w", padx=(4, 0))
        # Zoom to: fit the shown quadrants, or one quadrant, into the view
        zoom_to = tk.Frame(display)
        self._zoom_to_btns: dict[str, tk.Button] = {}
        for target in ["All", *_QUADRANTS]:
            btn = tk.Button(zoom_to, text=target, width=3, padx=0,
                            command=lambda t=target: self._zoom_to(t))
            btn.pack(side=tk.LEFT, padx=(0, 2))
            self._zoom_to_btns[target] = btn
        field(display, 4, "Zoom to", zoom_to)
        self._agg_method_var = tk.StringVar(value="Geometric Mean")
        field(display, 5, "Aggregate", ttk.Combobox(
            display, textvariable=self._agg_method_var, width=10, state="readonly",
            values=list(_AGG_METHODS),
        ))
        self._agg_method_var.trace_add("write", lambda *_: self._on_agg_method_change())
        self._fold_var = tk.BooleanVar(value=False)
        self._fold_var.trace_add("write", lambda *_: self._on_fold_change())
        field(display, 6, "Fold", tk.Checkbutton(
            display, text="Fold to triangle", variable=self._fold_var, indicatoron=False))

        # Quadrants: show or hide target/decoy links (row = Protein1, column =
        # Protein2), laid out like the quadrants in the matrix
        quadrants = group("Quadrants")
        self._quad_vars: dict[str, tk.BooleanVar] = {}
        self._quad_shown: dict[str, bool] = {name: True for name in _QUADRANTS}   # copy of the checkboxes, safe to read from the export thread
        for k, name in enumerate(_QUADRANTS):
            var = tk.BooleanVar(value=True)
            var.trace_add("write", lambda *_: self._on_quadrants_change())
            self._quad_vars[name] = var
            tk.Checkbutton(quadrants, text=name, variable=var).grid(
                row=k // 2, column=k % 2, sticky="w", padx=(0, 12))

        # Score cutoff: lower and upper limit on one slider
        cutoff = group("Score cutoff")
        self._cutoff_slider = RangeSlider(cutoff, command=self._on_cutoff,
                                          release_command=self._on_cutoff_release)
        self._cutoff_slider.grid(row=0, column=0, columnspan=3, sticky="ew")
        self._cutoff_var = tk.StringVar(value="–")
        tk.Label(cutoff, textvariable=self._cutoff_var, justify=tk.LEFT).grid(
            row=1, column=0, columnspan=2, sticky="w")
        tk.Button(cutoff, text="Reset", command=self._reset_cutoff).grid(row=1, column=2, sticky="e")

        # Selection: details about the clicked cell
        selection = group("Selection")
        self._sel_info_var = tk.StringVar(value="")
        tk.Label(selection, textvariable=self._sel_info_var, anchor="w", justify=tk.LEFT,
                 wraplength=wrap, fg="#1a6fb5", font=("TkFixedFont", 9)).grid(
            row=0, column=0, columnspan=3, sticky="w")

        # Footer: how many items are loaded and which renderer is used
        renderer = "pygame SDL2 (GPU)" if _pygame is not None else "numpy (CPU)"
        tk.Label(footer, text=f"renderer: {renderer}", fg="#aaaaaa").pack(side=tk.BOTTOM, anchor="w")
        self._info_var = tk.StringVar(value="")
        tk.Label(footer, textvariable=self._info_var, fg="gray").pack(side=tk.BOTTOM, anchor="w")

    def _scrollable(self, parent: tk.Frame) -> tk.Frame:
        """Return a frame that sits in a scrollable canvas, so the panel
        still works when the window is shorter than its content."""
        canvas = tk.Canvas(parent, highlightthickness=0)
        bar = tk.Scrollbar(parent, orient=tk.VERTICAL, command=canvas.yview)
        canvas.configure(yscrollcommand=bar.set)
        bar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        inner = tk.Frame(canvas)
        window = canvas.create_window(0, 0, window=inner, anchor="nw")
        inner.bind("<Configure>", lambda _e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(window, width=e.width))

        def is_inside(widget) -> bool:
            # Walk up the widget tree. Comparing the widget names as text doesn't
            # work, because ".!frame2" (the matrix area) starts with ".!frame" (this panel).
            while widget is not None:
                if widget is parent:
                    return True
                widget = widget.master
            return False

        def wheel(event: tk.Event) -> None:
            # this is bound for the whole app, so only react while the mouse is over the panel
            if is_inside(self.winfo_containing(event.x_root, event.y_root)):
                canvas.yview_scroll(-1 if event.delta > 0 else 1, "units")

        self.bind_all("<MouseWheel>", wheel, add="+")
        return inner

    def _browse_csv(self) -> None:
        p = filedialog.askopenfilename(filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if p:
            self._csv_var.set(p)

    def _browse_fasta(self) -> None:
        p = filedialog.askopenfilename(filetypes=[("FASTA", "*.fasta *.fa *.faa"), ("All", "*.*")])
        if p:
            self._fasta_var.set(p)

    # -------------------------------------------------------------------------
    # The matrix area
    #
    # Four canvases in a grid: column headers on top, row labels on the left,
    # the matrix itself, and an empty corner. Plus the fonts and sizes they use.
    # -------------------------------------------------------------------------

    def _build_matrix_area(self) -> None:
        outer = tk.Frame(self)
        outer.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=4, pady=4)
        outer.rowconfigure(1, weight=1)
        outer.columnconfigure(1, weight=1)

        # corner      | column headers
        # row labels  | matrix
        self._corner = tk.Frame(outer, bg=BG)
        self._corner.grid(row=0, column=0, sticky="nsew")
        self._ch_canvas = tk.Canvas(outer, bg=BG, highlightthickness=0)
        self._ch_canvas.grid(row=0, column=1, sticky="nsew")
        self._rl_canvas = tk.Canvas(outer, bg=BG, highlightthickness=0)
        self._rl_canvas.grid(row=1, column=0, sticky="nsew")
        self._mx = tk.Canvas(outer, bg=BG, highlightthickness=0)
        self._mx.grid(row=1, column=1, sticky="nsew")

        vbar = tk.Scrollbar(outer, orient=tk.VERTICAL, command=self._yscroll)
        vbar.grid(row=1, column=2, sticky="ns")
        hbar = tk.Scrollbar(outer, orient=tk.HORIZONTAL, command=self._xscroll)
        hbar.grid(row=2, column=1, sticky="ew")
        self._mx.configure(yscrollcommand=vbar.set, xscrollcommand=hbar.set)

        self._mx.bind("<Configure>", lambda _e: self._schedule())
        self._mx.bind("<MouseWheel>", self._wheel)
        self._mx.bind("<Button-4>", self._wheel)
        self._mx.bind("<Button-5>", self._wheel)
        self._mx.bind("<ButtonPress-1>", lambda e: self._press(e.x, e.y))
        self._mx.bind("<B1-Motion>", lambda e: self._drag_to(e.x, e.y))
        self._mx.bind("<ButtonRelease-1>", lambda e: self._release(e.x, e.y))
        self._mx.bind("<Up>",    lambda e: self._scroll(0, -1))
        self._mx.bind("<Down>",  lambda e: self._scroll(0,  1))
        self._mx.bind("<Left>",  lambda e: self._scroll(-1, 0))
        self._mx.bind("<Right>", lambda e: self._scroll( 1, 0))
        self._mx.focus_set()

    def _init_fonts(self) -> None:
        self._lf = tkfont.Font(family="TkFixedFont", size=LABEL_FONT_SIZE)
        self._cw = self._lf.measure("A")
        self._ch = self._lf.metrics("linespace")
        self._pil_lf = _find_pil_font(LABEL_FONT_SIZE + 2)
        self._pil_cf = _find_pil_font(max(6, CELL_PX_DEFAULT - 2))
        self._measure_col_labels(["W" * COL_W])   # until data is loaded
        self._sync_sizes()

    def _update_cell_font(self) -> None:
        self._pil_cf = _find_pil_font(max(6, self._cell_px - 2))
        self._text_cache.clear()
        if self._pg_screen is not None:
            self._pg_font = _pygame.font.Font(None, max(8, self._cell_px - 2))
            self._pg_cache.clear()

    @property
    def _label_w(self) -> int:
        return ROW_W * self._cw + 6

    @property
    def _label_base(self) -> int:
        """The height of the line under the column names, counted from the
        top of the header."""
        return (self._ch if self._section is not None else 0) + self._slant_h

    @property
    def _header_h(self) -> int:
        return self._label_base + self._ch // 2 + 2

    def _measure_col_labels(self, labels: list[str]) -> None:
        """Find out how much room the slanted column names need. The
        longest name needs the most."""
        widest = max(labels, key=self._pil_lf.getlength, default="")
        img, x, y = self._slanted_label(widest)
        self._slant_h = y + LABEL_GAP_PX
        self._slant_w = img.width - x

    def _sync_sizes(self) -> None:
        self._corner.configure(width=self._label_w, height=self._header_h)
        self._rl_canvas.configure(width=self._label_w)
        self._ch_canvas.configure(height=self._header_h)

    # -------------------------------------------------------------------------
    # The pygame renderer
    #
    # pygame draws into its own window, which lies on top of the matrix
    # canvas and covers it completely. A new picture is only drawn when
    # something changes. A short loop puts the finished picture back on
    # screen and picks up what the mouse does.
    # -------------------------------------------------------------------------

    def _init_pygame(self) -> None:
        # A frame without a background is a window that Tk never paints into.
        # Drawing straight into the canvas doesn't work well: Tk repaints the
        # canvas after every scroll step, and the picture flashes white.
        self._pg_host = tk.Frame(self._mx, bg="")
        self._pg_host.place(x=0, y=0, relwidth=1, relheight=1)
        for wheel in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            self._pg_host.bind(wheel, self._wheel)
        self.update_idletasks()
        os.environ["SDL_WINDOWID"] = str(self._pg_host.winfo_id())
        _pygame.display.quit()
        _pygame.display.init()
        _pygame.font.init()
        self._pg_sync_size()
        self._pg_font = _pygame.font.Font(None, max(8, self._cell_px - 2))
        self._pg_loop()

    def _pg_sync_size(self) -> bool:
        """Make the pygame display and the off-screen picture as big as the
        window pygame draws into. Returns True if the size changed (the
        picture is then blank)."""
        # Ask the pygame window itself, not the canvas: after a resize the
        # canvas has its new size a moment earlier, and pygame would then
        # still pick up the old size of its window.
        size = (max(1, self._pg_host.winfo_width()), max(1, self._pg_host.winfo_height()))
        if size == self._pg_size:
            return False
        self._pg_screen = _pygame.display.set_mode(size, 0, 32)
        self._pg_size = size
        self._pg_frame = _pygame.Surface(size)
        self._pg_frame.fill(_BG)
        return True

    def _pg_loop(self) -> None:
        """Runs about every 16 ms. It only handles the mouse and shows the
        last picture again; drawing a new one is _redraw's job."""
        try:
            self._poll_pygame_events()
            self._present()
        except Exception as exc:
            self._report_render_error(exc)
        self.after(16, self._pg_loop)

    def _poll_pygame_events(self) -> None:
        # inside the pygame window, the mouse buttons go to pygame instead of Tk
        for ev in _pygame.event.get():
            if ev.type == _pygame.MOUSEBUTTONDOWN and ev.button == 1:
                self._press(*ev.pos)
            elif ev.type == _pygame.MOUSEMOTION:
                if ev.buttons[0]:
                    self._drag_to(*ev.pos)
                else:
                    # the button was let go somewhere we couldn't see it
                    self._release(*ev.pos)
            elif ev.type == _pygame.MOUSEBUTTONUP and ev.button == 1:
                self._release(*ev.pos)

    def _present(self) -> None:
        """Copy the finished picture into the pygame window (takes about 0.4 ms)."""
        if self._pg_sync_size():
            self._schedule()   # the canvas changed size, so draw a picture that fits
        self._pg_screen.blit(self._pg_frame, (0, 0))
        _pygame.display.flip()

    def _report_render_error(self, exc: Exception) -> None:
        """Show a drawing error once instead of hiding it or repeating it
        every frame."""
        message = f"{type(exc).__name__}: {exc}"
        if message != self._last_render_error:
            self._last_render_error = message
            traceback.print_exception(exc)
            self._status_var.set(f"render error: {message}")

    # -------------------------------------------------------------------------
    # Scrolling
    #
    # The matrix, column headers and row labels scroll together. Redraws are
    # bundled, so many scroll events in a row only cause one redraw.
    # -------------------------------------------------------------------------

    @property
    def _pad(self) -> int:
        return OVERVIEW_PAD_PX if self._overview else 0

    def _update_scrollregion(self) -> None:
        size, pad = len(self._blocks) * self._cell_px, self._pad
        self._mx.configure(scrollregion=(-pad, -pad, size + pad, size + pad))
        self._ch_canvas.configure(scrollregion=(-pad, 0, size + pad, self._header_h))
        self._rl_canvas.configure(scrollregion=(0, -pad, self._label_w, size + pad))

    def _yscroll(self, *args) -> None:
        self._mx.yview(*args)
        self._rl_canvas.yview(*args)
        self._schedule()

    def _xscroll(self, *args) -> None:
        self._mx.xview(*args)
        self._ch_canvas.xview(*args)
        self._schedule()

    def _scroll(self, dx: int, dy: int) -> None:
        if dy:
            self._mx.yview_scroll(dy, "units")
            self._rl_canvas.yview_scroll(dy, "units")
        if dx:
            self._mx.xview_scroll(dx, "units")
            self._ch_canvas.xview_scroll(dx, "units")
        self._schedule()

    def _scroll_to(self, ox: float, oy: float) -> None:
        """Move the view so that pixel (ox, oy) of the whole matrix sits in
        the top-left corner. Tk keeps the view inside the matrix by itself."""
        size, pad = len(self._blocks) * self._cell_px, self._pad
        if size:
            fx, fy = (ox + pad) / (size + 2 * pad), (oy + pad) / (size + 2 * pad)
            self._mx.xview_moveto(fx)
            self._ch_canvas.xview_moveto(fx)
            self._mx.yview_moveto(fy)
            self._rl_canvas.yview_moveto(fy)
        self._schedule()

    def _visible_blocks(self) -> tuple[range, range]:
        """Which row blocks and column blocks are currently on screen."""
        n, cp, c = len(self._blocks), self._cell_px, self._mx
        rows = range(max(0, int(c.canvasy(0) // cp)), min(n, int(c.canvasy(c.winfo_height()) // cp) + 1))
        cols = range(max(0, int(c.canvasx(0) // cp)), min(n, int(c.canvasx(c.winfo_width()) // cp) + 1))
        return rows, cols

    def _schedule(self) -> None:
        """Ask for a redraw. Several requests before the next idle moment
        are combined into one."""
        if not self._pending:
            self._pending = True
            self.after(0, self._redraw)

    # -------------------------------------------------------------------------
    # Moving around with the mouse
    #
    # Like on a map: the wheel zooms in and out around the mouse pointer, and
    # holding the left button down drags the matrix. A press that is let go
    # without moving is a normal click and selects the cell.
    # -------------------------------------------------------------------------

    def _wheel(self, event: tk.Event) -> None:
        v = self._zoom_value
        # Above the smallest cell size a notch changes the size by a share of
        # itself, so it feels the same at every size. Below it, each slider
        # step already merges 1.3 times as many items.
        step = max(1, round(v * WHEEL_ZOOM_STEP))
        if event.num == 5 or event.delta < 0:   # zoom out
            v = v - 1 if v <= CELL_PX_MIN else max(CELL_PX_MIN, v - step)
        else:
            v = v + 1 if v < CELL_PX_MIN else v + step
        v = min(max(v, self._zoom_min), CELL_PX_MAX)
        # ask Tk where the pointer is; that also works above the pygame window
        anchor = (self.winfo_pointerx() - self._mx.winfo_rootx(),
                  self.winfo_pointery() - self._mx.winfo_rooty())
        self._on_zoom(str(v), anchor)
        self._zoom.set(v)

    def _press(self, x: int, y: int) -> None:
        self._mx.focus_set()   # so the arrow keys work again after typing in the control panel
        self._drag_start = (x, y, self._mx.canvasx(0), self._mx.canvasy(0))
        self._dragging = False

    def _drag_to(self, x: int, y: int) -> None:
        if self._drag_start is None:
            return
        x0, y0, ox, oy = self._drag_start
        if not self._dragging and max(abs(x - x0), abs(y - y0)) > DRAG_THRESHOLD_PX:
            self._dragging = True
            self._show_drag_cursor(True)
        if self._dragging:
            self._scroll_to(ox - (x - x0), oy - (y - y0))

    def _release(self, x: int, y: int) -> None:
        if self._drag_start is None:
            return
        if self._dragging:
            self._show_drag_cursor(False)
        else:
            self._select_at(x, y)
        self._drag_start = None
        self._dragging = False

    def _show_drag_cursor(self, dragging: bool) -> None:
        """Show the four-arrow cursor while the matrix is being dragged."""
        if self._pg_screen is not None:
            # Changing the cursor of the Tk canvas would make Tk paint the
            # canvas again, right over the pygame picture, which shows as a
            # flash. So the pygame window gets the cursor instead.
            _pygame.mouse.set_cursor(_pygame.SYSTEM_CURSOR_SIZEALL if dragging else _pygame.SYSTEM_CURSOR_ARROW)
        else:
            self._mx.configure(cursor="fleur" if dragging else "")

    # -------------------------------------------------------------------------
    # Loading data
    #
    # The CSV is read in a background thread so the window doesn't freeze;
    # the result is handed back to the Tk thread when it's ready.
    # -------------------------------------------------------------------------

    def _load(self) -> None:
        csv = self._csv_var.get().strip()
        level = self._level_var.get()
        order = self._order_var.get()
        fasta = self._fasta_var.get().strip() or None
        if not csv:
            messagebox.showerror("Error", "Please select a CSV file.")
            return
        if order == "sequence" and not fasta:
            messagebox.showerror("Error", "Order 'sequence' requires a FASTA file.")
            return
        try:
            species = int(self._species_var.get())
        except ValueError:
            messagebox.showerror("Error", "Species must be an integer NCBI taxid.")
            return

        self._load_btn.configure(state=tk.DISABLED)
        self._status_var.set("Loading…")

        def worker() -> None:
            try:
                seqs = read_fasta(fasta) if fasta else None
                if level == "residue":
                    items, sparse, decoy = build_residue_matrix(csv, n=None)   # always all residues
                    if seqs:
                        # with a FASTA, every residue of the loaded proteins is shown,
                        # not only the crosslinked ones
                        items, decoy = add_unlinked_residues(items, decoy, seqs)
                    items, section = sort_residues(items, decoy, order, seqs, species)
                else:
                    items, sparse, decoy = build_matrix(csv, n=None)   # always all proteins
                    items, section = sort_proteins(items, decoy, order, seqs, species)
                self.after(0, lambda: self._on_loaded(items, sparse, decoy, section, seqs, level))
            except Exception as exc:
                msg = str(exc)
                self.after(0, lambda: self._on_error(msg))

        threading.Thread(target=worker, daemon=True).start()

    def _on_loaded(self, items: list, sparse: dict, decoy: dict, section: dict | None,
                   seqs: dict[str, str] | None, level: str) -> None:
        self._level = level
        self._proteins = items
        self._sparse = sparse
        self._decoy = decoy
        self._section = section
        self._seqs = seqs or {}
        if level == "residue":
            self._col_labels = [str(r.pos)[:COL_W] for r in items]
            self._row_labels = [f"{short_accession(r.protein)}:{r.pos}"[:ROW_W] for r in items]
        else:
            self._col_labels = [p[:COL_W] for p in items]
            self._row_labels = [p[:ROW_W] for p in items]
        self._split = next((i for i, p in enumerate(items) if decoy.get(p, False)), len(items))
        scores = list(sparse.values())
        self._min_score = min(scores) if scores else 0.0
        self._max_score = max(scores) if scores else 1.0
        self._score_sorted = np.sort(scores) if scores else np.array([])
        self._item_index = {p: i for i, p in enumerate(items)}
        index = self._item_index
        self._link_rows = np.array([index[a] for a, _ in sparse], dtype=np.int64)
        self._link_cols = np.array([index[b] for _, b in sparse], dtype=np.int64)
        self._link_scores = np.array(scores, dtype=float)
        self._label_cache.clear()
        self._slant_cache.clear()
        self._measure_col_labels(self._col_labels)
        self._clear_selection()   # it pointed into the old data
        self._cut_lo, self._cut_hi = self._min_score, self._max_score
        self._cutoff_slider.set_bounds(self._cut_lo, self._cut_hi)
        self._update_cutoff_label()
        self._runs = self._find_runs()
        self._recompute_blocks()
        self._zoom_min = self._find_zoom_min()
        self._zoom.configure(from_=self._zoom_min)
        if self._zoom_value < self._zoom_min:   # the new matrix needs fewer block steps
            self._on_zoom(str(self._zoom_min))
            self._zoom.set(self._zoom_min)
        for target, btn in self._zoom_to_btns.items():
            btn.configure(state=tk.NORMAL if self._zoom_region(target) else tk.DISABLED)

        unit = "residues" if level == "residue" else "proteins"
        self._status_var.set(
            f"Loaded — {self._split} targets · {len(items) - self._split} decoys · {len(sparse)} links"
        )
        self._info_var.set(f"{len(items)} {unit}")
        self._load_btn.configure(state=tk.NORMAL)
        self._sync_sizes()
        self._update_scrollregion()
        for view in (self._mx.xview_moveto, self._mx.yview_moveto,
                     self._ch_canvas.xview_moveto, self._rl_canvas.yview_moveto):
            view(0)
        self._schedule()

    def _on_error(self, msg: str) -> None:
        self._load_btn.configure(state=tk.NORMAL)
        self._status_var.set(f"Error: {msg}")
        messagebox.showerror("Load error", msg)

    # -------------------------------------------------------------------------
    # Colours
    #
    # Without a colormap, low scores are light grey and high scores dark grey.
    # With one, the colour comes from the colormap and follows the rank of
    # the score among all scores (quantile).
    # -------------------------------------------------------------------------

    def _on_colormap_change(self) -> None:
        name = self._cmap_var.get()
        if name == "(none)":
            self._cmap_lut = None
        else:
            self._cmap_name = name
            self._cmap_lut = _make_lut(name)
        self._cell_rgb = self._colours_for(self._cell_scores)
        self._pg_cache.clear()
        self._schedule()

    def _cell_colour(self, score: float, ref: np.ndarray | None = None) -> tuple[int, int, int]:
        """The colour of one cell. ref are the sorted values the colour is
        measured against; without it, all scores."""
        if ref is None:
            ref, (mn, mx) = self._score_sorted, (self._min_score, self._max_score)
        else:
            mn, mx = (ref[0], ref[-1]) if len(ref) else (0.0, 1.0)
        if self._cmap_lut is not None:
            return _score_to_bg(score, self._cmap_lut, ref)
        grey = int(200 - (score - mn) / (mx - mn) * 160) if mx > mn else 120
        grey = max(40, min(200, grey))
        return grey, grey, grey

    def _colours_for(self, scores: np.ndarray) -> np.ndarray:
        """The colours of the cells, one RGB row per value. For scores this
        is the same as _cell_colour. Counts and root sums of squares get
        their own scale, from the smallest to the largest value on screen
        (see _reduce_block_scores)."""
        ref = self._colour_ref
        mn, mx = (ref[0], ref[-1]) if len(ref) else (0.0, 1.0)
        if self._cmap_lut is None:
            if mx > mn:
                grey = (200 - (scores - mn) / (mx - mn) * 160).astype(np.int64)
            else:
                grey = np.full(len(scores), 120, dtype=np.int64)
            grey = np.clip(grey, 40, 200).astype(np.uint8)
            return np.stack([grey, grey, grey], axis=1)
        ranked = ref
        if len(ranked) > 0:
            idx = (np.searchsorted(ranked, scores) / len(ranked) * 255).astype(np.int64)
        else:
            idx = np.full(len(scores), 128, dtype=np.int64)
        return self._cmap_lut[np.clip(idx, 0, 255)]

    @staticmethod
    def _cell_text(score: float, cell_px: int) -> str | None:
        """The score written inside a cell, but only once cells are big enough."""
        return f"{score:.2f}" if cell_px >= SCORE_THRESHOLD_PX else None

    # -------------------------------------------------------------------------
    # Score cutoff
    #
    # Cells whose score is outside the chosen range are hidden. The colours
    # don't change, so a score always has the same colour.
    # -------------------------------------------------------------------------

    def _in_cutoff(self, score: float) -> bool:
        return self._cut_lo <= score <= self._cut_hi

    def _quadrant_of(self, row: int, col: int) -> str:
        """TT, TD, DT or DD for the cell of item `row` and item `col`."""
        return _QUADRANTS[2 * (row >= self._split) + (col >= self._split)]

    def _quadrant_shown(self, row: int, col: int) -> bool:
        return self._quad_shown[self._quadrant_of(row, col)]

    def _shown_links(self) -> np.ndarray:
        """Which crosslinks are visible: inside the cutoff and in a quadrant
        that is switched on. One bool per link."""
        scores = self._link_scores
        shown = np.array([self._quad_shown[q] for q in _QUADRANTS])
        quadrant = 2 * (self._link_rows >= self._split) + (self._link_cols >= self._split)
        return (scores >= self._cut_lo) & (scores <= self._cut_hi) & shown[quadrant]

    def _on_quadrants_change(self) -> None:
        self._quad_shown = {name: var.get() for name, var in self._quad_vars.items()}
        self._group_links()
        self._reduce_block_scores()
        self._update_cutoff_label()
        self._pg_cache.clear()
        self._schedule()
        if self._selected is not None:
            self._select_block(*self._selected)   # update the text for the selected cell

    def _on_cutoff(self, lo: float, hi: float) -> None:
        self._cut_lo, self._cut_hi = lo, hi
        self._group_links()
        self._reduce_block_scores()
        self._update_cutoff_label()
        self._schedule()

    def _on_cutoff_release(self, lo: float, hi: float) -> None:
        if self._selected is not None:
            self._select_block(*self._selected)   # update the text for the selected cell

    def _reset_cutoff(self) -> None:
        if not self._sparse:
            return
        self._cutoff_slider.set_bounds(self._min_score, self._max_score)
        self._on_cutoff(self._min_score, self._max_score)
        self._on_cutoff_release(self._min_score, self._max_score)

    def _update_cutoff_label(self) -> None:
        shown = int(self._shown_links().sum())
        self._cutoff_var.set(f"{self._cut_lo:.4g} – {self._cut_hi:.4g}\n"
                             f"{shown} of {len(self._link_scores)} links shown")

    # -------------------------------------------------------------------------
    # Zooming out into blocks
    #
    # Below the smallest cell size, zooming out further merges neighbouring
    # items into blocks. A block never mixes targets with decoys; small
    # sections are packed into one block whole, big ones are cut. A block's
    # value is combined as chosen under
    # Aggregate. Big matrices get extra zoom steps, so the whole matrix can
    # always be zoomed out to fit on screen. "Zoom to" picks the zoom that
    # fits a quadrant (or all shown quadrants) into the view.
    # -------------------------------------------------------------------------

    def _on_zoom(self, value: str, anchor: tuple[float, float] | None = None) -> None:
        """Called by the zoom slider and by the mouse wheel. The spot under
        `anchor` (in canvas pixels; the middle of the view if not given)
        stays where it is while everything around it grows or shrinks."""
        v = int(float(value))
        if v == self._zoom_value:
            return   # the slider only reports what the wheel has already done
        self._zoom_value = v
        ax, ay = anchor or (self._mx.winfo_width() / 2, self._mx.winfo_height() / 2)
        item_x = self._item_at(self._mx.canvasx(0) + ax)
        item_y = self._item_at(self._mx.canvasy(0) + ay)

        self._cell_px = max(v, CELL_PX_MIN)   # cells don't get smaller; items get merged instead
        block_size = self._block_size_for(v)
        self._zoom_label_var.set(f"{v}px" if block_size == 1 else f"≤{block_size} → {CELL_PX_MIN}px")
        # hiding or showing the labels moves the matrix on screen; keep the
        # anchor on the same spot of the screen
        rx, ry = self._mx.winfo_rootx(), self._mx.winfo_rooty()
        self._set_overview(block_size > 1)
        ax += rx - self._mx.winfo_rootx()
        ay += ry - self._mx.winfo_rooty()
        self._update_cell_font()
        if block_size != self._block_size:
            self._block_size = block_size
            self._recompute_blocks()
            self._pg_cache.clear()
            if self._selected_items is not None:
                # the selected items now sit in other blocks
                self._select_block(*(self._item_to_block[i] for i in self._selected_items))
        self._update_scrollregion()
        self._scroll_to(self._pixel_of(item_x) - ax, self._pixel_of(item_y) - ay)

    def _set_overview(self, on: bool) -> None:
        """Hide the column headers and row labels while zoomed out into
        blocks, and show them again when zooming back in."""
        if on == self._overview:
            return
        self._overview = on
        for widget in (self._corner, self._ch_canvas, self._rl_canvas):
            if on:
                widget.grid_remove()
            else:
                widget.grid()
        self.update_idletasks()   # so the matrix canvas has its new size right away

    @staticmethod
    def _block_size_for(v: int) -> int:
        """How many items one block holds at zoom slider value v."""
        if v >= CELL_PX_MIN:
            return 1
        return max(2, round(BLOCK_ZOOM_GROWTH ** (CELL_PX_MIN - v)))

    def _blocks_for(self, bs: int) -> list[tuple[int, int]]:
        """Cut the items into blocks of at most bs items. A block never
        mixes targets with decoys. A section (e.g. one protein in residue
        mode) bigger than a block is cut into pieces of bs items; smaller
        sections are packed together whole, so that matrices with many small
        sections can still be zoomed out far."""
        blocks: list[tuple[int, int]] = []
        cur: tuple[int, int] | None = None   # the block being filled with whole sections
        for s, e in self._runs:
            if cur is not None and (s == self._split or e - cur[0] > bs or e - s >= bs):
                blocks.append(cur)
                cur = None
            if e - s >= bs:
                blocks += [(i, min(i + bs, e)) for i in range(s, e, bs)]
            else:
                cur = (cur[0] if cur else s, e)
        if cur is not None:
            blocks.append(cur)
        return blocks

    def _extent_px(self, v: int, start: int, end: int) -> int:
        """How many pixels the items start:end take up at zoom slider value
        v. start and end must lie on block borders (0, the target/decoy line
        or the number of items)."""
        blocks = sum(1 for s, e in self._blocks_for(self._block_size_for(v)) if start <= s and e <= end)
        return blocks * max(v, CELL_PX_MIN)

    def _find_runs(self) -> list[tuple[int, int]]:
        """The sections in axis order, also cut at the target/decoy line.
        Blocks are built from these (see _blocks_for)."""
        items, section, split = self._proteins, self._section, self._split
        runs, start = [], 0
        for j in range(1, len(items) + 1):
            if (j == len(items) or j == split
                    or (section is not None
                        and section.get(items[j], "") != section.get(items[start], ""))):
                runs.append((start, j))
                start = j
        return runs

    def _find_zoom_min(self) -> int:
        """The lowest zoom slider value: at least BLOCK_ZOOM_STEPS below
        CELL_PX_MIN, and further for big matrices, until the whole matrix
        fits into ZOOM_OUT_FIT_PX or blocks can't get any bigger (one block
        for the targets and one for the decoys)."""
        n = len(self._proteins)
        longest = max(self._split, n - self._split, 1)
        v = CELL_PX_MIN - BLOCK_ZOOM_STEPS
        while self._extent_px(v, 0, n) > ZOOM_OUT_FIT_PX and self._block_size_for(v) < longest:
            v -= 1
        return v

    def _zoom_region(self, target: str) -> tuple[int, int, int, int] | None:
        """The items (row start, row end, column start, column end) that
        `target` covers: one quadrant, or for "All" the box around all
        quadrants that are switched on. None if it holds no items."""
        n, split = len(self._proteins), self._split
        sides = {"T": (0, split), "D": (split, n)}
        names = [target] if target != "All" else [q for q in _QUADRANTS if self._quad_shown[q]]
        boxes = [(*sides[q[0]], *sides[q[1]]) for q in names]
        boxes = [b for b in boxes if b[0] < b[1] and b[2] < b[3]]
        if not boxes:
            return (0, n, 0, n) if target == "All" and n else None
        return (min(b[0] for b in boxes), max(b[1] for b in boxes),
                min(b[2] for b in boxes), max(b[3] for b in boxes))

    def _zoom_to(self, target: str) -> None:
        """Zoom so that the region fills the view as far as possible, and
        center it. A hidden quadrant is switched on first."""
        region = self._zoom_region(target)
        if region is None:
            return
        if target != "All" and not self._quad_shown[target]:
            self._quad_vars[target].set(True)   # the trace redraws with it shown
        r0, r1, c0, c1 = region
        # the whole area, with the labels if they are shown right now
        full_w, full_h = self._mx.winfo_width(), self._mx.winfo_height()
        if not self._overview:
            full_w, full_h = full_w + self._label_w, full_h + self._header_h
        # with cells the labels take their space; zoomed out into blocks they
        # are hidden and only the padding around the matrix is left free
        with_labels = (full_w - self._label_w, full_h - self._header_h)
        overview = (full_w - 2 * OVERVIEW_PAD_PX, full_h - 2 * OVERVIEW_PAD_PX)

        def fits(v: int) -> bool:
            fw, fh = with_labels if v >= CELL_PX_MIN else overview
            return self._extent_px(v, c0, c1) <= fw and self._extent_px(v, r0, r1) <= fh

        # the biggest zoom at which the region still fits
        v = next((v for v in range(CELL_PX_MAX, self._zoom_min - 1, -1) if fits(v)), self._zoom_min)
        self._on_zoom(str(v))
        self._zoom.set(v)
        w, h = self._mx.winfo_width(), self._mx.winfo_height()
        x0, x1 = self._pixel_of(c0), self._pixel_of(c1)
        y0, y1 = self._pixel_of(r0), self._pixel_of(r1)
        self._scroll_to(x0 - (w - (x1 - x0)) / 2, y0 - (h - (y1 - y0)) / 2)

    def _item_at(self, pixel: float) -> float:
        """Which item lies at this pixel of the whole matrix. The answer has
        decimals (2.5 = the middle of item 2), so that a spot can be found
        again exactly after zooming."""
        n = len(self._blocks)
        if n == 0:
            return 0.0
        pos = min(max(pixel / self._cell_px, 0), n)
        block = min(int(pos), n - 1)
        start, end = self._blocks[block]
        return start + (pos - block) * (end - start)

    def _pixel_of(self, item: float) -> float:
        """The opposite of _item_at."""
        if not self._blocks:
            return 0.0
        block = self._item_to_block[min(int(item), len(self._proteins) - 1)]
        start, end = self._blocks[block]
        return (block + (item - start) / (end - start)) * self._cell_px

    def _on_agg_method_change(self) -> None:
        self._agg_method = _AGG_METHODS.get(self._agg_method_var.get(), "geomean")
        self._reduce_block_scores()
        self._pg_cache.clear()
        self._schedule()
        if self._selected is not None:
            self._select_block(*self._selected)   # update the text for the selected cell

    def _on_fold_change(self) -> None:
        """Fold the matrix along its diagonal, or unfold it again. A selected
        cell above the diagonal moves to its mirror cell below it."""
        self._folded = self._fold_var.get()
        self._group_links()
        self._reduce_block_scores()
        self._pg_cache.clear()
        if self._selected is not None:
            rb, cb = self._selected
            if self._folded and rb < cb:
                rb, cb = cb, rb
                self._selected_items = self._selected_items and self._selected_items[::-1]
            self._select_block(rb, cb)
        self._schedule()

    def _recompute_blocks(self) -> None:
        n, split, bs = len(self._proteins), self._split, self._block_size
        blocks = self._blocks_for(bs)

        self._blocks = blocks
        self._item_to_block = [bi for bi, (s, e) in enumerate(blocks) for _ in range(s, e)]
        self._split_block = self._item_to_block[split] if 0 <= split < n else len(blocks)
        self._group_links()
        self._reduce_block_scores()

    def _group_links(self) -> None:
        """Sort every visible crosslink into its cell (row block of Protein1,
        column block of Protein2). Links hidden by the cutoff or by a
        quadrant checkbox are left out, so the block values only describe
        what you can see. A block never crosses the target/decoy line, so
        each cell lies in one quadrant. Folded, a crosslink above the
        diagonal goes to its mirror cell below it (TD then joins DT)."""
        scores = self._link_scores
        keep = self._shown_links()
        to_block = np.array(self._item_to_block, dtype=np.int64)
        n = max(1, len(self._blocks))
        rb, cb = to_block[self._link_rows[keep]], to_block[self._link_cols[keep]]
        if self._folded:
            rb, cb = np.maximum(rb, cb), np.minimum(rb, cb)
        cell = rb * n + cb
        # a stable sort keeps the crosslinks of one cell in their original order
        order = np.argsort(cell, kind="stable")
        cell = cell[order]
        first = np.flatnonzero(np.diff(cell, prepend=-1))   # where the next cell begins
        self._cell_rows, self._cell_cols = np.divmod(cell[first], n)
        self._group_scores = scores[keep][order]
        self._group_starts = first.tolist() + [len(cell)]

    def _reduce_block_scores(self) -> None:
        """Turn the crosslinks of every cell into the one value it shows: a
        score, or with "Count" the number of crosslinks."""
        scores, starts = self._group_scores, self._group_starts
        # folded, even a single cell can hold two crosslinks (both directions)
        combined = self._block_size > 1 or self._folded
        if combined:
            values = _aggregate_groups(scores, starts, self._agg_method)
        else:
            values = scores   # one item per block, so a cell holds exactly one crosslink
        self._cell_scores = np.array(values, dtype=float)
        if combined and self._agg_method in ("count", "rss"):
            # counts aren't scores, and a root sum of squares grows with the
            # number of crosslinks past every single score, so both are
            # coloured by their rank among the values on screen
            self._colour_ref = np.sort(self._cell_scores)
        else:
            self._colour_ref = self._score_sorted
        self._cell_rgb = self._colours_for(self._cell_scores)

    def _cell_at(self, rb: int, cb: int) -> tuple[float, int] | None:
        """The score shown in the cell at (row block, column block) and how
        many crosslinks went into it, or None if the cell is empty. Only
        crosslinks with Protein1 in the row and Protein2 in the column count."""
        n = max(1, len(self._blocks))
        # the cells are sorted row by row, so the right one can be looked up
        k = int(np.searchsorted(self._cell_rows * n + self._cell_cols, rb * n + cb))
        if k == len(self._cell_rows) or (self._cell_rows[k], self._cell_cols[k]) != (rb, cb):
            return None
        return float(self._cell_scores[k]), self._group_starts[k + 1] - self._group_starts[k]

    # -------------------------------------------------------------------------
    # Clicking on a cell
    #
    # Shows both proteins, their scores in both directions, and fetches the
    # full protein names from UniProt in the background.
    # -------------------------------------------------------------------------

    def _select_at(self, x: int, y: int) -> None:
        """Select the block under the mouse, or clear the selection if the
        click was outside the matrix."""
        if not self._blocks:
            return
        cp, n = self._cell_px, len(self._blocks)
        cb = (x + int(self._mx.canvasx(0))) // cp
        rb = (y + int(self._mx.canvasy(0))) // cp
        if self._folded and rb < cb:
            rb, cb = cb, rb   # the empty half of a folded matrix: take the mirror cell
        if 0 <= rb < n and 0 <= cb < n:
            self._selected_items = (self._blocks[rb][0], self._blocks[cb][0])
            self._select_block(rb, cb)
        else:
            self._clear_selection()
        self._schedule()

    def _clear_selection(self) -> None:
        self._selected = None
        self._selected_items = None
        self._sel_token += 1
        self._sel_info_var.set("")

    def _select_block(self, rb: int, cb: int) -> None:
        self._selected = (rb, cb)
        self._sel_token += 1
        token = self._sel_token
        row_items = self._proteins[slice(*self._blocks[rb])]
        col_items = self._proteins[slice(*self._blocks[cb])]
        score_line = self._score_line(rb, cb)

        def show(fetched: tuple | None = None) -> None:
            sides = [self._describe_side(items, fetched[k] if fetched else None, fetched is not None)
                     for k, items in enumerate((row_items, col_items))]
            self._sel_info_var.set("\n".join(sides + [score_line]))

        show()
        # UniProt only makes sense when each side is a single protein
        if len({protein_of(it) for it in row_items}) == 1 and len({protein_of(it) for it in col_items}) == 1:
            def fetch() -> None:
                fetched = tuple(self._fetch_uniprot(short_accession(protein_of(items[0])))
                                for items in (row_items, col_items))
                self.after(0, lambda: self._sel_token == token and show(fetched))
            threading.Thread(target=fetch, daemon=True).start()

    def _describe_side(self, items: list, uniprot: tuple[str, int] | None, fetched: bool) -> str:
        """One line of text about the row or column side of the selection."""
        decoys = {self._decoy.get(it, False) for it in items}
        tag = "target" if decoys == {False} else "decoy" if decoys == {True} else "mixed"
        sections = {self._section.get(it, "") for it in items} if self._section else set()
        section = next(iter(sections)) if len(sections) == 1 else ""
        sec = f"  [{section}]" if section else ""

        if uniprot:
            name, length = uniprot
            acc = short_accession(protein_of(items[0]))
            return f"{name}  ({length} aa, {tag}){sec}{self._pos_range(items)}  [{acc}]"
        if len(items) > 1:
            return f"{len(items)} items  (?, {tag}){sec}" if fetched else f"{len(items)} items ({tag}){sec}"
        protein = protein_of(items[0])
        if fetched:
            size = "?"
        elif protein in self._seqs:
            size = f"{len(self._seqs[protein])} aa"
        else:
            size = "fetching…"
        return f"{items[0]}  ({size}, {tag}){sec}"

    def _pos_range(self, items: list) -> str:
        if self._level != "residue":
            return ""
        postions = sorted(it.pos for it in items)
        if len(postions) == 1:
            return f"  pos {postions[0]}"
        return f"  pos {postions[0]}-{postions[-1]}"

    def _score_line(self, rb: int, cb: int) -> str:
        """The scores of the clicked cell (row→col) and of its mirror cell
        (col→row). On the diagonal they're the same cell, so it's shown once."""
        if self._folded:
            return self._folded_score_line(rb, cb)
        if rb == cb:
            return self._direction_score(rb, cb)
        return (f"row→col {self._direction_score(rb, cb)}   |   "
                f"col→row {self._direction_score(cb, rb)}")

    def _folded_score_line(self, rb: int, cb: int) -> str:
        """Folded, a cell holds both directions together. Between two single
        items the two directions are listed, and the combined value only if
        both directions are shown."""
        single = all(e - s == 1 for s, e in (self._blocks[rb], self._blocks[cb]))
        if single and rb == cb:
            return self._direction_score(rb, cb)
        lines = []
        if single:
            lines.append(f"row→col {self._direction_score(rb, cb)}   |   "
                         f"col→row {self._direction_score(cb, rb)}")
        cell = self._cell_at(rb, cb)
        if single and (cell is None or cell[1] < 2):
            return lines[0]   # nothing to combine
        if cell is None:
            lines.append("(no crosslinks)")
        elif self._agg_method == "count":
            lines.append(f"both directions: {cell[1]} crosslink(s)")
        else:
            score, n_links = cell
            lines.append(f"both directions: {_AGG_NAMES[self._agg_method]} of "
                         f"{n_links} crosslink(s): {score:.4f}")
        return "\n".join(lines)

    def _direction_score(self, rb: int, cb: int) -> str:
        (r_lo, r_hi), (c_lo, c_hi) = self._blocks[rb], self._blocks[cb]
        if r_hi - r_lo == 1 and c_hi - c_lo == 1:
            score = self._sparse.get((self._proteins[r_lo], self._proteins[c_lo]))
            if score is None:
                return "(no crosslink)"
            if not self._quadrant_shown(r_lo, c_lo):
                hidden = f" (hidden: {self._quadrant_of(r_lo, c_lo)} quadrant)"
            elif not self._in_cutoff(score):
                hidden = " (hidden by cutoff)"
            else:
                hidden = ""
            return f"score: {score:.4f}{hidden}"
        if not self._quadrant_shown(r_lo, c_lo):
            return f"(hidden: {self._quadrant_of(r_lo, c_lo)} quadrant)"
        cell = self._cell_at(rb, cb)
        if cell is None:
            return "(no crosslinks)"
        score, n_links = cell
        if self._agg_method == "count":
            return f"{n_links} crosslink(s)"
        return f"{_AGG_NAMES[self._agg_method]} of {n_links} crosslink(s): {score:.4f}"

    def _fetch_uniprot(self, accession: str) -> tuple[str, int] | None:
        """The full name and sequence length from UniProt, or None if that fails."""
        if accession in self._uniprot_cache:
            return self._uniprot_cache[accession]
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "xlms-matrix/1.0"})
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode())
        except Exception:
            return None
        name = (data.get("proteinDescription", {})
                    .get("recommendedName", {})
                    .get("fullName", {})
                    .get("value", accession))
        result = (name, int(data.get("sequence", {}).get("length", 0)))
        self._uniprot_cache[accession] = result
        return result

    # -------------------------------------------------------------------------
    # Drawing what's on screen
    #
    # Only the visible part is drawn: the matrix (through pygame or as a
    # numpy image), the column headers and the row labels.
    # -------------------------------------------------------------------------

    def _redraw(self) -> None:
        self._pending = False
        if not self._proteins:
            return
        try:
            if not self._overview:   # the labels are hidden while zoomed out into blocks
                self._draw_col_headers()
                self._draw_row_labels()
            self._draw_matrix()
        except Exception as exc:
            self._report_render_error(exc)
        else:
            self._last_render_error = None   # so the next error is reported again

    def _draw_matrix(self) -> None:
        if self._pg_screen is not None:
            self._draw_matrix_pygame()
        else:
            self._draw_matrix_numpy()

    def _visible_cells(self):
        """Every visible cell that has a value: the visible row and column
        blocks, then one array each for the cells' row block, column block,
        score and colour."""
        rows, cols = self._visible_blocks()
        r, c = self._cell_rows, self._cell_cols
        seen = (r >= rows.start) & (r < rows.stop) & (c >= cols.start) & (c < cols.stop)
        return rows, cols, r[seen], c[seen], self._cell_scores[seen], self._cell_rgb[seen]

    def _matrix_pixels(self, w: int, h: int) -> np.ndarray:
        """The visible part of the matrix as an RGB array of w x h pixels.
        The cells are first painted with a single pixel each and then blown
        up to the cell size. That's much faster than painting them one by one."""
        cp = self._cell_px
        ox, oy = int(self._mx.canvasx(0)), int(self._mx.canvasy(0))
        rows, cols, r, c, _, rgb = self._visible_cells()
        frame = np.full((h, w, 3), 255, dtype=np.uint8)
        if len(rows) and len(cols):
            small = np.full((len(rows), len(cols), 3), 255, dtype=np.uint8)
            small[r - rows.start, c - cols.start] = rgb
            big = small.repeat(cp, axis=0).repeat(cp, axis=1)
            # the first visible cell usually starts a little outside the canvas
            x, y = cols.start * cp - ox, rows.start * cp - oy
            part = big[max(0, -y):, max(0, -x):]
            target = frame[max(0, y):, max(0, x):]
            ph, pw = min(len(part), len(target)), min(part.shape[1], target.shape[1])
            target[:ph, :pw] = part[:ph, :pw]
        return frame

    def _score_texts(self):
        """The scores that are written into the visible cells, as (x, y,
        colour, text) with the top-left corner of the cell. Gives nothing
        while the cells are too small for text."""
        cp = self._cell_px
        if cp < SCORE_THRESHOLD_PX:
            return
        ox, oy = int(self._mx.canvasx(0)), int(self._mx.canvasy(0))
        _, _, r, c, scores, _ = self._visible_cells()
        for rb, cb, score in zip(r.tolist(), c.tolist(), scores.tolist()):
            yield (cb * cp - ox, rb * cp - oy, self._cell_colour(score, self._colour_ref),
                   self._cell_text(score, cp))

    def _overlay_positions(self) -> tuple[int, int, tuple[int, int] | None]:
        """Where the target/decoy lines go (None if there are no decoys) and
        the top-left corner of the selected cell."""
        cp, n = self._cell_px, len(self._blocks)
        ox, oy = int(self._mx.canvasx(0)), int(self._mx.canvasy(0))
        split = self._split_block if 0 < self._split_block < n else None
        sx = split * cp - ox if split is not None else None
        sy = split * cp - oy if split is not None else None
        sel = None
        if self._selected is not None:
            sel = (self._selected[1] * cp - ox, self._selected[0] * cp - oy)
        return sx, sy, sel

    def _matrix_box(self) -> tuple[int, int, int, int]:
        """Where the whole matrix lies in the view: left, top, right, bottom.
        Zoomed out it can be smaller than the view."""
        size = len(self._blocks) * self._cell_px
        ox, oy = int(self._mx.canvasx(0)), int(self._mx.canvasy(0))
        return -ox, -oy, size - ox, size - oy

    def _draw_matrix_pygame(self) -> None:
        """Draw the visible matrix into the off-screen picture, then show it."""
        self._pg_sync_size()
        w, h = self._pg_size
        screen, cp = self._pg_frame, self._cell_px
        _pygame.surfarray.blit_array(screen, self._matrix_pixels(w, h).swapaxes(0, 1))

        # cells with a score written into them are drawn on top, one by one
        for x, y, colour, text in self._score_texts():
            surf = self._pg_cache.get((colour, text))
            if surf is None:
                surf = _pygame.Surface((cp, cp))
                surf.fill(colour)
                ts = self._pg_font.render(text, True, _auto_fg(colour), colour)
                surf.blit(ts, ts.get_rect(center=(cp // 2, cp // 2)))
                self._pg_cache[(colour, text)] = surf
            screen.blit(surf, (x, y))

        sx, sy, sel = self._overlay_positions()
        bx0, by0, bx1, by1 = self._matrix_box()
        if sx is not None:
            _pygame.draw.line(screen, _SEP, (sx, max(0, by0)), (sx, min(h, by1)))
            if 0 <= sy <= h:
                _pygame.draw.line(screen, _SEP, (max(0, bx0), sy), (min(w, bx1), sy))
        if self._overview:
            _pygame.draw.rect(screen, _SEP, (bx0 - 1, by0 - 1, bx1 - bx0 + 2, by1 - by0 + 2), 1)
        if sel is not None:
            _pygame.draw.rect(screen, _HL, (*sel, cp, cp), max(2, cp // 8))

        self._present()

    def _draw_matrix_numpy(self) -> None:
        c = self._mx
        w, h = max(1, c.winfo_width()), max(1, c.winfo_height())
        cp = self._cell_px
        img = Image.fromarray(self._matrix_pixels(w, h))
        draw = ImageDraw.Draw(img)
        for x, y, colour, text in self._score_texts():
            mask, mx, my = self._text_mask(text)
            img.paste(_auto_fg(colour), (x + cp // 2 - mx, y + cp // 2 - my), mask)

        sx, sy, sel = self._overlay_positions()
        bx0, by0, bx1, by1 = self._matrix_box()
        if sx is not None:
            draw.line([(sx, max(0, by0)), (sx, min(h, by1))], fill=_SEP, width=1)
            if 0 <= sy <= h:
                draw.line([(max(0, bx0), sy), (min(w, bx1), sy)], fill=_SEP, width=1)
        if self._overview:
            draw.rectangle([bx0 - 1, by0 - 1, bx1, by1], outline=_SEP, width=1)
        if sel is not None:
            hx, hy = sel
            draw.rectangle([hx, hy, hx + cp - 1, hy + cp - 1], outline=_HL, width=max(2, cp // 8))

        self._mx_photo = ImageTk.PhotoImage(img)
        c.delete("all")
        c.create_image(int(c.canvasx(0)), int(c.canvasy(0)), image=self._mx_photo, anchor="nw")

    def _text_mask(self, text: str) -> tuple[Image.Image, int, int]:
        """A score written once and kept: a greyscale picture of the text,
        and where its centre is inside that picture."""
        cached = self._text_cache.get(text)
        if cached is None:
            w, h = 10 * self._cell_px, 2 * self._cell_px
            img = Image.new("L", (w, h), 0)
            ImageDraw.Draw(img).text((w // 2, h // 2), text, font=self._pil_cf, anchor="mm", fill=255)
            box = img.getbbox() or (0, 0, 1, 1)
            cached = self._text_cache[text] = (img.crop(box), w // 2 - box[0], h // 2 - box[1])
        return cached

    def _label_image(self, text: str) -> Image.Image:
        """A row name written once and kept, as a greyscale picture that is
        pasted wherever the name is needed."""
        img = self._label_cache.get(text)
        if img is None:
            lh = self._ch
            img = Image.new("L", (self._label_w, 2 * lh), 0)
            ImageDraw.Draw(img).text((2, lh), text, font=self._pil_lf, anchor="lm", fill=255)
            self._label_cache[text] = img
        return img

    def _slanted_label(self, text: str) -> tuple[Image.Image, int, int]:
        """A column name written once at LABEL_ANGLE, going up to the right:
        the greyscale picture, and the point (x, y) inside it where the text
        starts. That point is put under the middle of the column."""
        cached = self._slant_cache.get(text)
        if cached is None:
            lh, font = self._ch, self._pil_lf
            w = int(font.getlength(text)) + 4
            flat = Image.new("L", (w, lh), 0)
            ImageDraw.Draw(flat).text((2, lh // 2), text, font=font, anchor="lm", fill=255)
            img = flat.rotate(LABEL_ANGLE, resample=Image.Resampling.BICUBIC, expand=True)
            # rotate() turns around the middle of the picture, so the start
            # of the text ends up half the text length away from the middle
            a, d = math.radians(LABEL_ANGLE), w / 2 - 2
            x = img.width / 2 - d * math.cos(a)
            y = img.height / 2 + d * math.sin(a)
            cached = self._slant_cache[text] = (img, round(x), round(y))
        return cached

    def _section_runs(self, first_items: list, split: int):
        """Runs of neighbouring slots with the same section name, as
        (start, end, name). A run never crosses `split`. first_items[k] is
        the item that stands for slot k."""
        n = len(first_items)
        has_sep = 0 < split < n
        i = 0
        while i < n:
            label = self._section.get(first_items[i], "")
            j = i + 1
            while j < n and not (has_sep and j == split) and self._section.get(first_items[j], "") == label:
                j += 1
            yield i, j, label
            i = j

    def _draw_col_headers(self) -> None:
        c = self._ch_canvas
        w, h = max(1, c.winfo_width()), max(1, self._header_h)
        ox = int(c.canvasx(0))
        cp, lh, font = self._cell_px, self._ch, self._pil_lf
        n = len(self._blocks)
        split = self._split_block
        has_sep = 0 < split < n
        _, cols = self._visible_blocks()

        img = Image.new("RGB", (w, h), _BG)
        draw = ImageDraw.Draw(img)
        y = 0

        # section names on top (a block never spans two sections)
        if self._section is not None:
            first_items = [self._proteins[s] for s, _ in self._blocks]
            for i, j, label in self._section_runs(first_items, split):
                x0, x1 = i * cp - ox, j * cp - ox
                if x1 > 0 and x0 < w and label:
                    draw.text(((x0 + x1) // 2, y + lh // 2), label, font=font, anchor="mm", fill=_FG)
                    draw.rectangle([x0, y, x1 - 1, y + lh - 1], outline=_SEP)
        y = self._label_base   # the line under the names

        if has_sep:
            draw.line([(split * cp - ox, 0), (split * cp - ox, h)], fill=_SEP, width=1)

        # Names going up to the right, but only above single items. On a merged
        # block a name would look as if the whole block were that one protein.
        # Names of columns left of the view can still reach into it.
        labels = self._col_labels
        reach = math.ceil(self._slant_w / cp)
        for cb in range(max(0, cols.start - reach), cols.stop):
            s, e = self._blocks[cb]
            if e - s == 1 and labels[s]:
                label, lx, ly = self._slanted_label(labels[s])
                img.paste(_FG, (cb * cp + cp // 2 - ox - lx, y - LABEL_GAP_PX - ly), label)

        draw.line([(0, y), (w, y)], fill=_SEP, width=1)
        if has_sep:
            draw.text((split * cp - ox, y), "+", font=font, anchor="mm", fill=_SEP)

        self._ch_photo = ImageTk.PhotoImage(img)
        c.delete("all")
        c.create_image(ox, 0, image=self._ch_photo, anchor="nw")

    def _draw_row_labels(self) -> None:
        c = self._rl_canvas
        w, h = max(1, self._label_w), max(1, c.winfo_height())
        oy = int(c.canvasy(0))
        cp = self._cell_px
        n = len(self._blocks)
        split = self._split_block
        rows, _ = self._visible_blocks()

        img = Image.new("RGB", (w, h), _BG)
        draw = ImageDraw.Draw(img)
        for rb in rows:
            s, e = self._blocks[rb]
            if e - s == 1:
                label = self._row_labels[s]
            else:  # a merged block shows the name of its section instead
                label = self._section.get(self._proteins[s], "") if self._section is not None else ""
            if label:
                img.paste(_FG, (0, rb * cp + cp // 2 - oy - self._ch), self._label_image(label))

        if 0 < split < n and 0 <= split * cp - oy <= h:
            draw.line([(0, split * cp - oy), (w, split * cp - oy)], fill=_SEP, width=1)

        self._rl_photo = ImageTk.PhotoImage(img)
        c.delete("all")
        c.create_image(0, oy, image=self._rl_photo, anchor="nw")

    # -------------------------------------------------------------------------
    # Exporting the whole matrix as a PNG
    #
    # Unlike the screen, the export always shows every item in its own cell
    # (no blocks). The score cutoff still applies.
    # -------------------------------------------------------------------------

    def _export(self) -> None:
        if not self._proteins:
            messagebox.showerror("Export", "No data loaded.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".png", filetypes=[("PNG image", "*.png")],
            title="Export Matrix as PNG",
        )
        if not path:
            return
        self._status_var.set("Exporting…")
        threading.Thread(target=self._export_worker, args=(path,), daemon=True).start()

    def _export_worker(self, path: str) -> None:
        try:
            self._render_full_image().save(path)
            name = os.path.basename(path)
            self.after(0, lambda: self._status_var.set(f"Exported → {name}"))
        except Exception as exc:
            msg = str(exc)
            self.after(0, lambda: self._status_var.set(f"Export failed: {msg}"))

    def _render_full_image(self) -> Image.Image:
        cp = min(self._cell_px, 16)
        items = self._proteins
        n = len(items)
        split = self._split
        has_sep = 0 < split < n

        font_lbl = _find_pil_font(LABEL_FONT_SIZE + 2)
        font_cell = _find_pil_font(max(6, cp - 2))
        lh, lw = self._ch, self._label_w
        hh = self._header_h
        # the slanted names of the last columns reach past the matrix on the right
        total_w, total_h = lw + n * cp + self._slant_w, hh + n * cp

        img = Image.new("RGB", (total_w, total_h), _BG)
        draw = ImageDraw.Draw(img)

        # column headers: section names, then the item names going up to the right
        if self._section is not None:
            for i, j, label in self._section_runs(items, split):
                x0, x1 = lw + i * cp, lw + j * cp
                if label:
                    draw.text(((x0 + x1) // 2, lh // 2), label, font=font_lbl, anchor="mm", fill=_FG)
                    draw.rectangle([x0, 0, x1 - 1, lh - 1], outline=_SEP)
        y = self._label_base
        for ci, lb in enumerate(self._col_labels):
            if lb:
                label, lx, ly = self._slanted_label(lb)
                img.paste(_FG, (lw + ci * cp + cp // 2 - lx, y - LABEL_GAP_PX - ly), label)
        draw.line([(lw, y), (total_w, y)], fill=_SEP, width=1)
        if has_sep:
            draw.text((lw + split * cp, y), "+", font=font_lbl, anchor="mm", fill=_SEP)

        # row labels
        for ri in range(n):
            draw.text((2, hh + ri * cp + cp // 2), self._row_labels[ri], font=font_lbl, anchor="lm", fill=_FG)

        # lines between targets and decoys
        if has_sep:
            sx, sy = lw + split * cp, hh + split * cp
            draw.line([(0, sy), (lw, sy)], fill=_SEP, width=1)
            draw.line([(sx, 0), (sx, total_h)], fill=_SEP, width=1)
            draw.line([(lw, sy), (total_w, sy)], fill=_SEP, width=1)

        # the cells themselves: every crosslink that is shown, or folded, both
        # directions of an item pair combined in the cell below the diagonal
        keep = self._shown_links()
        rows, cols, scores = self._link_rows[keep], self._link_cols[keep], self._link_scores[keep]
        ref = None   # what the colours are measured against; None means all scores
        if self._folded:
            cell = np.maximum(rows, cols) * n + np.minimum(rows, cols)
            order = np.argsort(cell, kind="stable")
            cell = cell[order]
            first = np.flatnonzero(np.diff(cell, prepend=-1))   # where the next cell begins
            rows, cols = np.divmod(cell[first], n)
            scores = _aggregate_groups(scores[order], first.tolist() + [len(cell)], self._agg_method)
            if self._agg_method in ("count", "rss"):
                ref = np.sort(scores)   # like on screen, see _reduce_block_scores
        for ri, ci, score in zip(rows.tolist(), cols.tolist(), np.asarray(scores).tolist()):
            colour = self._cell_colour(score, ref)
            x0, y0 = lw + ci * cp, hh + ri * cp
            draw.rectangle([x0, y0, x0 + cp - 1, y0 + cp - 1], fill=colour)
            text = self._cell_text(score, cp)
            if text:
                draw.text((x0 + cp // 2, y0 + cp // 2), text, font=font_cell, anchor="mm",
                          fill=_auto_fg(colour))
        return img


# -----------------------------------------------------------------------------
# Starting the app
# -----------------------------------------------------------------------------

def main() -> None:
    MatrixApp(sys.argv[1] if len(sys.argv) > 1 else None).mainloop()


if __name__ == "__main__":
    main()

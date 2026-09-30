"""GUI matrix viewer for XL-MS crosslink data. See docs/usage.md for details.

Data comes from the same pipeline as the CLI (scripts/print_matrix.py). The
matrix is drawn either through an embedded pygame/SDL2 surface or, without
pygame, as a numpy/PIL image. Zooming out below CELL_PX_MIN merges
neighbouring items into blocks whose scores are aggregated."""
from __future__ import annotations

import json
import math
import os
import sys
import threading
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
    DEFAULT_N,
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

CELL_PX_MIN = 4
CELL_PX_MAX = 32
CELL_PX_DEFAULT = 10
SCORE_THRESHOLD_PX = 18   # cell_px at which the numeric score is drawn in each cell

BLOCK_ZOOM_STEPS = 20     # extra slider ticks below CELL_PX_MIN, used for block aggregation
BLOCK_ZOOM_GROWTH = 1.3   # per-tick multiplicative growth of block size (log-zoom feel)
BLOCK_SIZE_MAX = 500      # hard safety cap regardless of the formula above
_AGG_METHODS = {"Mean": "mean", "Geometric Mean": "geomean", "Max Score": "max"}   # dropdown -> key
_AGG_NAMES = {"mean": "Mean", "geomean": "Geometric mean", "max": "Max"}          # key -> info text

LABEL_FONT_SIZE = 9
SIDEBAR_W = 270           # width of the left control panel, in px
BG = "white"
_BG = (255, 255, 255)
_FG = (0, 0, 0)
_SEP = (153, 153, 153)
_HL = (0, 120, 215)  # selection highlight

_MONO_FONTS = [
    "C:/Windows/Fonts/consola.ttf",
    "C:/Windows/Fonts/cour.ttf",
    "C:/Windows/Fonts/lucon.ttf",
]


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def _find_pil_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    for path in _MONO_FONTS:
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def _make_lut(cmap_name: str) -> np.ndarray:
    """(256, 3) uint8 array mapping 0..255 -> RGB for a matplotlib colormap."""
    if not _MPL:
        return np.zeros((256, 3), dtype=np.uint8)
    cmap = matplotlib.colormaps[cmap_name]
    return (cmap(np.linspace(0, 1, 256))[:, :3] * 255).astype(np.uint8)


def _score_to_bg(
    score: float,
    min_s: float,
    max_s: float,
    lut: np.ndarray,
    score_sorted: np.ndarray | None = None,
) -> tuple[int, int, int]:
    """Colormap colour for a score: by rank if score_sorted is given
    (quantile normalisation), else linear between min_s and max_s."""
    if score_sorted is not None and len(score_sorted) > 0:
        idx = int(np.searchsorted(score_sorted, score) / len(score_sorted) * 255)
    elif max_s <= min_s:
        idx = 128
    else:
        idx = int((score - min_s) / (max_s - min_s) * 255)
    r, g, b = lut[max(0, min(255, idx))]
    return int(r), int(g), int(b)


def _auto_fg(bg_rgb: tuple[int, int, int]) -> tuple[int, int, int]:
    """Black or white, whichever contrasts better with bg_rgb (BT.601 luminance)."""
    lum = 0.299 * bg_rgb[0] + 0.587 * bg_rgb[1] + 0.114 * bg_rgb[2]
    return (0, 0, 0) if lum > 128 else (255, 255, 255)


def _aggregate(values: list[float], method: str) -> float:
    if method == "max":
        return max(values)
    if method == "geomean":
        positive = [v for v in values if v > 0]
        if positive:
            return math.exp(sum(math.log(v) for v in positive) / len(positive))
        # Geometric mean is undefined without positive values; fall back to the
        # arithmetic mean rather than dropping the pair (which would wrongly
        # render as "no crosslink" and hide real data).
    return sum(values) / len(values)


# --------------------------------------------------------------------------
# Two-handle range slider
# --------------------------------------------------------------------------

class RangeSlider(tk.Canvas):
    """Horizontal slider with two handles that select a [lo, hi] sub-range
    of [min, max]. command(lo, hi) fires while dragging, release_command(lo, hi)
    when the mouse button is released."""
    PAD = 9       # px between canvas edge and track end (room for a handle)
    HANDLE = 5    # handle half-width in px

    def __init__(self, master, command, release_command, **kw) -> None:
        kw.setdefault("height", 26)
        kw.setdefault("highlightthickness", 0)
        super().__init__(master, **kw)
        self._command = command
        self._release_command = release_command
        self._min, self._max = 0.0, 1.0
        self._lo, self._hi = 0.0, 1.0
        self._drag: str | None = None   # "lo", "hi" or None
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
            return self._min   # exact ends, so the extreme scores stay inclusive
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
        elif self._hi >= self._max:      # handles on top of each other at the right end
            self._drag = "lo"
        elif self._lo <= self._min:      # ... at the left end
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


# --------------------------------------------------------------------------
# Application
# --------------------------------------------------------------------------

class MatrixApp(tk.Tk):
    def __init__(self, initial_csv: str | None = None) -> None:
        super().__init__()
        self.title("XL-MS Matrix Viewer")
        self.geometry("1200x800")
        self.minsize(700, 450)

        # Data — axis items are protein names (str) or ResidueIds, see self._level
        self._level: str = "protein"
        self._proteins: list = []
        self._sparse: dict[tuple, float] = {}
        self._decoy: dict = {}
        self._section: dict | None = None
        self._seqs: dict[str, str] = {}
        self._split: int = 0              # index of the first decoy item
        self._col_labels: list[str] = []
        self._row_labels: list[str] = []
        self._item_index: dict = {}

        # Blocks — self._blocks replaces self._proteins in all geometry code.
        # At block_size 1 it is [(0,1), (1,2), ...], i.e. one block per item.
        self._cell_px: int = CELL_PX_DEFAULT
        self._block_size: int = 1
        self._blocks: list[tuple[int, int]] = []    # (start, end) item ranges
        self._item_to_block: list[int] = []
        self._split_block: int = 0
        self._agg_method: str = "mean"
        self._block_pair_raw: dict[tuple[int, int], list[float]] = {}
        self._block_scores: dict[tuple[int, int], float] = {}
        self._block_pair_n: dict[tuple[int, int], int] = {}

        # Colours
        self._cmap_name: str = "coolwarm"
        self._cmap_lut: np.ndarray | None = None   # None = greyscale
        self._norm_mode: str = "linear"
        self._min_score: float = 0.0
        self._max_score: float = 1.0
        self._score_sorted: np.ndarray = np.array([])

        # Score cutoff — cells with scores outside [lo, hi] are hidden
        self._cut_lo: float = float("-inf")
        self._cut_hi: float = float("inf")

        # Redraw / selection state
        self._pending: bool = False
        self._selected: tuple[int, int] | None = None   # (row block, col block)
        self._sel_token: int = 0                        # guards against stale UniProt replies
        self._uniprot_cache: dict[str, tuple[str, int]] = {}

        # Fonts (tk font for geometry only, PIL fonts for drawing)
        self._lf: tkfont.Font | None = None
        self._cw: int = 8
        self._ch: int = 13
        self._pil_lf: ImageFont.FreeTypeFont | ImageFont.ImageFont | None = None
        self._pil_cf: ImageFont.FreeTypeFont | ImageFont.ImageFont | None = None

        # PhotoImage references (must stay alive to prevent GC)
        self._mx_photo: ImageTk.PhotoImage | None = None
        self._ch_photo: ImageTk.PhotoImage | None = None
        self._rl_photo: ImageTk.PhotoImage | None = None

        # pygame renderer state (only used when pygame is available)
        self._pg_screen = None
        self._pg_size: tuple[int, int] = (0, 0)
        self._pg_font = None
        self._pg_cache: dict[tuple, object] = {}   # (bg colour, text) -> cell surface

        self._build_sidebar()
        self._build_matrix_area()

        self.update_idletasks()
        self._init_fonts()
        if _pygame is not None:
            self._init_pygame()

        if initial_csv:
            self._csv_var.set(initial_csv)
            self._load()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_sidebar(self) -> None:
        """All controls live in a fixed-width panel on the left, so the
        matrix gets the full window height."""
        outer = tk.Frame(self, width=SIDEBAR_W)
        outer.pack(side=tk.LEFT, fill=tk.Y, padx=(4, 0), pady=4)
        outer.pack_propagate(False)
        footer = tk.Frame(outer)
        footer.pack(side=tk.BOTTOM, fill=tk.X)
        side = self._scrollable(outer)
        wrap = SIDEBAR_W - 50   # wraplength for multi-line labels

        def group(title: str) -> tk.LabelFrame:
            frame = tk.LabelFrame(side, text=title, padx=6, pady=4)
            frame.pack(fill=tk.X, pady=(0, 6))
            frame.columnconfigure(1, weight=1)
            return frame

        def field(frame, row: int, label, widget, button=None) -> None:
            """label | widget [| button] on one grid row; label may be a StringVar."""
            if isinstance(label, tk.Variable):
                tk.Label(frame, textvariable=label).grid(row=row, column=0, sticky="w")
            else:
                tk.Label(frame, text=label).grid(row=row, column=0, sticky="w")
            widget.grid(row=row, column=1, columnspan=1 if button else 2, sticky="ew",
                        padx=(4, 0), pady=1)
            if button:
                button.grid(row=row, column=2, padx=(2, 0))

        # --- data: what to load
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
        self._only_linked_var = tk.BooleanVar(value=True)
        self._only_linked_check = tk.Checkbutton(
            data, text="Show only linked residues", variable=self._only_linked_var,
            state=tk.DISABLED,
        )
        self._only_linked_check.grid(row=3, column=0, columnspan=3, sticky="w")
        self._level_var.trace_add("write", lambda *_: self._on_level_change())
        self._n_label_var = tk.StringVar(value="N (proteins)")
        self._n_var = tk.IntVar(value=DEFAULT_N)
        field(data, 4, self._n_label_var,
              tk.Spinbox(data, textvariable=self._n_var, from_=10, to=5000, width=8))
        self._order_var = tk.StringVar(value="confidence")
        field(data, 5, "Order", ttk.Combobox(data, textvariable=self._order_var, width=10,
                                             state="readonly", values=ORDER_MODES))
        self._species_var = tk.StringVar(value="9606")
        field(data, 6, "Species", tk.Entry(data, textvariable=self._species_var, width=10))

        buttons = tk.Frame(data)
        buttons.grid(row=7, column=0, columnspan=3, sticky="ew", pady=(4, 0))
        self._load_btn = tk.Button(buttons, text="Load", command=self._load, width=8)
        self._load_btn.pack(side=tk.LEFT)
        tk.Button(buttons, text="Export PNG…", command=self._export).pack(side=tk.LEFT, padx=4)
        self._status_var = tk.StringVar(value="No data loaded.")
        tk.Label(data, textvariable=self._status_var, fg="gray", wraplength=wrap,
                 justify=tk.LEFT).grid(row=8, column=0, columnspan=3, sticky="w")

        # --- display: colours, zoom, aggregation
        display = group("Display")
        self._cmap_var = tk.StringVar(value="(none)")
        field(display, 0, "Colormap", ttk.Combobox(
            display, textvariable=self._cmap_var, width=10, state="readonly",
            values=["(none)", "coolwarm"] if _MPL else ["(none)"],
        ))
        self._cmap_var.trace_add("write", lambda *_: self._on_colormap_change())
        self._norm_var = tk.StringVar(value="linear")
        norm = tk.Frame(display)
        for label, value in [("Linear", "linear"), ("Quantile", "quantile")]:
            tk.Radiobutton(norm, text=label, variable=self._norm_var, value=value,
                           command=self._on_norm_change).pack(side=tk.LEFT)
        field(display, 1, "Norm", norm)
        if not _MPL:
            tk.Label(display, text="(matplotlib not available)", fg="red").grid(
                row=2, column=0, columnspan=3, sticky="w")

        self._zoom_label_var = tk.StringVar(value=f"{CELL_PX_DEFAULT}px")
        self._zoom = tk.Scale(
            display, from_=CELL_PX_MIN - BLOCK_ZOOM_STEPS, to=CELL_PX_MAX,
            orient=tk.HORIZONTAL, showvalue=False, command=self._on_zoom,
        )
        self._zoom.set(CELL_PX_DEFAULT)
        field(display, 3, "Zoom", self._zoom)
        tk.Label(display, textvariable=self._zoom_label_var, fg="gray").grid(
            row=4, column=1, columnspan=2, sticky="w", padx=(4, 0))
        self._agg_method_var = tk.StringVar(value="Mean")
        field(display, 5, "Aggregate", ttk.Combobox(
            display, textvariable=self._agg_method_var, width=10, state="readonly",
            values=list(_AGG_METHODS),
        ))
        self._agg_method_var.trace_add("write", lambda *_: self._on_agg_method_change())

        # --- score cutoff: bottom and top limit on one slider
        cutoff = group("Score cutoff")
        self._cutoff_slider = RangeSlider(cutoff, command=self._on_cutoff,
                                          release_command=self._on_cutoff_release)
        self._cutoff_slider.grid(row=0, column=0, columnspan=3, sticky="ew")
        self._cutoff_var = tk.StringVar(value="–")
        tk.Label(cutoff, textvariable=self._cutoff_var, justify=tk.LEFT).grid(
            row=1, column=0, columnspan=2, sticky="w")
        tk.Button(cutoff, text="Reset", command=self._reset_cutoff).grid(row=1, column=2, sticky="e")

        # --- selection info
        selection = group("Selection")
        self._sel_info_var = tk.StringVar(value="")
        tk.Label(selection, textvariable=self._sel_info_var, anchor="w", justify=tk.LEFT,
                 wraplength=wrap, fg="#1a6fb5", font=("TkFixedFont", 9)).grid(
            row=0, column=0, columnspan=3, sticky="w")

        # --- footer
        renderer = "pygame SDL2 (GPU)" if _pygame is not None else "numpy (CPU)"
        tk.Label(footer, text=f"renderer: {renderer}", fg="#aaaaaa").pack(side=tk.BOTTOM, anchor="w")
        self._info_var = tk.StringVar(value="")
        tk.Label(footer, textvariable=self._info_var, fg="gray").pack(side=tk.BOTTOM, anchor="w")

    def _scrollable(self, parent: tk.Frame) -> tk.Frame:
        """Return a frame inside a vertically scrollable canvas, so the
        sidebar stays usable when the window is shorter than its content."""
        canvas = tk.Canvas(parent, highlightthickness=0)
        bar = tk.Scrollbar(parent, orient=tk.VERTICAL, command=canvas.yview)
        canvas.configure(yscrollcommand=bar.set)
        bar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        inner = tk.Frame(canvas)
        window = canvas.create_window(0, 0, window=inner, anchor="nw")
        inner.bind("<Configure>", lambda _e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(window, width=e.width))

        def wheel(event: tk.Event) -> None:
            # global binding, but only scroll while the pointer is over the sidebar
            over = self.winfo_containing(event.x_root, event.y_root)
            if over is not None and str(over).startswith(str(parent)):
                canvas.yview_scroll(-1 if event.delta > 0 else 1, "units")

        self.bind_all("<MouseWheel>", wheel, add="+")
        return inner

    def _build_matrix_area(self) -> None:
        outer = tk.Frame(self)
        outer.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=4, pady=4)
        outer.rowconfigure(1, weight=1)
        outer.columnconfigure(1, weight=1)

        # corner | column headers
        # row labels | matrix
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
        self._mx.bind("<Button-1>", self._on_click)
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
        self._sync_sizes()

    def _update_cell_font(self) -> None:
        self._pil_cf = _find_pil_font(max(6, self._cell_px - 2))
        if self._pg_screen is not None:
            self._pg_font = _pygame.font.Font(None, max(8, self._cell_px - 2))
            self._pg_cache.clear()

    def _init_pygame(self) -> None:
        """Embed an SDL window into the matrix canvas and start its redraw loop."""
        self.update_idletasks()
        os.environ["SDL_WINDOWID"] = str(self._mx.winfo_id())
        _pygame.display.quit()
        _pygame.display.init()
        _pygame.font.init()
        size = (max(1, self._mx.winfo_width()), max(1, self._mx.winfo_height()))
        self._pg_screen = _pygame.display.set_mode(size, 0, 32)
        self._pg_size = size
        self._pg_font = _pygame.font.Font(None, max(8, self._cell_px - 2))
        self._pg_loop()

    def _pg_loop(self) -> None:
        """Continuous ~60 fps loop that keeps the pygame surface alive over tkinter repaints."""
        try:
            if self._proteins:
                self._draw_matrix()
            else:
                self._pg_screen.fill(_BG)
                _pygame.display.flip()
        except Exception:
            pass
        self.after(16, self._pg_loop)

    # ------------------------------------------------------------------
    # Geometry & scrolling
    # ------------------------------------------------------------------

    @property
    def _label_w(self) -> int:
        return ROW_W * self._cw + 6

    @property
    def _header_h(self) -> int:
        lines = (1 if self._section is not None else 0) + COL_W + 1
        return lines * self._ch + 2

    def _sync_sizes(self) -> None:
        self._corner.configure(width=self._label_w, height=self._header_h)
        self._rl_canvas.configure(width=self._label_w)
        self._ch_canvas.configure(height=self._header_h)

    def _update_scrollregion(self) -> None:
        size = len(self._blocks) * self._cell_px
        self._mx.configure(scrollregion=(0, 0, size, size))
        self._ch_canvas.configure(scrollregion=(0, 0, size, self._header_h))
        self._rl_canvas.configure(scrollregion=(0, 0, self._label_w, size))

    def _yscroll(self, *args) -> None:
        self._mx.yview(*args)
        self._rl_canvas.yview(*args)
        self._schedule()

    def _xscroll(self, *args) -> None:
        self._mx.xview(*args)
        self._ch_canvas.xview(*args)
        self._schedule()

    def _wheel(self, event: tk.Event) -> None:
        self._scroll(0, 3 if (event.num == 5 or event.delta < 0) else -3)

    def _scroll(self, dx: int, dy: int) -> None:
        if dy:
            self._mx.yview_scroll(dy, "units")
            self._rl_canvas.yview_scroll(dy, "units")
        if dx:
            self._mx.xview_scroll(dx, "units")
            self._ch_canvas.xview_scroll(dx, "units")
        self._schedule()

    def _visible_blocks(self) -> tuple[range, range]:
        """(visible row blocks, visible column blocks) of the matrix canvas."""
        n, cp, c = len(self._blocks), self._cell_px, self._mx
        rows = range(max(0, int(c.canvasy(0) // cp)), min(n, int(c.canvasy(c.winfo_height()) // cp) + 1))
        cols = range(max(0, int(c.canvasx(0) // cp)), min(n, int(c.canvasx(c.winfo_width()) // cp) + 1))
        return rows, cols

    def _schedule(self) -> None:
        """Coalesce redraw requests into one redraw on the next idle tick."""
        if not self._pending:
            self._pending = True
            self.after(0, self._redraw)

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def _on_level_change(self) -> None:
        level = self._level_var.get()
        self._only_linked_check.configure(state=tk.NORMAL if level == "residue" else tk.DISABLED)
        self._n_label_var.set(f"N ({level}s)")

    def _browse_csv(self) -> None:
        p = filedialog.askopenfilename(filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if p:
            self._csv_var.set(p)

    def _browse_fasta(self) -> None:
        p = filedialog.askopenfilename(filetypes=[("FASTA", "*.fasta *.fa *.faa"), ("All", "*.*")])
        if p:
            self._fasta_var.set(p)

    def _load(self) -> None:
        csv = self._csv_var.get().strip()
        level = self._level_var.get()
        order = self._order_var.get()
        fasta = self._fasta_var.get().strip() or None
        only_linked = self._only_linked_var.get()
        if not csv:
            messagebox.showerror("Error", "Please select a CSV file.")
            return
        if order == "sequence" and not fasta:
            messagebox.showerror("Error", "Order 'sequence' requires a FASTA file.")
            return
        if level == "residue" and not only_linked and not fasta:
            messagebox.showerror("Error", "Showing unlinked residues requires a FASTA file.")
            return
        try:
            species = int(self._species_var.get())
        except ValueError:
            messagebox.showerror("Error", "Species must be an integer NCBI taxid.")
            return
        n = self._n_var.get()

        self._load_btn.configure(state=tk.DISABLED)
        self._status_var.set("Loading…")

        def worker() -> None:
            try:
                seqs = read_fasta(fasta) if fasta else None
                if level == "residue":
                    items, sparse, decoy = build_residue_matrix(csv, n=n)
                    if not only_linked:
                        items, decoy = add_unlinked_residues(items, decoy, seqs or {})
                    items, section = sort_residues(items, decoy, order, seqs, species)
                else:
                    items, sparse, decoy = build_matrix(csv, n=n)
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
        self._cut_lo, self._cut_hi = self._min_score, self._max_score
        self._cutoff_slider.set_bounds(self._cut_lo, self._cut_hi)
        self._update_cutoff_label()
        self._recompute_blocks()

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

    # ------------------------------------------------------------------
    # Colours
    # ------------------------------------------------------------------

    def _on_colormap_change(self) -> None:
        name = self._cmap_var.get()
        if name == "(none)":
            self._cmap_lut = None
        else:
            self._cmap_name = name
            self._cmap_lut = _make_lut(name)
        self._pg_cache.clear()
        self._schedule()

    def _on_norm_change(self) -> None:
        self._norm_mode = self._norm_var.get()
        self._pg_cache.clear()
        self._schedule()

    def _cell_colour(self, score: float) -> tuple[int, int, int]:
        """Background colour of a cell: colormap if one is selected, else
        greyscale from light (low score) to dark (high score)."""
        if self._cmap_lut is not None:
            quantile = self._score_sorted if self._norm_mode == "quantile" else None
            return _score_to_bg(score, self._min_score, self._max_score, self._cmap_lut, quantile)
        mn, mx = self._min_score, self._max_score
        grey = int(200 - (score - mn) / (mx - mn) * 160) if mx > mn else 120
        grey = max(40, min(200, grey))
        return grey, grey, grey

    @staticmethod
    def _cell_text(score: float, cell_px: int) -> str | None:
        """The score printed inside a cell, once cells are big enough."""
        return f"{score:.2f}" if cell_px >= SCORE_THRESHOLD_PX else None

    # ------------------------------------------------------------------
    # Score cutoff
    # ------------------------------------------------------------------

    def _in_cutoff(self, score: float) -> bool:
        return self._cut_lo <= score <= self._cut_hi

    def _on_cutoff(self, lo: float, hi: float) -> None:
        """Slider moved: hide cells outside [lo, hi]. Colours keep the
        full-range scale, so only visibility changes."""
        self._cut_lo, self._cut_hi = lo, hi
        if self._block_size > 1:
            self._recompute_block_buckets()
            self._reduce_block_scores()
        self._update_cutoff_label()
        self._schedule()

    def _on_cutoff_release(self, lo: float, hi: float) -> None:
        if self._selected is not None:
            self._select_block(*self._selected)   # refresh the info text

    def _reset_cutoff(self) -> None:
        if not self._sparse:
            return
        self._cutoff_slider.set_bounds(self._min_score, self._max_score)
        self._on_cutoff(self._min_score, self._max_score)
        self._on_cutoff_release(self._min_score, self._max_score)

    def _update_cutoff_label(self) -> None:
        scores = self._score_sorted
        shown = int(np.searchsorted(scores, self._cut_hi, side="right")
                    - np.searchsorted(scores, self._cut_lo, side="left"))
        self._cutoff_var.set(f"{self._cut_lo:.4g} – {self._cut_hi:.4g}\n"
                             f"{shown} of {len(scores)} links shown")

    # ------------------------------------------------------------------
    # Zoom & block aggregation
    # ------------------------------------------------------------------

    def _on_zoom(self, value: str) -> None:
        v = int(float(value))
        if v >= CELL_PX_MIN:
            self._cell_px = v
            block_size = 1
            self._zoom_label_var.set(f"{v}px")
        else:
            # below the minimum cell size, zooming out merges items into blocks
            self._cell_px = CELL_PX_MIN
            steps = CELL_PX_MIN - v
            block_size = min(BLOCK_SIZE_MAX, max(2, round(BLOCK_ZOOM_GROWTH ** steps)))
            self._zoom_label_var.set(f"≤{block_size} → {CELL_PX_MIN}px")
        self._update_cell_font()
        if block_size != self._block_size:
            self._block_size = block_size
            self._recompute_blocks()
            self._pg_cache.clear()
        self._update_scrollregion()
        self._schedule()

    def _on_agg_method_change(self) -> None:
        self._agg_method = _AGG_METHODS.get(self._agg_method_var.get(), "mean")
        self._reduce_block_scores()
        self._pg_cache.clear()
        self._schedule()

    def _recompute_blocks(self) -> None:
        """Group consecutive items into blocks of at most block_size items.
        A block never crosses the target/decoy split or a section boundary."""
        items, section, split, bs = self._proteins, self._section, self._split, self._block_size
        n = len(items)
        has_sep = 0 < split < n

        def same_block(i: int, j: int) -> bool:
            if has_sep and j == split:
                return False
            return section is None or section.get(items[j], "") == section.get(items[i], "")

        blocks = []
        i = 0
        while i < n:
            j = i + 1
            if bs > 1:
                while j < min(n, i + bs) and same_block(i, j):
                    j += 1
            blocks.append((i, j))
            i = j

        self._blocks = blocks
        self._item_to_block = [bi for bi, (s, e) in enumerate(blocks) for _ in range(s, e)]
        self._split_block = self._item_to_block[split] if 0 <= split < n else len(blocks)
        self._recompute_block_buckets()
        self._reduce_block_scores()

    def _recompute_block_buckets(self) -> None:
        """Collect the raw scores falling into each directed block pair
        (row block of Protein1, column block of Protein2). Scores outside
        the cutoff are left out, so aggregates reflect only visible links."""
        buckets: dict[tuple[int, int], list[float]] = {}
        if self._block_size > 1:
            for (pi, pj), score in self._sparse.items():
                if not self._in_cutoff(score):
                    continue
                bi = self._item_to_block[self._item_index[pi]]
                bj = self._item_to_block[self._item_index[pj]]
                buckets.setdefault((bi, bj), []).append(score)
        self._block_pair_raw = buckets

    def _reduce_block_scores(self) -> None:
        self._block_scores = {k: _aggregate(v, self._agg_method) for k, v in self._block_pair_raw.items()}
        self._block_pair_n = {k: len(v) for k, v in self._block_pair_raw.items()}

    def _block_score(self, rb: int, cb: int) -> float | None:
        """Score shown in the cell at (row block, column block). Directional:
        only crosslinks with Protein1 in the row and Protein2 in the column."""
        if self._block_size <= 1:
            score = self._sparse.get((self._proteins[rb], self._proteins[cb]))
            return score if score is not None and self._in_cutoff(score) else None
        return self._block_scores.get((rb, cb))

    # ------------------------------------------------------------------
    # Selection
    # ------------------------------------------------------------------

    def _on_click(self, event: tk.Event) -> None:
        self._select_at(event.x, event.y)

    def _select_at(self, x: int, y: int) -> None:
        """Select the block under canvas-window pixel (x, y), or clear the selection."""
        if not self._blocks:
            return
        cp, n = self._cell_px, len(self._blocks)
        cb = (x + int(self._mx.canvasx(0))) // cp
        rb = (y + int(self._mx.canvasy(0))) // cp
        if 0 <= rb < n and 0 <= cb < n:
            self._select_block(rb, cb)
        else:
            self._clear_selection()
        self._schedule()

    def _clear_selection(self) -> None:
        self._selected = None
        self._sel_token += 1
        self._sel_info_var.set("")

    def _select_block(self, rb: int, cb: int) -> None:
        """Show info for a cell (a 1×1 block) or block, then look the proteins
        up on UniProt in the background if each side is a single protein."""
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
        if len({protein_of(it) for it in row_items}) == 1 and len({protein_of(it) for it in col_items}) == 1:
            def fetch() -> None:
                fetched = tuple(self._fetch_uniprot(short_accession(protein_of(items[0])))
                                for items in (row_items, col_items))
                self.after(0, lambda: self._sel_token == token and show(fetched))
            threading.Thread(target=fetch, daemon=True).start()

    def _describe_side(self, items: list, uniprot: tuple[str, int] | None, fetched: bool) -> str:
        """One info line for the row or column side of the selection."""
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
        positions = sorted(it.pos for it in items)
        if len(positions) == 1:
            return f"  pos {positions[0]}"
        return f"  pos {positions[0]}-{positions[-1]}"

    def _score_line(self, rb: int, cb: int) -> str:
        """Scores of both directions: row→col (the clicked cell) and col→row
        (its mirror cell). On the diagonal both are the same cell."""
        if rb == cb:
            return self._direction_score(rb, cb)
        return (f"row→col {self._direction_score(rb, cb)}   |   "
                f"col→row {self._direction_score(cb, rb)}")

    def _direction_score(self, rb: int, cb: int) -> str:
        (r_lo, r_hi), (c_lo, c_hi) = self._blocks[rb], self._blocks[cb]
        if r_hi - r_lo == 1 and c_hi - c_lo == 1:
            score = self._sparse.get((self._proteins[r_lo], self._proteins[c_lo]))
            if score is None:
                return "(no crosslink)"
            hidden = "" if self._in_cutoff(score) else " (hidden by cutoff)"
            return f"score: {score:.4f}{hidden}"
        score = self._block_scores.get((rb, cb))
        if score is None:
            return "(no crosslinks)"
        n_links = self._block_pair_n.get((rb, cb), 0)
        return f"{_AGG_NAMES[self._agg_method]} of {n_links} crosslink(s): {score:.4f}"

    def _fetch_uniprot(self, accession: str) -> tuple[str, int] | None:
        """(recommended full name, sequence length) from UniProt, or None."""
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

    # ------------------------------------------------------------------
    # Rendering (visible part only)
    # ------------------------------------------------------------------

    def _redraw(self) -> None:
        self._pending = False
        if not self._proteins:
            return
        self._draw_col_headers()
        self._draw_row_labels()
        self._draw_matrix()

    def _draw_matrix(self) -> None:
        if self._pg_screen is not None:
            self._draw_matrix_pygame()
        else:
            self._draw_matrix_numpy()

    def _visible_cells(self):
        """Yield (x, y, colour, text) for every non-empty visible cell,
        in canvas-window pixel coordinates."""
        cp = self._cell_px
        ox, oy = int(self._mx.canvasx(0)), int(self._mx.canvasy(0))
        rows, cols = self._visible_blocks()
        for rb in rows:
            for cb in cols:
                score = self._block_score(rb, cb)
                if score is not None:
                    yield cb * cp - ox, rb * cp - oy, self._cell_colour(score), self._cell_text(score, cp)

    def _overlay_positions(self) -> tuple[int, int, tuple[int, int] | None]:
        """Window-pixel x/y of the target/decoy separator (or None when
        absent) and the top-left of the selected cell."""
        cp, n = self._cell_px, len(self._blocks)
        ox, oy = int(self._mx.canvasx(0)), int(self._mx.canvasy(0))
        split = self._split_block if 0 < self._split_block < n else None
        sx = split * cp - ox if split is not None else None
        sy = split * cp - oy if split is not None else None
        sel = None
        if self._selected is not None:
            sel = (self._selected[1] * cp - ox, self._selected[0] * cp - oy)
        return sx, sy, sel

    def _draw_matrix_pygame(self) -> None:
        w, h = max(1, self._mx.winfo_width()), max(1, self._mx.winfo_height())
        if (w, h) != self._pg_size:
            self._pg_screen = _pygame.display.set_mode((w, h), 0, 32)
            self._pg_size = (w, h)
        screen, cp = self._pg_screen, self._cell_px
        screen.fill(_BG)

        for x, y, colour, text in self._visible_cells():
            surf = self._pg_cache.get((colour, text))
            if surf is None:
                surf = _pygame.Surface((cp, cp))
                surf.fill(colour)
                if text:
                    ts = self._pg_font.render(text, True, _auto_fg(colour), colour)
                    surf.blit(ts, ts.get_rect(center=(cp // 2, cp // 2)))
                self._pg_cache[(colour, text)] = surf
            screen.blit(surf, (x, y))

        sx, sy, sel = self._overlay_positions()
        if sx is not None:
            _pygame.draw.line(screen, _SEP, (sx, 0), (sx, h))
            if 0 <= sy <= h:
                _pygame.draw.line(screen, _SEP, (0, sy), (w, sy))
        if sel is not None:
            _pygame.draw.rect(screen, _HL, (*sel, cp, cp), max(2, cp // 8))

        # SDL owns the mouse inside its window, so clicks arrive here
        for ev in _pygame.event.get():
            if ev.type == _pygame.MOUSEBUTTONDOWN and ev.button == 1:
                self.after(0, lambda p=ev.pos: self._select_at(*p))

        _pygame.display.flip()

    def _draw_matrix_numpy(self) -> None:
        c = self._mx
        w, h = max(1, c.winfo_width()), max(1, c.winfo_height())
        cp = self._cell_px
        arr = np.full((h, w, 3), 255, dtype=np.uint8)
        texts: list[tuple[int, int, str, tuple[int, int, int]]] = []

        for x, y, colour, text in self._visible_cells():
            y0, y1 = max(0, y), min(h, y + cp)
            x0, x1 = max(0, x), min(w, x + cp)
            if y1 > y0 and x1 > x0:
                arr[y0:y1, x0:x1] = colour
                if text:
                    texts.append((x + cp // 2, y + cp // 2, text, _auto_fg(colour)))

        img = Image.fromarray(arr)
        draw = ImageDraw.Draw(img)
        for cx, cy, text, fg in texts:
            draw.text((cx, cy), text, font=self._pil_cf, anchor="mm", fill=fg)

        sx, sy, sel = self._overlay_positions()
        if sx is not None:
            draw.line([(sx, 0), (sx, h)], fill=_SEP, width=1)
            if 0 <= sy <= h:
                draw.line([(0, sy), (w, sy)], fill=_SEP, width=1)
        if sel is not None:
            hx, hy = sel
            draw.rectangle([hx, hy, hx + cp - 1, hy + cp - 1], outline=_HL, width=max(2, cp // 8))

        self._mx_photo = ImageTk.PhotoImage(img)
        c.delete("all")
        c.create_image(int(c.canvasx(0)), int(c.canvasy(0)), image=self._mx_photo, anchor="nw")

    def _section_runs(self, first_items: list, split: int):
        """Yield (start, end, label) for runs of equal section label; runs
        never cross `split`. first_items[k] is the item that represents slot k."""
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

        # section label row (every block is single-section by construction)
        if self._section is not None:
            first_items = [self._proteins[s] for s, _ in self._blocks]
            for i, j, label in self._section_runs(first_items, split):
                x0, x1 = i * cp - ox, j * cp - ox
                if x1 > 0 and x0 < w and label:
                    draw.text(((x0 + x1) // 2, y + lh // 2), label, font=font, anchor="mm", fill=_FG)
                    draw.rectangle([x0, y, x1 - 1, y + lh - 1], outline=_SEP)
            y += lh

        if has_sep:
            draw.line([(split * cp - ox, 0), (split * cp - ox, h)], fill=_SEP, width=1)

        # vertical labels — only for single-item blocks; a character label on a
        # multi-item block would suggest the whole block is one item
        labels = self._col_labels
        for char_idx in range(max((len(lb) for lb in labels), default=0)):
            for cb in cols:
                s, e = self._blocks[cb]
                if e - s == 1 and char_idx < len(labels[s]) and labels[s][char_idx] != " ":
                    draw.text((cb * cp + cp // 2 - ox, y + lh // 2), labels[s][char_idx],
                              font=font, anchor="mm", fill=_FG)
            y += lh

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
            else:  # multi-item block: show its section name, if any
                label = self._section.get(self._proteins[s], "") if self._section is not None else ""
            if label:
                draw.text((2, rb * cp + cp // 2 - oy), label, font=self._pil_lf, anchor="lm", fill=_FG)

        if 0 < split < n and 0 <= split * cp - oy <= h:
            draw.line([(0, split * cp - oy), (w, split * cp - oy)], fill=_SEP, width=1)

        self._rl_photo = ImageTk.PhotoImage(img)
        c.delete("all")
        c.create_image(0, oy, image=self._rl_photo, anchor="nw")

    # ------------------------------------------------------------------
    # PNG export (full matrix, one cell per item, no aggregation, cutoff applied)
    # ------------------------------------------------------------------

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
        hh = ((1 if self._section else 0) + COL_W + 1) * lh + 2
        total_w, total_h = lw + n * cp, hh + n * cp

        img = Image.new("RGB", (total_w, total_h), _BG)
        draw = ImageDraw.Draw(img)

        # column headers: section row, then vertical labels
        y = 0
        if self._section is not None:
            for i, j, label in self._section_runs(items, split):
                x0, x1 = lw + i * cp, lw + j * cp
                if label:
                    draw.text(((x0 + x1) // 2, y + lh // 2), label, font=font_lbl, anchor="mm", fill=_FG)
                    draw.rectangle([x0, y, x1 - 1, y + lh - 1], outline=_SEP)
            y += lh
        labels = self._col_labels
        for char_idx in range(max((len(lb) for lb in labels), default=0)):
            for ci, lb in enumerate(labels):
                if char_idx < len(lb) and lb[char_idx] != " ":
                    draw.text((lw + ci * cp + cp // 2, y + lh // 2), lb[char_idx],
                              font=font_lbl, anchor="mm", fill=_FG)
            y += lh
        draw.line([(lw, y), (total_w, y)], fill=_SEP, width=1)
        if has_sep:
            draw.text((lw + split * cp, y), "+", font=font_lbl, anchor="mm", fill=_SEP)

        # row labels
        for ri in range(n):
            draw.text((2, hh + ri * cp + cp // 2), self._row_labels[ri], font=font_lbl, anchor="lm", fill=_FG)

        # target/decoy separators
        if has_sep:
            sx, sy = lw + split * cp, hh + split * cp
            draw.line([(0, sy), (lw, sy)], fill=_SEP, width=1)
            draw.line([(sx, 0), (sx, total_h)], fill=_SEP, width=1)
            draw.line([(lw, sy), (total_w, sy)], fill=_SEP, width=1)

        # cells
        for ri, pi in enumerate(items):
            for ci, pj in enumerate(items):
                score = self._sparse.get((pi, pj))
                if score is None or not self._in_cutoff(score):
                    continue
                colour = self._cell_colour(score)
                x0, y0 = lw + ci * cp, hh + ri * cp
                draw.rectangle([x0, y0, x0 + cp - 1, y0 + cp - 1], fill=colour)
                text = self._cell_text(score, cp)
                if text:
                    draw.text((x0 + cp // 2, y0 + cp // 2), text, font=font_cell, anchor="mm",
                              fill=_auto_fg(colour))
        return img


def main() -> None:
    MatrixApp(sys.argv[1] if len(sys.argv) > 1 else None).mainloop()


if __name__ == "__main__":
    main()

"""
GUI matrix viewer for XL-MS crosslink data.

Usage:
    python gui/matrix_app.py [csv_path]
"""
from __future__ import annotations

import os
import sys
import threading
import tkinter as tk
import tkinter.filedialog as filedialog
import tkinter.font as tkfont
import tkinter.messagebox as messagebox
import tkinter.ttk as ttk
from typing import Optional

import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageTk

try:
    import pygame as _pygame
except Exception:
    _pygame = None  # type: ignore[assignment]

_RENDERER_LABEL = "pygame SDL2 (GPU)" if _pygame is not None else "numpy (CPU)"

try:
    import matplotlib.cm as _mpl_cm
    import matplotlib.pyplot as _mpl_plt
    _MPL = True
except Exception:
    _MPL = False

_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_DIR, "..", "scripts"))
sys.path.insert(0, os.path.join(_DIR, "..", "src"))

from print_matrix import (  # type: ignore[import]
    build_matrix,
    _sort_proteins,
    _build_thresholds,
    _score_to_dot,
    DEFAULT_N,
    COL_W,
    ROW_W,
)
from xlms.io import read_fasta  # type: ignore[import]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CELL_PX_MIN = 4
CELL_PX_MAX = 32
CELL_PX_DEFAULT = 10
SCORE_THRESHOLD_PX = 18  # cell_px at which scores replace dot symbols

LABEL_FONT_SIZE = 9
SEP_COLOR = "#999999"
BG = "white"

_BG = (255, 255, 255)
_FG = (0, 0, 0)
_SEP = (153, 153, 153)
_BG_PG  = (255, 255, 255)
_FG_PG  = (0, 0, 0)
_SEP_PG = (153, 153, 153)
_HL            = (0, 120, 215)   # highlight — PIL and pygame share the same RGB
_CONTOUR_COLOR = (200, 0, 0)     # dark-red contour isoline color

_MONO_FONTS = [
    "C:/Windows/Fonts/consola.ttf",
    "C:/Windows/Fonts/cour.ttf",
    "C:/Windows/Fonts/lucon.ttf",
]

# Characters that can appear inside a matrix cell
_GLYPH_CHARS = [*"·○●⬤", *"0123456789"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_pil_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    for path in _MONO_FONTS:
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def _build_glyph_atlas(
    font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
    cell_px: int,
) -> dict[str, np.ndarray]:
    """Pre-render each glyph to a (cell_px, cell_px, 3) uint8 numpy array."""
    atlas: dict[str, np.ndarray] = {}
    for ch in _GLYPH_CHARS:
        g = Image.new("RGB", (cell_px, cell_px), _BG)
        ImageDraw.Draw(g).text((cell_px // 2, cell_px // 2), ch,
                               font=font, anchor="mm", fill=_FG)
        atlas[ch] = np.array(g, dtype=np.uint8)
    return atlas


def _cell_label(score: float, cell_px: int, thresholds: list[float]) -> str:
    if cell_px >= SCORE_THRESHOLD_PX:
        return f"{score:.0f}"
    return _score_to_dot(score, thresholds)


def _make_lut(cmap_name: str) -> np.ndarray:
    """Return a (256, 3) uint8 array mapping [0..255] → RGB for the named colormap."""
    if not _MPL:
        return np.zeros((256, 3), dtype=np.uint8)
    try:
        import matplotlib
        cmap = matplotlib.colormaps[cmap_name]
    except (AttributeError, KeyError):
        cmap = _mpl_cm.get_cmap(cmap_name)  # type: ignore[attr-defined]
    return (cmap(np.linspace(0, 1, 256))[:, :3] * 255).astype(np.uint8)


def _score_to_bg(
    score: float, min_s: float, max_s: float, lut: np.ndarray
) -> tuple[int, int, int]:
    """Map a score to an RGB background colour via LUT lookup."""
    if max_s <= min_s:
        idx = 128
    else:
        idx = int((score - min_s) / (max_s - min_s) * 255)
        idx = max(0, min(255, idx))
    r, g, b = lut[idx]
    return int(r), int(g), int(b)


def _auto_fg(bg_rgb: tuple[int, int, int]) -> tuple[int, int, int]:
    """Return black or white for best contrast against bg_rgb (BT.601 luminance)."""
    lum = 0.299 * bg_rgb[0] + 0.587 * bg_rgb[1] + 0.114 * bg_rgb[2]
    return (0, 0, 0) if lum > 128 else (255, 255, 255)


def _compute_contour_segs(
    proteins: list[str],
    sparse: dict[tuple[str, str], float],
    split: int,
    n_levels: int,
) -> list[list[tuple[float, float]]]:
    """
    Return isoline polylines in fractional cell-index space (0..n-1).
    Coordinates are zoom/scroll-independent — multiply by cell_px at draw time.
    """
    if not _MPL or not proteins:
        return []
    n = len(proteins)
    grid = np.zeros((n, n), dtype=np.float32)
    for ri, pi in enumerate(proteins):
        for ci, pj in enumerate(proteins):
            v = sparse.get((min(pi, pj), max(pi, pj)))
            if v is not None:
                grid[ri, ci] = v
    try:
        fig, ax = _mpl_plt.subplots()
        cs = ax.contour(np.arange(n), np.arange(n), grid, levels=n_levels)
        segs: list[list[tuple[float, float]]] = []
        for collection in cs.collections:
            for path in collection.get_paths():
                verts = path.vertices
                segs.append([(float(v[0]), float(v[1])) for v in verts])
        _mpl_plt.close(fig)
        return segs
    except Exception:
        return []


# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------

class MatrixApp(tk.Tk):
    def __init__(self, initial_csv: str | None = None) -> None:
        super().__init__()
        self.title("XL-MS Matrix Viewer")
        self.geometry("1200x800")
        self.minsize(600, 400)

        # Data
        self._proteins: list[str] = []
        self._sparse: dict[tuple[str, str], float] = {}
        self._decoy: dict[str, bool] = {}
        self._section: dict[str, str] | None = None
        self._thresholds: list[float] = []
        self._split: int = 0

        # Render state
        self._cell_px: int = CELL_PX_DEFAULT
        self._pending: bool = False
        self._selected: tuple[int, int] | None = None

        # Heatmap colormap
        self._cmap_name: str = "viridis"
        self._cmap_lut: np.ndarray | None = None  # None = disabled
        self._min_score: float = 0.0
        self._max_score: float = 1.0

        # Contour overlay
        self._contour_segs: list[list[tuple[float, float]]] = []
        # _contour_enabled and _contour_levels are tk.Var, created in _build_controls

        # Tkinter fonts (for geometry only — not used for drawing)
        self._lf: tkfont.Font | None = None
        self._cw: int = 8
        self._ch: int = 13

        # PIL fonts
        self._pil_lf: ImageFont.FreeTypeFont | ImageFont.ImageFont | None = None
        self._pil_cf: ImageFont.FreeTypeFont | ImageFont.ImageFont | None = None

        # Glyph atlas: char -> (cell_px, cell_px, 3) uint8 array
        self._atlas: dict[str, np.ndarray] = {}

        # PhotoImage references (must stay alive to prevent GC)
        self._mx_photo: ImageTk.PhotoImage | None = None
        self._ch_photo: ImageTk.PhotoImage | None = None
        self._rl_photo: ImageTk.PhotoImage | None = None

        self._build_controls()
        self._build_matrix_area()
        self._build_statusbar()

        self.update_idletasks()
        self._init_fonts()
        if _pygame is not None:
            self._init_pygame_embed()

        if initial_csv:
            self._csv_var.set(initial_csv)
            self._load()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_controls(self) -> None:
        bar = tk.Frame(self, bd=1, relief=tk.GROOVE, pady=4)
        bar.pack(side=tk.TOP, fill=tk.X, padx=4, pady=2)

        row1 = tk.Frame(bar)
        row1.pack(fill=tk.X, padx=4)
        tk.Label(row1, text="CSV:").pack(side=tk.LEFT)
        self._csv_var = tk.StringVar()
        tk.Entry(row1, textvariable=self._csv_var, width=45).pack(side=tk.LEFT, padx=2)
        tk.Button(row1, text="Browse…", command=self._browse_csv).pack(side=tk.LEFT)
        tk.Label(row1, text="  FASTA:").pack(side=tk.LEFT)
        self._fasta_var = tk.StringVar()
        tk.Entry(row1, textvariable=self._fasta_var, width=35).pack(side=tk.LEFT, padx=2)
        tk.Button(row1, text="Browse…", command=self._browse_fasta).pack(side=tk.LEFT)

        row2 = tk.Frame(bar)
        row2.pack(fill=tk.X, padx=4, pady=(2, 0))
        tk.Label(row2, text="N:").pack(side=tk.LEFT)
        self._n_var = tk.IntVar(value=DEFAULT_N)
        tk.Spinbox(row2, textvariable=self._n_var, from_=10, to=5000, width=6).pack(side=tk.LEFT, padx=2)
        tk.Label(row2, text="  Order:").pack(side=tk.LEFT)
        self._order_var = tk.StringVar(value="confidence")
        ttk.Combobox(
            row2, textvariable=self._order_var, width=12, state="readonly",
            values=["confidence", "alpha", "cluster", "sequence", "pathway", "complex"],
        ).pack(side=tk.LEFT, padx=2)
        tk.Label(row2, text="  Species:").pack(side=tk.LEFT)
        self._species_var = tk.StringVar(value="9606")
        tk.Entry(row2, textvariable=self._species_var, width=8).pack(side=tk.LEFT, padx=2)
        self._load_btn = tk.Button(row2, text="Load", command=self._load, width=8)
        self._load_btn.pack(side=tk.LEFT, padx=8)
        tk.Button(row2, text="Export PNG…", command=self._export, width=12).pack(side=tk.LEFT, padx=4)
        self._status_var = tk.StringVar(value="No data loaded.")
        tk.Label(row2, textvariable=self._status_var, fg="gray").pack(side=tk.LEFT, padx=4)

        row3 = tk.Frame(bar)
        row3.pack(fill=tk.X, padx=4, pady=(2, 0))
        tk.Label(row3, text="Colormap:").pack(side=tk.LEFT)
        self._cmap_var = tk.StringVar(value="(none)")
        _cmap_choices = (
            ["(none)", "viridis", "plasma", "inferno", "magma", "cividis", "coolwarm", "RdYlGn"]
            if _MPL else ["(none)"]
        )
        ttk.Combobox(
            row3, textvariable=self._cmap_var, width=10,
            state="readonly", values=_cmap_choices,
        ).pack(side=tk.LEFT, padx=2)
        self._cmap_var.trace_add("write", lambda *_: self._on_colormap_change())

        tk.Label(row3, text="   Contours:").pack(side=tk.LEFT)
        self._contour_enabled = tk.BooleanVar(value=False)
        self._contour_levels = tk.IntVar(value=5)
        tk.Checkbutton(
            row3, text="Show", variable=self._contour_enabled,
            command=self._on_contour_change,
        ).pack(side=tk.LEFT)
        tk.Label(row3, text="Levels:").pack(side=tk.LEFT, padx=(6, 0))
        tk.Spinbox(
            row3, textvariable=self._contour_levels, from_=2, to=20,
            width=4, command=self._on_contour_change,
        ).pack(side=tk.LEFT, padx=2)
        if not _MPL:
            tk.Label(row3, text="(matplotlib not available)", fg="red").pack(side=tk.LEFT, padx=4)

    def _build_matrix_area(self) -> None:
        outer = tk.Frame(self)
        outer.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=4, pady=2)
        outer.rowconfigure(1, weight=1)
        outer.columnconfigure(1, weight=1)

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
        self._mx.bind("<Up>",    lambda e: self._arrow(0, -1))
        self._mx.bind("<Down>",  lambda e: self._arrow(0,  1))
        self._mx.bind("<Left>",  lambda e: self._arrow(-1, 0))
        self._mx.bind("<Right>", lambda e: self._arrow( 1, 0))
        self._mx.focus_set()

    def _build_statusbar(self) -> None:
        bar = tk.Frame(self, bd=1, relief=tk.GROOVE, pady=3)
        bar.pack(side=tk.BOTTOM, fill=tk.X, padx=4, pady=2)

        row1 = tk.Frame(bar)
        row1.pack(fill=tk.X)
        tk.Label(row1, text="Zoom:").pack(side=tk.LEFT, padx=(4, 0))
        self._zoom = tk.Scale(
            row1, from_=CELL_PX_MIN, to=CELL_PX_MAX,
            orient=tk.HORIZONTAL, length=180, showvalue=True,
            command=self._on_zoom,
        )
        self._zoom.set(CELL_PX_DEFAULT)
        self._zoom.pack(side=tk.LEFT, padx=4)
        self._info_var = tk.StringVar(value="")
        tk.Label(row1, textvariable=self._info_var, fg="gray").pack(side=tk.LEFT, padx=8)
        tk.Label(row1, text=f"renderer: {_RENDERER_LABEL}", fg="#aaaaaa").pack(side=tk.RIGHT, padx=8)

        self._sel_info_var = tk.StringVar(value="")
        tk.Label(bar, textvariable=self._sel_info_var, anchor="w",
                 fg="#1a6fb5", font=("TkFixedFont", 9)).pack(fill=tk.X, padx=8, pady=(0, 2))

    # ------------------------------------------------------------------
    # Font / geometry
    # ------------------------------------------------------------------

    def _init_fonts(self) -> None:
        self._lf = tkfont.Font(family="TkFixedFont", size=LABEL_FONT_SIZE)
        self._cw = self._lf.measure("A")
        self._ch = self._lf.metrics("linespace")
        self._pil_lf = _find_pil_font(LABEL_FONT_SIZE + 2)
        self._pil_cf = _find_pil_font(max(6, CELL_PX_DEFAULT - 2))
        self._atlas = _build_glyph_atlas(self._pil_cf, CELL_PX_DEFAULT)
        self._sync_sizes()

    def _update_cell_font(self) -> None:
        self._pil_cf = _find_pil_font(max(6, self._cell_px - 2))
        self._atlas = _build_glyph_atlas(self._pil_cf, self._cell_px)
        if _pygame is not None and hasattr(self, "_pg_screen"):
            self._pg_font = _pygame.font.Font(None, max(8, self._cell_px - 2))
            self._pg_score_cache = {}
            self._build_pg_atlas()

    def _init_pygame_embed(self) -> None:
        self.update_idletasks()
        os.environ["SDL_WINDOWID"] = str(self._mx.winfo_id())
        _pygame.display.quit()
        _pygame.display.init()
        _pygame.font.init()
        w = max(1, self._mx.winfo_width())
        h = max(1, self._mx.winfo_height())
        self._pg_screen: _pygame.Surface = _pygame.display.set_mode((w, h), 0, 32)  # type: ignore[name-defined]
        self._pg_size: tuple[int, int] = (w, h)
        self._pg_atlas: dict[str, _pygame.Surface] = {}  # type: ignore[name-defined]
        self._pg_font: _pygame.font.Font = _pygame.font.Font(None, max(8, self._cell_px - 2))  # type: ignore[name-defined]
        self._pg_score_cache: dict[str, _pygame.Surface] = {}  # type: ignore[name-defined]
        self._build_pg_atlas()
        self._pg_loop()

    def _pg_loop(self) -> None:
        """Continuous 60 fps loop that keeps the pygame surface alive over tkinter repaints."""
        try:
            if _pygame is not None and hasattr(self, "_pg_screen"):
                if self._proteins:
                    self._draw_matrix()
                else:
                    self._pg_screen.fill(_BG_PG)
                    _pygame.display.flip()
        except Exception:
            pass
        self.after(16, self._pg_loop)

    def _build_pg_atlas(self) -> None:
        cp = self._cell_px
        font = self._pil_cf
        self._pg_atlas = {}
        for ch in _GLYPH_CHARS:
            g = Image.new("RGB", (cp, cp), _BG)
            ImageDraw.Draw(g).text((cp // 2, cp // 2), ch, font=font, anchor="mm", fill=_FG)
            arr = np.array(g, dtype=np.uint8)
            surf = _pygame.surfarray.make_surface(np.swapaxes(arr, 0, 1))
            self._pg_atlas[ch] = surf.convert()

    @property
    def _label_w(self) -> int:
        return ROW_W * self._cw + 6

    @property
    def _header_h(self) -> int:
        lines = (1 if self._section is not None else 0) + COL_W + 1
        return lines * self._ch + 2

    def _sync_sizes(self) -> None:
        lw, hh = self._label_w, self._header_h
        self._corner.configure(width=lw, height=hh)
        self._rl_canvas.configure(width=lw)
        self._ch_canvas.configure(height=hh)

    def _virtual(self) -> tuple[int, int]:
        s = len(self._proteins) * self._cell_px
        return s, s

    def _update_scrollregion(self) -> None:
        vw, vh = self._virtual()
        self._mx.configure(scrollregion=(0, 0, vw, vh))
        self._ch_canvas.configure(scrollregion=(0, 0, vw, self._header_h))
        self._rl_canvas.configure(scrollregion=(0, 0, self._label_w, vh))

    # ------------------------------------------------------------------
    # Scroll / zoom
    # ------------------------------------------------------------------

    def _yscroll(self, *args) -> None:
        self._mx.yview(*args)
        self._rl_canvas.yview(*args)
        if _pygame is not None and hasattr(self, "_pg_screen") and self._proteins:
            self._draw_matrix()
        self._schedule(fast=True)

    def _xscroll(self, *args) -> None:
        self._mx.xview(*args)
        self._ch_canvas.xview(*args)
        if _pygame is not None and hasattr(self, "_pg_screen") and self._proteins:
            self._draw_matrix()
        self._schedule(fast=True)

    def _wheel(self, event: tk.Event) -> None:
        delta = 3 if (event.num == 5 or event.delta < 0) else -3
        self._mx.yview_scroll(delta, "units")
        self._rl_canvas.yview_scroll(delta, "units")
        if _pygame is not None and hasattr(self, "_pg_screen") and self._proteins:
            self._draw_matrix()
        self._schedule(fast=True)

    def _arrow(self, dx: int, dy: int) -> None:
        if dy:
            self._mx.yview_scroll(dy, "units")
            self._rl_canvas.yview_scroll(dy, "units")
        if dx:
            self._mx.xview_scroll(dx, "units")
            self._ch_canvas.xview_scroll(dx, "units")
        if _pygame is not None and hasattr(self, "_pg_screen") and self._proteins:
            self._draw_matrix()
        self._schedule(fast=True)

    def _on_zoom(self, value: str) -> None:
        self._cell_px = int(float(value))
        self._update_cell_font()
        self._update_scrollregion()
        self._schedule()

    def _on_click(self, event: tk.Event) -> None:
        """Handle left-click on the matrix canvas (numpy path fallback)."""
        if not self._proteins:
            return
        ox = int(self._mx.canvasx(0))
        oy = int(self._mx.canvasy(0))
        col = (event.x + ox) // self._cell_px
        row = (event.y + oy) // self._cell_px
        n = len(self._proteins)
        if 0 <= row < n and 0 <= col < n:
            self._on_cell_select(row, col)
        else:
            self._selected = None
            self._sel_info_var.set("")
        self._schedule()

    def _on_cell_select(self, row: int, col: int) -> None:
        self._selected = (row, col)
        pi = self._proteins[row]
        pj = self._proteins[col]
        key = (min(pi, pj), max(pi, pj))
        score = self._sparse.get(key)
        di = "decoy" if self._decoy.get(pi) else "target"
        dj = "decoy" if self._decoy.get(pj) else "target"
        if score is not None:
            self._sel_info_var.set(
                f"{pi}  ×  {pj}    score: {score:.2f}    {di} – {dj}"
            )
        else:
            self._sel_info_var.set(
                f"{pi}  ×  {pj}    no link    {di} – {dj}"
            )

    def _clear_selection(self) -> None:
        self._selected = None
        self._sel_info_var.set("")

    def _schedule(self, fast: bool = False) -> None:
        if not self._pending:
            self._pending = True
            self.after(0, self._redraw)

    # ------------------------------------------------------------------
    # File loading
    # ------------------------------------------------------------------

    def _browse_csv(self) -> None:
        p = filedialog.askopenfilename(filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if p:
            self._csv_var.set(p)

    def _browse_fasta(self) -> None:
        p = filedialog.askopenfilename(
            filetypes=[("FASTA", "*.fasta *.fa *.faa"), ("All", "*.*")]
        )
        if p:
            self._fasta_var.set(p)

    def _load(self) -> None:
        csv = self._csv_var.get().strip()
        if not csv:
            messagebox.showerror("Error", "Please select a CSV file.")
            return
        order = self._order_var.get()
        fasta = self._fasta_var.get().strip() or None
        if order == "sequence" and not fasta:
            messagebox.showerror("Error", "Order 'sequence' requires a FASTA file.")
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
                proteins, sparse, decoy = build_matrix(csv, n=n)
                seqs = read_fasta(fasta) if fasta else None
                proteins, section = _sort_proteins(proteins, sparse, decoy, order, seqs, species)
                self.after(0, lambda: self._on_loaded(proteins, sparse, decoy, section))
            except Exception as exc:
                msg = str(exc)
                self.after(0, lambda: self._on_error(msg))

        threading.Thread(target=worker, daemon=True).start()

    def _on_colormap_change(self) -> None:
        name = self._cmap_var.get()
        if name == "(none)":
            self._cmap_lut = None
        else:
            self._cmap_name = name
            self._cmap_lut = _make_lut(name)
        if hasattr(self, "_pg_score_cache"):
            self._pg_score_cache = {}
        self._schedule()

    def _on_contour_change(self) -> None:
        if self._contour_enabled.get() and self._proteins:
            self._contour_segs = _compute_contour_segs(
                self._proteins, self._sparse, self._split, self._contour_levels.get()
            )
        else:
            self._contour_segs = []
        self._schedule()

    def _on_loaded(
        self,
        proteins: list[str],
        sparse: dict[tuple[str, str], float],
        decoy: dict[str, bool],
        section: dict[str, str] | None,
    ) -> None:
        self._proteins = proteins
        self._sparse = sparse
        self._decoy = decoy
        self._section = section
        self._thresholds = _build_thresholds(sparse) if sparse else []
        self._split = next(
            (i for i, p in enumerate(proteins) if decoy.get(p, False)), len(proteins)
        )
        self._min_score = min(sparse.values()) if sparse else 0.0
        self._max_score = max(sparse.values()) if sparse else 1.0
        if self._cmap_lut is not None:
            self._cmap_lut = _make_lut(self._cmap_name)
        if self._contour_enabled.get():
            self._contour_segs = _compute_contour_segs(
                proteins, sparse, self._split, self._contour_levels.get()
            )
        else:
            self._contour_segs = []
        n_t = self._split
        n_d = len(proteins) - n_t
        self._status_var.set(f"Loaded — {n_t} targets · {n_d} decoys · {len(sparse)} links")
        self._info_var.set(f"{len(proteins)} proteins")
        self._load_btn.configure(state=tk.NORMAL)
        self._sync_sizes()
        self._update_scrollregion()
        self._mx.xview_moveto(0)
        self._mx.yview_moveto(0)
        self._ch_canvas.xview_moveto(0)
        self._rl_canvas.yview_moveto(0)
        self._schedule()

    def _on_error(self, msg: str) -> None:
        self._load_btn.configure(state=tk.NORMAL)
        self._status_var.set(f"Error: {msg}")
        messagebox.showerror("Load error", msg)

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def _export(self) -> None:
        if not self._proteins:
            messagebox.showerror("Export", "No data loaded.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG image", "*.png")],
            title="Export Matrix as PNG",
        )
        if not path:
            return
        self._status_var.set("Exporting…")
        threading.Thread(target=self._export_worker, args=(path,), daemon=True).start()

    def _export_worker(self, path: str) -> None:
        try:
            img = self._render_full_image()
            img.save(path)
            name = os.path.basename(path)
            self.after(0, lambda: self._status_var.set(f"Exported → {name}"))
        except Exception as exc:
            msg = str(exc)
            self.after(0, lambda: self._status_var.set(f"Export failed: {msg}"))

    def _render_full_image(self) -> Image.Image:
        cp = min(self._cell_px, 16)
        n = len(self._proteins)
        proteins = self._proteins
        sparse = self._sparse
        thresholds = self._thresholds
        split = self._split
        has_sep = 0 < split < n

        font_lbl = _find_pil_font(LABEL_FONT_SIZE + 2)
        font_cell = _find_pil_font(max(6, cp - 2))
        lh = self._ch
        lw = self._label_w

        hdr_lines = (1 if self._section else 0) + COL_W + 1
        hh = hdr_lines * lh + 2
        total_w = lw + n * cp
        total_h = hh + n * cp

        img = Image.new("RGB", (total_w, total_h), _BG)
        draw = ImageDraw.Draw(img)

        # ── Column headers ───────────────────────────────────────────────
        y = 0
        if self._section is not None:
            i = 0
            while i < n:
                sec = self._section.get(proteins[i], "")
                j = i + 1
                while j < n and not (has_sep and j == split) \
                        and self._section.get(proteins[j], "") == sec:
                    j += 1
                x0s = lw + i * cp;  x1s = lw + j * cp
                if sec:
                    draw.text(((x0s + x1s) // 2, y + lh // 2), sec,
                              font=font_lbl, anchor="mm", fill=_FG)
                    draw.rectangle([x0s, y, x1s - 1, y + lh - 1], outline=_SEP)
                i = j
            y += lh

        labels_col = [p[:COL_W] for p in proteins]
        max_char = max((len(lb) for lb in labels_col), default=0)
        for char_idx in range(max_char):
            for ci in range(n):
                lb = labels_col[ci]
                ch = lb[char_idx] if char_idx < len(lb) else " "
                if ch != " ":
                    draw.text((lw + ci * cp + cp // 2, y + lh // 2), ch,
                              font=font_lbl, anchor="mm", fill=_FG)
            y += lh

        draw.line([(lw, y), (total_w, y)], fill=_SEP, width=1)
        if has_sep:
            draw.text((lw + split * cp, y), "+", font=font_lbl, anchor="mm", fill=_SEP)

        # ── Row labels ───────────────────────────────────────────────────
        for ri in range(n):
            cy = hh + ri * cp + cp // 2
            draw.text((2, cy), proteins[ri][:ROW_W], font=font_lbl, anchor="lm", fill=_FG)
        if has_sep:
            draw.line([(0, hh + split * cp), (lw, hh + split * cp)], fill=_SEP, width=1)

        # ── Separator lines through full matrix ──────────────────────────
        if has_sep:
            sx = lw + split * cp
            draw.line([(sx, 0), (sx, total_h)], fill=_SEP, width=1)
            sy = hh + split * cp
            draw.line([(lw, sy), (total_w, sy)], fill=_SEP, width=1)

        # ── Matrix cells ─────────────────────────────────────────────────
        cmap_lut = self._cmap_lut
        for ri in range(n):
            pi = proteins[ri]
            for ci in range(n):
                pj = proteins[ci]
                score = sparse.get((min(pi, pj), max(pi, pj)))
                if score is None:
                    continue
                label = _cell_label(score, cp, thresholds)
                if cmap_lut is not None:
                    bg = _score_to_bg(score, self._min_score, self._max_score, cmap_lut)
                    x0c, y0c = lw + ci * cp, hh + ri * cp
                    draw.rectangle([x0c, y0c, x0c + cp - 1, y0c + cp - 1], fill=bg)
                    draw.text(
                        (x0c + cp // 2, y0c + cp // 2),
                        label, font=font_cell, anchor="mm", fill=_auto_fg(bg),
                    )
                else:
                    draw.text(
                        (lw + ci * cp + cp // 2, hh + ri * cp + cp // 2),
                        label, font=font_cell, anchor="mm", fill=_FG,
                    )

        for polyline in self._contour_segs:
            pts = [(lw + int(x * cp), hh + int(y * cp)) for x, y in polyline]
            if len(pts) >= 2:
                draw.line(pts, fill=_CONTOUR_COLOR, width=1)

        return img

    # ------------------------------------------------------------------
    # Viewport helpers
    # ------------------------------------------------------------------

    def _visible_cols(self) -> tuple[int, int]:
        n = len(self._proteins)
        cp = self._cell_px
        x0 = self._mx.canvasx(0)
        x1 = self._mx.canvasx(self._mx.winfo_width())
        return max(0, int(x0 // cp)), min(n, int(x1 // cp) + 1)

    def _visible_rows(self) -> tuple[int, int]:
        n = len(self._proteins)
        cp = self._cell_px
        y0 = self._mx.canvasy(0)
        y1 = self._mx.canvasy(self._mx.winfo_height())
        return max(0, int(y0 // cp)), min(n, int(y1 // cp) + 1)

    # ------------------------------------------------------------------
    # Drawing (PIL-based — one image per canvas per frame)
    # ------------------------------------------------------------------

    def _redraw(self) -> None:
        self._pending = False
        if not self._proteins:
            return
        self._draw_col_headers()
        self._draw_row_labels()
        self._draw_matrix()

    def _draw_matrix(self) -> None:
        c = self._mx
        w = max(1, c.winfo_width())
        h = max(1, c.winfo_height())
        ox = int(c.canvasx(0))
        oy = int(c.canvasy(0))

        cp = self._cell_px
        n = len(self._proteins)
        split = self._split
        has_sep = 0 < split < n
        c0, c1 = self._visible_cols()
        r0, r1 = self._visible_rows()
        proteins = self._proteins
        sparse = self._sparse
        thresholds = self._thresholds

        # ------------------------------------------------------------------
        # pygame SDL2 path — GPU blitting, no PIL readback
        # ------------------------------------------------------------------
        if _pygame is not None and hasattr(self, "_pg_screen"):
            if (w, h) != self._pg_size:
                self._pg_screen = _pygame.display.set_mode((w, h), 0, 32)
                self._pg_size = (w, h)

            screen = self._pg_screen
            atlas = self._pg_atlas
            atlas_np_pg = self._atlas   # numpy atlas for colorization of Unicode glyphs
            score_cache = self._pg_score_cache
            pg_font = self._pg_font
            screen.fill(_BG_PG)

            cmap_lut = self._cmap_lut
            for ri in range(r0, r1):
                pi = proteins[ri]
                row_y = ri * cp - oy
                for ci in range(c0, c1):
                    pj = proteins[ci]
                    score: Optional[float] = sparse.get((min(pi, pj), max(pi, pj)))
                    if score is None:
                        continue
                    label = _cell_label(score, cp, thresholds)
                    if cmap_lut is not None:
                        bg = _score_to_bg(score, self._min_score, self._max_score, cmap_lut)
                        fg = _auto_fg(bg)
                        cache_key = (label, fg, bg)
                        glyph = score_cache.get(cache_key)
                        if glyph is None:
                            np_glyph = atlas_np_pg.get(label)
                            if np_glyph is not None:
                                # Colorize numpy glyph: white→bg, black→fg
                                mask = (np_glyph[:, :, 0] > 200)[:, :, np.newaxis]
                                colored = np.where(mask, bg, fg).astype(np.uint8)
                                cell_surf = _pygame.surfarray.make_surface(np.swapaxes(colored, 0, 1))
                            else:
                                # ASCII score string — safe to render with pg_font
                                text_surf = pg_font.render(label, True, fg, bg)
                                cell_surf = _pygame.Surface((cp, cp))
                                cell_surf.fill(bg)
                                cell_surf.blit(text_surf, text_surf.get_rect(center=(cp // 2, cp // 2)))
                            score_cache[cache_key] = glyph = cell_surf
                        screen.blit(glyph, (ci * cp - ox, row_y))
                    else:
                        glyph = atlas.get(label)
                        if glyph is None:
                            glyph = score_cache.get(label)
                            if glyph is None:
                                text_surf = pg_font.render(label, True, _FG_PG, _BG_PG)
                                cell_surf = _pygame.Surface((cp, cp))
                                cell_surf.fill(_BG_PG)
                                r = text_surf.get_rect(center=(cp // 2, cp // 2))
                                cell_surf.blit(text_surf, r)
                                glyph = cell_surf
                                score_cache[label] = glyph
                        screen.blit(glyph, (ci * cp - ox, row_y))

            if has_sep:
                sx = split * cp - ox
                _pygame.draw.line(screen, _SEP_PG, (sx, 0), (sx, h))
                sy = split * cp - oy
                if 0 <= sy <= h:
                    _pygame.draw.line(screen, _SEP_PG, (0, sy), (w, sy))

            # Highlight selected cell
            if self._selected is not None:
                sel_r, sel_c = self._selected
                hx = sel_c * cp - ox
                hy = sel_r * cp - oy
                lw = max(2, cp // 8)
                _pygame.draw.rect(screen, _HL, (hx, hy, cp, cp), lw)

            # Contour overlay
            for polyline in self._contour_segs:
                pts = [(int(x * cp - ox), int(y * cp - oy)) for x, y in polyline]
                if len(pts) >= 2:
                    _pygame.draw.lines(screen, _CONTOUR_COLOR, False, pts, 1)

            # Process click events before flip
            for ev in _pygame.event.get():
                if ev.type == _pygame.MOUSEBUTTONDOWN and ev.button == 1:
                    col_c = (ev.pos[0] + ox) // cp
                    row_c = (ev.pos[1] + oy) // cp
                    if 0 <= row_c < n and 0 <= col_c < n:
                        self.after(0, lambda r=row_c, c_=col_c: self._on_cell_select(r, c_))
                    else:
                        self.after(0, self._clear_selection)

            _pygame.display.flip()
            return

        # ------------------------------------------------------------------
        # numpy/PIL fallback
        # ------------------------------------------------------------------
        atlas_np = self._atlas
        cmap_lut = self._cmap_lut
        arr = np.full((h, w, 3), 255, dtype=np.uint8)
        score_texts: list[tuple] = []  # (cx, cy, label[, fg])

        for ri in range(r0, r1):
            pi = proteins[ri]
            row_y = ri * cp - oy
            for ci in range(c0, c1):
                pj = proteins[ci]
                score = sparse.get((min(pi, pj), max(pi, pj)))
                if score is None:
                    continue
                label = _cell_label(score, cp, thresholds)
                col_x = ci * cp - ox
                if cmap_lut is not None:
                    bg = _score_to_bg(score, self._min_score, self._max_score, cmap_lut)
                    fg = _auto_fg(bg)
                    glyph = atlas_np.get(label)
                    iy0 = max(0, row_y); iy1 = min(h, row_y + cp)
                    ix0 = max(0, col_x); ix1 = min(w, col_x + cp)
                    if iy1 > iy0 and ix1 > ix0:
                        if glyph is not None:
                            # Colorize numpy glyph: white→bg, black→fg
                            gy0 = max(0, -row_y); gx0 = max(0, -col_x)
                            gy1 = gy0 + (iy1 - iy0); gx1 = gx0 + (ix1 - ix0)
                            sl = glyph[gy0:gy1, gx0:gx1]
                            mask = (sl[:, :, 0] > 200)[:, :, np.newaxis]
                            arr[iy0:iy1, ix0:ix1] = np.where(mask, bg, fg).astype(np.uint8)
                        else:
                            arr[iy0:iy1, ix0:ix1] = bg
                            score_texts.append((col_x + cp // 2, row_y + cp // 2, label, fg))
                else:
                    glyph = atlas_np.get(label)
                    if glyph is not None:
                        gy0 = max(0, -row_y);  iy0 = max(0, row_y)
                        gx0 = max(0, -col_x);  ix0 = max(0, col_x)
                        gy1 = min(cp, h - row_y); iy1 = iy0 + (gy1 - gy0)
                        gx1 = min(cp, w - col_x); ix1 = ix0 + (gx1 - gx0)
                        if gy1 > gy0 and gx1 > gx0:
                            arr[iy0:iy1, ix0:ix1] = glyph[gy0:gy1, gx0:gx1]
                    else:
                        cx = col_x + cp // 2
                        cy = row_y + cp // 2
                        if 0 <= cx <= w and 0 <= cy <= h:
                            score_texts.append((cx, cy, label))

        img = Image.fromarray(arr)
        draw = ImageDraw.Draw(img)

        for entry in score_texts:
            fg_color = entry[3] if len(entry) > 3 else _FG
            draw.text((entry[0], entry[1]), entry[2], font=self._pil_cf, anchor="mm", fill=fg_color)

        if has_sep:
            sx = split * cp - ox
            draw.line([(sx, 0), (sx, h)], fill=_SEP, width=1)
            sy = split * cp - oy
            if 0 <= sy <= h:
                draw.line([(0, sy), (w, sy)], fill=_SEP, width=1)

        if self._selected is not None:
            sel_r, sel_c = self._selected
            hx = sel_c * cp - ox
            hy = sel_r * cp - oy
            lw = max(2, cp // 8)
            draw.rectangle([hx, hy, hx + cp - 1, hy + cp - 1], outline=_HL, width=lw)

        for polyline in self._contour_segs:
            pts = [(int(x * cp - ox), int(y * cp - oy)) for x, y in polyline]
            if len(pts) >= 2:
                draw.line(pts, fill=_CONTOUR_COLOR, width=1)

        self._mx_photo = ImageTk.PhotoImage(img)
        c.delete("all")
        c.create_image(ox, oy, image=self._mx_photo, anchor="nw")

    def _draw_col_headers(self) -> None:
        c = self._ch_canvas
        w = max(1, c.winfo_width())
        h = max(1, self._header_h)
        ox = int(c.canvasx(0))

        cp = self._cell_px
        lh = self._ch
        n = len(self._proteins)
        split = self._split
        has_sep = 0 < split < n
        c0, c1 = self._visible_cols()
        pil_lf = self._pil_lf

        img = Image.new("RGB", (w, h), _BG)
        draw = ImageDraw.Draw(img)
        y = 0

        # Section label row
        if self._section is not None:
            i = 0
            while i < n:
                label = self._section.get(self._proteins[i], "")
                j = i + 1
                while j < n:
                    if has_sep and j == split:
                        break
                    if self._section.get(self._proteins[j], "") != label:
                        break
                    j += 1
                x0s = i * cp - ox
                x1s = j * cp - ox
                if x1s > 0 and x0s < w and label:
                    mid = (x0s + x1s) // 2
                    draw.text((mid, y + lh // 2), label, font=pil_lf, anchor="mm", fill=_FG)
                    draw.rectangle([x0s, y, x1s - 1, y + lh - 1], outline=_SEP)
                i = j
            y += lh

        # Vertical separator in header
        if has_sep:
            sx = split * cp - ox
            draw.line([(sx, 0), (sx, h)], fill=_SEP, width=1)

        # Stacked chars
        labels = [p[:COL_W] for p in self._proteins]
        max_len = max((len(lb) for lb in labels), default=0)
        for char_idx in range(max_len):
            for ci in range(c0, c1):
                lb = labels[ci]
                ch = lb[char_idx] if char_idx < len(lb) else " "
                if ch != " ":
                    cx = ci * cp + cp // 2 - ox
                    draw.text((cx, y + lh // 2), ch, font=pil_lf, anchor="mm", fill=_FG)
            y += lh

        # Horizontal rule
        draw.line([(0, y), (w, y)], fill=_SEP, width=1)
        if has_sep:
            sx = split * cp - ox
            draw.text((sx, y), "+", font=pil_lf, anchor="mm", fill=_SEP)

        self._ch_photo = ImageTk.PhotoImage(img)
        c.delete("all")
        c.create_image(ox, 0, image=self._ch_photo, anchor="nw")

    def _draw_row_labels(self) -> None:
        c = self._rl_canvas
        w = max(1, self._label_w)
        h = max(1, c.winfo_height())
        oy = int(c.canvasy(0))

        cp = self._cell_px
        n = len(self._proteins)
        split = self._split
        has_sep = 0 < split < n
        r0, r1 = self._visible_rows()
        pil_lf = self._pil_lf
        lh = self._ch

        img = Image.new("RGB", (w, h), _BG)
        draw = ImageDraw.Draw(img)

        for ri in range(r0, r1):
            label = self._proteins[ri][:ROW_W]
            cy = ri * cp + cp // 2 - oy
            draw.text((2, cy), label, font=pil_lf, anchor="lm", fill=_FG)

        if has_sep:
            sy = split * cp - oy
            if 0 <= sy <= h:
                draw.line([(0, sy), (w, sy)], fill=_SEP, width=1)

        self._rl_photo = ImageTk.PhotoImage(img)
        c.delete("all")
        c.create_image(0, oy, image=self._rl_photo, anchor="nw")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    initial = sys.argv[1] if len(sys.argv) > 1 else None
    MatrixApp(initial).mainloop()


if __name__ == "__main__":
    main()

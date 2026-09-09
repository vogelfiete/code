"""GUI matrix viewer for XL-MS crosslink data. See docs/usage.md for details."""
from __future__ import annotations

import math
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
    _MPL = True
except Exception:
    _MPL = False

_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_DIR, "..", "scripts"))
sys.path.insert(0, os.path.join(_DIR, "..", "src"))

from print_matrix import (  # type: ignore[import]
    build_matrix,
    build_residue_matrix,
    ResidueId,
    _sort_proteins,
    _sort_residues,
    _add_unlinked_residues,
    _short_accession,
    _build_thresholds,
    _score_to_dot,
    DEFAULT_N,
    COL_W,
    ROW_W,
)
from xlms.io import read_fasta  # type: ignore[import]

CELL_PX_MIN = 4
CELL_PX_MAX = 32
CELL_PX_DEFAULT = 10
SCORE_THRESHOLD_PX = 18  # cell_px at which scores replace dot symbols

BLOCK_ZOOM_STEPS = 20     # extra slider ticks reserved for block aggregation, below CELL_PX_MIN
BLOCK_ZOOM_GROWTH = 1.3   # per-tick multiplicative growth of block size (log-zoom feel)
BLOCK_SIZE_MAX = 500      # hard safety cap regardless of the formula above
_AGG_METHOD_LABELS = ["Mean", "Geometric Mean", "Max Score"]
_AGG_METHOD_MAP = {"Mean": "mean", "Geometric Mean": "geomean", "Max Score": "max"}

LABEL_FONT_SIZE = 9
SEP_COLOR = "#999999"
BG = "white"

_BG = (255, 255, 255)
_FG = (0, 0, 0)
_SEP = (153, 153, 153)
_BG_PG  = (255, 255, 255)
_FG_PG  = (0, 0, 0)
_SEP_PG = (153, 153, 153)
_HL = (0, 120, 215)  # selection highlight, shared between the PIL and pygame paths

_MONO_FONTS = [
    "C:/Windows/Fonts/consola.ttf",
    "C:/Windows/Fonts/cour.ttf",
    "C:/Windows/Fonts/lucon.ttf",
]

_GLYPH_CHARS = [*"·○●⬤", *"0123456789"]


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
    score: float,
    min_s: float,
    max_s: float,
    lut: np.ndarray,
    score_sorted: np.ndarray | None = None,
) -> tuple[int, int, int]:
    """Map a score to an RGB background colour via LUT lookup."""
    if score_sorted is not None and len(score_sorted) > 0:
        idx = int(np.searchsorted(score_sorted, score) / len(score_sorted) * 255)
    elif max_s <= min_s:
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


class MatrixApp(tk.Tk):
    def __init__(self, initial_csv: str | None = None) -> None:
        super().__init__()
        self.title("XL-MS Matrix Viewer")
        self.geometry("1200x800")
        self.minsize(600, 400)

        # Data — self._proteins holds either protein names (str) or
        # ResidueId(protein, pos), depending on self._level
        self._level: str = "protein"
        self._proteins: list = []
        self._sparse: dict[tuple, float] = {}
        self._decoy: dict = {}
        self._section: dict | None = None
        self._thresholds: list[float] = []
        self._split: int = 0
        self._col_labels: list[str] = []
        self._row_labels: list[str] = []

        # Render state
        self._cell_px: int = CELL_PX_DEFAULT
        self._pending: bool = False
        self._selected: tuple[int, int] | None = None

        # Aggregate-zoom / blocks — self._blocks is a drop-in replacement for
        # self._proteins in geometry code; at block_size==1 it's just
        # [(0,1),(1,2),...], so geometry math stays pixel-identical to before
        self._block_size: int = 1
        self._blocks: list[tuple[int, int]] = []
        self._item_to_block: list[int] = []
        self._split_block: int = 0
        self._item_index: dict = {}
        self._agg_method: str = "mean"
        self._block_pair_raw: dict[tuple[int, int], list[float]] = {}
        self._block_scores: dict[tuple[int, int], float] = {}
        self._block_pair_n: dict[tuple[int, int], int] = {}

        # Heatmap colormap
        self._cmap_name: str = "viridis"
        self._cmap_lut: np.ndarray | None = None  # None = disabled
        self._min_score: float = 0.0
        self._max_score: float = 1.0
        self._score_sorted: np.ndarray = np.array([])
        self._norm_mode: str = "linear"

        # UniProt lookup cache and stale-click guard
        self._uniprot_cache: dict[str, tuple[str, int]] = {}
        self._sel_token: int = 0

        # _norm_var is tk.Var, created in _build_controls

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
        tk.Label(row2, text="Level:").pack(side=tk.LEFT)
        self._level_var = tk.StringVar(value="protein")
        ttk.Combobox(
            row2, textvariable=self._level_var, width=8, state="readonly",
            values=["protein", "residue"],
        ).pack(side=tk.LEFT, padx=2)
        self._only_linked_var = tk.BooleanVar(value=True)
        self._only_linked_check = tk.Checkbutton(
            row2, text="Show only linked residues", variable=self._only_linked_var,
        )
        self._only_linked_check.pack(side=tk.LEFT, padx=(4, 0))
        self._level_var.trace_add("write", lambda *_: self._only_linked_check.configure(
            state=tk.NORMAL if self._level_var.get() == "residue" else tk.DISABLED
        ))
        self._only_linked_check.configure(state=tk.DISABLED)
        self._n_label_var = tk.StringVar(value="  N (proteins):")
        self._level_var.trace_add("write", lambda *_: self._n_label_var.set(
            f"  N ({self._level_var.get()}s):"
        ))
        tk.Label(row2, textvariable=self._n_label_var).pack(side=tk.LEFT)
        self._n_var = tk.IntVar(value=DEFAULT_N)
        tk.Spinbox(row2, textvariable=self._n_var, from_=10, to=5000, width=6).pack(side=tk.LEFT, padx=2)
        tk.Label(row2, text="  Order:").pack(side=tk.LEFT)
        self._order_var = tk.StringVar(value="confidence")
        ttk.Combobox(
            row2, textvariable=self._order_var, width=12, state="readonly",
            values=["confidence", "alpha", "cluster", "sequence", "pathway", "complex", "size"],
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
        _cmap_choices = (["(none)", "coolwarm"] if _MPL else ["(none)"])
        ttk.Combobox(
            row3, textvariable=self._cmap_var, width=10,
            state="readonly", values=_cmap_choices,
        ).pack(side=tk.LEFT, padx=2)
        self._cmap_var.trace_add("write", lambda *_: self._on_colormap_change())

        tk.Label(row3, text="   Norm:").pack(side=tk.LEFT)
        self._norm_var = tk.StringVar(value="linear")
        for _lbl, _val in [("Linear", "linear"), ("Quantile", "quantile")]:
            tk.Radiobutton(
                row3, text=_lbl, variable=self._norm_var, value=_val,
                command=self._on_norm_change,
            ).pack(side=tk.LEFT)

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
            row1, from_=CELL_PX_MIN - BLOCK_ZOOM_STEPS, to=CELL_PX_MAX,
            orient=tk.HORIZONTAL, length=200, showvalue=False,
            command=self._on_zoom,
        )
        self._zoom.set(CELL_PX_DEFAULT)
        self._zoom.pack(side=tk.LEFT, padx=4)
        self._zoom_label_var = tk.StringVar(value=f"{CELL_PX_DEFAULT}px")
        tk.Label(row1, textvariable=self._zoom_label_var, width=14, anchor="w").pack(side=tk.LEFT)

        tk.Label(row1, text="  Aggregate:").pack(side=tk.LEFT, padx=(8, 0))
        self._agg_method_var = tk.StringVar(value=_AGG_METHOD_LABELS[0])
        ttk.Combobox(
            row1, textvariable=self._agg_method_var, width=13, state="readonly",
            values=_AGG_METHOD_LABELS,
        ).pack(side=tk.LEFT, padx=2)
        self._agg_method_var.trace_add("write", lambda *_: self._on_agg_method_change())

        self._info_var = tk.StringVar(value="")
        tk.Label(row1, textvariable=self._info_var, fg="gray").pack(side=tk.LEFT, padx=8)
        tk.Label(row1, text=f"renderer: {_RENDERER_LABEL}", fg="#aaaaaa").pack(side=tk.RIGHT, padx=8)

        self._sel_info_var = tk.StringVar(value="")
        tk.Label(bar, textvariable=self._sel_info_var, anchor="w",
                 fg="#1a6fb5", font=("TkFixedFont", 9)).pack(fill=tk.X, padx=8, pady=(0, 2))

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
        s = len(self._blocks) * self._cell_px
        return s, s

    def _update_scrollregion(self) -> None:
        vw, vh = self._virtual()
        self._mx.configure(scrollregion=(0, 0, vw, vh))
        self._ch_canvas.configure(scrollregion=(0, 0, vw, self._header_h))
        self._rl_canvas.configure(scrollregion=(0, 0, self._label_w, vh))

    def _yscroll(self, *args) -> None:
        self._mx.yview(*args)
        self._rl_canvas.yview(*args)
        self._schedule(fast=True)

    def _xscroll(self, *args) -> None:
        self._mx.xview(*args)
        self._ch_canvas.xview(*args)
        self._schedule(fast=True)

    def _wheel(self, event: tk.Event) -> None:
        delta = 3 if (event.num == 5 or event.delta < 0) else -3
        self._mx.yview_scroll(delta, "units")
        self._rl_canvas.yview_scroll(delta, "units")
        self._schedule(fast=True)

    def _arrow(self, dx: int, dy: int) -> None:
        if dy:
            self._mx.yview_scroll(dy, "units")
            self._rl_canvas.yview_scroll(dy, "units")
        if dx:
            self._mx.xview_scroll(dx, "units")
            self._ch_canvas.xview_scroll(dx, "units")
        self._schedule(fast=True)

    def _on_zoom(self, value: str) -> None:
        v = int(float(value))
        if v >= CELL_PX_MIN:
            self._cell_px = v
            new_block_size = 1
            self._zoom_label_var.set(f"{v}px")
        else:
            self._cell_px = CELL_PX_MIN
            steps = CELL_PX_MIN - v
            new_block_size = min(BLOCK_SIZE_MAX, max(2, round(BLOCK_ZOOM_GROWTH ** steps)))
            self._zoom_label_var.set(f"≤{new_block_size} → {CELL_PX_MIN}px")
        self._update_cell_font()
        if new_block_size != self._block_size:
            self._block_size = new_block_size
            self._recompute_blocks()
            if hasattr(self, "_pg_score_cache"):
                self._pg_score_cache = {}
        self._update_scrollregion()
        self._schedule()

    def _on_agg_method_change(self) -> None:
        self._agg_method = _AGG_METHOD_MAP.get(self._agg_method_var.get(), "mean")
        self._reduce_block_scores()
        if hasattr(self, "_pg_score_cache"):
            self._pg_score_cache = {}
        self._schedule()

    def _recompute_blocks(self) -> None:
        n = len(self._proteins)
        proteins, section, split, bs = self._proteins, self._section, self._split, self._block_size
        if n == 0:
            self._blocks, self._item_to_block, self._split_block = [], [], 0
        else:
            if bs <= 1:
                blocks = [(i, i + 1) for i in range(n)]
            else:
                has_sep = 0 < split < n
                blocks = []
                i = 0
                while i < n:
                    sec0 = section.get(proteins[i], "") if section is not None else None
                    limit = min(n, i + bs)
                    j = i + 1
                    while j < limit:
                        if has_sep and j == split:
                            break
                        if section is not None and section.get(proteins[j], "") != sec0:
                            break
                        j += 1
                    blocks.append((i, j))
                    i = j
            item_to_block = [0] * n
            for bi, (s, e) in enumerate(blocks):
                item_to_block[s:e] = [bi] * (e - s)
            self._blocks = blocks
            self._item_to_block = item_to_block
            self._split_block = item_to_block[split] if 0 <= split < n else len(blocks)
        self._recompute_block_buckets()
        self._reduce_block_scores()

    def _recompute_block_buckets(self) -> None:
        buckets: dict[tuple[int, int], list[float]] = {}
        if self._block_size > 1 and self._blocks:
            item_to_block, idx = self._item_to_block, self._item_index
            for (pi, pj), score in self._sparse.items():
                bi, bj = idx.get(pi), idx.get(pj)
                if bi is None or bj is None:
                    continue  # defensive; sparse keys are always in self._proteins
                bi, bj = item_to_block[bi], item_to_block[bj]
                key = (bi, bj) if bi <= bj else (bj, bi)
                buckets.setdefault(key, []).append(score)
        self._block_pair_raw = buckets

    def _reduce_block_scores(self) -> None:
        method = self._agg_method
        scores: dict[tuple[int, int], float] = {}
        counts: dict[tuple[int, int], int] = {}
        for key, vals in self._block_pair_raw.items():
            counts[key] = len(vals)
            if method == "max":
                scores[key] = max(vals)
            elif method == "geomean":
                pos = [v for v in vals if v > 0]
                if pos:
                    scores[key] = math.exp(sum(math.log(v) for v in pos) / len(pos))
                else:
                    # No positive scores contribute to this block-pair. Geometric
                    # mean is undefined for non-positive inputs; fall back to the
                    # arithmetic mean rather than dropping the pair (which would
                    # wrongly render as "no crosslink" and hide real data).
                    scores[key] = sum(vals) / len(vals)
            else:  # "mean"
                scores[key] = sum(vals) / len(vals)
        self._block_scores = scores
        self._block_pair_n = counts

    def _on_click(self, event: tk.Event) -> None:
        """Handle left-click on the matrix canvas (numpy path fallback)."""
        if not self._proteins or not self._blocks:
            return
        ox = int(self._mx.canvasx(0))
        oy = int(self._mx.canvasy(0))
        cp = self._cell_px
        cb = (event.x + ox) // cp
        rb = (event.y + oy) // cp
        n = len(self._blocks)
        if 0 <= rb < n and 0 <= cb < n:
            self._on_block_select(rb, cb)
        else:
            self._selected = None
            self._sel_info_var.set("")
        self._schedule()

    def _item_protein(self, item) -> str:
        """The parent protein name for a protein-mode or residue-mode axis item."""
        return item.protein if self._level == "residue" else item

    def _on_cell_select(self, row: int, col: int) -> None:
        self._selected = (row, col)
        self._sel_token += 1
        token = self._sel_token
        pi = self._proteins[row]
        pj = self._proteins[col]
        key = (min(pi, pj), max(pi, pj))
        score = self._sparse.get(key)
        score_line = f"score: {score:.4f}" if score is not None else "(no crosslink)"

        def _fmt_local(item) -> str:
            protein = self._item_protein(item)
            tag = "decoy" if self._decoy.get(item) else "target"
            aa = f"{len(self._seqs[protein])} aa" if protein in self._seqs else "fetching…"
            sec = (self._section.get(item, "") if self._section else "")
            sec_part = f"  [{sec}]" if sec else ""
            return f"{item}  ({aa}, {tag}){sec_part}"

        self._sel_info_var.set(f"{_fmt_local(pi)}\n{_fmt_local(pj)}\n{score_line}")
        threading.Thread(
            target=self._fetch_uniprot_info, args=(pi, pj, token), daemon=True
        ).start()

    def _fetch_uniprot_accession(self, accession: str) -> tuple[str, int] | None:
        import urllib.request
        import json as _json

        if accession in self._uniprot_cache:
            return self._uniprot_cache[accession]
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "xlms-matrix/1.0"})
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = _json.loads(resp.read().decode())
            name = (data.get("proteinDescription", {})
                        .get("recommendedName", {})
                        .get("fullName", {})
                        .get("value", accession))
            length = int(data.get("sequence", {}).get("length", 0))
            result: tuple[str, int] = (name, length)
            self._uniprot_cache[accession] = result
            return result
        except Exception:
            return None

    def _fetch_uniprot_info(self, pi, pj, token: int) -> None:
        acc_i = _short_accession(self._item_protein(pi))
        acc_j = _short_accession(self._item_protein(pj))
        ri = self._fetch_uniprot_accession(acc_i)
        rj = self._fetch_uniprot_accession(acc_j)

        def _update() -> None:
            if self._sel_token != token:
                return
            key = (min(pi, pj), max(pi, pj))
            score = self._sparse.get(key)
            score_line = f"score: {score:.4f}" if score is not None else "(no crosslink)"

            def _line(item, acc: str, r: tuple[str, int] | None) -> str:
                tag = "decoy" if self._decoy.get(item) else "target"
                sec = (self._section.get(item, "") if self._section else "")
                sec_part = f"  [{sec}]" if sec else ""
                pos_part = f"  pos {item.pos}" if self._level == "residue" else ""
                if r:
                    return f"{r[0]}  ({r[1]} aa, {tag}){sec_part}{pos_part}  [{acc}]"
                return f"{item}  (?, {tag}){sec_part}"

            self._sel_info_var.set(
                f"{_line(pi, acc_i, ri)}\n{_line(pj, acc_j, rj)}\n{score_line}"
            )

        self.after(0, _update)

    def _agg_method_label(self) -> str:
        return {"mean": "Mean", "geomean": "Geometric mean", "max": "Max"}.get(self._agg_method, "Mean")

    def _block_score_line(self, rb: int, cb: int) -> str:
        key = (rb, cb) if rb <= cb else (cb, rb)
        score = self._block_scores.get(key)
        if score is None:
            return "(no crosslinks between these blocks)"
        n_pairs = self._block_pair_n.get(key, 0)
        return f"{self._agg_method_label()} of {n_pairs} crosslink(s): {score:.4f}"

    def _on_block_select(self, rb: int, cb: int) -> None:
        self._selected = (rb, cb)
        r_lo, r_hi = self._blocks[rb]
        c_lo, c_hi = self._blocks[cb]
        if r_hi - r_lo == 1 and c_hi - c_lo == 1:
            self._on_cell_select(r_lo, c_lo)
            return

        self._sel_token += 1
        token = self._sel_token
        row_items = self._proteins[r_lo:r_hi]
        col_items = self._proteins[c_lo:c_hi]
        score_line = self._block_score_line(rb, cb)

        def _fmt_block(items) -> str:
            decoys = {self._decoy.get(it, False) for it in items}
            tag = "target" if decoys == {False} else "decoy" if decoys == {True} else "mixed"
            secs = {self._section.get(it, "") for it in items} if self._section else set()
            sec = next(iter(secs)) if len(secs) == 1 else None
            sec_part = f"  [{sec}]" if sec else ""
            return f"{len(items)} items ({tag}){sec_part}"

        self._sel_info_var.set(f"{_fmt_block(row_items)}\n{_fmt_block(col_items)}\n{score_line}")

        row_proteins = {self._item_protein(it) for it in row_items}
        col_proteins = {self._item_protein(it) for it in col_items}
        if len(row_proteins) == 1 and len(col_proteins) == 1:
            threading.Thread(
                target=self._fetch_uniprot_info_block,
                args=(row_items, col_items, token),
                daemon=True,
            ).start()

    def _fetch_uniprot_info_block(self, row_items: list, col_items: list, token: int) -> None:
        acc_i = _short_accession(self._item_protein(row_items[0]))
        acc_j = _short_accession(self._item_protein(col_items[0]))
        ri = self._fetch_uniprot_accession(acc_i)
        rj = self._fetch_uniprot_accession(acc_j)

        def _pos_range(items: list) -> str:
            if self._level != "residue":
                return ""
            positions = sorted(it.pos for it in items)
            if len(positions) == 1:
                return f"  pos {positions[0]}"
            return f"  pos {positions[0]}-{positions[-1]}"

        def _update() -> None:
            if self._sel_token != token or self._selected is None:
                return
            rb, cb = self._selected
            score_line = self._block_score_line(rb, cb)

            def _line(items: list, acc: str, r: tuple[str, int] | None) -> str:
                decoys = {self._decoy.get(it, False) for it in items}
                tag = "target" if decoys == {False} else "decoy" if decoys == {True} else "mixed"
                sec = self._section.get(items[0], "") if self._section else ""
                sec_part = f"  [{sec}]" if sec else ""
                pos_part = _pos_range(items)
                if r:
                    return f"{r[0]}  ({r[1]} aa, {tag}){sec_part}{pos_part}  [{acc}]"
                return f"{len(items)} items  (?, {tag}){sec_part}"

            self._sel_info_var.set(
                f"{_line(row_items, acc_i, ri)}\n{_line(col_items, acc_j, rj)}\n{score_line}"
            )

        self.after(0, _update)

    def _clear_selection(self) -> None:
        self._selected = None
        self._sel_info_var.set("")

    def _schedule(self, fast: bool = False) -> None:
        if not self._pending:
            self._pending = True
            self.after(0, self._redraw)

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
        level = self._level_var.get()
        order = self._order_var.get()
        fasta = self._fasta_var.get().strip() or None
        only_linked = self._only_linked_var.get()
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
                        items, decoy = _add_unlinked_residues(items, decoy, seqs or {})
                    items, section = _sort_residues(items, sparse, decoy, order, seqs, species)
                else:
                    items, sparse, decoy = build_matrix(csv, n=n)
                    items, section = _sort_proteins(items, sparse, decoy, order, seqs, species)
                self.after(0, lambda: self._on_loaded(items, sparse, decoy, section, seqs, level))
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

    def _on_norm_change(self) -> None:
        self._norm_mode = self._norm_var.get()
        if hasattr(self, "_pg_score_cache"):
            self._pg_score_cache = {}
        self._schedule()

    def _on_loaded(
        self,
        proteins: list,
        sparse: dict[tuple, float],
        decoy: dict,
        section: dict | None,
        seqs: dict[str, str] | None = None,
        level: str = "protein",
    ) -> None:
        self._level = level
        self._proteins = proteins
        self._sparse = sparse
        self._decoy = decoy
        self._section = section
        self._seqs: dict[str, str] = seqs or {}
        if level == "residue":
            self._col_labels = [str(r.pos)[:COL_W] for r in proteins]
            self._row_labels = [f"{_short_accession(r.protein)}:{r.pos}"[:ROW_W] for r in proteins]
        else:
            self._col_labels = [p[:COL_W] for p in proteins]
            self._row_labels = [p[:ROW_W] for p in proteins]
        self._thresholds = _build_thresholds(sparse) if sparse else []
        self._split = next(
            (i for i, p in enumerate(proteins) if decoy.get(p, False)), len(proteins)
        )
        self._min_score = min(sparse.values()) if sparse else 0.0
        self._max_score = max(sparse.values()) if sparse else 1.0
        self._score_sorted = np.sort(list(sparse.values())) if sparse else np.array([])
        if self._cmap_lut is not None:
            self._cmap_lut = _make_lut(self._cmap_name)
        self._item_index = {p: i for i, p in enumerate(proteins)}
        self._recompute_blocks()
        n_t = self._split
        n_d = len(proteins) - n_t
        unit = "residues" if level == "residue" else "proteins"
        self._status_var.set(f"Loaded — {n_t} targets · {n_d} decoys · {len(sparse)} links")
        self._info_var.set(f"{len(proteins)} {unit}")
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

        # column headers
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

        labels_col = self._col_labels
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

        # row labels
        for ri in range(n):
            cy = hh + ri * cp + cp // 2
            draw.text((2, cy), self._row_labels[ri], font=font_lbl, anchor="lm", fill=_FG)
        if has_sep:
            draw.line([(0, hh + split * cp), (lw, hh + split * cp)], fill=_SEP, width=1)

        if has_sep:
            sx = lw + split * cp
            draw.line([(sx, 0), (sx, total_h)], fill=_SEP, width=1)
            sy = hh + split * cp
            draw.line([(lw, sy), (total_w, sy)], fill=_SEP, width=1)

        cmap_lut = self._cmap_lut
        for ri in range(n):
            pi = proteins[ri]
            for ci in range(n):
                pj = proteins[ci]
                score = sparse.get((min(pi, pj), max(pi, pj)))
                if score is None:
                    continue
                if cmap_lut is not None:
                    bg = _score_to_bg(score, self._min_score, self._max_score, cmap_lut,
                                      score_sorted=self._score_sorted if self._norm_mode == "quantile" else None)
                    x0c, y0c = lw + ci * cp, hh + ri * cp
                    draw.rectangle([x0c, y0c, x0c + cp - 1, y0c + cp - 1], fill=bg)
                    if cp >= SCORE_THRESHOLD_PX:
                        draw.text(
                            (x0c + cp // 2, y0c + cp // 2),
                            f"{score:.2f}", font=font_cell, anchor="mm", fill=_auto_fg(bg),
                        )
                else:  # greyscale default — no dot symbols
                    mn, mx = self._min_score, self._max_score
                    grey_val = int(200 - (score - mn) / (mx - mn) * 160) if mx > mn else 120
                    grey_val = max(40, min(200, grey_val))
                    bg = (grey_val, grey_val, grey_val)
                    x0c, y0c = lw + ci * cp, hh + ri * cp
                    draw.rectangle([x0c, y0c, x0c + cp - 1, y0c + cp - 1], fill=bg)
                    if cp >= SCORE_THRESHOLD_PX:
                        draw.text(
                            (x0c + cp // 2, y0c + cp // 2),
                            f"{score:.2f}", font=font_cell, anchor="mm", fill=_auto_fg(bg),
                        )

        return img

    def _visible_col_blocks(self) -> tuple[int, int]:
        n = len(self._blocks)
        cp = self._cell_px
        x0 = self._mx.canvasx(0)
        x1 = self._mx.canvasx(self._mx.winfo_width())
        return max(0, int(x0 // cp)), min(n, int(x1 // cp) + 1)

    def _visible_row_blocks(self) -> tuple[int, int]:
        n = len(self._blocks)
        cp = self._cell_px
        y0 = self._mx.canvasy(0)
        y1 = self._mx.canvasy(self._mx.winfo_height())
        return max(0, int(y0 // cp)), min(n, int(y1 // cp) + 1)

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
        n = len(self._blocks)
        split = self._split_block
        has_sep = 0 < split < n
        c0, c1 = self._visible_col_blocks()
        r0, r1 = self._visible_row_blocks()
        proteins = self._proteins
        sparse = self._sparse
        thresholds = self._thresholds

        if self._block_size <= 1:
            def _score_at(rb: int, cb: int) -> Optional[float]:
                pi, pj = proteins[rb], proteins[cb]
                return sparse.get((min(pi, pj), max(pi, pj)))
        else:
            block_scores = self._block_scores

            def _score_at(rb: int, cb: int) -> Optional[float]:
                key = (rb, cb) if rb <= cb else (cb, rb)
                return block_scores.get(key)

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
            for rb in range(r0, r1):
                row_y = rb * cp - oy
                for cb in range(c0, c1):
                    score: Optional[float] = _score_at(rb, cb)
                    if score is None:
                        continue
                    if cmap_lut is not None:
                        bg = _score_to_bg(score, self._min_score, self._max_score, cmap_lut,
                                      score_sorted=self._score_sorted if self._norm_mode == "quantile" else None)
                        bg_col = (int(bg[0]), int(bg[1]), int(bg[2]))
                        show_text = cp >= SCORE_THRESHOLD_PX
                        score_str = f"{score:.2f}" if show_text else ""
                        cache_key = ("cm", bg_col, score_str)
                        cell_surf = score_cache.get(cache_key)
                        if cell_surf is None:
                            cell_surf = _pygame.Surface((cp, cp))
                            cell_surf.fill(bg_col)
                            if show_text:
                                fg_col = _auto_fg(bg_col)
                                ts = pg_font.render(score_str, True, fg_col, bg_col)
                                cell_surf.blit(ts, ts.get_rect(center=(cp // 2, cp // 2)))
                            score_cache[cache_key] = cell_surf
                        screen.blit(cell_surf, (cb * cp - ox, row_y))
                    else:  # greyscale default — no dot symbols
                        mn, mx = self._min_score, self._max_score
                        grey_val = int(200 - (score - mn) / (mx - mn) * 160) if mx > mn else 120
                        grey_val = max(40, min(200, grey_val))
                        bg_col = (grey_val, grey_val, grey_val)
                        show_text = cp >= SCORE_THRESHOLD_PX
                        score_str = f"{score:.2f}" if show_text else ""
                        cache_key = ("gs", grey_val, score_str)
                        cell_surf = score_cache.get(cache_key)
                        if cell_surf is None:
                            cell_surf = _pygame.Surface((cp, cp))
                            cell_surf.fill(bg_col)
                            if show_text:
                                fg_col = _auto_fg(bg_col)
                                ts = pg_font.render(score_str, True, fg_col, bg_col)
                                cell_surf.blit(ts, ts.get_rect(center=(cp // 2, cp // 2)))
                            score_cache[cache_key] = cell_surf
                        screen.blit(cell_surf, (cb * cp - ox, row_y))

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

            # Process click events before flip
            for ev in _pygame.event.get():
                if ev.type == _pygame.MOUSEBUTTONDOWN and ev.button == 1:
                    col_c = (ev.pos[0] + ox) // cp
                    row_c = (ev.pos[1] + oy) // cp
                    if 0 <= row_c < n and 0 <= col_c < n:
                        self.after(0, lambda r=row_c, c_=col_c: self._on_block_select(r, c_))
                    else:
                        self.after(0, self._clear_selection)

            _pygame.display.flip()
            return

        atlas_np = self._atlas
        cmap_lut = self._cmap_lut
        arr = np.full((h, w, 3), 255, dtype=np.uint8)
        score_texts: list[tuple] = []  # (cx, cy, label[, fg])

        for rb in range(r0, r1):
            row_y = rb * cp - oy
            for cb in range(c0, c1):
                score = _score_at(rb, cb)
                if score is None:
                    continue
                col_x = cb * cp - ox
                if cmap_lut is not None:
                    bg = _score_to_bg(score, self._min_score, self._max_score, cmap_lut,
                                      score_sorted=self._score_sorted if self._norm_mode == "quantile" else None)
                    iy0 = max(0, row_y); iy1 = min(h, row_y + cp)
                    ix0 = max(0, col_x); ix1 = min(w, col_x + cp)
                    if iy1 > iy0 and ix1 > ix0:
                        arr[iy0:iy1, ix0:ix1] = bg
                        if cp >= SCORE_THRESHOLD_PX:
                            score_texts.append((col_x + cp // 2, row_y + cp // 2,
                                                f"{score:.2f}", _auto_fg(bg)))
                else:  # greyscale default — no dot symbols
                    mn, mx = self._min_score, self._max_score
                    grey_val = int(200 - (score - mn) / (mx - mn) * 160) if mx > mn else 120
                    grey_val = max(40, min(200, grey_val))
                    bg = (grey_val, grey_val, grey_val)
                    iy0 = max(0, row_y); iy1 = min(h, row_y + cp)
                    ix0 = max(0, col_x); ix1 = min(w, col_x + cp)
                    if iy1 > iy0 and ix1 > ix0:
                        arr[iy0:iy1, ix0:ix1] = bg
                        if cp >= SCORE_THRESHOLD_PX:
                            score_texts.append((col_x + cp // 2, row_y + cp // 2,
                                                f"{score:.2f}", _auto_fg(bg)))

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
        n = len(self._blocks)
        split = self._split_block
        has_sep = 0 < split < n
        c0, c1 = self._visible_col_blocks()
        pil_lf = self._pil_lf
        proteins = self._proteins
        blocks = self._blocks

        img = Image.new("RGB", (w, h), _BG)
        draw = ImageDraw.Draw(img)
        y = 0

        # Section label row — walked over blocks; each block is guaranteed
        # single-section by construction, so adjacent same-label blocks
        # still merge into one wide, legible span exactly as before
        if self._section is not None:
            i = 0
            while i < n:
                label = self._section.get(proteins[blocks[i][0]], "")
                j = i + 1
                while j < n:
                    if has_sep and j == split:
                        break
                    if self._section.get(proteins[blocks[j][0]], "") != label:
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

        # Stacked chars — only for single-item blocks; a multi-item block's
        # per-char label would misleadingly suggest the whole block is one item
        labels = self._col_labels
        max_len = max((len(lb) for lb in labels), default=0)
        for char_idx in range(max_len):
            for cb in range(c0, c1):
                s, e = blocks[cb]
                if e - s != 1:
                    continue
                lb = labels[s]
                ch = lb[char_idx] if char_idx < len(lb) else " "
                if ch != " ":
                    cx = cb * cp + cp // 2 - ox
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
        n = len(self._blocks)
        split = self._split_block
        has_sep = 0 < split < n
        r0, r1 = self._visible_row_blocks()
        pil_lf = self._pil_lf
        lh = self._ch
        proteins = self._proteins
        blocks = self._blocks

        img = Image.new("RGB", (w, h), _BG)
        draw = ImageDraw.Draw(img)

        for rb in range(r0, r1):
            s, e = blocks[rb]
            cy = rb * cp + cp // 2 - oy
            if e - s == 1:
                label = self._row_labels[s]
            elif self._section is not None:
                label = self._section.get(proteins[s], "")
            else:
                label = ""
            if label:
                draw.text((2, cy), label, font=pil_lf, anchor="lm", fill=_FG)

        if has_sep:
            sy = split * cp - oy
            if 0 <= sy <= h:
                draw.line([(0, sy), (w, sy)], fill=_SEP, width=1)

        self._rl_photo = ImageTk.PhotoImage(img)
        c.delete("all")
        c.create_image(0, oy, image=self._rl_photo, anchor="nw")


def main() -> None:
    initial = sys.argv[1] if len(sys.argv) > 1 else None
    MatrixApp(initial).mainloop()


if __name__ == "__main__":
    main()

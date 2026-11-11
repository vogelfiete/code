"""A window for looking at XL-MS crosslink data as a matrix.
See docs/usage.md for how to use it.

The data comes from the same code the command line tool uses
(scripts/print_matrix.py). The matrix is drawn either through an embedded
pygame window, which is faster, or without pygame as a plain numpy/PIL image.
The axes are a tree: complexes hold proteins, proteins hold residues, and
clicking the axis labels opens or closes them. When you zoom out far
enough, neighbouring rows and columns are merged into blocks and their
scores are combined."""
from __future__ import annotations

import json
import math
import os
import sys
import threading
import traceback
import warnings
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

from fetch_features import FEATURE_TYPES, features_for  # type: ignore[import]
from print_matrix import (  # type: ignore[import]
    COL_W,
    ORDER_MODES,
    ROW_W,
    add_unlinked_residues,
    build_matrix,
    build_residue_matrix,
    complex_hierarchy,
    find_sequence,
    seq_index,
    short_accession,
    sort_proteins,
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
BLOCK_ZOOM_GROWTH = 1.12  # each step makes blocks this much bigger (feels like a log zoom)
ZOOM_OUT_FIT_PX = 200     # big matrices get more block steps, until the whole matrix fits into this many pixels
OVERVIEW_PAD_PX = 12      # zoomed out into blocks: empty space around the matrix
MIN_VISIBLE_PX = 40       # the view can be moved past the matrix edges until only this much of it is left
AXIS_TICK_PX = 4          # length of a position tick mark on the residue axis
AXIS_TICK_GAP_PX = 6      # minimum empty space kept between two tick labels
_AGG_METHODS = {"Geometric Mean": "geomean", "Max Score": "max",
                "Root Sum of Squares": "rss", "Count": "count"}                    # dropdown -> key
_QUADRANTS = ["TT", "TD", "DT", "DD"]   # row target/decoy + column target/decoy; index = 2*row_decoy + col_decoy
# Which side of a crosslink the rows and the columns are: the matrix has a
# direction, the CSV row A,B is in row A, column B
_AXIS_TITLES = ["row: Protein1 ↓", "col: Protein2 →"]
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

# Annotation tracks: one thin coloured strip per feature type along the
# residue axes, in a fixed colour per type
TRACK_PX = 4              # thickness of one track
TRACK_GAP_PX = 1          # space between two tracks
_FEATURE_COLOURS = dict(zip(FEATURE_TYPES, [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]))

# What the "Show" buttons open: every group and protein is set to fit
_SHOW_LEVELS = ["Top", "Complex", "Protein", "Residue"]

# Score distribution above the cutoff slider
HIST_H = 56               # height of the plot in pixels
HIST_BINS = 40
_TARGET_COLOUR = "#2a78d6"
_DECOY_COLOUR = "#eb6834"
_OUTSIDE_COLOUR = "#d4d4d4"   # bars outside the cutoff

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


def _wrap_label(text: str, font, width: float, max_lines: int = 2) -> list[str]:
    """Break a section name into at most max_lines lines no wider than
    `width` pixels, between words where possible. A word longer than a
    whole line is broken between letters, and what still doesn't fit is
    cut off with "…". Empty if not even one letter fits."""
    def fits(s: str) -> bool:
        return font.getlength(s) <= width

    lines: list[str] = []
    cur = ""
    for word in text.split():
        cand = f"{cur} {word}" if cur else word
        if fits(cand):
            cur = cand
            continue
        if cur:
            lines.append(cur)
        while not fits(word):
            k = next((k for k in range(len(word) - 1, 0, -1) if fits(word[:k])), 0)
            if k == 0:
                return []
            lines.append(word[:k])
            word = word[k:]
        cur = word
    if cur:
        lines.append(cur)
    if len(lines) > max_lines:
        last = lines[max_lines - 1]
        while last and not fits(last + "…"):
            last = last[:-1]
        lines = lines[:max_lines - 1] + ([last.rstrip() + "…"] if last else [])
    return lines


def _make_lut(cmap_name: str) -> np.ndarray:
    """A table of 256 RGB colours taken from a matplotlib colormap."""
    if not _MPL:
        return np.zeros((256, 3), dtype=np.uint8)
    cmap = matplotlib.colormaps[cmap_name]
    return (cmap(np.linspace(0, 1, 256))[:, :3] * 255).astype(np.uint8)


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


def _group_cells(rows: np.ndarray, cols: np.ndarray, scores: np.ndarray, n: int, folded: bool):
    """Sort crosslinks (row, column, score) into the cells of an n x n grid.
    Folded, a crosslink above the diagonal goes to its mirror cell below it.
    Returns the row and column of every cell that has crosslinks (sorted row
    by row), the scores ordered cell by cell, and where each cell's part of
    them starts (plus the end), ready for _aggregate_groups."""
    if folded:
        rows, cols = np.maximum(rows, cols), np.minimum(rows, cols)
    cell = rows * n + cols
    # a stable sort keeps the crosslinks of one cell in their original order
    order = np.argsort(cell, kind="stable")
    cell = cell[order]
    first = np.flatnonzero(np.diff(cell, prepend=-1))   # where the next cell begins
    cell_rows, cell_cols = np.divmod(cell[first], n)
    return cell_rows, cell_cols, scores[order], first.tolist() + [len(cell)]


# -----------------------------------------------------------------------------
# A slider with two handles
#
# Tkinter doesn't have one, so this draws its own on a canvas. The left
# handle sets the lower limit and the right handle the upper limit; the two
# can't pass each other.
# -----------------------------------------------------------------------------

class RangeSlider(tk.Canvas):
    """Pick a range [lo, hi] inside [min, max]. command(lo, hi) is called
    while you drag, release_command(lo, hi) when you let go. Above the
    track it shows how the scores are spread (set_distribution), on the
    same scale as the handles, with the part outside the range greyed out."""
    PAD = 9       # space between the canvas edge and the end of the track
    HANDLE = 5    # half the width of a handle
    TRACK_H = 26  # height of the part with the handles

    def __init__(self, master, command, release_command, **kw) -> None:
        kw.setdefault("height", HIST_H + self.TRACK_H)
        kw.setdefault("highlightthickness", 0)
        super().__init__(master, **kw)
        self._command = command
        self._release_command = release_command
        self._min, self._max = 0.0, 1.0
        self._lo, self._hi = 0.0, 1.0
        self._drag: str | None = None   # which handle is being dragged: "lo", "hi" or None
        # number of target and decoy scores per bin, see set_distribution
        self._counts: tuple[np.ndarray, np.ndarray] | None = None
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

    def set_distribution(self, targets: np.ndarray, decoys: np.ndarray) -> None:
        """Show how the target and decoy scores are spread over [min, max]:
        targets as bars, decoys as a line, both counted in HIST_BINS bins."""
        edges = np.linspace(self._min, self._max, HIST_BINS + 1) if self._max > self._min else None
        if edges is None:
            self._counts = None
        else:
            self._counts = (np.histogram(targets, edges)[0], np.histogram(decoys, edges)[0])
        self._redraw()

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

    def _draw_distribution(self) -> None:
        """Bars for the targets, a line for the decoys, on one common
        logarithmic count scale (the many low-scoring decoys would flatten
        everything else). Bins outside [lo, hi] are grey. A small legend
        sits in the top left."""
        if self._counts is None:
            return
        targets, decoys = self._counts
        top, bottom = 14, HIST_H - 2   # room for the legend above the bars
        highest = max(int(targets.max(initial=0)), int(decoys.max(initial=0)), 1)
        edges = np.linspace(self._min, self._max, HIST_BINS + 1)

        def y(count: int) -> float:
            return bottom - math.log1p(count) / math.log1p(highest) * (bottom - top)

        decoy_line: list[float] = []
        for k in range(HIST_BINS):
            x0, x1 = self._x(edges[k]), self._x(edges[k + 1])
            inside = edges[k + 1] >= self._lo and edges[k] <= self._hi
            if targets[k]:
                self.create_rectangle(x0 + 0.5, y(targets[k]), x1 - 0.5, bottom,
                                      fill=_TARGET_COLOUR if inside else _OUTSIDE_COLOUR, width=0)
            decoy_line += [x0, y(decoys[k]), x1, y(decoys[k])]
        if decoys.any():
            self.create_line(*decoy_line, fill=_DECOY_COLOUR, width=2)
        self.create_line(self.PAD, bottom, self.winfo_width() - self.PAD, bottom, fill="#999999")
        # the cutoff as two thin lines over the plot
        for value in (self._lo, self._hi):
            self.create_line(self._x(value), top, self._x(value), bottom, fill="#1a6fb5", dash=(2, 2))
        # legend: colour mark + name in normal text colour
        x = self.PAD
        for colour, name, line in ((_TARGET_COLOUR, "target", False), (_DECOY_COLOUR, "decoy", True)):
            if line:
                self.create_line(x, 6, x + 10, 6, fill=colour, width=2)
            else:
                self.create_rectangle(x, 2, x + 10, 10, fill=colour, width=0)
            item = self.create_text(x + 14, 6, text=name, anchor="w", font=("TkDefaultFont", 8))
            x = self.bbox(item)[2] + 10
        self.create_text(self.winfo_width() - self.PAD, 6, text=f"max {highest}, log", anchor="e",
                         font=("TkDefaultFont", 8), fill="gray")

    def _redraw(self) -> None:
        self.delete("all")
        self._draw_distribution()
        y = HIST_H + self.TRACK_H // 2
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
# The axis tree
#
# What stands on the axes is a tree: groups of proteins (complexes, or
# places with "Group rest by location") can sit inside bigger groups, a
# group holds proteins, and a protein holds its residues. An open node shows
# its children on the axes, a closed one stands there as a single row and
# column. Rows and columns always show the same nodes.
# -----------------------------------------------------------------------------

class AxisNode:
    """A group, a protein or a residue. key is the group name, the protein
    name or the ResidueId."""
    __slots__ = ("kind", "key", "name", "decoy", "parent", "children", "open", "depth")

    def __init__(self, kind: str, key, name: str, decoy: bool, parent: AxisNode | None) -> None:
        self.kind = kind          # "group", "protein" or "residue"
        self.key = key
        self.name = name
        self.decoy = decoy
        self.parent = parent
        self.children: list[AxisNode] = []
        self.open = kind == "group"   # proteins start closed
        self.depth = 0 if parent is None else parent.depth + 1
        if parent is not None:
            parent.children.append(self)

    def path(self) -> list[AxisNode]:
        """The nodes above this one, the biggest first."""
        path, node = [], self.parent
        while node is not None:
            path.insert(0, node)
            node = node.parent
        return path

    def walk(self):
        """This node and everything under it."""
        yield self
        for child in self.children:
            yield from child.walk()

    def protein(self) -> str | None:
        """The protein this node is or belongs to (None for a group)."""
        if self.kind == "protein":
            return self.key
        return self.key.protein if self.kind == "residue" else None

    def group(self) -> str:
        """The name of the nearest group above this node, or ""."""
        node = self.parent
        while node is not None and node.kind != "group":
            node = node.parent
        return node.name if node is not None else ""


def _build_tree(proteins: list[str], decoy: dict[str, bool], sections: dict[str, str] | None,
                parents: dict[str, str], residues: dict[str, list]) -> list[AxisNode]:
    """The top nodes of the tree. proteins are in axis order (sections stay
    together), sections names the group of each protein, parents the bigger
    group of a group, residues the ResidueIds of each protein sorted by
    position. Targets and decoys never share a group."""
    top: list[AxisNode] = []
    groups: dict[tuple[bool, str], AxisNode] = {}

    def group(name: str, is_decoy: bool) -> AxisNode:
        if (is_decoy, name) not in groups:
            bigger = parents.get(name)
            parent = group(bigger, is_decoy) if bigger else None
            node = groups[(is_decoy, name)] = AxisNode("group", name, name, is_decoy, parent)
            if parent is None:
                top.append(node)
        return groups[(is_decoy, name)]

    for p in proteins:
        is_decoy = decoy.get(p, False)
        parent = group(sections[p], is_decoy) if sections and p in sections else None
        node = AxisNode("protein", p, p, is_decoy, parent)
        if parent is None:
            top.append(node)
        for r in residues.get(p, ()):
            AxisNode("residue", r, str(r.pos), is_decoy, node)
    return top


def _slot_labels(node: AxisNode) -> tuple[str, str]:
    """The column and the row name of a node. A closed group gets a "+" in
    front, so it can be told apart from a protein."""
    if node.kind == "residue":
        return str(node.key.pos)[:COL_W], f"{short_accession(node.key.protein)}:{node.key.pos}"[:ROW_W]
    name = "+" + node.name if node.kind == "group" else node.name
    return name[:COL_W], name[:ROW_W]


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

        # The loaded data: the best score of every protein pair and of every
        # residue pair (empty without residue positions in the CSV), each
        # once more as arrays (row, column, score) of protein or residue
        # numbers, because numpy is much faster on these
        self._prot_sparse: dict[tuple, float] = {}
        self._res_sparse: dict[tuple, float] = {}
        self._prot_index: dict[str, int] = {}
        self._res_index: dict = {}
        self._res_protein: np.ndarray = np.zeros(0, dtype=np.int64)   # protein number of each residue
        self._prot_links = (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0))
        self._res_links = (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0))
        self._seqs: dict[str, str] = {}
        self._seq_index: dict[str, str] = {}   # self._seqs, also keyed by bare accession

        # The axis tree and what of it is on the axes right now: the slots
        # are the nodes that stand on the axes, in order (see _rebuild_axis)
        self._tree: list[AxisNode] = []
        self._slots: list[AxisNode] = []
        self._slot_paths: list[list[AxisNode]] = []   # the open nodes above each slot
        self._band_levels: int = 0        # rows of group/protein names above the slot names
        self._uniform: bool = True        # all slots are proteins, or all are residues
        self._has_residue_slots: bool = False
        self._split: int = 0              # index of the first decoy slot
        self._col_labels: list[str] = []
        self._row_labels: list[str] = []
        # Every crosslink that the slots show, as (row slot, column slot,
        # score): protein pairs between closed proteins, residue pairs where
        # at least one protein is open
        self._link_rows: np.ndarray = np.zeros(0, dtype=np.int64)
        self._link_cols: np.ndarray = np.zeros(0, dtype=np.int64)
        self._link_scores: np.ndarray = np.zeros(0)

        # Annotations from UniProt: {protein: [[type, start, end, text], ...]},
        # None until they are first switched on
        self._features: dict[str, list] | None = None
        self._features_loading: bool = False
        self._shown_features: list[str] = []
        self._slot_tracks: np.ndarray = np.zeros((0, 0), dtype=bool)    # (feature type, slot)
        self._block_tracks: np.ndarray = np.zeros((0, 0), dtype=bool)   # (feature type, block)

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
        self._legend_photo: ImageTk.PhotoImage | None = None
        # Writing text is slow, so every name and every score is only written once
        self._label_cache: dict[str, Image.Image] = {}
        self._slant_cache: dict[str, tuple[Image.Image, int, int]] = {}
        self._vname_cache: dict[tuple[str, int], Image.Image] = {}
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
        self._order_var = tk.StringVar(value="confidence")
        field(data, 3, "Order", ttk.Combobox(data, textvariable=self._order_var, width=10,
                                             state="readonly", values=ORDER_MODES))
        # only used by the "complex" order: group the rest by cellular location
        self._locations_var = tk.BooleanVar(value=False)
        locations = tk.Checkbutton(data, text="Group rest by location",
                                   variable=self._locations_var)
        field(data, 5, "", locations)
        locations.grid_remove()   # grid() brings it back with the same options
        self._order_var.trace_add("write", lambda *_: (
            locations.grid() if self._order_var.get() == "complex" else locations.grid_remove()))

        buttons = tk.Frame(data)
        buttons.grid(row=6, column=0, columnspan=3, sticky="ew", pady=(4, 0))
        self._load_btn = tk.Button(buttons, text="Load", command=self._load, width=8)
        self._load_btn.pack(side=tk.LEFT)
        tk.Button(buttons, text="Export PNG…", command=self._export).pack(side=tk.LEFT, padx=4)
        self._status_var = tk.StringVar(value="No data loaded.")
        tk.Label(data, textvariable=self._status_var, fg="gray", wraplength=wrap,
                 justify=tk.LEFT).grid(row=7, column=0, columnspan=3, sticky="w")

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
        # Show: open or close every group and protein at once (single ones
        # are opened by clicking their name, closed by clicking their band)
        show = tk.Frame(display)
        for k, level in enumerate(_SHOW_LEVELS):
            tk.Button(show, text=level, padx=1,
                      command=lambda lv=level: self._show_level(lv)).grid(row=0, column=k, sticky="ew")
            show.columnconfigure(k, weight=1)
        tk.Label(display, text="Show").grid(row=8, column=0, sticky="w", pady=(4, 0))
        show.grid(row=9, column=0, columnspan=3, sticky="ew", pady=1)
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
        # what the colours mean, redrawn whenever that changes (_update_legend)
        self._legend = tk.Label(display, bg=BG, bd=0)
        self._legend.grid(row=7, column=0, columnspan=3, sticky="w", pady=(4, 0))

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

        # Annotations: UniProt features as coloured tracks along residue axes
        annotations = group("Annotations")
        self._feature_vars: dict[str, tk.BooleanVar] = {}
        for k, name in enumerate(FEATURE_TYPES):
            var = tk.BooleanVar(value=False)
            var.trace_add("write", lambda *_: self._on_features_change())
            self._feature_vars[name] = var
            cell = tk.Frame(annotations)
            cell.grid(row=k // 2, column=k % 2, sticky="w", padx=(0, 6))
            tk.Label(cell, bg=_FEATURE_COLOURS[name], width=1, height=1, bd=0).pack(side=tk.LEFT)
            tk.Checkbutton(cell, text=name, variable=var, padx=0).pack(side=tk.LEFT)
        self._features_var = tk.StringVar(value="shown along residue rows and columns")
        tk.Label(annotations, textvariable=self._features_var, fg="gray", wraplength=wrap,
                 justify=tk.LEFT).grid(row=4, column=0, columnspan=2, sticky="w")

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
        self._corner = tk.Canvas(outer, bg=BG, highlightthickness=0)
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

        # clicking a name opens it, clicking a band above/left of names closes
        # it; the wheel zooms here too, around the pointer
        self._ch_canvas.bind("<ButtonRelease-1>", lambda e: self._on_label_click(e, "col"))
        self._rl_canvas.bind("<ButtonRelease-1>", lambda e: self._on_label_click(e, "row"))
        for canvas in (self._ch_canvas, self._rl_canvas):
            for wheel in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                canvas.bind(wheel, self._wheel)
        self._mx.bind("<Configure>", lambda _e: (self._update_scrollregion(), self._schedule()))
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
    def _band_h(self) -> int:
        """The height of one row of group/protein names on top (and the
        width of one on the left): room for two lines of text."""
        return 2 * self._ch

    @property
    def _section_band(self) -> int:
        """All rows of group/protein names together, one per level of open
        nodes above the slots."""
        return self._band_levels * self._band_h

    @property
    def _name_w(self) -> int:
        """The width of the row names alone."""
        return ROW_W * self._cw + 6

    @property
    def _tracks_px(self) -> int:
        """The room the annotation tracks take, next to the names. Only
        there while residues are on the axes."""
        if not self._shown_features or not self._has_residue_slots:
            return 0
        return len(self._shown_features) * (TRACK_PX + TRACK_GAP_PX) + TRACK_GAP_PX

    def _names_shown(self, overview: bool) -> bool:
        """Whether the single slot names are shown. Zoomed out into blocks
        they can't be read; residues keep a sparse position axis instead."""
        return not overview or self._has_residue_slots

    def _headers_shown(self, overview: bool) -> bool:
        """Whether there are column headers and row labels at all. Zoomed out
        into blocks, group names (complexes, …) are still worth showing:
        that's where links between groups are looked at."""
        return self._names_shown(overview) or self._band_levels > 0

    def _label_sizes(self, overview: bool) -> tuple[int, int, int]:
        """(width of everything left of the matrix, height of the line under
        the column names counted from the top, height of the whole header)
        while zoomed out into blocks or not. The group/protein names come
        first, then the slot names if they are shown, then the annotation
        tracks next to the matrix."""
        names = self._names_shown(overview)
        tracks = self._tracks_px
        width = self._section_band + (self._name_w if names else 0) + tracks
        base = self._section_band + (self._slant_h if names else 0)
        return width, base, base + self._ch // 2 + 2 + tracks

    @property
    def _label_w(self) -> int:
        return self._label_sizes(self._overview)[0]

    @property
    def _label_base(self) -> int:
        return self._label_sizes(self._overview)[1]

    @property
    def _header_h(self) -> int:
        return self._label_sizes(self._overview)[2]

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
        self._draw_corner()

    def _draw_corner(self) -> None:
        """The axis titles in the empty corner above the row names, if they fit."""
        c = self._corner
        c.delete("all")
        w, base = self._label_w, self._label_base
        lines = _AXIS_TITLES
        if self._lf is None or max(self._lf.measure(t) for t in lines) > w - 4 or base < 2 * self._ch:
            return
        for k, text in enumerate(reversed(lines)):
            c.create_text(w - 4, base - 2 - k * self._ch, text=text, anchor="se",
                          font=self._lf, fill="gray")

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

    def _margins(self) -> tuple[int, int]:
        """How far the view can be moved past the matrix edges, left/right
        and up/down. Almost a whole view, like on a map: that way zooming
        keeps the spot under the mouse in place even near an edge."""
        return (max(self._pad, self._mx.winfo_width() - MIN_VISIBLE_PX),
                max(self._pad, self._mx.winfo_height() - MIN_VISIBLE_PX))

    def _update_scrollregion(self) -> None:
        """Fit the scroll region to the matrix size and the view size. The
        view keeps showing the same part of the matrix."""
        ox, oy = self._mx.canvasx(0), self._mx.canvasy(0)
        size = len(self._blocks) * self._cell_px
        mx, my = self._margins()
        self._mx.configure(scrollregion=(-mx, -my, size + mx, size + my))
        self._ch_canvas.configure(scrollregion=(-mx, 0, size + mx, self._header_h))
        self._rl_canvas.configure(scrollregion=(0, -my, self._label_w, size + my))
        self._scroll_to(ox, oy)

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
        the top-left corner. Tk keeps the view inside the scroll region
        (the matrix plus _margins) by itself."""
        size = len(self._blocks) * self._cell_px
        mx, my = self._margins()
        if size:
            fx, fy = (ox + mx) / (size + 2 * mx), (oy + my) / (size + 2 * my)
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
        order = self._order_var.get()
        fasta = self._fasta_var.get().strip() or None
        locations = self._locations_var.get()
        if not csv:
            messagebox.showerror("Error", "Please select a CSV file.")
            return
        if order == "sequence" and not fasta:
            messagebox.showerror("Error", "Order 'sequence' requires a FASTA file.")
            return

        self._load_btn.configure(state=tk.DISABLED)
        self._status_var.set("Loading…")

        def worker() -> None:
            try:
                # Mismatched protein names between the CSV and the FASTA (or
                # other data quirks) only raise a warning, not an error, and
                # that warning would otherwise go to a console the GUI user
                # doesn't see — so it's collected here and shown in Status.
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    seqs = read_fasta(fasta) if fasta else None
                    proteins, prot_sparse, decoy = build_matrix(csv, n=None)   # always all proteins
                    proteins, sections = sort_proteins(proteins, decoy, order, seqs, locations)
                    parents = complex_hierarchy(set(sections.values())) if order == "complex" and sections else {}
                    # the residues, so that proteins can be opened
                    try:
                        residues, res_sparse, res_decoy = build_residue_matrix(csv, n=None)
                        if seqs:
                            # with a FASTA, every residue of the loaded proteins is shown,
                            # not only the crosslinked ones
                            residues, res_decoy = add_unlinked_residues(residues, res_decoy, seqs)
                    except ValueError:
                        residues, res_sparse = [], {}
                        warnings.warn("no residue positions in the CSV, proteins can't be opened")
                    notes = [str(w.message) for w in caught]
                data = dict(proteins=proteins, prot_sparse=prot_sparse, decoy=decoy, sections=sections,
                            parents=parents, residues=residues, res_sparse=res_sparse, seqs=seqs)
                self.after(0, lambda: self._on_loaded(data, notes))
            except Exception as exc:
                msg = str(exc)
                self.after(0, lambda: self._on_error(msg))

        threading.Thread(target=worker, daemon=True).start()

    def _on_loaded(self, data: dict, notes: list[str] = ()) -> None:
        """Take over what the load worker read (see _load): build the axis
        tree and the link arrays, then show the proteins."""
        proteins, decoy = data["proteins"], data["decoy"]
        self._prot_sparse, self._res_sparse = data["prot_sparse"], data["res_sparse"]
        self._seqs = data["seqs"] or {}
        self._seq_index = seq_index(self._seqs)   # also reachable by bare accession

        by_protein: dict[str, list] = {}
        for r in data["residues"]:
            by_protein.setdefault(r.protein, []).append(r)
        for rs in by_protein.values():
            rs.sort(key=lambda r: r.pos)
        self._tree = _build_tree(proteins, decoy, data["sections"], data["parents"], by_protein)

        # numbers for proteins and residues, and the links between them as arrays
        self._prot_index = {p: i for i, p in enumerate(proteins)}
        residues = [n.key for top in self._tree for n in top.walk() if n.kind == "residue"]
        self._res_index = {r: i for i, r in enumerate(residues)}
        self._res_protein = np.array([self._prot_index[r.protein] for r in residues], dtype=np.int64)

        def arrays(sparse: dict, index: dict):
            return (np.array([index[a] for a, _ in sparse], dtype=np.int64),
                    np.array([index[b] for _, b in sparse], dtype=np.int64),
                    np.array(list(sparse.values()), dtype=float))
        self._prot_links = arrays(self._prot_sparse, self._prot_index)
        self._res_links = arrays(self._res_sparse, self._res_index)

        scores = np.concatenate([self._prot_links[2], self._res_links[2]])
        self._min_score = float(scores.min()) if len(scores) else 0.0
        self._max_score = float(scores.max()) if len(scores) else 1.0
        # colours are measured against the scores of the level on screen (see _rebuild_axis)
        self._score_sorted = np.sort(self._prot_links[2])
        self._label_cache.clear()
        self._slant_cache.clear()
        self._vname_cache.clear()
        # measured once over every name that can come up, so the header
        # doesn't change size when something is opened
        self._measure_col_labels([_slot_labels(n)[0] for top in self._tree for n in top.walk()])
        self._cut_lo, self._cut_hi = self._min_score, self._max_score
        self._cutoff_slider.set_bounds(self._cut_lo, self._cut_hi)
        if self._features is not None:
            self._features = None   # the new proteins need their own
            self._load_features()

        self._show_level("Protein")   # also moves the view to the top-left corner

        n_targets = sum(1 for p in proteins if not decoy.get(p, False))
        status = (f"Loaded — {n_targets} target · {len(proteins) - n_targets} decoy proteins · "
                  f"{len(self._prot_sparse)} protein links · {len(self._res_sparse)} residue links")
        if notes:
            # e.g. CSV protein names the FASTA didn't recognise — shortened so
            # one long list of names can't push the rest of the panel off screen
            shortened = [n if len(n) <= 220 else n[:220] + "…" for n in notes]
            status += "\n⚠ " + "\n⚠ ".join(shortened)
        self._status_var.set(status)
        self._load_btn.configure(state=tk.NORMAL)

    # -------------------------------------------------------------------------
    # Opening and closing parts of the axis tree
    #
    # Clicking a name on an axis opens that group or protein, clicking one of
    # the bands above (or left of) the names closes it again. The "Show"
    # buttons set the whole tree at once. Afterwards the slots, the links
    # they show and everything that depends on them are worked out again.
    # -------------------------------------------------------------------------

    def _show_level(self, level: str) -> None:
        """Open or close every node: "Top" closes the biggest groups,
        "Complex" the smallest ones (those without groups inside), "Protein"
        shows every protein and "Residue" every residue."""
        for top in self._tree:
            for node in top.walk():
                if node.kind == "protein":
                    node.open = level == "Residue"
                elif node.kind == "group":
                    if level == "Top":
                        node.open = node.depth > 0
                    elif level == "Complex":
                        node.open = any(c.kind == "group" for c in node.children)
                    else:
                        node.open = True
        self._apply_axis()
        self._scroll_to(-self._pad, -self._pad)

    def _on_label_click(self, event: tk.Event, axis: str) -> None:
        """A click on the column headers or the row labels: on a band, close
        that group or protein; on a name, open it."""
        if not self._blocks or self._dragging:
            return
        if axis == "col":
            along, across = self._ch_canvas.canvasx(event.x), event.y
        else:
            along, across = self._rl_canvas.canvasy(event.y), event.x
        block = int(along // self._cell_px)
        if not 0 <= block < len(self._blocks):
            return
        start, end = self._blocks[block]
        if across < self._section_band:
            path = self._slot_paths[start]
            level = int(across // self._band_h)
            if level < len(path):
                self._set_open(path[level], False)
        elif end - start == 1 and self._names_shown(self._overview):
            node = self._slots[start]
            if node.children:
                self._set_open(node, True)
            elif node.kind == "protein":
                self._status_var.set(f"{node.name}: no residue positions to open")

    def _set_open(self, node: AxisNode, is_open: bool) -> None:
        """Open or close one node. Its first row/column stays where it is on
        screen."""
        def first_slot() -> int:
            return next(k for k, s in enumerate(self._slots) if s is node or node in self._slot_paths[k])
        k = first_slot()
        screen_x = self._pixel_of(k) - self._mx.canvasx(0)
        screen_y = self._pixel_of(k) - self._mx.canvasy(0)
        node.open = is_open
        self._apply_axis()
        k = first_slot()
        self._scroll_to(self._pixel_of(k) - screen_x, self._pixel_of(k) - screen_y)

    def _rebuild_axis(self) -> None:
        """Work out the slots (the nodes standing on the axes), their names
        and bands, and which crosslinks they show."""
        slots: list[AxisNode] = []

        def walk(nodes: list[AxisNode]) -> None:
            for node in nodes:
                if node.open and node.children:
                    walk(node.children)
                else:
                    slots.append(node)
        walk(self._tree)
        self._slots = slots
        self._slot_paths = [node.path() for node in slots]
        self._band_levels = max((node.depth for node in slots), default=0)
        kinds = {node.kind for node in slots}
        self._uniform = kinds in ({"protein"}, {"residue"})
        self._has_residue_slots = "residue" in kinds
        self._split = next((k for k, node in enumerate(slots) if node.decoy), len(slots))
        labels = [_slot_labels(node) for node in slots]
        self._col_labels = [c for c, _ in labels]
        self._row_labels = [r for _, r in labels]

        # Which slot every protein and every residue is in. An open protein
        # has no slot of its own (-1): its residues have one each.
        prot_slot = np.full(len(self._prot_index), -1, dtype=np.int64)
        res_slot = np.full(len(self._res_index), -1, dtype=np.int64)
        for k, node in enumerate(slots):
            if node.kind == "residue":
                res_slot[self._res_index[node.key]] = k
            else:
                for n in node.walk():
                    if n.kind == "protein":
                        prot_slot[self._prot_index[n.key]] = k
        prot_open = prot_slot < 0
        res_slot = np.where(res_slot >= 0, res_slot, prot_slot[self._res_protein])

        # protein links between closed proteins, residue links where at
        # least one of the two proteins is open
        pr, pc, ps = self._prot_links
        rr, rc, rs = self._res_links
        use_p = ~(prot_open[pr] | prot_open[pc])
        use_r = prot_open[self._res_protein[rr]] | prot_open[self._res_protein[rc]]
        self._link_rows = np.concatenate([prot_slot[pr[use_p]], res_slot[rr[use_r]]])
        self._link_cols = np.concatenate([prot_slot[pc[use_p]], res_slot[rc[use_r]]])
        self._link_scores = np.concatenate([ps[use_p], rs[use_r]])
        # the colours follow the scores of the level on screen
        self._score_sorted = np.sort(rs if kinds == {"residue"} else ps if kinds == {"protein"}
                                     else self._link_scores)
        self._slot_tracks = self._tracks_for(slots)

    def _apply_axis(self) -> None:
        """After the slots changed: blocks, zoom range, layout and the score
        distribution follow. The selection is cleared."""
        self._rebuild_axis()
        self._clear_selection()
        self._runs = self._find_runs()
        self._recompute_blocks()
        self._zoom_min = self._find_zoom_min()
        self._zoom.configure(from_=self._zoom_min)
        if self._zoom_value < self._zoom_min:   # the new axes need fewer block steps
            self._on_zoom(str(self._zoom_min))
            self._zoom.set(self._zoom_min)
        for target, btn in self._zoom_to_btns.items():
            btn.configure(state=tk.NORMAL if self._zoom_region(target) else tk.DISABLED)
        self._update_distribution()
        self._update_cutoff_label()
        n_res = sum(1 for node in self._slots if node.kind == "residue")
        self._info_var.set(f"{len(self._slots)} rows ({n_res} residues)")
        # the bands and names can take a different amount of room now
        self._apply_header_layout()
        self._update_scrollregion()
        self._schedule()

    # -------------------------------------------------------------------------
    # Annotation tracks
    #
    # UniProt features (helices, binding sites, ...) of the chosen types are
    # drawn as thin coloured strips along every residue row and column. They
    # are fetched in the background the first time a type is switched on.
    # -------------------------------------------------------------------------

    def _on_features_change(self) -> None:
        self._shown_features = [name for name, var in self._feature_vars.items() if var.get()]
        if self._shown_features and self._features is None:
            self._load_features()
            return
        self._update_tracks()

    def _load_features(self) -> None:
        """Fetch the features of every loaded protein in the background."""
        if self._features_loading or not self._prot_index:
            return
        self._features_loading = True
        self._features_var.set("fetching from UniProt…")
        proteins = list(self._prot_index)

        def worker() -> None:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                features = features_for(proteins)
            note = str(caught[0].message) if caught else None
            self.after(0, lambda: self._on_features_loaded(features, note))
        threading.Thread(target=worker, daemon=True).start()

    def _on_features_loaded(self, features: dict[str, list], note: str | None) -> None:
        self._features_loading = False
        self._features = features
        n = sum(1 for f in features.values() if f)
        self._features_var.set(note or f"{n} of {len(features)} proteins annotated; "
                                        "shown along residue rows and columns")
        self._update_tracks()

    def _update_tracks(self) -> None:
        if not self._slots:
            return
        self._slot_tracks = self._tracks_for(self._slots)
        self._block_tracks = self._tracks_per_block()
        self._apply_header_layout()   # the tracks take room next to the names
        self._update_scrollregion()
        self._schedule()

    def _tracks_for(self, slots: list[AxisNode]) -> np.ndarray:
        """For each shown feature type and each slot: whether that residue
        lies inside such a feature. Proteins and groups never do."""
        tracks = np.zeros((len(self._shown_features), len(slots)), dtype=bool)
        if self._features is None:
            return tracks
        covered: dict[tuple[str, str], np.ndarray] = {}   # (protein, type) -> one bool per position

        def coverage(protein: str, kind: str) -> np.ndarray:
            if (protein, kind) not in covered:
                spans = [(f[1], f[2]) for f in self._features.get(protein, ()) if f[0] == kind]
                mask = np.zeros(max((e for _, e in spans), default=0) + 1, dtype=bool)
                for start, end in spans:
                    mask[start:end + 1] = True
                covered[(protein, kind)] = mask
            return covered[(protein, kind)]

        for t, kind in enumerate(self._shown_features):
            for k, node in enumerate(slots):
                if node.kind == "residue":
                    mask = coverage(node.key.protein, kind)
                    tracks[t, k] = node.key.pos < len(mask) and mask[node.key.pos]
        return tracks

    def _tracks_per_block(self) -> np.ndarray:
        """_slot_tracks per block: a block has a feature if any of its slots has it."""
        if not self._blocks or not self._slot_tracks.size:
            return np.zeros((len(self._shown_features), len(self._blocks)), dtype=bool)
        starts = [s for s, _ in self._blocks]
        return np.logical_or.reduceat(self._slot_tracks, starts, axis=1)

    def _features_at(self, node: AxisNode) -> list[str]:
        """The features a residue lies in, as short texts."""
        if self._features is None or node.kind != "residue":
            return []
        pos = node.key.pos
        out = []
        for kind, start, end, text in self._features.get(node.key.protein, ()):
            if start <= pos <= end:
                where = f"{start}" if start == end else f"{start}-{end}"
                out.append(f"{kind} {where}" + (f" ({text})" if text else ""))
        return out

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
        self._update_legend()
        self._schedule()

    def _colours_for(self, scores: np.ndarray, ref: np.ndarray | None = None) -> np.ndarray:
        """The colours of the values, one RGB row each. ref are the sorted
        values the colours are measured against; without it, the ones the
        cells on screen use (all scores, or for counts and root sums of
        squares the values on screen, see _reduce_block_scores). Without a
        colormap: grey from light (smallest) to dark (largest), linear in
        between. With one: the colour follows the value's rank in ref
        (quantile), so the median gets the middle colour."""
        ref = self._colour_ref if ref is None else ref
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

    def _colour_caption(self, folded: bool | None = None, merged: bool | None = None) -> str:
        """What the cell colours stand for, for the colour legend. Screen
        settings unless folded/merged are given (the export has no blocks)."""
        folded = self._folded if folded is None else folded
        merged = (self._block_size > 1 or not self._uniform) if merged is None else merged
        if not (folded or merged):
            what = "score"
        elif self._agg_method == "count":
            what = "crosslinks per cell"
        else:
            what = f"{_AGG_NAMES[self._agg_method]} of scores"
        return f"{what} ({'by rank' if self._cmap_lut is not None else 'linear'})"

    def _legend_image(self, width: int, ref: np.ndarray, caption: str) -> Image.Image:
        """A colour bar for the values in ref (sorted), as wide as `width`:
        the caption, the bar, and the smallest, middle and largest value
        under it. The bar uses the same colouring as the cells, so with a
        colormap the middle is the median and without one the midpoint."""
        lh, font, pad = self._ch, self._pil_lf, 4
        bar_h = 10
        img = Image.new("RGB", (max(1, width), 2 * lh + bar_h + 4), _BG)
        draw = ImageDraw.Draw(img)
        draw.text((pad, 0), caption, font=font, fill=_FG)
        bar_w = width - 2 * pad
        if len(ref) == 0 or bar_w < 2:
            return img
        frac = np.linspace(0, 1, bar_w)
        if self._cmap_lut is not None:
            values = ref[np.round(frac * (len(ref) - 1)).astype(np.int64)]
        else:
            values = ref[0] + frac * (ref[-1] - ref[0])
        bar = np.broadcast_to(self._colours_for(values, ref), (bar_h, bar_w, 3))
        img.paste(Image.fromarray(np.ascontiguousarray(bar)), (pad, lh + 1))
        draw.rectangle([pad, lh + 1, pad + bar_w - 1, lh + bar_h], outline=_SEP)
        y = lh + bar_h + 3
        for k, anchor in ((0, "la"), (len(values) // 2, "ma"), (len(values) - 1, "ra")):
            x = pad + (k if anchor != "ma" else bar_w // 2)
            draw.text((x, y), f"{values[k]:.3g}", font=font, anchor=anchor, fill=_FG)
        return img

    def _update_legend(self) -> None:
        """Redraw the colour legend in the control panel."""
        if not self._slots or self._pil_lf is None:
            return
        img = self._legend_image(SIDEBAR_W - 40, self._colour_ref, self._colour_caption())
        self._legend_photo = ImageTk.PhotoImage(img)
        self._legend.configure(image=self._legend_photo, height=img.height)

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
        self._update_distribution()
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
        if not len(self._score_sorted) and not self._prot_index:
            return
        self._cutoff_slider.set_bounds(self._min_score, self._max_score)
        self._on_cutoff(self._min_score, self._max_score)
        self._on_cutoff_release(self._min_score, self._max_score)

    def _update_distribution(self) -> None:
        """The score distribution over the cutoff slider: the crosslinks on
        the axes in the quadrants that are switched on, TT as targets and
        the rest (with at least one decoy side) as decoys."""
        rows, cols, scores = self._link_rows, self._link_cols, self._link_scores
        quadrant = 2 * (rows >= self._split) + (cols >= self._split)
        shown = np.array([self._quad_shown[q] for q in _QUADRANTS])[quadrant]
        self._cutoff_slider.set_distribution(scores[shown & (quadrant == 0)],
                                             scores[shown & (quadrant > 0)])

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
        stays where it is while everything around it grows or shrinks. The
        wheel may move the view past the matrix edges for that (see
        _margins); the slider stays inside them where it can."""
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
        ox, oy = self._pixel_of(item_x) - ax, self._pixel_of(item_y) - ay
        if anchor is None:
            # the slider: like before the view could go past the edges,
            # don't show empty space next to the matrix where it can be helped
            size, pad = len(self._blocks) * self._cell_px, self._pad
            w, h = self._mx.winfo_width(), self._mx.winfo_height()
            ox = max(-pad, min(ox, size + pad - w))
            oy = max(-pad, min(oy, size + pad - h))
        self._scroll_to(ox, oy)

    @property
    def _headers_hidden(self) -> bool:
        """Whether the column headers and row labels are removed from the
        layout right now (see _headers_shown)."""
        return not self._headers_shown(self._overview)

    def _set_overview(self, on: bool) -> None:
        """Zoomed out into blocks, the item names go: the headers and row
        labels shrink to the section names, or disappear when there are no
        sections. Residues keep a sparse position axis instead. Zooming
        back in brings everything back."""
        if on == self._overview:
            return
        self._overview = on
        self._apply_header_layout()

    def _apply_header_layout(self) -> None:
        """Show or hide the headers and row labels and give them the right
        size, for the current zoom and data (see _set_overview)."""
        hide = self._headers_hidden
        for widget in (self._corner, self._ch_canvas, self._rl_canvas):
            if hide:
                widget.grid_remove()
            else:
                widget.grid()
        self._sync_sizes()        # the labels take a different amount of room
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
        """Stretches of slots under the same node (the residues of one
        protein, the proteins of one group), also cut at the target/decoy
        line. Blocks are built from these (see _blocks_for)."""
        slots, split = self._slots, self._split
        runs, start = [], 0
        for j in range(1, len(slots) + 1):
            if j == len(slots) or j == split or slots[j].parent is not slots[start].parent:
                runs.append((start, j))
                start = j
        return runs

    def _find_zoom_min(self) -> int:
        """The lowest zoom slider value: at least BLOCK_ZOOM_STEPS below
        CELL_PX_MIN, and further for big matrices, until the whole matrix
        fits into ZOOM_OUT_FIT_PX or blocks can't get any bigger (one block
        for the targets and one for the decoys)."""
        n = len(self._slots)
        longest = max(self._split, n - self._split, 1)
        v = CELL_PX_MIN - BLOCK_ZOOM_STEPS
        while self._extent_px(v, 0, n) > ZOOM_OUT_FIT_PX and self._block_size_for(v) < longest:
            v -= 1
        return v

    def _zoom_region(self, target: str) -> tuple[int, int, int, int] | None:
        """The items (row start, row end, column start, column end) that
        `target` covers: one quadrant, or for "All" the box around all
        quadrants that are switched on. None if it holds no items."""
        n, split = len(self._slots), self._split
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
        if not self._headers_hidden:
            full_w, full_h = full_w + self._label_w, full_h + self._header_h

        def fits(v: int) -> bool:
            # what's left for the matrix at zoom v: the labels of that zoom
            # take their room, and zoomed out into blocks also the padding
            overview = v < CELL_PX_MIN
            label_w, _, header_h = self._label_sizes(overview)
            if not self._headers_shown(overview):
                label_w = header_h = 0
            pad = 2 * OVERVIEW_PAD_PX if overview else 0
            fw, fh = full_w - label_w - pad, full_h - header_h - pad
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
        block = self._item_to_block[min(int(item), len(self._slots) - 1)]
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
        n, split, bs = len(self._slots), self._split, self._block_size
        blocks = self._blocks_for(bs)

        self._blocks = blocks
        self._item_to_block = [bi for bi, (s, e) in enumerate(blocks) for _ in range(s, e)]
        self._split_block = self._item_to_block[split] if 0 <= split < n else len(blocks)
        self._block_tracks = self._tracks_per_block()
        self._group_links()
        self._reduce_block_scores()

    def _group_links(self) -> None:
        """Sort every visible crosslink into its cell (row block of Protein1,
        column block of Protein2). Links hidden by the cutoff or by a
        quadrant checkbox are left out, so the block values only describe
        what you can see. A block never crosses the target/decoy line, so
        each cell lies in one quadrant. Folded, a crosslink above the
        diagonal goes to its mirror cell below it (TD then joins DT)."""
        keep = self._shown_links()
        to_block = np.array(self._item_to_block, dtype=np.int64)
        rb, cb = to_block[self._link_rows[keep]], to_block[self._link_cols[keep]]
        (self._cell_rows, self._cell_cols, self._group_scores,
         self._group_starts) = _group_cells(rb, cb, self._link_scores[keep],
                                            max(1, len(self._blocks)), self._folded)

    def _reduce_block_scores(self) -> None:
        """Turn the crosslinks of every cell into the one value it shows: a
        score, or with "Count" the number of crosslinks."""
        scores, starts = self._group_scores, self._group_starts
        # folded, even a single cell can hold two crosslinks (both
        # directions), and so can a closed group or a closed protein next to
        # residues
        combined = self._merges_links
        if combined:
            values = _aggregate_groups(scores, starts, self._agg_method)
        else:
            values = scores   # one slot of one kind per block, so a cell holds exactly one crosslink
        self._cell_scores = np.array(values, dtype=float)
        if combined and self._agg_method in ("count", "rss"):
            # counts aren't scores, and a root sum of squares grows with the
            # number of crosslinks past every single score, so both are
            # coloured by their rank among the values on screen
            self._colour_ref = np.sort(self._cell_scores)
        else:
            self._colour_ref = self._score_sorted
        self._cell_rgb = self._colours_for(self._cell_scores)
        self._update_legend()

    @property
    def _merges_links(self) -> bool:
        """Whether a cell on screen can hold more than one crosslink."""
        return self._block_size > 1 or self._folded or not self._uniform

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
        row_nodes = self._slots[slice(*self._blocks[rb])]
        col_nodes = self._slots[slice(*self._blocks[cb])]
        score_line = self._score_line(rb, cb)
        # the annotations of a single residue on either side
        notes = [f"{side}: {'; '.join(found)}"
                 for side, nodes in (("row", row_nodes), ("col", col_nodes))
                 if len(nodes) == 1 and (found := self._features_at(nodes[0]))]

        def show(fetched: tuple | None = None) -> None:
            sides = [self._describe_side(nodes, fetched[k] if fetched else None, fetched is not None)
                     for k, nodes in enumerate((row_nodes, col_nodes))]
            self._sel_info_var.set("\n".join(sides + [score_line] + notes))

        show()
        # UniProt only makes sense when each side is (part of) a single protein
        row_proteins = {n.protein() for n in row_nodes}
        col_proteins = {n.protein() for n in col_nodes}
        if len(row_proteins) == 1 and len(col_proteins) == 1 and None not in row_proteins | col_proteins:
            def fetch() -> None:
                fetched = tuple(self._fetch_uniprot(short_accession(next(iter(proteins))))
                                for proteins in (row_proteins, col_proteins))
                self.after(0, lambda: self._sel_token == token and show(fetched))
            threading.Thread(target=fetch, daemon=True).start()

    def _describe_side(self, nodes: list[AxisNode], uniprot: tuple[str, int] | None,
                       fetched: bool) -> str:
        """One line of text about the row or column side of the selection."""
        decoys = {n.decoy for n in nodes}
        tag = "target" if decoys == {False} else "decoy" if decoys == {True} else "mixed"
        groups = {n.group() for n in nodes}
        group = next(iter(groups)) if len(groups) == 1 else ""
        sec = f"  [{group}]" if group else ""

        if uniprot:
            name, length = uniprot
            acc = short_accession(nodes[0].protein())
            return f"{name}  ({length} aa, {tag}){sec}{self._pos_range(nodes)}  [{acc}]"
        if len(nodes) > 1:
            return f"{len(nodes)} items  (?, {tag}){sec}" if fetched else f"{len(nodes)} items ({tag}){sec}"
        node = nodes[0]
        if node.kind == "group":
            n_proteins = sum(1 for n in node.walk() if n.kind == "protein")
            return f"{node.name}  ({n_proteins} proteins, {tag}){sec}"
        seq = find_sequence(node.protein(), self._seq_index)
        if fetched:
            size = "?"
        elif seq is not None:
            size = f"{len(seq)} aa"
        else:
            size = "fetching…"
        return f"{node.key}  ({size}, {tag}){sec}"

    def _pos_range(self, nodes: list[AxisNode]) -> str:
        postions = sorted(n.key.pos for n in nodes if n.kind == "residue")
        if not postions:
            return ""
        if len(postions) == 1:
            return f"  pos {postions[0]}"
        return f"  pos {postions[0]}-{postions[-1]}"

    def _single_pair(self, rb: int, cb: int) -> tuple[float | None, int, int] | None:
        """For a cell between one protein and one protein, or one residue
        and one residue: (its score or None, row slot, column slot). None
        for every other cell, whose value is combined from several links."""
        (r_lo, r_hi), (c_lo, c_hi) = self._blocks[rb], self._blocks[cb]
        if r_hi - r_lo != 1 or c_hi - c_lo != 1:
            return None
        a, b = self._slots[r_lo], self._slots[c_lo]
        if a.kind != b.kind or a.kind == "group":
            return None
        sparse = self._res_sparse if a.kind == "residue" else self._prot_sparse
        return sparse.get((a.key, b.key)), r_lo, c_lo

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
        proteins or residues the two directions are listed, and the
        combined value only if both directions are shown."""
        single = self._single_pair(rb, cb) is not None
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
        pair = self._single_pair(rb, cb)
        if pair is not None:
            score, row, col = pair
            if score is None:
                return "(no crosslink)"
            if not self._quadrant_shown(row, col):
                hidden = f" (hidden: {self._quadrant_of(row, col)} quadrant)"
            elif not self._in_cutoff(score):
                hidden = " (hidden by cutoff)"
            else:
                hidden = ""
            return f"score: {score:.4f}{hidden}"
        r_lo, c_lo = self._blocks[rb][0], self._blocks[cb][0]
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
        if not self._slots:
            return
        try:
            if not self._headers_hidden:   # hidden while zoomed out into blocks (not for residues)
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
        _, _, r, c, scores, rgb = self._visible_cells()
        for rb, cb, score, colour in zip(r.tolist(), c.tolist(), scores.tolist(), rgb.tolist()):
            yield cb * cp - ox, rb * cp - oy, tuple(colour), self._cell_text(score, cp)

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
            img = Image.new("L", (self._name_w, 2 * lh), 0)
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

    def _draw_name_lines(self, draw: ImageDraw.ImageDraw, lines: list[str],
                         cx: float, cy: float, fill) -> None:
        """Write the lines of a section name centred on (cx, cy)."""
        for k, line in enumerate(lines):
            y = cy + (k - (len(lines) - 1) / 2) * self._ch
            draw.text((cx, y), line, font=self._pil_lf, anchor="mm", fill=fill)

    def _section_name_top(self, draw: ImageDraw.ImageDraw, label: str, x0: int, x1: int,
                          level: int = 0) -> None:
        """A group or protein name above the columns x0:x1, in up to two
        lines, with its box, in band row `level`. The name is centred on the
        part that is in the picture, so a band wider than the view still
        shows its name."""
        band = self._band_h
        y0 = level * band
        draw.rectangle([x0, y0, x1 - 1, y0 + band - 1], outline=_SEP)
        x0, x1 = max(x0, 0), min(x1, draw.im.size[0])
        lines = _wrap_label(label, self._pil_lf, x1 - x0 - 4)
        self._draw_name_lines(draw, lines, (x0 + x1) / 2, y0 + band / 2, _FG)

    def _section_name_left(self, img: Image.Image, label: str, y0: int, y1: int,
                           level: int = 0) -> None:
        """A group or protein name left of the rows y0:y1, in up to two
        lines read from bottom to top, with its box, in band column `level`.
        Like on top, the name is centred on the part that is in the picture."""
        band = self._band_h
        x0 = level * band
        ImageDraw.Draw(img).rectangle([x0, y0, x0 + band - 1, y1 - 1], outline=_SEP)
        y0, y1 = max(y0, 0), min(y1, img.height)
        length = y1 - y0
        if length < 1:
            return
        key = (label, length)
        mask = self._vname_cache.get(key)
        if mask is None:
            if len(self._vname_cache) > 512:   # every zoom level makes new ones
                self._vname_cache.clear()
            flat = Image.new("L", (length, band), 0)
            lines = _wrap_label(label, self._pil_lf, length - 4)
            self._draw_name_lines(ImageDraw.Draw(flat), lines, length / 2, band / 2, 255)
            mask = self._vname_cache[key] = flat.transpose(Image.Transpose.ROTATE_90)
        img.paste(_FG, (x0, y0), mask)

    def _band_runs(self, first_slots: list[int], level: int, split: int):
        """Runs of neighbouring positions whose slots have the same open node
        at band row `level`, as (start, end, node); node is None where there
        is no node that deep. A run never crosses `split`. first_slots[k] is
        the slot that stands for position k (the first slot of a block)."""
        n = len(first_slots)
        has_sep = 0 < split < n

        def node_at(k: int) -> AxisNode | None:
            path = self._slot_paths[first_slots[k]]
            return path[level] if level < len(path) else None

        i = 0
        while i < n:
            node = node_at(i)
            j = i + 1
            while j < n and not (has_sep and j == split) and node_at(j) is node:
                j += 1
            yield i, j, node
            i = j

    def _draw_bands(self, img: Image.Image, first_slots: list[int], split: int,
                    cp: int, offset: int, axis: str) -> None:
        """The names of the open groups and proteins, one band row per
        level, above the columns (axis "col") or left of the rows ("row").
        Position k starts at pixel k * cp + offset along the axis."""
        draw = ImageDraw.Draw(img)
        limit = img.width if axis == "col" else img.height
        for level in range(self._band_levels):
            for i, j, node in self._band_runs(first_slots, level, split):
                a, b = i * cp + offset, j * cp + offset
                if node is None or b <= 0 or a >= limit:
                    continue
                if axis == "col":
                    self._section_name_top(draw, node.name, a, b, level)
                else:
                    self._section_name_left(img, node.name, a, b, level)

    def _draw_tracks(self, draw: ImageDraw.ImageDraw, masks: np.ndarray, positions,
                     cp: int, offset: int, across: int, axis: str) -> None:
        """The annotation tracks: for every shown feature type a strip of
        TRACK_PX, starting `across` pixels from the top (columns) or the
        left (rows), with the positions that have the feature coloured."""
        for t, kind in enumerate(self._shown_features):
            a0 = across + TRACK_GAP_PX + t * (TRACK_PX + TRACK_GAP_PX)
            colour = _FEATURE_COLOURS[kind]
            for k in positions:
                if masks[t, k]:
                    p0 = k * cp + offset
                    if axis == "col":
                        draw.rectangle([p0, a0, p0 + cp - 1, a0 + TRACK_PX - 1], fill=colour)
                    else:
                        draw.rectangle([a0, p0, a0 + TRACK_PX - 1, p0 + cp - 1], fill=colour)

    def _draw_col_headers(self) -> None:
        c = self._ch_canvas
        w, h = max(1, c.winfo_width()), max(1, self._header_h)
        ox = int(c.canvasx(0))
        cp, font = self._cell_px, self._pil_lf
        n = len(self._blocks)
        split = self._split_block
        has_sep = 0 < split < n
        _, cols = self._visible_blocks()

        img = Image.new("RGB", (w, h), _BG)
        draw = ImageDraw.Draw(img)

        # group and protein names on top, one row per level. A name too
        # narrow even for two lines gets shortened, or none, but the box
        # stays so the grouping is still visible.
        self._draw_bands(img, [s for s, _ in self._blocks], split, cp, -ox, "col")
        y = self._label_base   # the line under the names

        if has_sep:
            draw.line([(split * cp - ox, 0), (split * cp - ox, h)], fill=_SEP, width=1)

        if not self._names_shown(self._overview):
            pass   # zoomed out into blocks without residues: only the band names
        elif self._overview:
            # Zoomed out into blocks: too many slots to name one by one, so
            # only the name of one slot per block is shown, as sparse
            # ticks spaced out enough not to overlap.
            last_x = -math.inf
            for cb in cols:
                s, _ = self._blocks[cb]
                label = self._col_labels[s]
                if not label:
                    continue
                x = cb * cp + cp // 2 - ox
                half = font.getlength(label) / 2
                if x - half < last_x:
                    continue
                draw.line([(x, y - AXIS_TICK_PX), (x, y)], fill=_SEP, width=1)
                draw.text((x, y - AXIS_TICK_PX - 1), label, font=font, anchor="mb", fill=_FG)
                last_x = x + half + AXIS_TICK_GAP_PX
        else:
            # Names going up to the right (every block is one slot here).
            # Names of columns left of the view can still reach into it.
            reach = math.ceil(self._slant_w / cp)
            for cb in range(max(0, cols.start - reach), cols.stop):
                s, _ = self._blocks[cb]
                if self._col_labels[s]:
                    label, lx, ly = self._slanted_label(self._col_labels[s])
                    img.paste(_FG, (cb * cp + cp // 2 - ox - lx, y - LABEL_GAP_PX - ly), label)

        draw.line([(0, y), (w, y)], fill=_SEP, width=1)
        if has_sep:
            draw.text((split * cp - ox, y), "+", font=font, anchor="mm", fill=_SEP)
        if self._tracks_px:
            self._draw_tracks(draw, self._block_tracks, cols, cp, -ox, h - self._tracks_px, "col")

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
        right = w - self._tracks_px   # where the names end; the tracks come after

        img = Image.new("RGB", (w, h), _BG)
        draw = ImageDraw.Draw(img)
        # the same band names as on top, so the groups of the rows can be
        # told apart as well
        self._draw_bands(img, [s for s, _ in self._blocks], split, cp, -oy, "row")
        if not self._names_shown(self._overview):
            pass   # zoomed out into blocks without residues: only the band names
        elif self._overview:
            # Zoomed out into blocks: show the name of one slot per block,
            # as sparse ticks spaced out enough not to overlap.
            last_y = -math.inf
            gap = self._ch + 2
            for rb in rows:
                s, _ = self._blocks[rb]
                label = self._row_labels[s]
                y = rb * cp + cp // 2 - oy
                if not label or y - gap / 2 < last_y:
                    continue
                draw.text((right - AXIS_TICK_PX - 2, y), label, font=self._pil_lf, anchor="rm", fill=_FG)
                draw.line([(right - AXIS_TICK_PX, y), (right, y)], fill=_SEP, width=1)
                last_y = y + gap / 2
        else:
            for rb in rows:   # every block is one slot here
                label = self._row_labels[self._blocks[rb][0]]
                if label:
                    img.paste(_FG, (self._section_band, rb * cp + cp // 2 - oy - self._ch),
                              self._label_image(label))
        if self._tracks_px:
            self._draw_tracks(draw, self._block_tracks, rows, cp, -oy, right, "row")

        if 0 < split < n and 0 <= split * cp - oy <= h:
            draw.line([(0, split * cp - oy), (w, split * cp - oy)], fill=_SEP, width=1)

        self._rl_photo = ImageTk.PhotoImage(img)
        c.delete("all")
        c.create_image(0, oy, image=self._rl_photo, anchor="nw")

    # -------------------------------------------------------------------------
    # Exporting the whole matrix as a PNG
    #
    # Unlike the screen, the export always shows every slot in its own cell
    # (no blocks), with the axes opened as they are on screen. The score
    # cutoff and the quadrant checkboxes still apply.
    # -------------------------------------------------------------------------

    def _export(self) -> None:
        if not self._slots:
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
        n = len(self._slots)
        split = self._split
        has_sep = 0 < split < n

        font_lbl = _find_pil_font(LABEL_FONT_SIZE + 2)
        font_cell = _find_pil_font(max(6, cp - 2))
        # every slot has its own cell here, so all names are shown, even
        # while the screen is zoomed out into blocks
        lw, base, hh = self._label_sizes(overview=False)
        tracks = self._tracks_px
        # the slanted names of the last columns reach past the matrix on the right
        total_w, matrix_bottom = lw + n * cp + self._slant_w, hh + n * cp
        legend_w = min(max(total_w - 8, 1), 360)
        total_h = matrix_bottom + 2 * self._ch + 24   # room for the colour legend

        img = Image.new("RGB", (total_w, total_h), _BG)
        draw = ImageDraw.Draw(img)

        # which side of a crosslink rows and columns are, in the empty corner
        for k, text in enumerate(reversed(_AXIS_TITLES)):
            draw.text((lw - 4, base - 2 - k * self._ch), text, font=font_lbl, anchor="rs", fill=_SEP)

        # group and protein names on top and on the left, then the slot
        # names going up to the right
        every_slot = list(range(n))
        self._draw_bands(img, every_slot, split, cp, lw, "col")
        self._draw_bands(img, every_slot, split, cp, hh, "row")
        y = base
        for ci, lb in enumerate(self._col_labels):
            if lb:
                label, lx, ly = self._slanted_label(lb)
                img.paste(_FG, (lw + ci * cp + cp // 2 - lx, y - LABEL_GAP_PX - ly), label)
        draw.line([(lw, y), (total_w, y)], fill=_SEP, width=1)
        if has_sep:
            draw.text((lw + split * cp, y), "+", font=font_lbl, anchor="mm", fill=_SEP)

        # row labels
        for ri in range(n):
            draw.text((self._section_band + 2, hh + ri * cp + cp // 2), self._row_labels[ri],
                      font=font_lbl, anchor="lm", fill=_FG)

        # annotation tracks next to the matrix
        if tracks:
            self._draw_tracks(draw, self._slot_tracks, every_slot, cp, lw, hh - tracks, "col")
            self._draw_tracks(draw, self._slot_tracks, every_slot, cp, hh, lw - tracks, "row")

        # lines between targets and decoys
        if has_sep:
            sx, sy = lw + split * cp, hh + split * cp
            draw.line([(0, sy), (lw, sy)], fill=_SEP, width=1)
            draw.line([(sx, 0), (sx, matrix_bottom)], fill=_SEP, width=1)
            draw.line([(lw, sy), (total_w, sy)], fill=_SEP, width=1)

        # the cells themselves: every crosslink that is shown, or where a
        # cell holds several (folded, or a closed group or protein), their
        # combined value
        keep = self._shown_links()
        rows, cols, scores = self._link_rows[keep], self._link_cols[keep], self._link_scores[keep]
        ref = self._score_sorted   # what the colours are measured against
        if self._folded or not self._uniform:
            rows, cols, grouped, starts = _group_cells(rows, cols, scores, max(1, n), self._folded)
            scores = _aggregate_groups(grouped, starts, self._agg_method)
            if self._agg_method in ("count", "rss"):
                ref = np.sort(scores)   # like on screen, see _reduce_block_scores
        colours = self._colours_for(np.asarray(scores, dtype=float), ref)
        for ri, ci, score, colour in zip(rows.tolist(), cols.tolist(),
                                         np.asarray(scores).tolist(), colours.tolist()):
            colour = tuple(colour)
            x0, y0 = lw + ci * cp, hh + ri * cp
            draw.rectangle([x0, y0, x0 + cp - 1, y0 + cp - 1], fill=colour)
            text = self._cell_text(score, cp)
            if text:
                draw.text((x0 + cp // 2, y0 + cp // 2), text, font=font_cell, anchor="mm",
                          fill=_auto_fg(colour))

        # what the colours mean, under the matrix
        caption = self._colour_caption(folded=self._folded, merged=not self._uniform)
        img.paste(self._legend_image(legend_w, ref, caption), (4, matrix_bottom + 8))
        return img


# -----------------------------------------------------------------------------
# Starting the app
# -----------------------------------------------------------------------------

def main() -> None:
    MatrixApp(sys.argv[1] if len(sys.argv) > 1 else None).mainloop()


if __name__ == "__main__":
    main()

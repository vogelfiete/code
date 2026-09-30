"""Print a sparse adjacency matrix for the top-N most-confident unique proteins
(or residues) found in a crosslink CSV file. See docs/usage.md for the full CLI
reference.

Pipeline:  CSV -> build_matrix / build_residue_matrix -> sort_proteins /
sort_residues -> print_matrix.  The GUI (gui/matrix_app.py) reuses everything
up to the rendering step."""
from __future__ import annotations

import argparse
import io
import os
import sys
import tracemalloc
import warnings
from typing import NamedTuple

# Ensure Unicode characters render correctly on Windows consoles
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", line_buffering=True)

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from fetch_group_order import order_by_kegg, order_by_string, short_accession
from xlms.io import read_fasta

DEFAULT_N = 1000

COL_W = 10   # column header width, in chars
ROW_W = 12   # row label width, in chars
CELL_W = 2   # one dot/dash plus a space

N_LEVELS = 4
DOTS = ["·", "○", "●", "⬤"]   # small to large, keep in sync with N_LEVELS
EMPTY = " "

_DECOY_NAME_PREFIXES = ("decoy", "rev_", "contam_")


# --------------------------------------------------------------------------
# Axis items
# --------------------------------------------------------------------------

class ResidueId(NamedTuple):
    """One crosslinked residue: a protein plus a sequence position.

    Hashable and orderable by (protein, pos), so it's a drop-in for the
    protein-name strings used as the axis unit and sparse-dict key
    everywhere in this module.
    """
    protein: str
    pos: int

    def __str__(self) -> str:
        return f"{self.protein}:{self.pos}"


def protein_of(item) -> str:
    """Parent protein of an axis item (a protein name or a ResidueId)."""
    return item.protein if isinstance(item, ResidueId) else item


# --------------------------------------------------------------------------
# Small value parsers
# --------------------------------------------------------------------------

def _truncate(name: str, width: int) -> str:
    return name[:width]


def _is_ambiguous(value) -> bool:
    """Search tools report several candidates as e.g. "P1;P2" or "114;53"."""
    return ";" in str(value)


def _parse_position(val) -> int | None:
    """Residue position as int, or None if missing or ambiguous."""
    if pd.isna(val) or _is_ambiguous(val):
        return None
    try:
        return int(float(val))
    except (TypeError, ValueError):
        return None


def _combine_pep_link(pep_val, link_val) -> int | None:
    """Absolute position = peptide start + 1-based link offset - 1."""
    pep = _parse_position(pep_val)
    link = _parse_position(link_val)
    if pep is None or link is None:
        return None
    return pep + link - 1


def _to_bool(val) -> bool:
    if isinstance(val, bool):
        return val
    return str(val).strip().lower() in ("true", "1", "yes")


def _is_decoy_by_name(name: str) -> bool:
    return name.lower().startswith(_DECOY_NAME_PREFIXES)


# --------------------------------------------------------------------------
# Column detection (all names are matched case-insensitively)
# --------------------------------------------------------------------------

def _detect_decoy_cols(lower: dict[str, str]) -> dict[str, str]:
    """Format "B" = Decoy1/Decoy2 (per protein), format "A" = isTT/isDD
    (per link), or "none"."""
    if "decoy1" in lower and "decoy2" in lower:
        return {"format": "B", "d1": lower["decoy1"], "d2": lower["decoy2"]}
    if "istt" in lower and "isdd" in lower:
        return {"format": "A", "istt": lower["istt"], "isdd": lower["isdd"]}
    return {"format": "none"}


def _detect_position_cols(lower: dict[str, str]) -> dict[str, str] | None:
    """Format "direct": the column holds the absolute residue position
    (SeqPos1/2, Pos1/2 or Position1/2). Format "pep_link": PepPos1/2 (peptide
    start) plus LinkPos1/2 (1-based offset in the peptide). None if absent."""
    for a, b in (("seqpos1", "seqpos2"), ("pos1", "pos2"), ("position1", "position2")):
        if a in lower and b in lower:
            return {"format": "direct", "pos1": lower[a], "pos2": lower[b]}
    if all(k in lower for k in ("peppos1", "peppos2", "linkpos1", "linkpos2")):
        return {
            "format": "pep_link",
            "pep1": lower["peppos1"], "link1": lower["linkpos1"],
            "pep2": lower["peppos2"], "link2": lower["linkpos2"],
        }
    return None


def _positions(df: pd.DataFrame, pos_info: dict[str, str]) -> tuple[list, list]:
    """Per-row absolute residue positions (int or None) for both link ends."""
    if pos_info["format"] == "direct":
        return ([_parse_position(v) for v in df[pos_info["pos1"]]],
                [_parse_position(v) for v in df[pos_info["pos2"]]])
    return ([_combine_pep_link(p, l) for p, l in zip(df[pos_info["pep1"]], df[pos_info["link1"]])],
            [_combine_pep_link(p, l) for p, l in zip(df[pos_info["pep2"]], df[pos_info["link2"]])])


# --------------------------------------------------------------------------
# Decoy classification
# --------------------------------------------------------------------------

def _classify_proteins(
    df: pd.DataFrame,
    names1: list[str],
    names2: list[str],
    decoy_info: dict[str, str],
    proteins: set[str],
) -> dict[str, bool]:
    """Return {protein: is_decoy} for every protein in `proteins`.
    Rows are applied in order, so for conflicting rows the last one wins."""
    result: dict[str, bool] = {}
    fmt = decoy_info["format"]

    if fmt == "B":
        # each end of a link carries its own decoy flag
        for p1, p2, d1, d2 in zip(names1, names2, df[decoy_info["d1"]], df[decoy_info["d2"]]):
            if p1 in proteins:
                result[p1] = _to_bool(d1)
            if p2 in proteins:
                result[p2] = _to_bool(d2)

    elif fmt == "A":
        # TT links make both ends targets, DD links make both ends decoys
        for p1, p2, tt, dd in zip(names1, names2, df[decoy_info["istt"]], df[decoy_info["isdd"]]):
            if _to_bool(tt):
                label = False
            elif _to_bool(dd):
                label = True
            else:
                continue  # TD links say nothing about which end is the decoy
            for p in (p1, p2):
                if p in proteins:
                    result[p] = label
        # proteins only seen in TD links: guess from the name
        for p in proteins - result.keys():
            result[p] = _is_decoy_by_name(p)

    # unclassified proteins (or fmt == "none") default to target
    for p in proteins:
        result.setdefault(p, False)
    return result


# --------------------------------------------------------------------------
# Matrix construction
# --------------------------------------------------------------------------

def build_matrix(path: str, n: int = DEFAULT_N):
    """Protein level. Returns (proteins, sparse, is_decoy):
        proteins — N unique protein names, targets first then decoys,
                   confidence order preserved within each group
        sparse   — {(protein1, protein2): best_score}, directional: the
                   key is (Protein1, Protein2) exactly as in the CSV row,
                   so A-B and B-A are separate entries (row A/col B vs
                   row B/col A)
        is_decoy — {protein: bool}
    """
    return _build(path, n, level="protein")


def build_residue_matrix(path: str, n: int = DEFAULT_N):
    """Residue level — same as build_matrix, but the axis unit is a
    ResidueId(protein, position). Needs residue-position columns."""
    return _build(path, n, level="residue")


def _build(path: str, n: int, level: str):
    df = pd.read_csv(path)
    lower = {c.lower(): c for c in df.columns}
    missing = [c for c in ("protein1", "protein2", "score") if c not in lower]
    if missing:
        raise ValueError(f"CSV is missing required column(s): {', '.join(missing)}")
    decoy_info = _detect_decoy_cols(lower)
    pos_info = _detect_position_cols(lower) if level == "residue" else None
    if level == "residue" and pos_info is None:
        raise ValueError(
            "Residue-level mode requires residue-position columns "
            "(e.g. SeqPos1/SeqPos2, or PepPos1/PepPos2 + LinkPos1/LinkPos2). "
            "Found columns: " + ", ".join(df.columns)
        )

    df = df.sort_values(lower["score"], ascending=False)
    names1 = [str(v) for v in df[lower["protein1"]]]
    names2 = [str(v) for v in df[lower["protein2"]]]
    scores = [float(v) for v in df[lower["score"]]]

    # the two axis items of every row; None where a residue position is unusable
    if pos_info is None:
        items1, items2 = names1, names2
    else:
        pos1, pos2 = _positions(df, pos_info)
        items1 = [ResidueId(p, x) if x is not None else None for p, x in zip(names1, pos1)]
        items2 = [ResidueId(p, x) if x is not None else None for p, x in zip(names2, pos2)]

    # pass 1: walk rows highest-score-first until we have N unique items
    selected: dict = {}
    n_skipped = 0
    for name1, name2, a, b in zip(names1, names2, items1, items2):
        if _is_ambiguous(name1) or _is_ambiguous(name2):
            continue
        if a is None or b is None:
            n_skipped += 1
            continue
        selected.setdefault(a)
        selected.setdefault(b)
        if len(selected) >= n:
            break
    if n_skipped:
        warnings.warn(f"{n_skipped} row(s) skipped: missing/ambiguous residue position")
    if not selected:
        raise ValueError("No usable crosslinks were found in the CSV.")

    protein_is_decoy = _classify_proteins(
        df, names1, names2, decoy_info, {protein_of(x) for x in selected}
    )
    is_decoy = {x: protein_is_decoy[protein_of(x)] for x in selected}
    # targets first, then decoys; sorted() is stable, so confidence order is kept
    items = sorted(selected, key=lambda x: is_decoy[x])

    # pass 2: best score per directed item pair (row = Protein1, col = Protein2)
    wanted = set(items)
    sparse: dict[tuple, float] = {}
    for a, b, score in zip(items1, items2, scores):
        if a in wanted and b in wanted and score > sparse.get((a, b), float("-inf")):
            sparse[(a, b)] = score

    return items, sparse, is_decoy


def add_unlinked_residues(
    residues: list[ResidueId],
    residue_is_decoy: dict[ResidueId, bool],
    sequences: dict[str, str],
) -> tuple[list[ResidueId], dict[ResidueId, bool]]:
    """Also include every other residue position of each protein already
    present (needs the full sequence). The added residues have no score
    entries, so they render as empty cells."""
    proteins = list(dict.fromkeys(r.protein for r in residues))
    protein_is_decoy = {r.protein: residue_is_decoy[r] for r in residues}

    missing = [p for p in proteins if p not in sequences]
    if missing:
        warnings.warn(
            f"{len(missing)} protein(s) not found in FASTA — only linked "
            f"residues shown for them: {missing}"
        )

    expanded = list(residues)
    is_decoy = dict(residue_is_decoy)
    existing = set(residues)
    for protein in proteins:
        for pos in range(1, len(sequences.get(protein, "")) + 1):
            rid = ResidueId(protein, pos)
            if rid not in existing:
                existing.add(rid)
                expanded.append(rid)
                is_decoy[rid] = protein_is_decoy[protein]
    return expanded, is_decoy


# --------------------------------------------------------------------------
# Axis ordering
# --------------------------------------------------------------------------
# Each ordering takes (group, sequences, species) and returns
# (reordered group, {item: section_label} or None).

def _order_alpha(group, sequences, species):
    return sorted(group), None


def _order_sequence(group, sequences, species):
    """Cluster by 3-mer composition (crude but fast proxy for sequence similarity)."""
    import numpy as np
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import cdist, squareform

    sequences = sequences or {}
    present = [p for p in group if p in sequences]
    missing = [p for p in group if p not in sequences]
    if missing:
        warnings.warn(f"{len(missing)} protein(s) not found in FASTA, appended at end: {missing}")
    if len(present) < 2:
        return present + missing, None

    k = 3
    kmer_index: dict[str, int] = {}
    for p in present:
        seq = sequences[p]
        for i in range(len(seq) - k + 1):
            kmer_index.setdefault(seq[i:i + k], len(kmer_index))

    freq = np.zeros((len(present), len(kmer_index)))
    for row, p in enumerate(present):
        seq = sequences[p]
        for i in range(len(seq) - k + 1):
            freq[row, kmer_index[seq[i:i + k]]] += 1
        total = freq[row].sum()
        if total > 0:
            freq[row] /= total

    D = cdist(freq, freq, metric="cosine")
    np.fill_diagonal(D, 0.0)
    order = leaves_list(linkage(squareform(D), method="average"))
    return [present[i] for i in order] + missing, None


def _order_size(group, sequences, species):
    """Longest sequence first; proteins without a sequence go last."""
    def key(p: str) -> int:
        return -len(sequences[p]) if sequences and p in sequences else 1
    return sorted(group, key=key), None


def _order_pathway(group, sequences, species):
    return order_by_kegg(group, species)


def _order_complex(group, sequences, species):
    return order_by_string(group, species)


_ORDERINGS = {
    "alpha": _order_alpha,
    "sequence": _order_sequence,
    "pathway": _order_pathway,
    "complex": _order_complex,
    "size": _order_size,
}
ORDER_MODES = ["confidence", *_ORDERINGS]   # "confidence" = keep the input order


def _apply_order(group, mode, sequences, species=9606):
    ordering = _ORDERINGS.get(mode)
    if ordering is None or len(group) < 2:
        return group, None
    return ordering(group, sequences, species)


def sort_proteins(
    proteins: list[str],
    protein_is_decoy: dict[str, bool],
    mode: str,
    sequences: dict[str, str] | None = None,
    species: int = 9606,
) -> tuple[list[str], dict[str, str] | None]:
    """Reorder proteins by `mode`, keeping targets before decoys."""
    targets = [p for p in proteins if not protein_is_decoy.get(p, False)]
    decoys = [p for p in proteins if protein_is_decoy.get(p, False)]
    t_ordered, t_sections = _apply_order(targets, mode, sequences, species)
    d_ordered, d_sections = _apply_order(decoys, mode, sequences, species)
    sections = {**(t_sections or {}), **(d_sections or {})} or None
    return t_ordered + d_ordered, sections


def sort_residues(
    residues: list[ResidueId],
    residue_is_decoy: dict[ResidueId, bool],
    mode: str,
    sequences: dict[str, str] | None = None,
    species: int = 9606,
) -> tuple[list[ResidueId], dict[ResidueId, str]]:
    """Order the parent proteins with sort_proteins, then list each protein's
    residues by position. The parent protein is always the section label."""
    # first-occurrence order keeps mode="confidence" meaningful
    proteins = list(dict.fromkeys(r.protein for r in residues))
    protein_is_decoy = {r.protein: residue_is_decoy[r] for r in residues}

    ordered_proteins, _ = sort_proteins(proteins, protein_is_decoy, mode, sequences, species)

    by_protein: dict[str, list[ResidueId]] = {}
    for r in residues:
        by_protein.setdefault(r.protein, []).append(r)
    ordered = [r for p in ordered_proteins for r in sorted(by_protein[p], key=lambda r: r.pos)]
    return ordered, {r: r.protein for r in ordered}


# --------------------------------------------------------------------------
# Console rendering
# --------------------------------------------------------------------------

def _build_thresholds(sparse: dict[tuple, float]) -> list[float]:
    """Upper bounds of the first N_LEVELS-1 dot sizes (equal-width bins)."""
    max_score = max(sparse.values())
    return [max_score * (i + 1) / N_LEVELS for i in range(N_LEVELS - 1)]


def _score_to_dot(score: float, thresholds: list[float]) -> str:
    for dot, t in zip(DOTS, thresholds):
        if score <= t:
            return dot
    return DOTS[-1]


def _col_separator_line(n: int, split: int) -> str:
    """Horizontal rule with a '+' at the target/decoy boundary."""
    dashes = "-" * CELL_W
    if 0 < split < n:
        return dashes * split + "-+" + dashes * (n - split)
    return dashes * n


def _section_header_line(items: list, sections: dict, split: int, row_prefix: str) -> str:
    """`=label=` spans over each run of items that share a section."""
    n = len(items)
    has_sep = 0 < split < n
    line = row_prefix
    i = 0
    while i < n:
        if has_sep and i == split:
            line += "|"
        label = sections.get(items[i], "")
        j = i + 1
        while j < n and not (has_sep and j == split) and sections.get(items[j], "") == label:
            j += 1
        span = (j - i) * CELL_W
        line += ("=" + label[:span - 2].center(span - 2) + "=")[:span].ljust(span)
        i = j
    return line


def print_matrix(
    items: list,
    sparse: dict[tuple, float],
    is_decoy: dict,
    sections: dict | None = None,
    *,
    col_label_fn=None,
    row_label_fn=None,
    legend_source: str = "col",
) -> None:
    """Render the dot matrix to stdout. `items` are protein names or
    ResidueIds. legend_source ("col" or "row") picks which label set the
    truncated-name legend is built from."""
    col_label_fn = col_label_fn or (lambda p: _truncate(p, COL_W))
    row_label_fn = row_label_fn or (lambda p: _truncate(p, ROW_W))

    n = len(items)
    col_labels = [col_label_fn(p) for p in items]
    row_labels = [row_label_fn(p) for p in items]
    thresholds = _build_thresholds(sparse) if sparse else []
    split = next((i for i, p in enumerate(items) if is_decoy.get(p, False)), n)
    has_sep = 0 < split < n
    row_prefix = " " * (ROW_W + 2)

    if sections:
        print(_section_header_line(items, sections, split, row_prefix))

    # column labels are printed vertically, one character row at a time
    for char_idx in range(max(len(lb) for lb in col_labels)):
        line = row_prefix
        for i, lb in enumerate(col_labels):
            if has_sep and i == split:
                line += "|"
            line += (lb[char_idx] if char_idx < len(lb) else " ").center(CELL_W)
        print(line)
    print(row_prefix + _col_separator_line(n, split))

    for row, pi in enumerate(items):
        if has_sep and row == split:
            print("-" * (ROW_W + 2) + _col_separator_line(n, split))
        line = row_labels[row].ljust(ROW_W) + "  "
        for col, pj in enumerate(items):
            if has_sep and col == split:
                line += "|"
            score = sparse.get((pi, pj))
            line += (EMPTY if score is None else _score_to_dot(score, thresholds)).ljust(CELL_W)
        print(line)

    if sparse:
        print("\nScale (max score {:.2f}):".format(max(sparse.values())))
        lower = 0.0
        for dot, t in zip(DOTS, thresholds):
            print(f"  {dot}   {lower:.2f} - {t:.2f}")
            lower = t
        print(f"  {DOTS[-1]}   {lower:.2f}+")

    legend_labels = col_labels if legend_source == "col" else row_labels
    truncated = [(short, str(full)) for short, full in zip(legend_labels, items) if short != str(full)]
    if truncated:
        print("\nLegend (truncated -> full name):")
        for short, full in truncated:
            print(f"  {short}  ->  {full}")


def _format_bytes(n_bytes: int) -> str:
    if n_bytes < 1024:
        return f"{n_bytes} B"
    if n_bytes < 1024 ** 2:
        return f"{n_bytes / 1024:.2f} KB"
    return f"{n_bytes / 1024 ** 2:.2f} MB"


def _sparse_size_bytes(sparse: dict[tuple, float]) -> int:
    total = sys.getsizeof(sparse)
    for (a, b), v in sparse.items():
        total += sys.getsizeof(a) + sys.getsizeof(b) + sys.getsizeof(v) + sys.getsizeof((a, b))
    return total


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Print XL-MS adjacency matrix")
    parser.add_argument("csv", help="Path to the crosslink CSV file")
    parser.add_argument("--level", choices=["protein", "residue"], default="protein",
                        help="Axis granularity: protein (default) or residue "
                             "(requires residue-position columns in the CSV)")
    parser.add_argument("--n", type=int, default=DEFAULT_N,
                        help=f"Number of unique proteins/residues (default {DEFAULT_N})")
    parser.add_argument("--order", choices=ORDER_MODES, default="confidence",
                        help="Axis ordering: confidence (default), alpha, "
                             "sequence, pathway (KEGG), complex (STRING), "
                             "size (longest sequence first)")
    parser.add_argument("--fasta", default=None,
                        help="FASTA file path (required for --order sequence "
                             "and --include-unlinked)")
    parser.add_argument("--species", type=int, default=9606,
                        help="NCBI taxonomy ID for external DB queries (default 9606 = human)")
    parser.add_argument("--include-unlinked", action="store_true",
                        help="Also show residues with no crosslinks, not just linked "
                             "ones (residue mode only, requires --fasta)")
    args = parser.parse_args()

    if args.order == "sequence" and not args.fasta:
        parser.error("--fasta is required when --order sequence")
    if args.include_unlinked and not args.fasta:
        parser.error("--fasta is required when --include-unlinked")

    residue_level = args.level == "residue"
    tracemalloc.start()
    try:
        build = build_residue_matrix if residue_level else build_matrix
        items, sparse, is_decoy = build(args.csv, n=args.n)
    except ValueError as exc:
        parser.error(str(exc))

    sequences = read_fasta(args.fasta) if args.fasta else None
    if residue_level:
        if args.include_unlinked:
            items, is_decoy = add_unlinked_residues(items, is_decoy, sequences or {})
        items, sections = sort_residues(items, is_decoy, args.order, sequences, args.species)
        labels = dict(
            col_label_fn=lambda r: _truncate(str(r.pos), COL_W),
            row_label_fn=lambda r: _truncate(f"{short_accession(r.protein)}:{r.pos}", ROW_W),
            legend_source="row",
        )
    else:
        items, sections = sort_proteins(items, is_decoy, args.order, sequences, args.species)
        labels = {}

    _, peak_traced = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    print_matrix(items, sparse, is_decoy, sections, **labels)

    n_targets = sum(1 for p in items if not is_decoy.get(p, False))
    unit = "Residues" if residue_level else "Proteins"
    print(
        f"\n{unit}: {n_targets} targets, {len(items) - n_targets} decoys"
        f"  |  sparse store: {_format_bytes(_sparse_size_bytes(sparse))} ({len(sparse)} links)"
        f"  |  peak traced: {_format_bytes(peak_traced)}"
    )


if __name__ == "__main__":
    main()

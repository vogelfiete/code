"""
Print a sparse adjacency matrix for the top-N most-confident unique proteins
found in a crosslink CSV file.

Usage:
    python scripts/print_matrix.py <csv_path> [options]

Options:
    --n INT           Number of unique proteins to display (default 1000)
    --order MODE      Axis ordering — one of:
                        confidence  (default) highest-scoring proteins first
                        alpha       alphabetical
                        cluster     hierarchical clustering by crosslink scores
                        sequence    hierarchical clustering by k-mer similarity
                        pathway     group by KEGG pathway (requires --species)
                        complex     group by STRING complex/pathway (requires --species)
    --fasta PATH      FASTA file (required for --order sequence)
    --species INT     NCBI taxonomy ID for pathway/complex ordering (default 9606 = human)

Examples:
    python scripts/print_matrix.py data/links.csv
    python scripts/print_matrix.py data/links.csv --n 200 --order alpha
    python scripts/print_matrix.py data/links.csv --order sequence --fasta data/proteins.fasta
    python scripts/print_matrix.py data/links.csv --order pathway --species 9913
"""
from __future__ import annotations

import argparse
import io
import sys
import os
import tracemalloc
import warnings
from typing import Optional

# Ensure Unicode characters render correctly on Windows consoles
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", line_buffering=True)

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from xlms.io import read_fasta

DEFAULT_N = 1000 # number of unique proteins to display

COL_W = 10   # max chars for column header labels
ROW_W = 12   # max chars for row labels
CELL_W = 2   # chars per data cell: 1 dot/dash + 1 space

N_LEVELS = 4                    # number of distinct dot sizes
DOTS = ["·", "○", "●", "⬤"]   # small → large; change alongside N_LEVELS
EMPTY = " "                     # shown when no link exists between two proteins

_DECOY_NAME_PREFIXES = ("decoy", "rev_", "contam_")

# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------

def _truncate(name: str, width: int) -> str:
    return name[:width] if len(name) > width else name


def _is_ambiguous(name: str) -> bool:
    return ";" in str(name)


def _to_bool(val) -> bool:
    if isinstance(val, bool):
        return val
    return str(val).strip().lower() in ("true", "1", "yes")


def _is_decoy_by_name(name: str) -> bool:
    return name.lower().startswith(_DECOY_NAME_PREFIXES)


# ---------------------------------------------------------------------------
# Column detection
# ---------------------------------------------------------------------------

def _detect_protein_cols(columns: list[str]) -> tuple[str, str]:
    col_map = {c.lower(): c for c in columns}
    return col_map["protein1"], col_map["protein2"]


def _detect_score_col(columns: list[str]) -> str:
    col_map = {c.lower(): c for c in columns}
    return col_map["score"]


def _detect_decoy_cols(columns: list[str]) -> dict[str, str]:
    """
    Returns a dict describing which decoy format is present:
      format "B": {"format":"B", "d1": <Decoy1 col>, "d2": <Decoy2 col>}
      format "A": {"format":"A", "istt": ..., "isdd": ...}
      format "none": {"format":"none"}
    """
    col_map = {c.lower(): c for c in columns}
    if "decoy1" in col_map and "decoy2" in col_map:
        return {"format": "B", "d1": col_map["decoy1"], "d2": col_map["decoy2"]}
    if "istt" in col_map and "isdd" in col_map:
        return {"format": "A", "istt": col_map["istt"], "isdd": col_map["isdd"]}
    return {"format": "none"}


# ---------------------------------------------------------------------------
# Per-protein decoy classification
# ---------------------------------------------------------------------------

def _classify_proteins(
    df: "pd.DataFrame",
    p1_col: str,
    p2_col: str,
    decoy_info: dict[str, str],
    all_proteins: set[str],
) -> dict[str, bool]:
    """Return {protein_name: is_decoy} for every protein in all_proteins."""
    result: dict[str, bool] = {}
    fmt = decoy_info.get("format", "none")

    if fmt == "B":
        d1_col = decoy_info["d1"]
        d2_col = decoy_info["d2"]
        for _, row in df.iterrows():
            p1, p2 = str(row[p1_col]), str(row[p2_col])
            if p1 in all_proteins:
                result[p1] = _to_bool(row[d1_col])
            if p2 in all_proteins:
                result[p2] = _to_bool(row[d2_col])

    elif fmt == "A":
        istt_col = decoy_info["istt"]
        isdd_col = decoy_info["isdd"]
        for _, row in df.iterrows():
            p1, p2 = str(row[p1_col]), str(row[p2_col])
            istt = _to_bool(row[istt_col])
            isdd = _to_bool(row[isdd_col])
            if istt:
                if p1 in all_proteins:
                    result[p1] = False  # target
                if p2 in all_proteins:
                    result[p2] = False
            elif isdd:
                if p1 in all_proteins:
                    result[p1] = True   # decoy
                if p2 in all_proteins:
                    result[p2] = True
        # Proteins that only appear in TD links: name-based fallback
        for p in all_proteins:
            if p not in result:
                result[p] = _is_decoy_by_name(p)

    # fmt == "none" or any unclassified protein: default to target
    for p in all_proteins:
        result.setdefault(p, False)

    return result


# ---------------------------------------------------------------------------
# Score → dot mapping
# ---------------------------------------------------------------------------

def _build_thresholds(sparse: dict[tuple[str, str], float]) -> list[float]:
    max_score = max(sparse.values())
    return [max_score * (i + 1) / N_LEVELS for i in range(N_LEVELS - 1)]


def _score_to_dot(score: float, thresholds: list[float]) -> str:
    for i, t in enumerate(thresholds):
        if score <= t:
            return DOTS[i]
    return DOTS[-1]


# ---------------------------------------------------------------------------
# Axis ordering
# ---------------------------------------------------------------------------

def _apply_order(
    group: list[str],
    sparse: dict[tuple[str, str], float],
    mode: str,
    sequences: dict[str, str] | None,
    species: int = 9606,
) -> tuple[list[str], dict[str, str] | None]:
    """Return (reordered group, section_label_dict or None)."""
    if len(group) < 2 or mode == "confidence":
        return group, None

    if mode == "alpha":
        return sorted(group), None

    if mode == "cluster":
        import numpy as np
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import squareform

        n = len(group)
        max_score = max(sparse.values()) if sparse else 1.0
        S = np.zeros((n, n))
        for i, pi in enumerate(group):
            for j, pj in enumerate(group):
                S[i, j] = sparse.get((min(pi, pj), max(pi, pj)), 0.0)
        D = max_score - S
        np.fill_diagonal(D, 0.0)
        Z = linkage(squareform(D), method="average")
        order = leaves_list(Z)
        return [group[i] for i in order], None

    if mode == "sequence":
        import numpy as np
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import cdist, squareform

        # Separate proteins that have sequences from those that don't
        missing = [p for p in group if p not in (sequences or {})]
        present = [p for p in group if p in (sequences or {})]
        if missing:
            warnings.warn(
                f"{len(missing)} protein(s) not found in FASTA, appended at end: {missing}"
            )
        if len(present) < 2:
            return present + missing, None

        # Build 3-mer frequency matrix
        k = 3
        all_kmers: dict[str, int] = {}
        for p in present:
            seq = sequences[p]  # type: ignore[index]
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i + k]
                if kmer not in all_kmers:
                    all_kmers[kmer] = len(all_kmers)

        mat = np.zeros((len(present), len(all_kmers)))
        for row_idx, p in enumerate(present):
            seq = sequences[p]  # type: ignore[index]
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i + k]
                mat[row_idx, all_kmers[kmer]] += 1
            total = mat[row_idx].sum()
            if total > 0:
                mat[row_idx] /= total

        D = cdist(mat, mat, metric="cosine")
        np.fill_diagonal(D, 0.0)
        Z = linkage(squareform(D), method="average")
        order = leaves_list(Z)
        return [present[i] for i in order] + missing, None

    if mode == "pathway":
        from fetch_group_order import order_by_kegg
        return order_by_kegg(group, species)

    if mode == "complex":
        from fetch_group_order import order_by_string
        return order_by_string(group, species)

    if mode == "size":
        def _size_key(p: str) -> int:
            return -len(sequences[p]) if sequences and p in sequences else 1
        return sorted(group, key=_size_key), None

    return group, None  # fallback


def _sort_proteins(
    proteins: list[str],
    sparse: dict[tuple[str, str], float],
    protein_is_decoy: dict[str, bool],
    mode: str,
    sequences: dict[str, str] | None = None,
    species: int = 9606,
) -> tuple[list[str], dict[str, str] | None]:
    """Reorder proteins by `mode`, keeping targets before decoys."""
    targets = [p for p in proteins if not protein_is_decoy.get(p, False)]
    decoys  = [p for p in proteins if protein_is_decoy.get(p, False)]
    t_ordered, t_sec = _apply_order(targets, sparse, mode, sequences, species)
    d_ordered, d_sec = _apply_order(decoys,  sparse, mode, sequences, species)
    combined = {**(t_sec or {}), **(d_sec or {})} if (t_sec or d_sec) else None
    return t_ordered + d_ordered, combined


# ---------------------------------------------------------------------------
# Matrix builder
# ---------------------------------------------------------------------------

def build_matrix(
    path: str,
    n: int = DEFAULT_N,
) -> tuple[list[str], dict[tuple[str, str], float], dict[str, bool]]:
    """
    Returns:
        proteins        — N unique protein names, targets first then decoys,
                          confidence order preserved within each group
        sparse          — {(min_prot, max_prot): best_score}
        protein_is_decoy — {protein_name: bool}
    """
    df = pd.read_csv(path)
    columns = list(df.columns)
    p1_col, p2_col = _detect_protein_cols(columns)
    score_col = _detect_score_col(columns)
    decoy_info = _detect_decoy_cols(columns)

    # Keep decoy columns alongside the three core columns
    decoy_cols = [v for k, v in decoy_info.items() if k not in ("format",)]
    keep_cols = list(dict.fromkeys([p1_col, p2_col, score_col] + decoy_cols))
    df = df[keep_cols].copy()
    df.sort_values(score_col, ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)

    # Pass 1 — collect N unique proteins in confidence order, skip ambiguous
    selected: dict[str, None] = {}
    for _, row in df.iterrows():
        p1, p2 = str(row[p1_col]), str(row[p2_col])
        if _is_ambiguous(p1) or _is_ambiguous(p2):
            continue
        selected.setdefault(p1)
        selected.setdefault(p2)
        if len(selected) >= n:
            break
    proteins_unordered = list(selected.keys())

    # Classify each protein as target (False) or decoy (True)
    protein_is_decoy = _classify_proteins(
        df, p1_col, p2_col, decoy_info, set(proteins_unordered)
    )

    # Sort: targets first, then decoys; preserve confidence order within each group
    conf_rank = {p: i for i, p in enumerate(proteins_unordered)}
    proteins = sorted(proteins_unordered, key=lambda p: (int(protein_is_decoy[p]), conf_rank[p]))

    # Pass 2 — build sparse score dict for links between selected proteins
    protein_set = set(proteins)
    sparse: dict[tuple[str, str], float] = {}
    for _, row in df.iterrows():
        p1, p2 = str(row[p1_col]), str(row[p2_col])
        if p1 in protein_set and p2 in protein_set:
            key = (min(p1, p2), max(p1, p2))
            score = float(row[score_col])
            if score > sparse.get(key, float("-inf")):
                sparse[key] = score

    return proteins, sparse, protein_is_decoy


# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------

def _col_separator_line(n: int, split: int) -> str:
    """Horizontal rule with a '+' at the target/decoy boundary."""
    dashes = "-" * CELL_W
    if 0 < split < n:
        return (dashes * split) + "-+" + (dashes * (n - split))
    return dashes * n


def _section_header_line(
    proteins: list[str],
    protein_section: dict[str, str],
    split: int,
    row_prefix: str,
) -> str:
    """Build the section-label line printed above the stacked column headers."""
    n = len(proteins)
    has_sep = 0 < split < n
    line = row_prefix
    i = 0
    while i < n:
        if has_sep and i == split:
            line += "|"
        label = protein_section.get(proteins[i], "")
        j = i + 1
        while j < n:
            if has_sep and j == split:
                break
            if protein_section.get(proteins[j], "") != label:
                break
            j += 1
        count = j - i
        span  = count * CELL_W
        inner      = label[:span - 2] if len(label) > span - 2 else label
        cell_plain = (("=" + inner.center(span - 2, " ") + "=") if span >= 2 else "=")[:span].ljust(span)
        line += cell_plain
        i = j
    return line


def print_matrix(
    proteins: list[str],
    sparse: dict[tuple[str, str], float],
    protein_is_decoy: dict[str, bool],
    protein_section: dict[str, str] | None = None,
) -> None:
    n = len(proteins)
    labels = [_truncate(p, COL_W) for p in proteins]
    thresholds = _build_thresholds(sparse) if sparse else []

    split = next((i for i, p in enumerate(proteins) if protein_is_decoy.get(p, False)), n)
    has_sep = 0 < split < n

    row_prefix = " " * (ROW_W + 2)

    # --- Section label header (pathway / complex modes only) ---
    if protein_section:
        print(_section_header_line(proteins, protein_section, split, row_prefix))

    # --- Column headers (stacked vertically) ---
    max_label_len = max(len(lb) for lb in labels)
    for char_idx in range(max_label_len):
        line = row_prefix
        for i, lb in enumerate(labels):
            if has_sep and i == split:
                line += "|"
            char = lb[char_idx] if char_idx < len(lb) else " "
            line += char.center(CELL_W)
        print(line)

    # Horizontal rule under headers
    print(row_prefix + _col_separator_line(n, split))

    # --- Data rows ---
    for row_idx, pi in enumerate(proteins):
        # Row separator between target and decoy blocks
        if has_sep and row_idx == split:
            print("-" * (ROW_W + 2) + _col_separator_line(n, split))

        row_label = _truncate(pi, ROW_W).ljust(ROW_W)
        line = row_label + "  "
        for j, pj in enumerate(proteins):
            if has_sep and j == split:
                line += "|"
            key = (min(pi, pj), max(pi, pj))
            score: Optional[float] = sparse.get(key)
            dot_char = _score_to_dot(score, thresholds) if score is not None else EMPTY
            line += dot_char.ljust(CELL_W)
        print(line)

    # --- Scale legend ---
    if sparse:
        max_score = max(sparse.values())
        print("\nScale (max score {:.2f}):".format(max_score))
        lower = 0.0
        for i, t in enumerate(thresholds):
            print(f"  {DOTS[i]}   {lower:.2f} - {t:.2f}")
            lower = t
        print(f"  {DOTS[-1]}   {lower:.2f}+")

    # --- Name legend for truncated labels ---
    truncated = [(short, full) for short, full in zip(labels, proteins) if short != full]
    if truncated:
        print("\nLegend (truncated -> full name):")
        for short, full in truncated:
            print(f"  {short}  ->  {full}")


# ---------------------------------------------------------------------------
# Memory helpers
# ---------------------------------------------------------------------------

def _format_bytes(n_bytes: int) -> str:
    if n_bytes < 1024:
        return f"{n_bytes} B"
    if n_bytes < 1024 ** 2:
        return f"{n_bytes / 1024:.2f} KB"
    return f"{n_bytes / 1024 ** 2:.2f} MB"


def _sparse_size_bytes(sparse: dict[tuple[str, str], float]) -> int:
    import sys as _sys
    total = _sys.getsizeof(sparse)
    for (a, b), v in sparse.items():
        total += _sys.getsizeof(a) + _sys.getsizeof(b) + _sys.getsizeof(v) + _sys.getsizeof((a, b))
    return total


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Print XL-MS adjacency matrix")
    parser.add_argument("csv", help="Path to the crosslink CSV file")
    parser.add_argument("--n", type=int, default=DEFAULT_N,
                        help=f"Number of unique proteins (default {DEFAULT_N})")
    parser.add_argument("--order",
                        choices=["confidence", "alpha", "cluster", "sequence",
                                 "pathway", "complex", "size"],
                        default="confidence",
                        help="Axis ordering: confidence (default), alpha, cluster, "
                             "sequence, pathway (KEGG), complex (STRING), "
                             "size (longest sequence first)")
    parser.add_argument("--fasta", default=None,
                        help="FASTA file path (required for --order sequence)")
    parser.add_argument("--species", type=int, default=9606,
                        help="NCBI taxonomy ID for external DB queries (default 9606 = human)")
    args = parser.parse_args()

    if args.order == "sequence" and not args.fasta:
        parser.error("--fasta is required when --order sequence")

    tracemalloc.start()

    proteins, sparse, protein_is_decoy = build_matrix(args.csv, n=args.n)

    sequences = read_fasta(args.fasta) if args.fasta else None
    proteins, protein_section = _sort_proteins(
        proteins, sparse, protein_is_decoy, args.order, sequences, args.species
    )

    _, peak_traced = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    print_matrix(proteins, sparse, protein_is_decoy, protein_section)

    sparse_bytes = _sparse_size_bytes(sparse)
    n_targets = sum(1 for p in proteins if not protein_is_decoy.get(p, False))
    n_decoys = len(proteins) - n_targets
    print(
        f"\nProteins: {n_targets} targets, {n_decoys} decoys"
        f"  |  sparse store: {_format_bytes(sparse_bytes)} ({len(sparse)} links)"
        f"  |  peak traced: {_format_bytes(peak_traced)}"
    )


if __name__ == "__main__":
    main()

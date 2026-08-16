"""Print a sparse adjacency matrix for the top-N most-confident unique proteins
found in a crosslink CSV file. See docs/usage.md for the full CLI reference."""
from __future__ import annotations

import argparse
import io
import sys
import os
import tracemalloc
import warnings
from typing import NamedTuple, Optional

# Ensure Unicode characters render correctly on Windows consoles
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", line_buffering=True)

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from xlms.io import read_fasta

DEFAULT_N = 1000

COL_W = 10   # column header width, in chars
ROW_W = 12   # row label width, in chars
CELL_W = 2   # one dot/dash plus a space

N_LEVELS = 4
DOTS = ["·", "○", "●", "⬤"]   # small to large, keep in sync with N_LEVELS
EMPTY = " "

_DECOY_NAME_PREFIXES = ("decoy", "rev_", "contam_")


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


def _truncate(name: str, width: int) -> str:
    return name[:width] if len(name) > width else name


def _is_ambiguous(name: str) -> bool:
    return ";" in str(name)


def _parse_position(val) -> int | None:
    """Parse a residue-position value. Returns None for missing values and
    for ambiguous ones (some search tools report multiple candidate
    positions as e.g. "114;53" — same convention as ambiguous protein
    names — which can't be pinned to a single residue)."""
    if pd.isna(val):
        return None
    if _is_ambiguous(val):
        return None
    try:
        return int(float(val))
    except (TypeError, ValueError):
        return None


def _combine_pep_link(pep_val, link_val) -> int | None:
    """Resolve an absolute residue position from a peptide start position
    plus the crosslinked residue's 1-based offset within that peptide:
    pep_start + link_offset - 1. Returns None if either value is missing
    or ambiguous."""
    pep = _parse_position(pep_val)
    link = _parse_position(link_val)
    if pep is None or link is None:
        return None
    return pep + link - 1


def _short_accession(name: str) -> str:
    """sp|P12345|PROT_HUMAN desc -> P12345; anything without '|' is returned
    as its first whitespace-delimited token."""
    if "|" in name:
        parts = name.split("|")
        if len(parts) >= 2:
            return parts[1].strip()
    return name.split()[0].strip()


def _to_bool(val) -> bool:
    if isinstance(val, bool):
        return val
    return str(val).strip().lower() in ("true", "1", "yes")


def _is_decoy_by_name(name: str) -> bool:
    return name.lower().startswith(_DECOY_NAME_PREFIXES)


def _detect_protein_cols(columns: list[str]) -> tuple[str, str]:
    col_map = {c.lower(): c for c in columns}
    return col_map["protein1"], col_map["protein2"]


def _detect_score_col(columns: list[str]) -> str:
    col_map = {c.lower(): c for c in columns}
    return col_map["score"]


def _detect_decoy_cols(columns: list[str]) -> dict[str, str]:
    """Figure out which decoy column format the CSV uses and return a dict
    describing it: format "B" (Decoy1/Decoy2), format "A" (isTT/isDD), or
    "none" if neither is present."""
    col_map = {c.lower(): c for c in columns}
    if "decoy1" in col_map and "decoy2" in col_map:
        return {"format": "B", "d1": col_map["decoy1"], "d2": col_map["decoy2"]}
    if "istt" in col_map and "isdd" in col_map:
        return {"format": "A", "istt": col_map["istt"], "isdd": col_map["isdd"]}
    return {"format": "none"}


def _detect_position_cols(columns: list[str]) -> dict[str, str] | None:
    """Figure out how residue positions are expressed in this CSV and return
    a dict describing it:
      format "direct":   {"format":"direct", "pos1": <col>, "pos2": <col>}
        — the column already holds the residue's absolute position in the
        protein (SeqPos1/SeqPos2 or a common synonym).
      format "pep_link": {"format":"pep_link", "pep1":<col>, "link1":<col>,
                           "pep2":<col>, "link2":<col>}
        — PepPos1/PepPos2 (where the linked peptide starts in the protein)
        plus LinkPos1/LinkPos2 (the crosslinked residue's 1-based offset
        within that peptide); the absolute position is pep + link - 1.
    Returns None if neither is present."""
    col_map = {c.lower(): c for c in columns}
    for a, b in (("seqpos1", "seqpos2"), ("pos1", "pos2"), ("position1", "position2")):
        if a in col_map and b in col_map:
            return {"format": "direct", "pos1": col_map[a], "pos2": col_map[b]}
    if all(k in col_map for k in ("peppos1", "peppos2", "linkpos1", "linkpos2")):
        return {
            "format": "pep_link",
            "pep1": col_map["peppos1"], "link1": col_map["linkpos1"],
            "pep2": col_map["peppos2"], "link2": col_map["linkpos2"],
        }
    return None


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
        # a protein that only shows up in TD links never gets classified above,
        # so fall back to guessing from its name
        for p in all_proteins:
            if p not in result:
                result[p] = _is_decoy_by_name(p)

    # unclassified proteins (or fmt == "none") default to target
    for p in all_proteins:
        result.setdefault(p, False)

    return result


def _build_thresholds(sparse: dict[tuple[str, str], float]) -> list[float]:
    max_score = max(sparse.values())
    return [max_score * (i + 1) / N_LEVELS for i in range(N_LEVELS - 1)]


def _score_to_dot(score: float, thresholds: list[float]) -> str:
    for i, t in enumerate(thresholds):
        if score <= t:
            return DOTS[i]
    return DOTS[-1]


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

        missing = [p for p in group if p not in (sequences or {})]
        present = [p for p in group if p in (sequences or {})]
        if missing:
            warnings.warn(
                f"{len(missing)} protein(s) not found in FASTA, appended at end: {missing}"
            )
        if len(present) < 2:
            return present + missing, None

        # cluster by 3-mer composition (crude but fast proxy for sequence similarity)
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

    return group, None


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


def _sort_residues(
    residues: list[ResidueId],
    sparse: dict[tuple[ResidueId, ResidueId], float],
    residue_is_decoy: dict[ResidueId, bool],
    mode: str,
    sequences: dict[str, str] | None = None,
    species: int = 9606,
) -> tuple[list[ResidueId], dict[ResidueId, str]]:
    """Order residues by `mode` without reimplementing any ordering mode:
    order the distinct parent proteins with the existing protein-level
    _sort_proteins (via a score-aggregated protein sparse dict), then expand
    each protein into its selected residues, sorted by position. The parent
    protein is always returned as a per-residue section label, so column/row
    headers can group same-protein residues regardless of `mode`."""
    # first-occurrence order, not sorted — keeps e.g. mode="confidence" (a
    # no-op passthrough) meaningful instead of silently alphabetizing
    proteins = list(dict.fromkeys(r.protein for r in residues))
    protein_is_decoy = {r.protein: residue_is_decoy[r] for r in residues}

    protein_sparse: dict[tuple[str, str], float] = {}
    for (r1, r2), score in sparse.items():
        if r1.protein == r2.protein:
            continue
        key = (min(r1.protein, r2.protein), max(r1.protein, r2.protein))
        if score > protein_sparse.get(key, float("-inf")):
            protein_sparse[key] = score

    ordered_proteins, _ = _sort_proteins(
        proteins, protein_sparse, protein_is_decoy, mode, sequences, species
    )

    by_protein: dict[str, list[ResidueId]] = {}
    for r in residues:
        by_protein.setdefault(r.protein, []).append(r)

    ordered_residues: list[ResidueId] = []
    for protein in ordered_proteins:
        ordered_residues.extend(sorted(by_protein[protein], key=lambda r: r.pos))

    residue_section = {r: r.protein for r in ordered_residues}
    return ordered_residues, residue_section


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

    decoy_cols = [v for k, v in decoy_info.items() if k not in ("format",)]
    keep_cols = list(dict.fromkeys([p1_col, p2_col, score_col] + decoy_cols))
    df = df[keep_cols].copy()
    df.sort_values(score_col, ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)

    # walk rows highest-score-first and keep collecting until we have N unique proteins
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

    protein_is_decoy = _classify_proteins(
        df, p1_col, p2_col, decoy_info, set(proteins_unordered)
    )

    # targets first, then decoys, confidence order preserved within each group
    conf_rank = {p: i for i, p in enumerate(proteins_unordered)}
    proteins = sorted(proteins_unordered, key=lambda p: (int(protein_is_decoy[p]), conf_rank[p]))

    # second pass over the rows to fill in the sparse score matrix
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


def build_residue_matrix(
    path: str,
    n: int = DEFAULT_N,
) -> tuple[list[ResidueId], dict[tuple[ResidueId, ResidueId], float], dict[ResidueId, bool]]:
    """Residue-level counterpart of build_matrix — same two-pass algorithm,
    but the axis unit is a (protein, position) ResidueId instead of a bare
    protein name.

    Returns:
        residues          — N unique residues, targets first then decoys,
                            confidence order preserved within each group
        sparse             — {(min_res, max_res): best_score}
        residue_is_decoy   — {residue: bool}, inherited from the parent protein
    """
    df = pd.read_csv(path)
    columns = list(df.columns)
    p1_col, p2_col = _detect_protein_cols(columns)
    score_col = _detect_score_col(columns)
    decoy_info = _detect_decoy_cols(columns)
    pos_info = _detect_position_cols(columns)
    if pos_info is None:
        raise ValueError(
            "Residue-level mode requires residue-position columns "
            "(e.g. SeqPos1/SeqPos2, or PepPos1/PepPos2 + LinkPos1/LinkPos2). "
            "Found columns: " + ", ".join(columns)
        )

    if pos_info["format"] == "direct":
        pos1_col, pos2_col = pos_info["pos1"], pos_info["pos2"]
        pos_raw_cols = [pos1_col, pos2_col]
        def _get_pos1(row):
            return _parse_position(row[pos1_col])
        def _get_pos2(row):
            return _parse_position(row[pos2_col])
    else:
        pep1_col, link1_col = pos_info["pep1"], pos_info["link1"]
        pep2_col, link2_col = pos_info["pep2"], pos_info["link2"]
        pos_raw_cols = [pep1_col, link1_col, pep2_col, link2_col]
        def _get_pos1(row):
            return _combine_pep_link(row[pep1_col], row[link1_col])
        def _get_pos2(row):
            return _combine_pep_link(row[pep2_col], row[link2_col])

    decoy_cols = [v for k, v in decoy_info.items() if k not in ("format",)]
    keep_cols = list(dict.fromkeys([p1_col, p2_col, score_col] + pos_raw_cols + decoy_cols))
    df = df[keep_cols].copy()
    df.sort_values(score_col, ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)

    # walk rows highest-score-first and keep collecting until we have N unique residues
    selected: dict[ResidueId, None] = {}
    n_skipped = 0
    for _, row in df.iterrows():
        p1, p2 = str(row[p1_col]), str(row[p2_col])
        if _is_ambiguous(p1) or _is_ambiguous(p2):
            continue
        pos1 = _get_pos1(row)
        pos2 = _get_pos2(row)
        if pos1 is None or pos2 is None:
            n_skipped += 1
            continue
        selected.setdefault(ResidueId(p1, pos1))
        selected.setdefault(ResidueId(p2, pos2))
        if len(selected) >= n:
            break
    residues_unordered = list(selected.keys())
    if n_skipped:
        warnings.warn(f"{n_skipped} row(s) skipped: missing/ambiguous residue position")
    if not residues_unordered:
        raise ValueError("No residues with valid position data were found in the CSV.")

    protein_is_decoy = _classify_proteins(
        df, p1_col, p2_col, decoy_info, {r.protein for r in residues_unordered}
    )
    residue_is_decoy = {r: protein_is_decoy[r.protein] for r in residues_unordered}

    # targets first, then decoys, confidence order preserved within each group
    conf_rank = {r: i for i, r in enumerate(residues_unordered)}
    residues = sorted(
        residues_unordered, key=lambda r: (int(residue_is_decoy[r]), conf_rank[r])
    )

    # second pass over the rows to fill in the sparse score matrix
    residue_set = set(residues)
    sparse: dict[tuple[ResidueId, ResidueId], float] = {}
    for _, row in df.iterrows():
        pos1 = _get_pos1(row)
        pos2 = _get_pos2(row)
        if pos1 is None or pos2 is None:
            continue
        r1 = ResidueId(str(row[p1_col]), pos1)
        r2 = ResidueId(str(row[p2_col]), pos2)
        if r1 in residue_set and r2 in residue_set:
            key = (min(r1, r2), max(r1, r2))
            score = float(row[score_col])
            if score > sparse.get(key, float("-inf")):
                sparse[key] = score

    return residues, sparse, residue_is_decoy


def _add_unlinked_residues(
    residues: list[ResidueId],
    residue_is_decoy: dict[ResidueId, bool],
    sequences: dict[str, str],
) -> tuple[list[ResidueId], dict[ResidueId, bool]]:
    """Expand a linked-residue list to also include every other residue
    position of each protein already present, using full sequences from a
    FASTA file. `sparse` is left untouched — the added residues simply have
    no score entries, so they render as empty cells. Proteins missing from
    `sequences` are left as-is (can't be expanded without a sequence)."""
    linked_proteins = list(dict.fromkeys(r.protein for r in residues))
    protein_is_decoy = {r.protein: residue_is_decoy[r] for r in residues}
    existing = set(residues)

    missing = [p for p in linked_proteins if p not in sequences]
    if missing:
        warnings.warn(
            f"{len(missing)} protein(s) not found in FASTA — only linked "
            f"residues shown for them: {missing}"
        )

    expanded = list(residues)
    is_decoy = dict(residue_is_decoy)
    for protein in linked_proteins:
        seq = sequences.get(protein)
        if not seq:
            continue
        for pos in range(1, len(seq) + 1):
            rid = ResidueId(protein, pos)
            if rid not in existing:
                existing.add(rid)
                expanded.append(rid)
                is_decoy[rid] = protein_is_decoy[protein]

    return expanded, is_decoy


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
    proteins: list,
    sparse: dict[tuple, float],
    protein_is_decoy: dict,
    protein_section: dict | None = None,
    *,
    col_label_fn=None,
    row_label_fn=None,
    legend_source: str = "col",
) -> None:
    """Render the dot matrix to stdout. `proteins` is the ordered list of
    axis units (protein names in protein mode, ResidueIds in residue mode).
    col_label_fn/row_label_fn default to the plain truncated-name labels used
    in protein mode; pass overrides for residue mode. legend_source picks
    which label set ("col" or "row") the truncated-name legend is built from."""
    col_label_fn = col_label_fn or (lambda p: _truncate(p, COL_W))
    row_label_fn = row_label_fn or (lambda p: _truncate(p, ROW_W))

    n = len(proteins)
    labels = [col_label_fn(p) for p in proteins]
    row_labels = [row_label_fn(p) for p in proteins]
    thresholds = _build_thresholds(sparse) if sparse else []

    split = next((i for i, p in enumerate(proteins) if protein_is_decoy.get(p, False)), n)
    has_sep = 0 < split < n

    row_prefix = " " * (ROW_W + 2)

    if protein_section:
        print(_section_header_line(proteins, protein_section, split, row_prefix))

    # column headers are printed one character row at a time, stacked vertically
    max_label_len = max(len(lb) for lb in labels)
    for char_idx in range(max_label_len):
        line = row_prefix
        for i, lb in enumerate(labels):
            if has_sep and i == split:
                line += "|"
            char = lb[char_idx] if char_idx < len(lb) else " "
            line += char.center(CELL_W)
        print(line)

    print(row_prefix + _col_separator_line(n, split))

    for row_idx, pi in enumerate(proteins):
        if has_sep and row_idx == split:
            print("-" * (ROW_W + 2) + _col_separator_line(n, split))

        line = row_labels[row_idx].ljust(ROW_W) + "  "
        for j, pj in enumerate(proteins):
            if has_sep and j == split:
                line += "|"
            key = (min(pi, pj), max(pi, pj))
            score: Optional[float] = sparse.get(key)
            dot_char = _score_to_dot(score, thresholds) if score is not None else EMPTY
            line += dot_char.ljust(CELL_W)
        print(line)

    if sparse:
        max_score = max(sparse.values())
        print("\nScale (max score {:.2f}):".format(max_score))
        lower = 0.0
        for i, t in enumerate(thresholds):
            print(f"  {DOTS[i]}   {lower:.2f} - {t:.2f}")
            lower = t
        print(f"  {DOTS[-1]}   {lower:.2f}+")

    legend_labels = labels if legend_source == "col" else row_labels
    truncated = [(short, str(full)) for short, full in zip(legend_labels, proteins) if short != str(full)]
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


def _sparse_size_bytes(sparse: dict[tuple[str, str], float]) -> int:
    import sys as _sys
    total = _sys.getsizeof(sparse)
    for (a, b), v in sparse.items():
        total += _sys.getsizeof(a) + _sys.getsizeof(b) + _sys.getsizeof(v) + _sys.getsizeof((a, b))
    return total


def main() -> None:
    parser = argparse.ArgumentParser(description="Print XL-MS adjacency matrix")
    parser.add_argument("csv", help="Path to the crosslink CSV file")
    parser.add_argument("--level", choices=["protein", "residue"], default="protein",
                        help="Axis granularity: protein (default) or residue "
                             "(requires SeqPos1/SeqPos2 columns in the CSV)")
    parser.add_argument("--n", type=int, default=DEFAULT_N,
                        help=f"Number of unique proteins/residues (default {DEFAULT_N})")
    parser.add_argument("--order",
                        choices=["confidence", "alpha", "cluster", "sequence",
                                 "pathway", "complex", "size"],
                        default="confidence",
                        help="Axis ordering: confidence (default), alpha, cluster, "
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

    tracemalloc.start()

    try:
        if args.level == "residue":
            items, sparse, is_decoy = build_residue_matrix(args.csv, n=args.n)
        else:
            items, sparse, is_decoy = build_matrix(args.csv, n=args.n)
    except ValueError as exc:
        parser.error(str(exc))

    sequences = read_fasta(args.fasta) if args.fasta else None
    if args.level == "residue":
        if args.include_unlinked:
            items, is_decoy = _add_unlinked_residues(items, is_decoy, sequences or {})
        items, section = _sort_residues(
            items, sparse, is_decoy, args.order, sequences, args.species
        )
        col_label_fn = lambda r: _truncate(str(r.pos), COL_W)
        row_label_fn = lambda r: _truncate(f"{_short_accession(r.protein)}:{r.pos}", ROW_W)
        legend_source = "row"
    else:
        items, section = _sort_proteins(
            items, sparse, is_decoy, args.order, sequences, args.species
        )
        col_label_fn = row_label_fn = None
        legend_source = "col"

    _, peak_traced = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    print_matrix(items, sparse, is_decoy, section,
                 col_label_fn=col_label_fn, row_label_fn=row_label_fn,
                 legend_source=legend_source)

    sparse_bytes = _sparse_size_bytes(sparse)
    n_targets = sum(1 for p in items if not is_decoy.get(p, False))
    n_decoys = len(items) - n_targets
    unit = "Residues" if args.level == "residue" else "Proteins"
    print(
        f"\n{unit}: {n_targets} targets, {n_decoys} decoys"
        f"  |  sparse store: {_format_bytes(sparse_bytes)} ({len(sparse)} links)"
        f"  |  peak traced: {_format_bytes(peak_traced)}"
    )


if __name__ == "__main__":
    main()

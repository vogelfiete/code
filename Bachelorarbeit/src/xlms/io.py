"""Reading crosslink CSV files and FASTA sequence files into the types from
xlms.models."""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Union

import pandas as pd
from Bio import SeqIO

from xlms.models import CrossLink, CrossLinkDataset

_BOOL_COLS = ["isDecoy", "isTT", "isTD", "isDD"]

_CSV_DTYPE: dict = {
    "Id": str,
    "Protein1": str,
    "Protein2": str,
    "SeqPos1": "Int64",
    "SeqPos2": "Int64",
    "Score": float,
}


# -----------------------------------------------------------------------------
# Reading crosslinks
#
# The CSV can say whether a link is target or decoy in one of two ways:
# four flag columns (isDecoy/isTT/isTD/isDD), or one decoy flag per protein
# plus a DecoyType column. Both end up as the same CrossLink objects.
# -----------------------------------------------------------------------------

def _normalise_bool(series: pd.Series) -> pd.Series:
    lowered = series.astype(str).str.strip().str.lower()
    return lowered.map({"true": True, "1": True, "false": False, "0": False}).astype(bool)


def read_csv(path: Union[str, Path]) -> List[CrossLink]:
    df = pd.read_csv(
        path,
        dtype=_CSV_DTYPE,
        true_values=["True", "true", "1"],
        false_values=["False", "false", "0"],
    )
    columns = set(df.columns)
    if {"isDecoy", "isTT", "isTD", "isDD"}.issubset(columns):
        for col in _BOOL_COLS:
            df[col] = _normalise_bool(df[col])
    elif {"Decoy1", "Decoy2", "DecoyType"}.issubset(columns):
        df["isDecoy"] = _normalise_bool(df["Decoy1"]) | _normalise_bool(df["Decoy2"])
        decoy_type = df["DecoyType"].astype(str).str.strip().str.upper()
        df["isTT"] = decoy_type == "TT"
        df["isTD"] = decoy_type == "TD"
        df["isDD"] = decoy_type == "DD"
    else:
        raise ValueError(
            "CSV must contain either {isDecoy, isTT, isTD, isDD} "
            "or {Decoy1, Decoy2, DecoyType}. "
            f"Found columns: {sorted(columns)}"
        )

    return [
        CrossLink(
            id=row["Id"],
            protein1=row["Protein1"],
            protein2=row["Protein2"],
            seq_pos1=int(row["SeqPos1"]),
            seq_pos2=int(row["SeqPos2"]),
            score=float(row["Score"]),
            is_decoy=bool(row["isDecoy"]),
            is_tt=bool(row["isTT"]),
            is_td=bool(row["isTD"]),
            is_dd=bool(row["isDD"]),
        )
        for _, row in df.iterrows()
    ]


# -----------------------------------------------------------------------------
# Reading sequences
#
# The FASTA header (up to the first space) becomes the protein name, the
# sequence is stored in upper case.
# -----------------------------------------------------------------------------

def read_fasta(path: Union[str, Path]) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8-sig") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            sequences[record.id] = str(record.seq).upper()
    return sequences


# -----------------------------------------------------------------------------
# Loading both at once
#
# Also warns about proteins in the CSV that the FASTA file doesn't know,
# which usually means the two files name proteins differently.
# -----------------------------------------------------------------------------

def load_dataset(
    csv_path: Union[str, Path],
    fasta_path: Union[str, Path],
) -> CrossLinkDataset:
    crosslinks = read_csv(csv_path)
    sequences = read_fasta(fasta_path)

    fasta_ids = set(sequences.keys())
    csv_proteins = {xl.protein1 for xl in crosslinks} | {xl.protein2 for xl in crosslinks}
    missing = csv_proteins - fasta_ids
    if missing:
        warnings.warn(
            f"Protein names in CSV not found in FASTA (check header format): {sorted(missing)}"
        )

    return CrossLinkDataset(crosslinks=crosslinks, sequences=sequences)

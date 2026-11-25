"""Where the crosslinks come from: a CSV file (plus an optional FASTA), or
the crosslinking API. The command line tool and the GUI both use this, so
they take the same data options.

    --csv links.csv [--fasta proteins.fasta]
    --uuid <result set id> [<more ids> ...] [--api URL]

plus, for either source, --min-score, --sample and --seed to load fewer links."""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from xlms.api import DEFAULT_API_URL, fetch_matches, fetch_sequences
from xlms.io import read_fasta


# -----------------------------------------------------------------------------
# Command line options
#
# At most one of --csv or --uuid; with required=True (the default) exactly
# one. The rest is optional.
# -----------------------------------------------------------------------------

def add_data_args(parser: argparse.ArgumentParser, required: bool = True) -> None:
    source = parser.add_mutually_exclusive_group(required=required)
    source.add_argument("--csv", help="Path to a crosslink CSV file")
    source.add_argument("--uuid", nargs="+",
                        help="Load the crosslinks of these result sets from the API")
    parser.add_argument("--fasta", default=None,
                        help="FASTA file with the protein sequences (with --uuid the "
                             "sequences come from the API, unless a FASTA is given)")
    parser.add_argument("--api", default=DEFAULT_API_URL,
                        help=f"Address of the API's matrix endpoints (default {DEFAULT_API_URL})")
    parser.add_argument("--min-score", type=float, default=None,
                        help="Only load crosslinks with at least this score")
    parser.add_argument("--sample", type=int, default=None,
                        help="Only load this many crosslinks, picked at random "
                             "(after --min-score)")
    parser.add_argument("--seed", type=int, default=None,
                        help="With --sample: pick the same crosslinks every time")


def describe(args: argparse.Namespace) -> str:
    """A short name for the data, for window titles and messages."""
    return os.path.basename(args.csv) if args.csv else "result set " + ", ".join(args.uuid)


# -----------------------------------------------------------------------------
# Loading
#
# Gives back the crosslink table (what build_matrix reads) and the protein
# sequences, or None for the sequences if there aren't any.
# -----------------------------------------------------------------------------

def load_data(args: argparse.Namespace) -> tuple[pd.DataFrame, dict[str, str] | None]:
    if args.csv:
        table = _filter(pd.read_csv(args.csv), args.min_score, args.sample, args.seed)
        sequences = read_fasta(args.fasta) if args.fasta else None
    else:
        table = fetch_matches(args.api, args.uuid, args.min_score, args.sample, args.seed)
        sequences = read_fasta(args.fasta) if args.fasta else fetch_sequences(args.api, args.uuid)
    return table, sequences


def _filter(df: pd.DataFrame, min_score, sample, seed) -> pd.DataFrame:
    """The same score cutoff and random sample the API does, for CSV files."""
    if min_score is not None:
        score = next((c for c in df.columns if c.lower() == "score"), None)
        if score is None:
            raise ValueError("--min-score needs a Score column in the CSV")
        df = df[df[score] >= min_score]
    if sample is not None and sample < len(df):
        df = df.sample(n=sample, random_state=seed)
    return df

"""Reading crosslinks and protein sequences from the crosslinking API
(the xi2 matrix endpoints) instead of from CSV and FASTA files."""
from __future__ import annotations

import json
import urllib.parse
import urllib.request

import pandas as pd

DEFAULT_API_URL = "http://localhost:8080/xi2/matrix"
_USER_AGENT = {"User-Agent": "xlms-matrix/1.0"}


# -----------------------------------------------------------------------------
# Talking to the API
#
# Every request names one or more result sets by their UUID. The answer is
# always JSON.
# -----------------------------------------------------------------------------

def _get_json(api_url: str, endpoint: str, params: list[tuple[str, object]]):
    query = urllib.parse.urlencode([(k, v) for k, v in params if v is not None])
    url = f"{api_url.rstrip('/')}/{endpoint}?{query}"
    req = urllib.request.Request(url, headers=_USER_AGENT)
    with urllib.request.urlopen(req, timeout=600) as resp:
        return json.load(resp)


# -----------------------------------------------------------------------------
# Crosslinks
#
# The API hands back the same columns a crosslink CSV has (Protein1,
# Protein2, Score, Decoy1, Decoy2, SeqPos1, SeqPos2), so the table can go
# straight into build_matrix like a CSV would.
# -----------------------------------------------------------------------------

def fetch_matches(
    api_url: str,
    uuids: list[str],
    min_score: float | None = None,
    sample: int | None = None,
    seed: int | None = None,
) -> pd.DataFrame:
    """The crosslinks of the given result sets as a table. With min_score
    only links scoring at least that much are sent; with sample only that
    many links, picked at random (the same ones again for the same seed)."""
    params = [("uuid", u) for u in uuids]
    params += [("min_score", min_score), ("sample", sample), ("seed", seed)]
    rows = _get_json(api_url, "matches", params)
    return pd.DataFrame(rows, columns=["Protein1", "Protein2", "Score", "Decoy1", "Decoy2",
                                       "SeqPos1", "SeqPos2"])


# -----------------------------------------------------------------------------
# Sequences
#
# The same {protein: sequence} dict that read_fasta gives, for the proteins
# the crosslinks of the result sets touch.
# -----------------------------------------------------------------------------

def fetch_sequences(api_url: str, uuids: list[str]) -> dict[str, str]:
    return _get_json(api_url, "sequences", [("uuid", u) for u in uuids])

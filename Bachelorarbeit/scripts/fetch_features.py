"""Sequence annotations from UniProt (helices, binding sites, domains, ...)
for colouring the residue axes of the GUI.

For every protein UniProt lists "features": a type, the first and last
residue, and a short description. They are asked for in batches and cached
as JSON under ~/.cache/xlms_matrix/, so each protein is only looked up once."""
from __future__ import annotations

import json
import warnings

from fetch_group_order import _cache_path, _get, _load_cache, _save_cache, short_accession
from fetch_group_order import _UNIPROT_ACCESSION, _UNIPROT_BATCH

# The feature types that can be shown: UniProt's name -> its field name in
# the UniProt REST API
FEATURE_TYPES = {
    "Helix": "ft_helix",
    "Beta strand": "ft_strand",
    "Turn": "ft_turn",
    "Binding site": "ft_binding",
    "Active site": "ft_act_site",
    "Domain": "ft_domain",
    "Transmembrane": "ft_transmem",
    "Modified residue": "ft_mod_res",
}


def _describe(feature: dict) -> str:
    """The short text of a feature: its description, or for a binding site
    the ligand."""
    text = feature.get("description", "")
    ligand = feature.get("ligand", {}).get("name", "")
    return text or ligand


def _fetch(accessions: list[str]) -> dict[str, list[list]]:
    """The features of each accession as {accession: [[type, start, end,
    text], ...]}. Accessions UniProt doesn't know get an empty list."""
    result: dict[str, list[list]] = {acc: [] for acc in accessions}
    fields = ",".join(["accession", *FEATURE_TYPES.values()])
    for k in range(0, len(accessions), _UNIPROT_BATCH):
        batch = accessions[k:k + _UNIPROT_BATCH]
        url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions="
               + ",".join(batch) + f"&fields={fields}&format=json")
        for entry in json.loads(_get(url)).get("results", []):
            features = []
            for f in entry.get("features", []):
                loc = f.get("location", {})
                start = loc.get("start", {}).get("value")
                end = loc.get("end", {}).get("value")
                if f.get("type") in FEATURE_TYPES and start is not None and end is not None:
                    features.append([f["type"], int(start), int(end), _describe(f)])
            result[entry["primaryAccession"]] = features
    return result


def features_for(proteins: list[str]) -> dict[str, list[list]]:
    """{protein: [[type, start, end, text], ...]} for the proteins, as far
    as UniProt knows them (protein names like "sp|P12345|..." are fine).
    Positions count from 1, start and end both included."""
    cache_file = _cache_path("uniprot", "features.json")
    cache: dict[str, list[list]] = _load_cache(cache_file) or {}
    # "P12345-2" is an isoform of P12345
    bare = {p: short_accession(p).split("-")[0] for p in proteins}
    wanted = {acc for acc in bare.values() if _UNIPROT_ACCESSION.fullmatch(acc)}
    to_fetch = sorted(wanted - cache.keys())
    if to_fetch:
        print(f"[UniProt] Asking for the features of {len(to_fetch)} protein(s)…")
        try:
            cache.update(_fetch(to_fetch))
            _save_cache(cache_file, cache)
        except Exception as e:
            warnings.warn(f"UniProt: feature request failed — {e}")
    return {p: cache.get(bare[p], []) for p in proteins}

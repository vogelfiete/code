"""Fetch biological groupings from KEGG (pathways) or STRING DB (complexes) and
return a reordered protein list so proteins sharing a pathway or complex end up
next to each other, along with a per-protein section label. Used by print_matrix.py
via order_by_kegg() and order_by_string()."""
from __future__ import annotations

import json
import os
import time
import urllib.parse
import urllib.request
import warnings
from pathlib import Path
from typing import Optional

_CACHE_ROOT = Path.home() / ".cache" / "xlms_matrix"


def _cache_path(species: str, name: str) -> Path:
    p = _CACHE_ROOT / str(species)
    p.mkdir(parents=True, exist_ok=True)
    return p / name


def _load_cache(path: Path) -> Optional[dict]:
    if path.exists():
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return None


def _save_cache(path: Path, data: dict) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f)


def _get(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "xlms-matrix/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return resp.read().decode("utf-8")


def _post(url: str, data: dict) -> list:
    encoded = urllib.parse.urlencode(data).encode("utf-8")
    req = urllib.request.Request(url, data=encoded, headers={"User-Agent": "xlms-matrix/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _parse_identifier(name: str) -> str:
    """
    Strip XL-MS protein names to a bare identifier for database queries.
    sp|P12345|PROT_HUMAN desc  ->  P12345
    P12345                     ->  P12345
    PROT_HUMAN                 ->  PROT_HUMAN
    """
    # UniProt FASTA headers look like 'sp|ACC|ENTRY desc' or 'tr|ACC|ENTRY desc'
    if "|" in name:
        parts = name.split("|")
        if len(parts) >= 2:
            return parts[1].strip()
    return name.split()[0].strip()


def _cluster_by_memberships(
    group: list[str],
    memberships: dict[str, set[str]],
) -> list[str]:
    """
    Given a dict of {protein: {term_id, ...}}, cluster proteins by shared
    term co-membership and return them in dendrogram leaf order.
    Proteins with no memberships are appended at the end.
    """
    import numpy as np
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform

    present = [p for p in group if memberships.get(p)]
    absent  = [p for p in group if not memberships.get(p)]

    if len(present) < 2:
        return present + absent

    n = len(present)
    co = np.zeros((n, n))
    for i, pi in enumerate(present):
        for j, pj in enumerate(present):
            co[i, j] = len(memberships[pi] & memberships[pj])

    max_co = co.max() or 1.0
    D = max_co - co
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D), method="average")
    order = leaves_list(Z)
    return [present[i] for i in order] + absent


_NCBI_TO_KEGG: dict[int, str] = {
    9606:  "hsa",   # Homo sapiens
    10090: "mmu",   # Mus musculus
    10116: "rno",   # Rattus norvegicus
    9913:  "bta",   # Bos taurus
    9823:  "ssc",   # Sus scrofa
    7227:  "dme",   # Drosophila melanogaster
    6239:  "cel",   # C. elegans
    4932:  "sce",   # S. cerevisiae
}


def _ncbi_to_kegg_code(taxid: int) -> str:
    if taxid in _NCBI_TO_KEGG:
        return _NCBI_TO_KEGG[taxid]
    # Query KEGG for unknown species
    text = _get("https://rest.kegg.jp/list/organism")
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 4 and str(taxid) in parts[3]:
            return parts[1]
    raise ValueError(f"KEGG organism code not found for NCBI taxid {taxid}")


def _kegg_load_uniprot_map(org: str) -> dict[str, str]:
    """
    Download (and cache) the full UniProt accession → KEGG gene ID map for org.
    Uses https://rest.kegg.jp/conv/{org}/uniprot — a single request for all proteins.
    """
    cache_file = _cache_path(f"{org}_kegg", "uniprot_map.json")
    cached = _load_cache(cache_file)
    if cached is not None:
        return cached

    print(f"[KEGG] Downloading UniProt→gene map for {org} (one-time, will be cached)…")
    text = _get(f"https://rest.kegg.jp/conv/{org}/uniprot")
    result: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            # parts[0] = "up:P12345", parts[1] = "bta:280705"
            uniprot_acc = parts[0].strip().removeprefix("up:")
            kegg_gene   = parts[1].strip()
            result[uniprot_acc] = kegg_gene

    _save_cache(cache_file, result)
    return result


def _kegg_load_pathway_map(org: str) -> dict[str, list[str]]:
    """
    Download (and cache) the full KEGG gene ID → pathway IDs map for org.
    Uses https://rest.kegg.jp/link/pathway/{org} — a single request.
    """
    cache_file = _cache_path(f"{org}_kegg", "pathway_map.json")
    cached = _load_cache(cache_file)
    if cached is not None:
        return cached

    print(f"[KEGG] Downloading gene→pathway map for {org} (one-time, will be cached)…")
    text = _get(f"https://rest.kegg.jp/link/pathway/{org}")
    result: dict[str, list[str]] = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            gene_id    = parts[0].strip()
            pathway_id = parts[1].strip()
            result.setdefault(gene_id, []).append(pathway_id)

    _save_cache(cache_file, result)
    return result


def _kegg_load_pathway_names(org: str) -> dict[str, str]:
    """
    Download (and cache) human-readable pathway names for the organism.
    Uses https://rest.kegg.jp/list/pathway/{org} — a single request.
    Returns {pathway_id: display_name}, e.g. {"path:hsa00010": "Glycolysis / Gluconeogenesis"}.
    """
    cache_file = _cache_path(f"{org}_kegg", "pathway_names.json")
    cached = _load_cache(cache_file)
    if cached is not None:
        return cached

    print(f"[KEGG] Downloading pathway names for {org} (one-time, will be cached)…")
    text = _get(f"https://rest.kegg.jp/list/pathway/{org}")
    result: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            pathway_id   = parts[0].strip()
            # "Glycolysis / Gluconeogenesis - Homo sapiens (human)" → first part only
            display_name = parts[1].split(" - ")[0].strip()
            result[pathway_id] = display_name

    _save_cache(cache_file, result)
    return result


def _assign_primary_groups(
    proteins: list[str],
    memberships: dict[str, set[str]],
    name_map: dict[str, str],
    max_groups: int = 10,
) -> dict[str, str]:
    """
    Assign each protein to exactly one named section for display.
    Picks the top max_groups terms by number of matching selected proteins,
    then greedily assigns each protein to its highest-ranked matching term.
    Proteins with no match are assigned to "Unknown".
    Returns {protein_name: section_display_label}.
    """
    from collections import Counter
    term_counts = Counter(term for p in proteins for term in memberships.get(p, set()))
    top_terms = [term for term, _ in term_counts.most_common(max_groups)]

    assignment: dict[str, str] = {}
    for protein in proteins:
        pws = memberships.get(protein, set())
        for term in top_terms:
            if term in pws:
                assignment[protein] = name_map.get(term, term)
                break
        else:
            assignment[protein] = "Unknown"
    return assignment


def _group_and_sort(
    proteins: list[str],
    assignment: dict[str, str],
) -> list[str]:
    """
    Order proteins by their section (largest section first, Unknown last),
    preserving the original order within each section.
    """
    from collections import defaultdict
    groups: dict[str, list[str]] = defaultdict(list)
    for p in proteins:
        groups[assignment[p]].append(p)

    section_order = sorted(
        [s for s in groups if s != "Unknown"],
        key=lambda s: -len(groups[s]),
    )
    if "Unknown" in groups:
        section_order.append("Unknown")

    return [p for s in section_order for p in groups[s]]


def order_by_kegg(proteins: list[str], species: int = 9606) -> tuple[list[str], dict[str, str]]:
    """
    Order proteins so those sharing KEGG pathways appear adjacent.
    Downloads org-wide mapping tables once and caches them to disk.
    Returns (ordered_proteins, {protein: section_label}).
    """
    org = _ncbi_to_kegg_code(species)
    uniprot_map   = _kegg_load_uniprot_map(org)    # {uniprot_acc: kegg_gene_id}
    pathway_map   = _kegg_load_pathway_map(org)    # {kegg_gene_id: [pathway_ids]}
    pathway_names = _kegg_load_pathway_names(org)  # {pathway_id: display_name}

    memberships: dict[str, set[str]] = {}
    for protein in proteins:
        identifier = _parse_identifier(protein)
        gene_id = uniprot_map.get(identifier)
        if gene_id:
            memberships[protein] = set(pathway_map.get(gene_id, []))
        else:
            memberships[protein] = set()
            warnings.warn(
                f"KEGG: '{protein}' (id: '{identifier}') not in {org} UniProt map"
                " — appended at end"
            )

    assignment = _assign_primary_groups(proteins, memberships, pathway_names)
    ordered    = _group_and_sort(proteins, assignment)
    return ordered, assignment


_STRING_BASE = "https://string-db.org/api/json"
_STRING_COMPLEX_CATEGORIES = {"CORUM", "PPI_hub_proteins", "KEGG_Pathways"}


def _string_resolve_ids(identifiers: list[str], taxid: int) -> dict[str, str]:
    """Return {input_name: string_id} mapping."""
    try:
        results = _post(f"{_STRING_BASE}/get_string_ids", {
            "identifiers": "\r".join(identifiers),
            "species": taxid,
            "limit": 1,
            "echo_query": 1,
            "caller_identity": "xlms_matrix",
        })
        return {r["queryItem"]: r["stringId"] for r in results if "stringId" in r}
    except Exception as e:
        warnings.warn(f"STRING: ID resolution failed — {e}")
        return {}


def _string_fetch_enrichment(string_ids: list[str], taxid: int) -> list[dict]:
    """Return enrichment records from STRING."""
    try:
        return _post(f"{_STRING_BASE}/enrichment", {
            "identifiers": "\r".join(string_ids),
            "species": taxid,
            "caller_identity": "xlms_matrix",
        })
    except Exception as e:
        warnings.warn(f"STRING: enrichment fetch failed — {e}")
        return []


def order_by_string(proteins: list[str], species: int = 9606) -> tuple[list[str], dict[str, str]]:
    """
    Order proteins so those sharing STRING complex / pathway annotations appear
    adjacent. Uses two batch HTTP calls total, with local disk caching.
    Returns (ordered_proteins, {protein: section_label}).
    """
    cache_file      = _cache_path(f"{species}_string", "complexes.json")
    names_cache_file = _cache_path(f"{species}_string", "term_names.json")
    cache: dict[str, list[str]] = _load_cache(cache_file) or {}
    term_names: dict[str, str]  = _load_cache(names_cache_file) or {}

    identifiers = [_parse_identifier(p) for p in proteins]
    id_map = dict(zip(proteins, identifiers))  # protein -> bare identifier
    to_resolve = [p for p in proteins if id_map[p] not in cache]

    if to_resolve:
        string_ids_map = _string_resolve_ids([id_map[p] for p in to_resolve], species)
        time.sleep(1)

        if string_ids_map:
            enrichment = _string_fetch_enrichment(list(string_ids_map.values()), species)
            time.sleep(1)

            string_memberships: dict[str, set[str]] = {sid: set() for sid in string_ids_map.values()}
            for record in enrichment:
                if record.get("category") in _STRING_COMPLEX_CATEGORIES:
                    term = record["term"]
                    desc = record.get("description") or term
                    term_names[term] = desc
                    for gene in record.get("inputGenes", "").split(","):
                        gene = gene.strip()
                        for sid in string_ids_map.values():
                            if sid.endswith(gene) or gene in sid:
                                string_memberships[sid].add(term)

            for protein in to_resolve:
                bare = id_map[protein]
                if bare in string_ids_map:
                    sid = string_ids_map[bare]
                    cache[bare] = sorted(string_memberships.get(sid, set()))
                else:
                    cache[bare] = []
                    warnings.warn(f"STRING: no entry found for '{protein}' — appended at end")

        _save_cache(cache_file, cache)
        _save_cache(names_cache_file, term_names)

    memberships: dict[str, set[str]] = {
        p: set(cache.get(id_map[p], [])) for p in proteins
    }
    assignment = _assign_primary_groups(proteins, memberships, term_names)
    ordered    = _group_and_sort(proteins, assignment)
    return ordered, assignment

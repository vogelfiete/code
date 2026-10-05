"""Group proteins by biology for the "pathway" and "complex" orderings.

KEGG tells us which pathways a protein is part of, STRING which complexes.
Proteins that share a pathway or complex are put next to each other, and each
protein gets the name of its group as a section label for the matrix headers.
print_matrix.py uses this through order_by_kegg() and order_by_string().

Everything downloaded is cached as JSON under ~/.cache/xlms_matrix/, so each
lookup only hits the internet once."""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable, Optional

_CACHE_ROOT = Path.home() / ".cache" / "xlms_matrix"
_USER_AGENT = {"User-Agent": "xlms-matrix/1.0"}


# -----------------------------------------------------------------------------
# Turning protein names into database identifiers
#
# Crosslink files often use full UniProt names like "sp|P12345|PROT_HUMAN".
# The databases only want the accession in the middle.
# -----------------------------------------------------------------------------

def short_accession(name: str) -> str:
    """The bare identifier of a protein name, for database lookups.
    sp|P12345|PROT_HUMAN desc  ->  P12345
    P12345                     ->  P12345
    PROT_HUMAN                 ->  PROT_HUMAN
    """
    parts = name.split("|")
    if len(parts) >= 2:  # UniProt style 'sp|ACC|ENTRY desc' or 'tr|ACC|ENTRY desc'
        return parts[1].strip()
    return name.split()[0].strip()


# -----------------------------------------------------------------------------
# Downloading and caching
#
# Small wrappers around urllib plus reading and writing the JSON cache files.
# -----------------------------------------------------------------------------

def _cache_path(folder: str, name: str) -> Path:
    p = _CACHE_ROOT / folder
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
    req = urllib.request.Request(url, headers=_USER_AGENT)
    with urllib.request.urlopen(req, timeout=30) as resp:
        return resp.read().decode("utf-8")


def _post(url: str, data: dict) -> list:
    encoded = urllib.parse.urlencode(data).encode("utf-8")
    req = urllib.request.Request(url, data=encoded, headers=_USER_AGENT)
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


# -----------------------------------------------------------------------------
# Grouping proteins into named sections
#
# Used by both KEGG and STRING: once we know which terms (pathways or
# complexes) each protein belongs to, every protein is put into exactly one
# section, and the sections are sorted from largest to smallest.
# -----------------------------------------------------------------------------

def _assign_primary_groups(
    proteins: list[str],
    memberships: dict[str, set[str]],
    name_map: dict[str, str],
    max_groups: int = 10,
) -> dict[str, str]:
    """Give every protein one section. Only the max_groups terms that most
    proteins share become sections; a protein goes to the first of those it
    belongs to, or to "Unknown"."""
    term_counts = Counter(term for p in proteins for term in memberships.get(p, set()))
    top_terms = [term for term, _ in term_counts.most_common(max_groups)]

    assignment: dict[str, str] = {}
    for protein in proteins:
        terms = memberships.get(protein, set())
        match = next((t for t in top_terms if t in terms), None)
        assignment[protein] = "Unknown" if match is None else name_map.get(match, match)
    return assignment


def _group_and_sort(proteins: list[str], assignment: dict[str, str]) -> list[str]:
    """Biggest section first and "Unknown" last; inside a section the
    proteins keep the order they came in."""
    groups: dict[str, list[str]] = defaultdict(list)
    for p in proteins:
        groups[assignment[p]].append(p)
    sections = sorted((s for s in groups if s != "Unknown"), key=lambda s: -len(groups[s]))
    if "Unknown" in groups:
        sections.append("Unknown")
    return [p for s in sections for p in groups[s]]


def _order_by_memberships(proteins, memberships, term_names):
    assignment = _assign_primary_groups(proteins, memberships, term_names)
    return _group_and_sort(proteins, assignment), assignment


# -----------------------------------------------------------------------------
# KEGG pathways
#
# KEGG has three big tables per organism: UniProt accession to KEGG gene,
# gene to pathways, and pathway names. Each is downloaded once and cached.
# -----------------------------------------------------------------------------

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
    # not one of the usual species: look it up in KEGG's list of organisms
    for line in _get("https://rest.kegg.jp/list/organism").splitlines():
        parts = line.split("\t")
        if len(parts) >= 4 and str(taxid) in parts[3]:
            return parts[1]
    raise ValueError(f"KEGG organism code not found for NCBI taxid {taxid}")


def _cached_kegg_table(
    org: str,
    cache_name: str,
    endpoint: str,
    what: str,
    build: Callable[[list[list[str]]], dict],
) -> dict:
    """Download one of KEGG's tab-separated tables (once per organism), turn
    the rows into a dict with `build` and keep the result in the cache."""
    cache_file = _cache_path(f"{org}_kegg", cache_name)
    cached = _load_cache(cache_file)
    if cached is not None:
        return cached

    print(f"[KEGG] Downloading {what} for {org} (one-time, will be cached)…")
    text = _get(f"https://rest.kegg.jp/{endpoint}")
    rows = [parts for parts in (line.split("\t") for line in text.splitlines()) if len(parts) >= 2]
    result = build(rows)
    _save_cache(cache_file, result)
    return result


def _uniprot_to_gene(rows):
    # "up:P12345" -> "bta:280705"
    return {r[0].strip().removeprefix("up:"): r[1].strip() for r in rows}


def _gene_to_pathways(rows):
    result: dict[str, list[str]] = {}
    for r in rows:
        result.setdefault(r[0].strip(), []).append(r[1].strip())
    return result


def _pathway_names(rows):
    # "Glycolysis / Gluconeogenesis - Homo sapiens (human)" -> just the first part
    return {r[0].strip(): r[1].split(" - ")[0].strip() for r in rows}


def order_by_kegg(proteins: list[str], species: int = 9606) -> tuple[list[str], dict[str, str]]:
    """Put proteins that share KEGG pathways next to each other.
    Returns (ordered proteins, {protein: section name})."""
    org = _ncbi_to_kegg_code(species)
    uniprot_map = _cached_kegg_table(org, "uniprot_map.json", f"conv/{org}/uniprot",
                                     "UniProt→gene map", _uniprot_to_gene)
    pathway_map = _cached_kegg_table(org, "pathway_map.json", f"link/pathway/{org}",
                                     "gene→pathway map", _gene_to_pathways)
    pathway_names = _cached_kegg_table(org, "pathway_names.json", f"list/pathway/{org}",
                                       "pathway names", _pathway_names)

    memberships: dict[str, set[str]] = {}
    for protein in proteins:
        identifier = short_accession(protein)
        gene_id = uniprot_map.get(identifier)
        memberships[protein] = set(pathway_map.get(gene_id, [])) if gene_id else set()
        if not gene_id:
            warnings.warn(
                f"KEGG: '{protein}' (id: '{identifier}') not in {org} UniProt map"
                " — appended at end"
            )
    return _order_by_memberships(proteins, memberships, pathway_names)


# -----------------------------------------------------------------------------
# STRING complexes
#
# STRING first translates our identifiers into its own IDs, then an
# enrichment query tells us which complexes and pathways they belong to.
# Results are cached per protein, so later runs only ask about new ones.
# -----------------------------------------------------------------------------

_STRING_BASE = "https://string-db.org/api/json"
_STRING_COMPLEX_CATEGORIES = {"CORUM", "PPI_hub_proteins", "KEGG_Pathways"}


def _string_resolve_ids(identifiers: list[str], taxid: int) -> dict[str, str]:
    """Map our identifiers to STRING IDs: {input name: STRING id}."""
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
    """The enrichment records STRING finds for these proteins."""
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
    """Put proteins that share STRING complexes or pathways next to each other.
    Only proteins that aren't cached yet are sent to STRING (two requests).
    Returns (ordered proteins, {protein: section name})."""
    cache_file = _cache_path(f"{species}_string", "complexes.json")
    names_cache_file = _cache_path(f"{species}_string", "term_names.json")
    cache: dict[str, list[str]] = _load_cache(cache_file) or {}      # bare id -> terms
    term_names: dict[str, str] = _load_cache(names_cache_file) or {}

    bare = {p: short_accession(p) for p in proteins}
    to_resolve = [p for p in proteins if bare[p] not in cache]

    if to_resolve:
        string_ids = _string_resolve_ids([bare[p] for p in to_resolve], species)
        time.sleep(1)

        if string_ids:
            enrichment = _string_fetch_enrichment(list(string_ids.values()), species)
            time.sleep(1)

            terms_of: dict[str, set[str]] = {sid: set() for sid in string_ids.values()}
            for record in enrichment:
                if record.get("category") not in _STRING_COMPLEX_CATEGORIES:
                    continue
                term = record["term"]
                term_names[term] = record.get("description") or term
                for gene in record.get("inputGenes", "").split(","):
                    gene = gene.strip()
                    for sid in string_ids.values():
                        if sid.endswith(gene) or gene in sid:
                            terms_of[sid].add(term)

            for protein in to_resolve:
                sid = string_ids.get(bare[protein])
                cache[bare[protein]] = [] if sid is None else sorted(terms_of.get(sid, set()))
                if sid is None:
                    warnings.warn(f"STRING: no entry found for '{protein}' — appended at end")

        _save_cache(cache_file, cache)
        _save_cache(names_cache_file, term_names)

    memberships = {p: set(cache.get(bare[p], [])) for p in proteins}
    return _order_by_memberships(proteins, memberships, term_names)

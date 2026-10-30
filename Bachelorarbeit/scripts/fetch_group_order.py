"""Group proteins by biology for the "pathway" and "complex" orderings.

KEGG tells us which pathways a protein is part of, Gene Ontology (through
UniProt) which complexes. Proteins that share a pathway or complex are put
next to each other, and each protein gets the name of its group as a section
label for the matrix headers. print_matrix.py uses this through
order_by_kegg() and order_by_go_complex().

Everything downloaded is cached as JSON under ~/.cache/xlms_matrix/, so each
lookup only hits the internet once."""
from __future__ import annotations

import json
import re
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


# -----------------------------------------------------------------------------
# Grouping proteins into named sections
#
# Used by both KEGG and GO: once we know which terms (pathways or
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
# Complexes from Gene Ontology
#
# UniProt lists the GO "cellular component" terms of every protein. Only the
# terms that are a kind of "protein-containing complex" (GO:0032991) count as
# complexes; QuickGO tells us which terms those are. Both are cached, the
# UniProt answers per protein, so later runs only ask about new ones.
# -----------------------------------------------------------------------------

_GO_COMPLEX_ROOT = "GO:0032991"   # protein-containing complex
_UNIPROT_BATCH = 100              # accessions per UniProt request
# what a UniProt accession looks like (from the UniProt documentation)
_UNIPROT_ACCESSION = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")


def _go_complex_terms() -> set[str]:
    """All GO terms that are a kind of protein-containing complex (about 2000)."""
    cache_file = _cache_path("go", "complex_terms.json")
    cached = _load_cache(cache_file)
    if cached is None:
        print("[GO] Downloading the list of complex terms (one-time, will be cached)…")
        url = (f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
               f"{_GO_COMPLEX_ROOT}/descendants?relations=is_a")
        data = json.loads(_get(url))
        cached = {"terms": data["results"][0].get("descendants", [])}
        _save_cache(cache_file, cached)
    return set(cached["terms"])


def _uniprot_cellular_components(accessions: list[str]) -> dict[str, list[list[str]]]:
    """The GO cellular component terms of each accession, as
    {accession: [[GO id, name], ...]}. Accessions UniProt doesn't know get
    an empty list."""
    result: dict[str, list[list[str]]] = {acc: [] for acc in accessions}
    for k in range(0, len(accessions), _UNIPROT_BATCH):
        batch = accessions[k:k + _UNIPROT_BATCH]
        url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions="
               + ",".join(batch) + "&fields=accession,go_c&format=tsv")
        # each line: accession, then "name [GO:id]; name [GO:id]; ..."
        for line in _get(url).splitlines()[1:]:
            acc, _, terms = line.partition("\t")
            result[acc] = [[m.group(2), m.group(1)]
                           for m in re.finditer(r"\s*([^;]+?) \[(GO:\d+)\]", terms)]
    return result


def order_by_go_complex(proteins: list[str], species: int = 9606) -> tuple[list[str], dict[str, str]]:
    """Put proteins that share a GO complex next to each other. The species
    isn't needed, a UniProt accession already belongs to one species.
    Returns (ordered proteins, {protein: section name})."""
    cache_file = _cache_path("go", "uniprot_cellular_component.json")
    cache: dict[str, list[list[str]]] = _load_cache(cache_file) or {}

    # "P12345-2" is an isoform of P12345, which carries the annotation
    bare = {p: short_accession(p).split("-")[0] for p in proteins}
    no_accession = [p for p in proteins if not _UNIPROT_ACCESSION.fullmatch(bare[p])]
    to_fetch = sorted({bare[p] for p in proteins if p not in no_accession and bare[p] not in cache})
    if no_accession:
        warnings.warn(f"GO: {len(no_accession)} protein(s) without a UniProt accession — appended at end")

    if to_fetch:
        print(f"[GO] Asking UniProt about {len(to_fetch)} protein(s)…")
        try:
            cache.update(_uniprot_cellular_components(to_fetch))
            _save_cache(cache_file, cache)
        except Exception as e:
            warnings.warn(f"GO: UniProt request failed — {e}")

    complexes = _go_complex_terms() - {_GO_COMPLEX_ROOT}   # the root itself says nothing
    memberships: dict[str, set[str]] = {}
    term_names: dict[str, str] = {}
    for protein in proteins:
        terms = cache.get(bare[protein], [])
        memberships[protein] = {go for go, _ in terms if go in complexes}
        term_names.update((go, name) for go, name in terms)
    return _order_by_memberships(proteins, memberships, term_names)

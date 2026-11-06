"""Group proteins by biology for the "complex" ordering.

Gene Ontology (through UniProt) tells us which complexes a protein is part
of. Proteins that share a complex are put next to each other, and each
protein gets the name of its group as a section label for the matrix
headers. print_matrix.py uses this through order_by_go_complex().

Everything downloaded is cached as JSON under ~/.cache/xlms_matrix/, so each
lookup only hits the internet once."""
from __future__ import annotations

import json
import re
import urllib.request
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

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
# Once we know which GO terms each protein belongs to, every protein is
# put into exactly one section, and the sections are sorted from largest
# to smallest.
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


def order_by_go_complex(
    proteins: list[str], locations: bool = False,
) -> tuple[list[str], dict[str, str]]:
    """Put proteins that share a GO complex next to each other (no species
    needed, a UniProt accession already belongs to one species).
    With `locations`, proteins without one of the complex sections are
    grouped by their other cellular component terms ("cytoplasm",
    "plasma membrane", ...) instead of all going to "Unknown"; these
    location sections come after the complex sections.
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

    all_complexes = _go_complex_terms()
    complexes = all_complexes - {_GO_COMPLEX_ROOT}   # the root itself says nothing
    memberships: dict[str, set[str]] = {}
    term_names: dict[str, str] = {}
    for protein in proteins:
        terms = cache.get(bare[protein], [])
        memberships[protein] = {go for go, _ in terms if go in complexes}
        term_names.update((go, name) for go, name in terms)
    if not locations:
        return _order_by_memberships(proteins, memberships, term_names)

    assignment = _assign_primary_groups(proteins, memberships, term_names)
    in_complex = [p for p in proteins if assignment[p] != "Unknown"]
    rest = [p for p in proteins if assignment[p] == "Unknown"]
    # second round for the rest, on the terms that aren't complexes
    places = {p: {go for go, _ in cache.get(bare[p], []) if go not in all_complexes} for p in rest}
    assignment.update(_assign_primary_groups(rest, places, term_names))
    return _group_and_sort(in_complex, assignment) + _group_and_sort(rest, assignment), assignment

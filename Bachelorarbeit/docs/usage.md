# XL-MS Matrix Viewer — Usage

Tools for turning a crosslink CSV into a protein-protein (or, optionally,
residue-residue) adjacency matrix, either as a text dump on the console or
as an interactive GUI.

## Setup

The project expects Python 3.11+ with `pandas`, `biopython`, `scipy`, `numpy`
and `Pillow` installed. `matplotlib` and `pygame` are optional — without them
the GUI falls back to greyscale cells and CPU-only (numpy) rendering.

```
pip install pandas biopython scipy numpy Pillow matplotlib pygame
```

## Input data

- **Crosslink CSV** — needs `Protein1`, `Protein2` and `Score` columns
  (case-insensitive), plus one of two decoy-label formats:
  - `Decoy1` / `Decoy2` (per-protein booleans), or
  - `isTT` / `isDD` (per-link target/decoy classification)
  - For `--level residue` (see below), also needs residue-position columns,
    in one of two formats:
    - **direct**: `SeqPos1`/`SeqPos2`, `Pos1`/`Pos2`, or `Position1`/
      `Position2` — the residue's absolute position in the protein.
    - **peptide + link offset**: `PepPos1`/`PepPos2` (where the linked
      peptide starts in the protein) together with `LinkPos1`/`LinkPos2`
      (the crosslinked residue's 1-based offset within that peptide). The
      absolute position is computed as `PepPos + LinkPos - 1`.
    Both are case-insensitive; whichever pair is present is used
    automatically.
- **FASTA** (optional) — protein sequences, used by `--order sequence`,
  `--order size`, and the GUI's per-protein length lookup.

An example CSV is provided at `data/example_500.csv`. It only has
`Protein1`/`Protein2`/`Score` columns, so `--level residue` isn't available
on it — use it for protein-level mode only.

## Command line: `scripts/print_matrix.py`

```
python scripts/print_matrix.py <csv_path> [options]
```

| Option | Description |
|---|---|
| `--level {protein,residue}` | axis granularity (default `protein`) |
| `--n INT` | number of unique proteins/residues to display (default 1000) |
| `--order MODE` | axis ordering, see below (default `confidence`) |
| `--fasta PATH` | FASTA file, required for `--order sequence` and `--include-unlinked` |
| `--species INT` | NCBI taxonomy ID for pathway/complex lookups (default 9606, human) |
| `--include-unlinked` | residue mode: also show residues with no crosslinks, not just linked ones (needs `--fasta`) |

`--order` modes:

- `confidence` — highest-scoring proteins first (default)
- `alpha` — alphabetical
- `cluster` — hierarchical clustering by crosslink score similarity
- `sequence` — hierarchical clustering by k-mer similarity (needs `--fasta`)
- `pathway` — grouped by KEGG pathway (needs `--species`)
- `complex` — grouped by STRING complex/pathway (needs `--species`)
- `size` — longest sequence first (needs `--fasta` for real lengths)

In `--level residue` mode, every `--order` mode above is applied to each
residue's *parent protein* (so pathway/complex/cluster/etc. group residues
the same way they'd group whole proteins), and residues within a protein are
then sorted by position. Protein-name section headers are always shown above
the columns in residue mode, regardless of `--order`.

Examples:

```
python scripts/print_matrix.py data/links.csv
python scripts/print_matrix.py data/links.csv --n 200 --order alpha
python scripts/print_matrix.py data/links.csv --order sequence --fasta data/proteins.fasta
python scripts/print_matrix.py data/links.csv --order pathway --species 9913
python scripts/print_matrix.py data/links.csv --level residue --order cluster
```

Output is a dot-matrix printed to the console (targets before decoys, split
by a `+`/`-` separator), followed by a score legend, a truncated-name legend,
and a memory usage summary.

### Residue-level mode

Column labels are the bare residue position (e.g. `142`); the parent protein
is shown via the always-on section header instead. Row labels are
`accession:position` (e.g. `P00001:142`) — truncated protein names are
unreadable at the CLI's row-label width, so only the accession is used. The
legend maps truncated row labels to the full `protein:position` string.

A protein contributing only one or two residues gets a very narrow section
header block that can't fit its name (renders as `==`) — this is a limit of
the fixed 2-char cell width shared with protein mode, not a bug; the row
label and legend still identify the protein.

By default only residues that are an actual crosslink endpoint in the CSV
are shown (`--n` selects the top-N most confident of these). Pass
`--include-unlinked` (with `--fasta`) to also add every other residue
position of each protein that already has at least one linked residue, so
the matrix shows the linked residues' position within the full sequence.
This can add a lot of columns for long proteins — the `--n` cutoff still
limits how many *linked* residues (and therefore which proteins) are pulled
in, it just doesn't cap the unlinked expansion within those proteins.
Unlinked residues always render as empty cells (no crosslink score exists
for them).

Minimal example CSV with direct position columns:

```csv
Protein1,Protein2,Score,SeqPos1,SeqPos2,Decoy1,Decoy2
sp|P00001|PROT001_HUMAN,sp|P00002|PROT002_HUMAN,95.2,12,340,False,False
sp|P00001|PROT001_HUMAN,sp|P00002|PROT002_HUMAN,80.1,55,340,False,False
sp|P00001|PROT001_HUMAN,sp|P00003|PROT003_HUMAN,60.0,12,20,False,False
```

The same data using the peptide + link-offset format instead (resolves to
the identical residue positions, since `10+3-1=12` and `330+11-1=340`):

```csv
Protein1,Protein2,Score,PepPos1,LinkPos1,PepPos2,LinkPos2,Decoy1,Decoy2
sp|P00001|PROT001_HUMAN,sp|P00002|PROT002_HUMAN,95.2,10,3,330,11,False,False
sp|P00001|PROT001_HUMAN,sp|P00002|PROT002_HUMAN,80.1,50,6,330,11,False,False
sp|P00001|PROT001_HUMAN,sp|P00003|PROT003_HUMAN,60.0,10,3,15,6,False,False
```

```
python scripts/print_matrix.py residues.csv --level residue
```

## GUI: `gui/matrix_app.py`

```
python gui/matrix_app.py [csv_path]
```

Opens a Tkinter window with the same CSV/FASTA/`--level`/`--n`/`--order`/
`--species` controls as the CLI, plus:

- **Level** — `protein` (default) or `residue`; switches the axis unit and
  updates the `N` label accordingly. Requires the CSV to have residue-
  position columns (see above).
- **Show only linked residues** — checkbox next to Level, enabled only in
  residue mode. Checked (default) matches the CLI default: only crosslinked
  residues are shown. Unchecking it requires a FASTA file and adds every
  other residue of each protein already in the matrix, same as
  `--include-unlinked`.
- **Zoom** slider — cell size in pixels; past a threshold, cells show the
  numeric score instead of a dot.
- **Colormap** — `(none)` for greyscale, or a matplotlib colormap (e.g.
  `coolwarm`) if matplotlib is installed.
- **Norm** — linear or quantile score normalization for the colormap.
- Click a cell to see both proteins' names, decoy status, section (if
  pathway/complex ordering is active), and the crosslink score — in residue
  mode, also the residue position. This also triggers a background UniProt
  lookup for the full protein name and length.
- **Export PNG…** — renders the full (unscrolled) matrix to a PNG file.

If `pygame` is installed, the matrix canvas renders through an embedded SDL2
surface for smoother scrolling on large matrices; otherwise it falls back to
a pure numpy/PIL renderer.

## Pathway/complex lookups: `scripts/fetch_group_order.py`

Not run directly — it's the backend for `--order pathway` / `--order complex`
and exposes two functions:

```python
order_by_kegg(proteins: list[str], species: int) -> tuple[list[str], dict[str, str]]
order_by_string(proteins: list[str], species: int) -> tuple[list[str], dict[str, str]]
```

Both return the reordered protein list plus a `{protein: section_label}`
dict used for the section headers. Results are cached on disk under
`~/.cache/xlms_matrix/` so repeated lookups for the same species don't hit
KEGG/STRING again.

## `xlms` package (`src/xlms`)

Shared data loading used by the scripts above, also usable standalone:

```python
from xlms import read_csv, read_fasta, load_dataset

crosslinks = read_csv("data/links.csv")          # list[CrossLink]
sequences = read_fasta("data/proteins.fasta")    # dict[str, str]
dataset = load_dataset("data/links.csv", "data/proteins.fasta")

dataset.target_target()   # crosslinks classified TT
dataset.above_score(40.0) # crosslinks with score >= 40.0
```

## Tests

```
python tests/matrix.py
```

Runs a set of smoke tests against `xlms.io` covering both decoy CSV formats
and FASTA parsing.

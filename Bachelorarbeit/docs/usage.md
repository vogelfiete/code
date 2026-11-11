# XL-MS Matrix Viewer — Usage

Tools for turning a crosslink CSV into a protein-protein (or, optionally,
residue-residue) adjacency matrix, either as a text dump on the console or
as an interactive GUI.

## Setup

The project expects Python 3.11+ with `pandas`, `biopython`, `scipy`, `numpy`
and `Pillow` installed. `matplotlib` and `pygame` are optional — without them
the GUI falls back to greyscale cells and CPU-only (numpy) rendering. All of
them are listed in `requirements.txt`. From the `Bachelorarbeit` folder:

**Windows**

```
python -m venv venv
venv\Scripts\activate
python -m pip install -r requirements.txt
```

**Linux / macOS**

```
python3 -m venv venv
source venv/bin/activate
python -m pip install -r requirements.txt
```

The GUI needs Tkinter, which comes with Python on Windows and macOS. On Linux
it is a system package that pip can't install:

- Debian / Ubuntu: `sudo apt install python3-tk` (and `python3-venv` for the
  venv step)
- Fedora: `sudo dnf install python3-tkinter`
- Arch: `sudo pacman -S tk`

For the labels the GUI looks for a monospace font (Consolas on Windows,
DejaVu Sans Mono or Liberation Mono on Linux, Menlo on macOS) and otherwise
uses Pillow's built-in font. On Linux the pygame renderer is placed inside
the Tk window, which works under X11 and XWayland. If the matrix stays empty,
uninstall pygame (`python -m pip uninstall pygame`) to use the numpy renderer.

## Input data

- **Crosslink CSV** — needs `Protein1`, `Protein2` and `Score` columns
  (case-insensitive), plus one of two decoy-label formats:
  - `Decoy1` / `Decoy2` (per-protein booleans), or
  - `isTT` / `isDD` (per-link target/decoy classification)
  - For `--level residue` (see below) and for opening proteins in the GUI,
    also needs residue-position columns, in one of two formats:
    - **direct**: `SeqPos1`/`SeqPos2`, `Pos1`/`Pos2`, `Position1`/
      `Position2`, or `fromSite`/`ToSite` (xiFDR "Links" files) — the
      residue's absolute position in the protein.
    - **peptide + link offset**: `PepPos1`/`PepPos2` (where the linked
      peptide starts in the protein) together with `LinkPos1`/`LinkPos2`
      (the crosslinked residue's 1-based offset within that peptide). The
      absolute position is computed as `PepPos + LinkPos - 1`.
    Both are case-insensitive; whichever pair is present is used
    automatically.
- **FASTA** (optional) — protein sequences, used by `--order sequence`,
  `--order size`, and the GUI's per-protein length lookup.

## Matrix orientation

The matrix is **directional, not symmetric**: each CSV row is drawn in the
cell at **row = `Protein1`, column = `Protein2`**. A row `A,B,x` fills cell
(row A, column B); a row `B,A,y` fills cell (row B, column A). A–B and B–A
are never merged, so a pair that only appears one way round in the CSV
fills only one of the two mirrored cells. If the same direction appears in
several rows, the highest score is shown. In residue mode the same applies
to the two residues of a row.

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
| `--fasta PATH` | FASTA file, required for `--order sequence`; in residue mode it also adds the residues without crosslinks |
| `--locations` | with `--order complex`: group proteins without a complex by their other GO cellular component terms (cytoplasm, membrane, …) instead of putting them all in "Unknown" |

`--order` modes:

- `confidence` — highest-scoring proteins first (default)
- `alpha` — alphabetical
- `sequence` — hierarchical clustering by k-mer similarity (needs `--fasta`)
- `complex` — grouped by protein complex, taken from the Gene Ontology
  "cellular component" terms in UniProt (only terms that are a kind of
  *protein-containing complex*, GO:0032991); needs UniProt accessions.
  With `--locations`, proteins outside the complex sections get location
  sections, placed after the complexes. Complexes that are part of a bigger
  one (the large and the small ribosomal subunit are part of the ribosome)
  stay next to each other, and the bigger complex gets its own banner line
  above the section line
- `size` — longest sequence first (needs `--fasta` for real lengths)

In `--level residue` mode, every `--order` mode above is applied to each
residue's *parent protein* (so complex/etc. group residues
the same way they'd group whole proteins), and residues within a protein are
then sorted by position. Protein-name section headers are always shown above
the columns in residue mode, regardless of `--order`.

Examples:

```
python scripts/print_matrix.py data/links.csv
python scripts/print_matrix.py data/links.csv --n 200 --order alpha
python scripts/print_matrix.py data/links.csv --order sequence --fasta data/proteins.fasta
python scripts/print_matrix.py data/links.csv --order complex --locations
python scripts/print_matrix.py data/links.csv --level residue --order alpha
```

Output is a dot-matrix printed to the console (targets before decoys, split
by a `+`/`-` separator), followed by a score legend, a truncated-name legend,
and a memory usage summary. The first line says which side of a crosslink
the axes are (rows = Protein1, columns = Protein2). With sections
(`complex`), each section also gets a `[name]` line above its first row, so
the sections can be read off on the left as well as on top; a section inside
a bigger complex is written `[bigger > section]`.

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

`--n` selects the top-N most confident residues that are an actual
crosslink endpoint in the CSV. With `--fasta`, every other residue
position of each protein that already has at least one linked residue is
added too, so the matrix shows the linked residues' position within the
full sequence. Without a FASTA only the crosslinked residues are shown.
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

Opens a Tkinter window. All controls sit in a panel on the left, so the
matrix uses the full window height; the panel scrolls (mouse wheel over it)
if the window is too short to show everything. The panel has four groups:

**Data** — the same CSV/FASTA/`--order` options as the CLI, the **Group
rest by location** checkbox (`--locations`, only shown for the `complex`
order), plus **Load**, **Export PNG…** and the load status. Unlike the CLI
there is no `--n` and no `--level`: the GUI always loads all proteins and,
if the CSV has residue-position columns, all their residues too, and starts
at protein level (see *The axis tree* below). With a FASTA file, every
residue of the loaded proteins is there, also those without crosslinks
(same as the CLI); without a FASTA only the crosslinked residues are.

- **Export PNG…** — renders the full (unscrolled, unblocked) matrix to a PNG
  file, opened as it is on screen, with the score cutoff applied, the axis
  titles in the corner, the annotation tracks and the colour legend under
  the matrix.

The empty corner above the row names says which side of a crosslink the
axes are: `row: Protein1 ↓`, `col: Protein2 →` (the CSV row A,B is in row A,
column B).

**The axis tree** — the axes are a tree: with the `complex` order,
proteins sit in their complex, and complexes that are part of a bigger one
sit in that (e.g. *ribosome* → *cytosolic ribosome* → *cytosolic large
ribosomal subunit*); every protein holds its residues. Rows and columns
always show the same nodes. Each open group or protein gets a band with its
name above the columns and, turned to read bottom to top, left of the row
labels — one band row per level — so links between two groups can be read
off from both axes.

- **Click a name** (a column header or a row label) to open it: a protein
  opens into its residues, a closed group (written `+name`) into what it
  holds.
- **Click a band** to close that group or protein again: it becomes a single
  row and column, whose cells combine all its crosslinks with the
  **Aggregate** method.
- **Show** (in Display) sets the whole tree at once: **Top** closes the
  biggest groups, **Complex** the smallest ones (so every complex is one
  row), **Protein** shows every protein and **Residue** every residue.

Proteins and residues can be mixed: a closed protein next to an open one
gets the crosslinks between it and each of the open protein's residues.
Opening or closing keeps the clicked node where it was on screen. Proteins
can only be opened when the CSV has residue positions.

A band name that doesn't fit on one line is wrapped onto a second; if it
still doesn't fit, it's cut off with "…" (clicking a cell shows the full
name in Selection). A band bigger than the view keeps its name centred on
the visible part.

**Display**

- **Colormap** — `(none)` for greyscale, or a matplotlib colormap (e.g.
  `coolwarm`) if matplotlib is installed.
  With a colormap, the colour follows the rank of a score among all scores
  (quantile), so the colours are spread evenly even when most scores lie
  close together.
  The colour legend under **Fold** shows what the colours mean right now:
  what a cell holds (score, an Aggregate of scores, or crosslinks per cell),
  whether colours follow the rank or are linear, and the smallest, middle
  and largest value.
- **Zoom** — cell size in pixels; from 18 px on, cells show the numeric
  score. Zooming below 4 px merges neighbouring items into blocks. The
  slider zooms around the middle of the view. Big matrices (e.g. residue
  level) get extra zoom steps, so the whole matrix can be zoomed out to
  about 200 px. Blocks never cross the target/decoy line. Runs of nodes
  under the same group or protein smaller than a block are packed into one
  block whole; bigger ones are cut into several blocks. While zoomed out
  into blocks, the single item names are hidden (they can't be read at that
  size anyway), with a little padding and a grey border around the matrix.
  The band names of open groups and proteins stay on both axes, since
  that's the view where links between groups are looked at; without any
  bands the headers go completely and the matrix uses the whole area.
  While residues are on the axes the names stay instead,
  as a sparse position axis — one tick per block, spaced out so the numbers
  don't overlap — so there's still a sense of where you are even zoomed
  out past individual residues.
- **Zoom to** — **All** zooms and centers so that all quadrants that are
  switched on fill the view (counting the room the labels take at the
  resulting zoom); **TT**, **TD**, **DT** and **DD** do the same
  for one quadrant and switch it on first if it was hidden. Buttons for
  quadrants without items (e.g. no decoys loaded) are greyed out.
- **Aggregate** — what a block shows: the geometric mean of its scores
  (default), the maximum score, the root sum of squares of its scores
  (√(s₁² + s₂² + …), which grows with both the scores and the number of
  crosslinks), or the count of crosslinks in it. Root sums of squares and
  counts get their own colour scale, by their rank among all values on
  screen.
- **Fold to triangle** — a toggle button that folds the matrix along its
  diagonal: the crosslinks A→B and B→A of a pair end up together in one
  cell below the diagonal, combined with the **Aggregate** method (also for
  single cells, e.g. the geometric mean of both directions), and the upper
  half stays empty. TD links fold into the DT quadrant. Clicking the empty
  half selects the mirror cell; the selection shows both directions and the
  combined value. **Export PNG…** exports the folded matrix too.
- **Show** — **Top**, **Complex**, **Protein**, **Residue**, see *The axis
  tree* above.

Column names in the header are written at 45°, going up to the right from
the middle of their column.

**Quadrants** — four checkboxes, laid out like the matrix: **TT**
(target row, target column), **TD** (target row, decoy column), **DT** and
**DD**. Rows are Protein1, columns Protein2. Unchecking one hides all links
in that quadrant; like the cutoff, they are then left out of block
aggregates, the PNG export and the "links shown" count. A clicked cell in a
hidden quadrant is marked `(hidden: TD quadrant)` and so on.

**Annotations** — UniProt sequence features as coloured tracks along the
residue rows and columns, one thin strip per type next to the matrix:
**Helix**, **Beta strand**, **Turn**, **Binding site**, **Active site**,
**Domain**, **Transmembrane** and **Modified residue**, each with its own
colour (shown next to its checkbox). The features of all loaded proteins
are fetched in the background the first time a type is ticked, and cached
in `~/.cache/xlms_matrix/uniprot/features.json`. Tracks only appear while
residues are on the axes; zoomed out into blocks, a block is coloured if
any of its residues has the feature. Clicking a residue cell lists the
features it lies in under Selection.

**Score cutoff** — one slider with two handles: the left handle is the
bottom limit, the right handle the top limit (both inclusive). Cells whose
score is outside the range are hidden; they are also left out of block
aggregates and of the PNG export. Colours keep the full-range scale, so a
score has the same colour whatever the cutoff. Above the slider, on the same
score scale, a plot shows how the scores of the crosslinks on the axes are
spread: targets (TT) as blue bars, decoys (TD, DT, DD) as an orange line,
on a logarithmic count scale; only the quadrants that are switched on count,
and bins outside the cutoff turn grey. The label shows the range and how
many links are visible; **Reset** restores the full range. Loading new data
resets the cutoff.

**Moving around** — the matrix works like a map. Turn the mouse wheel over
it (or over the column headers or row labels) to zoom in and out; the spot
under the pointer stays where it is, also near the matrix edges and when
blocks start or stop — for that the view can be moved past the edges until
only a small part of the matrix is left. Hold the left mouse button and
drag to move the matrix (the column headers and row labels move along). The
scrollbars and the arrow keys scroll as well.

**Selection** — click a cell to see both sides' names (a protein, a
residue, or a closed group with its number of proteins), decoy status, the
group they are in, and the crosslink scores in both directions: `row→col`
(the clicked cell) and `col→row` (its mirror cell) — for residues also the
position and the annotations they lie in. A cell hidden by the
cutoff shows its score marked `(hidden by cutoff)`. When zoomed out into
blocks or on a closed group, the aggregated score of each direction is
shown. Clicking also
triggers a background UniProt lookup for the full protein name and length.
A click only counts if the mouse wasn't dragged in between, and the
selection stays on the same proteins when you zoom (opening or closing
clears it).

If `pygame` is installed, the matrix canvas renders through an embedded SDL2
surface for smoother scrolling on large matrices; otherwise it falls back to
a pure numpy/PIL renderer.

## Complex lookups: `scripts/fetch_group_order.py`

Not run directly — it's the backend for `--order complex` and exposes:

```python
order_by_go_complex(proteins: list[str], locations: bool = False) -> tuple[list[str], dict[str, str]]
```

It returns the reordered protein list plus a `{protein: section_label}`
dict used for the section headers.

```python
complex_hierarchy(sections) -> dict[str, str]
```

For those section names it returns the bigger complexes they are part of,
as `{name: name one level up}`. It follows the GO *part of* relations
(through QuickGO, cached in `go/ancestors_*.json`), keeps only assemblies
(kinds of protein-containing complex or of membraneless organelle — GO
files the ribosome there), not places like *cytosol*, and only those that
hold at least two of the sections; a general term (*ribosome*) that a more
specific one (*cytosolic ribosome*) already covers is dropped. The
`complex` ordering uses it to keep subunits together, the CLI for its
banner lines and the GUI for the axis tree. Results are cached on disk under
`~/.cache/xlms_matrix/` so repeated lookups don't hit UniProt or
QuickGO again. `order_by_go_complex` asks UniProt for the GO cellular
component terms of each accession (cached per protein in
`go/uniprot_cellular_component.json`) and keeps only the terms below
"protein-containing complex", which it gets from QuickGO once
(`go/complex_terms.json`). Isoform accessions (`P12345-2`) use the
annotation of the main entry; names that aren't UniProt accessions (e.g.
`decoy:P01267`) end up in "Unknown". Delete the cache files to fetch
fresh annotations.

## Annotations: `scripts/fetch_features.py`

Backend for the GUI's annotation tracks:

```python
features_for(proteins: list[str]) -> dict[str, list[list]]
```

returns `{protein: [[type, start, end, text], ...]}` for the feature types
in `FEATURE_TYPES`, from the UniProt REST API in batches of 100 (positions
from 1, both ends included; `text` is the description, or the ligand of a
binding site). Cached per accession in `uniprot/features.json`.

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

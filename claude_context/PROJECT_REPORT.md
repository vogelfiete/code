# Bachelorarbeit — Project Context Report (for Claude)

> First written 2026-09-30 from a full read of every source file, the git history, the real
> datasets in `Desktop/BA/`, and live runs of the tests + CLI (HEAD then: `65cc9af` 09092026_zooming).
> **Updated 2026-09-30** after the refactor + directional-matrix change (uncommitted at time of writing
> — check `git log`/`git status`). Line numbers refer to that refactored state; re-verify after new commits.
> Paths are relative to `Desktop/BA/code/` unless absolute.

---

## 0. TL;DR

- **What:** A bachelor thesis tool for **XL-MS (crosslinking mass spectrometry)** data. It turns a
  crosslink CSV into a **protein×protein (or residue×residue) adjacency matrix**. Two front-ends:
  a **console dot-matrix CLI** and a **Tkinter GUI** (zoom, block aggregation, colormaps,
  click-to-inspect, PNG export).
- **The matrix is DIRECTIONAL (user requirement):** cell (row A, col B) = best score of CSV rows
  with `Protein1=A, Protein2=B`; rows `B,A` go to (row B, col A). A–B and B–A are never merged.
- **Core thesis idea:** *axis ordering* so structure becomes visible — confidence, alphabet,
  sequence k-mer similarity, KEGG pathway, STRING complex, size. Targets always before decoys.
  (`cluster` by score was **removed** on user request, "for now" — may come back; it would need a
  symmetric distance, e.g. max/mean of both directions — ask the user.)
- **User decisions to respect:** no FDR filtering (intentional); keep the two CSV readers as they are
  (`xlms.read_csv` vs print_matrix reader, see §4.2).
- **Stack:** Python 3.13 (venv at `code/venv`), pandas, numpy, scipy, biopython, Pillow, optional
  matplotlib + pygame. No packaging; scripts use `sys.path.insert`.
- **User:** German student (git user `vogelnestfiete`), commits named `DDMMYYYY_topic`. Wants
  simple, readable code. Commits themself — only commit when asked.
- **Run anything with `venv/Scripts/python.exe`** — system Python lacks biopython.

---

## 1. Directory map

```
Desktop/BA/
├── 2SDA-rep1/2SDA-rep1/                 REAL DATA (xiFDR 2.3.11 output, ~95 MB Links CSV)
├── SCE_no_eDR_FBS_DSSO_NEW/SCE_.../     REAL DATA (xiFDR 2.3.11 output, ~7 MB Links CSV)
├── *.zip                                 zipped copies of the above
└── code/                                 ← git repo root
    ├── .gitignore        (/.venv /99eur /99euro /c /java /Python /websites — leftovers)
    ├── venv/             Python 3.13.2 venv (the one actually used; NOT .venv)
    ├── claude_context/   ← THIS REPORT
    └── Bachelorarbeit/
        ├── requirements.txt     pip packages with minimum versions (2026-10-13); Linux setup notes in docs/usage.md
        ├── pyrightconfig.json   points to ../.venv (WRONG — real venv is ../venv)
        ├── data/example_500.csv synthetic: 31,525 rows, 482 proteins, cols Protein1,Protein2,Score only
        ├── docs/usage.md        user-facing docs (CRLF line endings; keep in sync with features)
        ├── src/xlms/            small "library": models.py, io.py, __init__.py
        ├── scripts/print_matrix.py       CORE: CSV→matrix + orderings + CLI renderer (~615 lines)
        ├── scripts/fetch_group_order.py  KEGG / STRING ordering backends + short_accession (~280 lines)
        ├── gui/matrix_app.py             Tkinter GUI (~1080 lines, sectioned)
        └── tests/matrix.py               4 smoke tests for xlms.io only (plain asserts, no pytest)
```

`__pycache__/*.pyc` files ARE tracked in git (show as modified). Not ignored.

---

## 2. Domain primer

- **XL-MS:** a crosslinker (here DSSO or "2SDA") covalently links two spatially close residues
  (usually lysines). MS identifies both peptides → two proteins + positions. Inter-protein links ⇒
  protein–protein interaction evidence.
- **Target/Decoy:** searches include decoy proteins; links are **TT / TD / DD**. The tool shows decoys
  as a separate block after targets (separator `|`/`+` in CLI, grey lines in GUI).
- **Score:** search-engine score (higher = better). Scale differs per file (2SDA Links: median 1.75, 75 % below 2.7, max 58 — checked 2026-10-13; CSM ~24,
  ppi ~164–187, example CSV 0–1).
- **xiFDR levels:** CSM (per spectrum) → PeptidePairs → **Links** (unique residue pairs) → **ppi**
  (protein pairs). ppi files list each pair in one direction only; Links/CSM contain both
  directions for many pairs (2SDA Links: ~27.7k pairs both ways).
- **Residue position:** absolute = `PepPos + LinkPos − 1`.

---

## 3. Architecture & data flow

```
CSV ──► build_matrix / build_residue_matrix  (thin wrappers) ──► _build(path, n, level)   [print_matrix.py]
        ├─ column detection via lower-case map (protein1, protein2, score required → ValueError)
        ├─ sort rows by Score desc; per-row items: names (protein) or ResidueId (residue; None if pos bad)
        ├─ pass 1: collect N unique items (skip ambiguous ';' names; count rows w/o positions → warning)
        ├─ _classify_proteins → {protein: is_decoy};  items = targets first, decoys second (stable)
        └─ pass 2: sparse {(item_from_Protein1, item_from_Protein2): best_score}   ← DIRECTIONAL
     ──► [residue + include-unlinked] add_unlinked_residues(FASTA)
     ──► sort_proteins / sort_residues(items, is_decoy, mode, sequences, species)
          └─ _apply_order per target group and decoy group → _ORDERINGS[mode](group, sequences, species)
               pathway/complex → fetch_group_order.order_by_kegg / order_by_string (HTTP + disk cache)
     ──► CLI: print_matrix() to stdout      GUI: MatrixApp._on_loaded() → blocks → render
```

**Invariants**
- Axis item = `str` protein name or `ResidueId(protein, pos)` NamedTuple; `protein_of(item)` gives the
  parent protein for both. Same code serves both levels.
- **Sparse keys are directional `(row_item, col_item)` = (Protein1 side, Protein2 side)**; lookup for
  cell (i, j) is `sparse.get((items[i], items[j]))`. Duplicate same-direction rows keep the max.
  **Never** reintroduce `(min, max)` keys.
- Orderings no longer see scores (no `sparse` argument) — they order by name/sequence/annotation only.
- `split` = index of first decoy.
- GUI imports public names from print_matrix: `COL_W, DEFAULT_N, ORDER_MODES, ROW_W,
  add_unlinked_residues, build_matrix, build_residue_matrix, protein_of, short_accession,
  sort_proteins, sort_residues`. Changing them breaks the GUI.
- `src/xlms` is decoupled: CLI/GUI only use `xlms.io.read_fasta`.

---

## 4. File-by-file reference

### 4.1 `src/xlms/models.py`
- `CrossLink` frozen/slots dataclass: `id, protein1, protein2, seq_pos1, seq_pos2, score,
  is_decoy, is_tt, is_td, is_dd`. `CrossLinkDataset(crosslinks, sequences)` with filters
  `target_target()`, `target_decoy()`, `decoy_decoy()`, `above_score(t)`.

### 4.2 `src/xlms/io.py` (user: leave as is)
- `read_csv(path)`: strict, case-sensitive; needs `Id, Protein1, Protein2, SeqPos1, SeqPos2, Score`
  + either `{isDecoy,isTT,isTD,isDD}` or `{Decoy1,Decoy2,DecoyType}`. Per-link TT/TD/DD. Only used by
  tests; can't read xiFDR files (no `Id`, no `SeqPos*`).
- `read_fasta(path)` → `{record.id: SEQ_UPPER}` (keys = full header token, see gotcha §8.1).
- `load_dataset(csv, fasta)`.

### 4.3 `scripts/print_matrix.py`
Constants `DEFAULT_N=1000, COL_W=10, ROW_W=12, CELL_W=2, N_LEVELS=4, DOTS=["·","○","●","⬤"]`,
decoy name prefixes `("decoy","rev_","contam_")`. Rewraps stdout as UTF-8 at import.
Sections: axis items / value parsers / column detection / decoy classification / construction /
ordering / console rendering / CLI.

| Name (line) | Purpose |
|---|---|
| `ResidueId` (46), `protein_of` (60) | axis item type; parent protein of any item |
| `_parse_position` (78), `_combine_pep_link` (88) | int or None (NaN / ambiguous `"114;53"` → None) |
| `_detect_decoy_cols` (111) | format "B" Decoy1/2, "A" isTT/isDD, "none" (takes lower-case col map) |
| `_detect_position_cols` (121) | "direct" SeqPos/Pos/Position pairs, "pep_link" PepPos+LinkPos |
| `_positions` (137) | per-row position lists for both ends |
| `_classify_proteins` (150) | {protein: is_decoy}; last row wins; format A: TD-only → name guess |
| `build_matrix` (196) / `build_residue_matrix` (209) → `_build` (215) | unified two-pass build, directional sparse |
| `add_unlinked_residues` (278) | adds every position 1..len(seq) (needs FASTA key match) |
| `_order_*` (315–366), `_ORDERINGS` (370), `ORDER_MODES` (377) | ordering functions + registry |
| `sort_proteins` (387) / `sort_residues` (403) | targets/decoys ordered separately; residue sections = protein |
| `print_matrix` (469) | console render (vertical col labels, dots, scale legend, name legend) |
| `main` (552) | argparse `csv --level --n --order --fasta --species --include-unlinked` |

### 4.4 `scripts/fetch_group_order.py`
- `short_accession` (22): `sp|P12345|X desc` → `P12345`; else first token. Single shared definition.
- Cache `~/.cache/xlms_matrix/<folder>/*.json` (exists locally: `hsa_kegg/` full, `1423_string/` empty).
- KEGG: `_ncbi_to_kegg_code` (hardcoded 8 species else KEGG lookup), `_cached_kegg_table` (136) =
  generic download-TSV-parse-cache; builders `_uniprot_to_gene`, `_gene_to_pathways`, `_pathway_names`.
  `order_by_kegg` (175).
- STRING: `_string_resolve_ids`, `_string_fetch_enrichment`, `order_by_string` (236); categories
  `CORUM, PPI_hub_proteins, KEGG_Pathways`; gene matching `sid.endswith(gene) or gene in sid`.
- Shared: `_assign_primary_groups` (top 10 terms, first match, else "Unknown"), `_group_and_sort`
  (largest section first, Unknown last), `_order_by_memberships`.

### 4.5 `gui/matrix_app.py`
Module helpers: `_find_pil_font`, `_make_lut` (`matplotlib.colormaps`), `_score_to_bg`
(quantile via searchsorted only — the linear colormap scale and the Norm control were removed 2026-10-10 on user request; greyscale without a colormap is still linear), `_auto_fg`, `_aggregate(values, method)`.
`RangeSlider(tk.Canvas)`: custom two-handle slider (`set_bounds`, `values`, `command` while
dragging, `release_command` on release; handles can't cross; ends map exactly to min/max).
`MatrixApp(tk.Tk)`, sectioned with comment banners (use grep for line numbers):
- **UI construction:** `_build_sidebar()` = fixed-width (`SIDEBAR_W=270`) **left panel**, vertically
  scrollable via `_scrollable()` (global `<MouseWheel>` binding that only acts when the pointer is over
  the panel). Groups: **Data** (CSV/FASTA, Level, only-linked, N, Order, Species, Load, Export, status),
  **Display** (Colormap, Zoom, Aggregate), **Score cutoff** (RangeSlider, label, Reset),
  **Selection** (wrapped info text); footer = item count + renderer. No top bar / status bar any more
  (added 2026-09-30 on user request: bigger matrix area). Matrix grid packs to the right of it.
  `_init_pygame` embeds SDL via `SDL_WINDOWID` into `_pg_host`, a background-less `tk.Frame` placed over the whole matrix canvas (since 2026-10-10: drawing straight into the canvas flashed white, because Tk repaints the canvas after every scroll step; changing the canvas cursor did the same, so the drag cursor is set through pygame). Over that frame the mouse buttons only reach pygame; the wheel is bound on the frame, the arrow keys on `_mx`. pygame renders into an off-screen `_pg_frame` only from `_redraw` (real changes); `_pg_loop` (16 ms) just polls pygame events (`_poll_pygame_events`: left button down/motion/up → `_press`/`_drag_to`/`_release`) and re-shows the frame via `_present()` (~0.4 ms; `_pg_sync_size` handles resize). Render errors go to `_report_render_error` (stderr + status line, once per message). Tk `<MouseWheel>` works over the SDL window (verified with real OS input), no forwarding needed.
- **Score cutoff** (user request 2026-09-30): `_cut_lo/_cut_hi` (reset to min/max score on each load),
  `_in_cutoff(score)` inclusive. Applied in `_block_score` (single cells), `_recompute_block_buckets`
  (filtered *before* aggregation, so counts/aggregates only use visible links), `_render_full_image`,
  and `_direction_score` (shows `(hidden by cutoff)`). Colours deliberately keep the **full-range**
  scale (user decision). `_on_cutoff` re-buckets + redraws; `_on_cutoff_release` refreshes selection.
- **Geometry & scrolling**: `_visible_blocks()` → (row range, col range); `_schedule()` coalesces redraws;
  `_scroll_to(ox, oy)` = pixel-exact view position on all three canvases (Tk clamps to the scroll region).
- **Moving around with the mouse** (added 2026-10-10, user request "like Google Maps"): `_wheel` zooms
  around the pointer (calls `_on_zoom(v, anchor)` then moves the slider; step = 15 % of the cell size, 1 slider
  step in the block range; the wheel no longer scrolls). `_press`/`_drag_to`/`_release` are shared by the Tk
  bindings and the pygame events: moving more than `DRAG_THRESHOLD_PX` pans via `_scroll_to`, otherwise the
  release is a click → `_select_at`. A pygame motion event without the button ends a drag (lost button-up).
- **Loading** (459): `_load` validates, worker thread runs the CLI pipeline, `after(0)` → `_on_loaded`.
- **Colours** (561): `_cell_colour(score)` (colormap LUT or greyscale 200→40), `_cell_text(score, px)`
  (score text at ≥18 px). Used by both renderers and export.
- **Zoom & blocks**: slider −16..32; <4 → cell 4 px, `block_size = round(1.3^(4−v))` ∈ [2,500].
  `_on_zoom(value, anchor=None)` keeps the item under the anchor (default: view centre) in place via
  `_item_at`/`_pixel_of` (fractional item positions, because blocks change what a pixel covers); it ignores a
  value equal to `_zoom_value` (Tk's Scale fires its command late, after the wheel already applied it).
  `_recompute_blocks` (never crosses split/section) → `_group_links` (numpy: links → cells, cutoff applied,
  stable sort so sums keep the old order; **directional** key `row block * n + col block`) →
  `_reduce_block_scores` (`_aggregate_groups`; exact same floats as before). Aggregate options since 2026-10-13 (user request): Geometric Mean (default), Max Score, Count — Mean removed (still used inside `_aggregate` as the geomean fallback). Count is coloured against `_colour_ref` = sorted cell counts instead of `_score_sorted`. Result = cell arrays
  `_cell_rows/_cell_cols/_cell_scores/_cell_rgb` (sorted row by row) used for drawing; `_cell_at(rb, cb)` =
  (value, link count) for the selection text. The old dicts `_block_scores`/`_block_pair_n`/`_block_pair_raw`
  and `_block_score()` are gone.
- **Selection**: `_selected_items` remembers the clicked row/column item, so zooming re-selects the block
  that contains them; `_on_loaded` clears the selection. `_select_at(x, y)` (called from `_release`) → `_select_block`
  shows `_describe_side` for row/col + `_score_line`, then a UniProt fetch thread (only if each side
  is one protein; `_sel_token` guards stale replies; `_clear_selection` also bumps the token).
  `_score_line` shows **both directions**: `row→col … | col→row …` (single line on the diagonal).
- **Rendering** (reworked 2026-10-10 for speed, pictures verified identical): `_visible_cells()` masks the
  cell arrays to the viewport; `_matrix_pixels(w, h)` paints one pixel per cell and `repeat`s it up to the
  cell size (shared by both renderers; pygame gets it via `surfarray.blit_array`); `_score_texts()` yields the
  cells that get a score written in (≥18 px) — pygame cell-surface cache keyed (colour, text), numpy pastes
  cached `_text_mask`s; `_colours_for(scores)` = array version of `_cell_colour`. Headers/labels paste cached
  `_label_image(text, vertical)` masks instead of `draw.text` per character. `_overlay_positions` for
  separator + highlight; `_section_runs` shared by header + export.
- **Export** (990): `_render_full_image` — full matrix, one cell per item (no aggregation),
  `cp = min(cell_px, 16)` → never shows score text.

### 4.6 `tests/matrix.py`
`venv/Scripts/python.exe Bachelorarbeit/tests/matrix.py` → 4 tests pass. Covers only `xlms.io`.
No tests for print_matrix / GUI (regression harnesses were used ad hoc — see §9).

---

## 5. Ordering modes (`_ORDERINGS`, print_matrix.py)

| Mode | Method | Needs | Sections? |
|---|---|---|---|
| `confidence` | passthrough (order of first appearance in score-desc walk) | – | no |
| `alpha` | `sorted()` | – | no |
| `sequence` | 3-mer frequency vectors, cosine distance, average linkage; missing seqs appended | FASTA | no |
| `pathway` | KEGG | network/cache, `--species` | **yes** |
| `complex` | STRING enrichment | network/cache, `--species` | **yes** |
| `size` | longest sequence first; unknown length last | FASTA (else no-op) | no |
| ~~`cluster`~~ | removed 2026-09-30 (user: "for now") | | |

Residue mode applies the mode to parent proteins; section = protein name, always on. Groups of <2 unchanged.

---

## 6. The real datasets (in `Desktop/BA/`, outside the repo)

Both are **xiFDR 2.3.11** exports (+ `.mzid`): `CSM, Linear_PSM, Linear_Peptides, PeptidePairs, Links,
ppi, proteingroups, summary, xiVIEW`, `mzid.summary`.

| Dataset | Organism hint | Links rows | CSM rows | ppi rows | proteins |
|---|---|---|---|---|---|
| `SCE_no_eDR_FBS_DSSO_NEW` | Bovine (FBS; FETUA_BOVIN, ALBU_BOVIN; OX=9913), DSSO | 17,662 | 18,493 | 10,588 | ~1,133 |
| `2SDA-rep1` | ATP synthase sample (`AtpG`), 2SDA, fractions f5–f10; organism unverified | 217,227 | 236,379 | 61,241 | ~428 |

- Protein IDs are bare accessions (`P12763`); decoys `decoy:P01267`; ambiguous groups use `;`.
- Decoy columns `Decoy1/Decoy2` → format "B". Also `isTT/isTD/isDD, fdr, fdrGroup (Self/Between)`.
- **Links** positions are `fromSite`/`ToSite` → NOT recognized → residue mode fails on Links.
  CSM/xiVIEW have `PepPos1/2 + LinkPos1/2` → works; CSM also `ProteinLinkPos1/2` (not detected).
- Many decoys; only 151 SCE Links rows at fdr ≤ 5% — **no FDR filtering, by design**.
- Perf after refactor: 2SDA Links `--n 200` ≈ 5 s (was 61 s with `iterrows`).

---

## 7. How to run

```bash
cd Desktop/BA/code
PY=venv/Scripts/python.exe
$PY Bachelorarbeit/tests/matrix.py
$PY Bachelorarbeit/scripts/print_matrix.py Bachelorarbeit/data/example_500.csv --n 12
$PY Bachelorarbeit/scripts/print_matrix.py ../SCE_no_eDR_FBS_DSSO_NEW/SCE_no_eDR_FBS_DSSO_NEW/SCE_no_eDR_FBS_DSSO_NEW_CSM_xiFDR2.3.11.csv --level residue --n 50
$PY Bachelorarbeit/gui/matrix_app.py [csv]
```
venv: biopython 1.88, matplotlib 3.11.2, numpy 2.5.3, pandas 3.0.6, pillow 12.3.0, pygame 2.6.1,
scipy 1.18.1. Windows; Git Bash + PowerShell. When scripting the GUI headlessly: drive it from inside
`mainloop()` (worker threads use `after()`), patch `messagebox.showerror`, and write files with
`encoding="utf-8"` (cp1252 default breaks on `→`).

---

## 8. Gotchas & known limitations (check before "fixing")

1. **FASTA key mismatch:** `read_fasta` keys = full header token (`sp|P12763|FETUA_BOVIN`), xiFDR
   uses bare `P12763` → `sequence`, `size`, `--include-unlinked`, GUI "aa" silently find nothing.
2. **KEGG section labels never get names:** `pathway_names.json` keys are `hsa01100`, memberships use
   `path:hsa01100` → labels show raw ids. Pre-existing; told user, not fixed yet.
3. `pyrightconfig.json` points to `../.venv`; actual venv is `code/venv`.
4. Two CSV readers with different contracts (§4.2) — intentional for now.
5. `_build` may return **N+1** items (adds both ends of a row before checking `>= n`).
6. ppi-level files and `example_500.csv` give triangular matrices (each pair listed once) — expected.
7. (FIXED 2026-09-30) pygame idle loop used to rebuild the frame every 16 ms (34–72 % main thread busy while idle) and swallow errors; now idle cost ≈ 0 redraws, errors reported.
8. Export ignores block aggregation, never shows score text; very large N → huge PNG.
9. STRING gene matching is substring-based → possible false memberships.
10. Residue-mode CLI section headers render as `==` for proteins with 1–2 residues (by design).
11. `print_matrix.py` rewraps `sys.stdout` at import time (also when the GUI imports it).
13. (FIXED 2026-10-10, see Selection in §4.5) **Stale selection (reproduced 2026-09-30):** `_selected` (block indices) is not cleared on
    `_on_loaded` or when zoom changes the blocks → highlight on the wrong cell, and
    `_on_cutoff_release` → `_select_block` raises IndexError after zooming out / loading smaller data.
14. **N spinbox:** non-numeric N → `self._n_var.get()` raises TclError in `_load` (uncaught, silent).
15. (FIXED 2026-10-10) Clicking the matrix now calls `focus_set()` in `_press`.
16. (FIXED 2026-10-10, user target: redraw ≤ 50 ms) Redraw measured in a 1500×950 window (matrix canvas
    1111×758), median ms, before → after: pygame worst case 282 → 32 (2SDA Links N=1000, 4 px), numpy
    renderer 340 → 48; zoom step into the block range 82 → 36. Remaining cost at 4 px is mostly the header
    strip (~14 ms) and row labels (~6 ms). Tk canvas scrolling in the numpy renderer can still expose white
    strips for a frame.
12. `__pycache__` binaries are committed to git. Python files written by Claude use LF; git autocrlf converts.

---

## 9. History

| Date | Commit | What changed |
|---|---|---|
| 2026-06-28 | `b24c836` | xlms package, print_matrix CLI, fetch_group_order, tests |
| 2026-07-13 | `3ba323b` | GUI created |
| 2026-08-04 | `6b86079`, `19fa90a` | example CSV, GUI work, PyInstaller build (then exe removed) |
| 2026-08-16 | `b051747`, `c6e35bc` | build artefacts removed; residue level, pep+link positions, include-unlinked, docs |
| 2026-09-09 | `65cc9af` | GUI zoom-out block aggregation (mean/geomean/max) |
| 2026-09-30 | (uncommitted) | Claude refactor + directional matrix + cluster removed (details below) |

**2026-09-30 refactor (Claude), verified:**
- Phase A, behaviour-preserving: unified `_build`, `iterrows` → column lists, ordering registry,
  public `sort_*`/`add_unlinked_residues`/`protein_of`/`short_accession`, KEGG loaders merged, dead code
  removed (glyph atlases, `_cluster_by_memberships`, duplicate constants, unused params), GUI
  colour/selection/click logic unified. Verified: 18 golden CLI outputs byte-identical; old-vs-new
  KEGG (real cache) and STRING (mocked) identical; GUI harness 696 states (numpy + pygame image
  hashes, selection texts, exports) identical.
- Phase B: directional sparse keys + GUI both-direction info, `cluster` removed, docs updated.
  Verified: items/decoys unchanged on all datasets and `max(A→B, B→A)` == old symmetric value.
- Regression harness scripts lived in the session scratchpad (not kept in repo).

---

## 10. Working notes for future sessions

- When adding a CLI feature, mirror it in the GUI `_load` worker and in `docs/usage.md`.
- Keep the `str | ResidueId` duck-typing; use `protein_of()`.
- Keep keys directional everywhere (sparse, block buckets, lookups, export).
- Anything touching `_blocks` must keep "block_size 1 ⇒ identical to per-item behaviour" and must not
  let blocks cross the target/decoy split or section boundaries.
- User values readability/simplicity; prove refactors with before/after output comparisons.
- Possible next topics (unconfirmed): reinstating a cluster ordering for directional data,
  `fromSite/ToSite` & `ProteinLinkPos*` column support, FASTA accession normalization, KEGG
  name-key fix (§8.2), thesis figures via Export PNG.

---

## 11. Code style (user request, 2026-10-05)

- Every source file is split into **parts with plain-language headers** (no numbering), e.g.
  `# ---- / # Building the matrix / # <one or two sentences on what this part does> / # ----`.
  Comments and docstrings are written in a natural, non-terse voice. Keep this style when adding code;
  put new code into the matching part.
- **Two deliberate typos — do NOT fix them:** `n_skiped` (local in `print_matrix._build`) and
  `postions` (local in `MatrixApp._pos_range`). The user asked for them.
- Behaviour was verified unchanged after the restructure (18 CLI goldens identical; GUI harness 484
  pygame + 308 numpy states identical; cutoff/sidebar checks pass). Line numbers elsewhere in this
  report are outdated since — grep for function names.
- 2026-10-10 harness (session scratchpad, not in repo): frame hashes of 324 view states × 2 datasets × 2
  renderers were identical before/after the render rework; behaviour checks for anchor/drag/selection passed.
  pygame mouse handling can be tested by `pygame.event.post` of MOUSEBUTTONDOWN/MOUSEMOTION/MOUSEBUTTONUP.
- Test-harness lesson: GUI tests that use `winfo_containing` or real input need
  `app.attributes("-topmost", True)`, and runs can be disturbed by the user's own mouse while the
  window is open (seen as a scrolled first dataset) — repeat a run before trusting a diff.

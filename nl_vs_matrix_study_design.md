# Node-Link vs. Adjacency Matrix for Crosslink Network Interpretation: Study Design

*Brainstorm document — comparison framework, task battery, and axis-ordering evaluation for the CLMS adjacency matrix viewer project*

## 1. Why this comparison is a real open question

Node-link diagrams and adjacency matrices are the two dominant ways of drawing a graph, and the visualization literature has spent two decades trying to settle when each one wins. The result is not a clean verdict but a set of task- and size-dependent trade-offs, which is exactly the shape your project's research question has — so the study design below leans on that literature rather than re-deriving it from scratch.

The foundational controlled experiment is Ghoniem, Fekete & Castagliola (2004, 2005), who compared the two representations across seven generic graph tasks and found a crossover effect: for graphs below roughly twenty vertices, node-link diagrams were read faster and more accurately across almost every task, but above that threshold — and increasingly so as density rose — the matrix overtook node-link on most tasks except path-following, which stayed a node-link strength at every size (Ghoniem et al., 2004, IEEE InfoVis; Ghoniem, Fekete & Castagliola, 2005, *Information Visualization* 4(2)). Keller, Eckert & Clarkson (2006, *Information Visualization* 5(1)) reached a similar conclusion working with engineering connectivity models rather than social networks, reinforcing that the effect isn't domain-specific.

That said, the picture is not settled. A much larger, more recent online experiment by Ren et al. (2019, *Network Science*) — using both 20- and 50-node networks — found node-link diagrams winning on both accuracy and speed overall, with the matrix's disadvantage shrinking (but not disappearing) as networks grew, and with participants' matrix performance improving substantially with practice while node-link performance stayed flat. This is worth citing precisely because it complicates the older "matrices win at scale" story: it suggests the *ordering* of the matrix and the *familiarity* of the reader matter as much as raw network size, which is directly relevant to your third research question below.

Two further threads matter for your specific project:

**Axis ordering is not a cosmetic choice — it is close to the whole game for a matrix.** Behrisch, Bach, Riche, Schreck & Fekete's EuroVis STAR (2016, *Computer Graphics Forum* 35(3)) is the standard survey of matrix reordering algorithms (spectral, hierarchical-clustering-based, TSP/Barycenter, Reverse Cuthill-McKee, and others) and makes the case, with examples, that a badly-ordered matrix hides exactly the block/cluster structure a well-ordered one reveals. Ren et al. (2019) independently confirmed the perceptual stakes: they found that even the *choice between two reasonable orderings* (degree-sorted vs. attribute-sorted) measurably biased participants' impressions of network size and structure. For the specific reordering method your viewer already defaults to, Bar-Joseph, Gifford & Jaakkola's (2001, *Bioinformatics* 17(Suppl. 1)) fast optimal-leaf-ordering algorithm is the standard reference — it was built for exactly this kind of "make blocks visible" problem in hierarchical-clustering heatmaps, which is structurally the same problem as a protein-crosslink adjacency matrix.

**Task taxonomy gives you a principled way to design the interpretation tasks**, rather than inventing them ad hoc. Lee, Plaisant, Parr, Fekete & Henry (2006, BELIV workshop, ACM) proposed a taxonomy of graph-reading tasks organized into topology-based (adjacency, accessibility, common connections, connectivity/clusters), attribute-based (on nodes, on links), browsing (path-following, revisiting), and overview tasks, and explicitly note — consistent with Ghoniem et al. — that matrices are stronger for cluster/neighborhood identification while node-link diagrams are stronger for path tasks. Section 4 below builds your task battery directly from this taxonomy, adapted to crosslink data.

Finally, two lines of work speak directly to your domain and to the "why not just use a hybrid" question a reviewer might raise: Henry, Fekete & McGuffin's NodeTrix (2007, IEEE TVCG 13(6)) shows that stitching matrix "blobs" for dense clusters into an otherwise node-link layout is a recognized third option, useful context for discussing why you are testing pure representations rather than jumping straight to a hybrid. And in the crosslinking-MS world specifically, xiNET (Combe, Fischer & Rappsilber, 2015, *Molecular & Cellular Proteomics* 14(4)) and its successor xiVIEW (2024, *Journal of Molecular Biology*) are the community-standard visualization tools, and both are node-link-first — which is itself a useful framing for your project: if matrices demonstrably surface neighbourhoods that the field's current node-link-based tools miss, that is a concrete, citable gap your viewer would be filling.

## 2. Research questions

The task description you were given maps onto three testable questions:

- **RQ1 (representation comparison).** For the same crosslink dataset, do a node-link diagram and an adjacency matrix support equally accurate and equally fast interpretation, or does one outperform the other, and does this depend on task type (per Lee et al.'s taxonomy) and on dataset properties (size, density, modularity)?
- **RQ2 (matrix-only patterns).** Are there structural patterns — specifically, neighbourhoods of densely interlinked proteins ("candidate sub-complexes") and bridge proteins connecting two neighbourhoods — that are reliably identified in the matrix view but missed, or identified more slowly and less confidently, in the node-link view?
- **RQ3 (axis ordering).** Within the matrix condition, does the algorithm used to order the protein axis (arbitrary/alphabetical vs. hierarchical-clustering + optimal leaf order vs. an alternative reordering algorithm) change how effectively participants can detect neighbourhoods and bridges?

## 3. Comparison framework: what grounds the two representations can be compared on

To compare "Node-Link output" and "Matrix output" fairly, three families of attributes need to be fixed, varied, and measured. This is the standard structure for a graph-visualization user study (see Ghoniem et al., 2005 and Lee et al., 2006 for the template this follows).

### 3.1 Data-level attributes (properties of the crosslink dataset — held constant or systematically varied across trials)

| Attribute | Why it matters | How to vary it |
|---|---|---|
| Network size (number of proteins) | Ghoniem et al. found a crossover near ~20 nodes; Ren et al. tested 20 vs. 50 | Include a small (<20), medium (20–50) and large (>50 protein) dataset |
| Density (crosslinks / possible pairs) | Both Ghoniem and Behrisch et al. tie matrix advantage to density | Include a sparse and a dense dataset at comparable size |
| Modularity / community structure | This is what "neighbourhoods" literally means — need datasets with real or planted sub-complexes | Prefer datasets from known multi-subunit complexes (e.g., a ribosomal or proteasome CLMS dataset) where sub-complex membership is externally known, so "correct" neighbourhoods can be scored objectively |
| Degree distribution (hub proteins) | Hub proteins (e.g., scaffolding subunits) read very differently in each representation | Include at least one dataset with one or two clear high-degree hubs |
| Link confidence / support (spectral counts per crosslink) | Lets you design attribute-based tasks distinct from pure topology tasks | Use the score/PSM-count field already in your CLMS data |
| Self-links vs. inter-protein links | CLMS-specific distinction with no node-link/matrix literature equivalent — worth calling out explicitly | Keep both types present; note how each representation encodes a self-link (matrix diagonal cell vs. node-link self-loop) |

For datasets where no ground-truth sub-complex structure is available, consider generating one synthetic benchmark network with **planted** community structure using the Lancichinetti–Fortunato–Radicchi (LFR) benchmark generator (Lancichinetti, Fortunato & Radicchi, 2008, *Physical Review E* 78, 046110), relabeling generic nodes as "proteins." This buys you an objective accuracy metric for RQ2/RQ3 (precision/recall of identified neighbourhoods against the planted ones) that a real, ground-truth-free dataset cannot give you.

### 3.2 Representation-level attributes (the independent variable you are actually testing)

- **Representation type** (2 levels): node-link diagram (force-directed layout, as currently produced by standard tools like xiNET/xiVIEW) vs. adjacency matrix (your PyQtGraph viewer).
- **Matrix axis ordering** (only applies within the matrix arm; this is RQ3): (a) arbitrary/alphabetical order as a null baseline; (b) hierarchical clustering with optimal leaf ordering (Bar-Joseph et al., 2001) — your current default; (c) one alternative reordering algorithm from a different family, e.g. spectral ordering or Reverse Cuthill-McKee (both surveyed in Behrisch et al., 2016), to check whether HC+OLO's advantage (if any) is specific to that algorithm or general to "any good reordering"; (d) where ground truth exists, a domain/ground-truth ordering (proteins grouped by known complex membership) as a ceiling/upper-bound condition.
- **Visual encoding held constant** across conditions as confounds to control for: colour scale for link confidence, node/row labelling scheme, canvas size, and interaction affordances (zoom/pan/hover-highlight) available in both views, so that differences are attributable to layout logic rather than styling.

### 3.3 Task-level and outcome attributes (what you measure)

Outcomes should be collected per task, per representation/ordering condition, per dataset:

- **Accuracy** — correctness against ground truth (planted or literature-derived complex membership) or against a consensus/expert answer where no ground truth exists.
- **Response time** — time to first (or final) answer.
- **Confidence** — self-rated (e.g., 5-point Likert) confidence in the answer, which lets you distinguish "wrong but sure" errors from appropriately hedged ones.
- **Perceived workload** — a short NASA-TLX (Hart & Staveland, 1988) after each representation block, to capture effort/frustration differences that speed and accuracy alone can miss.
- **Preference** — forced-choice or ranking of representations per dataset size/density regime, since Ghoniem's and Ren's results suggest preference may not track objective performance.
- **Error type** — classify errors using Lee et al.'s (2006) task categories (e.g., an adjacency error vs. a connectivity/cluster error), so you can say *which kind* of reasoning breaks down in which representation, not just that accuracy differs.
- **Ordering-specific outcome (RQ3 only)** — for the neighbourhood-detection task specifically, precision and recall of identified neighbourhoods against the planted/known ones, and "time to first correctly identified neighbourhood," compared across ordering conditions.

## 4. Task battery for interpreting the crosslink network

Tasks are grouped using Lee et al.'s (2006) taxonomy and phrased for a crosslink-network context. For each, a literature-grounded hypothesis about which representation should win is given — these are hypotheses to test, not assumptions to bake in.

**Topology-based tasks**

1. *Adjacency* — "Which protein has the greatest number of distinct crosslink partners?" Hypothesis: matrix wins as network size grows (row/column-sum read is a single visual scan; Ghoniem et al., 2005).
2. *Common connections* — "List every pair of proteins that share at least three crosslink partners in common." A direct test of matrix block-reading; hypothesis: matrix wins, and more so under a good ordering (RQ3).
3. *Connectivity / shortest path* — "What is the shortest crosslink chain connecting Protein A to Protein Z?" Hypothesis: node-link wins at every size (the one consistent finding across Ghoniem, Keller, and Ren).
4. **Neighbourhood / cluster identification (core RQ2 task)** — "Identify groups of proteins that are densely crosslinked to each other but only sparsely linked to the rest of the network (candidate sub-complexes)." Run this task once per matrix ordering condition, and once in node-link, to answer both RQ2 and RQ3 directly.
5. *Bridge identification* — "Identify any protein(s) that connect two otherwise separate neighbourhoods." Tests whether a structurally important but low-degree "articulation" protein is more visible as an off-diagonal block in a well-ordered matrix than as a single edge in a dense node-link tangle.
6. *Outlier / isolate detection* — "Identify proteins with only one crosslink (weakly connected proteins)."

**Attribute-based tasks**

7. *On links* — "Identify the crosslinked pair(s) with the highest spectral/PSM support (highest-confidence link)."
8. *On nodes* — "Identify which proteins have only self-links (intra-protein crosslinks) versus proteins with only inter-protein links."

**Browsing and overview tasks**

9. *Overview / density estimation* — "Would you describe this crosslink network as sparse, moderately connected, or dense? Is it a single connected mass or several separate modules?" A fast, low-effort holistic judgement, useful as a first pass before the more effortful tasks above.
10. *Follow-path / revisit* (node-link-specific browsing task, included for completeness per Lee et al.'s taxonomy, though it has no natural matrix equivalent) — "Starting at Protein A, trace the crosslink path to Protein D, then confirm whether the same path can be retraced without error." Useful primarily as a manipulation check that participants can use the node-link tool correctly, and as an illustration of the one task type where node-link's advantage should be least contested.

## 5. RQ3 sub-study: isolating the effect of axis ordering

Because ordering is the more novel and more actionable part of your question (it is a design choice inside a tool you control, not just a fixed property of a representation), it deserves its own tightly controlled repeated-measures block, separate from the broader representation comparison:

- Take the same dataset(s) with known or planted neighbourhood structure and render the matrix under the ordering conditions in §3.2 (arbitrary, HC+optimal-leaf-order, one alternative algorithm, and — where available — ground-truth ordering).
- Within each ordering, run task 4 (neighbourhood identification) and task 5 (bridge identification) only — these are the two tasks whose whole premise is that blocks/patterns must be visually contiguous to be seen, which is exactly what reordering controls.
- Counterbalance the order in which participants see the four orderings (a Latin square), since Ren et al. (2019) showed matrix-reading performance improves with practice within a session — without counterbalancing, a later-seen ordering could look artificially better simply because the participant has had more practice reading matrices, not because that ordering is better.
- Report both accuracy/recall against the planted ground truth and participants' self-reported confidence, since Behrisch et al.'s survey emphasizes that a bad ordering doesn't just slow people down — it can make them confidently identify the *wrong* blocks as neighbourhoods.

## 6. Practical notes for implementing this in the viewer

A few implementation-facing suggestions that follow directly from the design above: expose the ordering algorithm as a switchable setting in the GUI (arbitrary / HC+OLO / alternative) rather than hard-coding HC+OLO, since the comparison in §5 requires switching between them on the same dataset; log task start/answer/confidence timestamps directly from the viewer if you run the study inside your own tool, which gives you response-time data for free; and consider building or sourcing one dataset with genuinely known sub-complex ground truth (e.g., a well-characterized multi-subunit complex) specifically so that RQ2/RQ3 accuracy can be scored objectively rather than relying only on inter-rater agreement.

## 7. Limitations to flag up front

Learning effects across representations are a known confound (Ren et al., 2019) and should be handled with counterbalancing and/or a between-subjects first exposure. "Correctness" for neighbourhood-detection tasks is only as good as your ground truth — real CLMS data rarely comes with a certified list of sub-complexes, which is the main argument for including at least one synthetic (LFR-based) dataset alongside real ones. Finally, most of the cited node-link/matrix literature used generic or social-network graphs, not protein interaction data specifically — part of the contribution of this project is testing whether their conclusions (matrix advantage at scale/density, node-link advantage for paths) actually transfer to crosslink networks, which have their own quirks (self-links, confidence-weighted edges, sparse-but-highly-structured sub-complexes) that generic benchmark graphs don't have.

## References

1. Ghoniem, M., Fekete, J.-D., & Castagliola, P. (2004). A comparison of the readability of graphs using node-link and matrix-based representations. *Proceedings of the IEEE Symposium on Information Visualization (InfoVis 2004)*. https://dl.acm.org/doi/10.5555/1038262.1038777
2. Ghoniem, M., Fekete, J.-D., & Castagliola, P. (2005). On the readability of graphs using node-link and matrix-based representations: A controlled experiment and statistical analysis. *Information Visualization*, 4(2), 114–135. https://doi.org/10.1057/palgrave.ivs.9500092
3. Keller, R., Eckert, C. M., & Clarkson, P. J. (2006). Matrices or node-link diagrams: Which visual representation is better for visualising connectivity models? *Information Visualization*, 5(1), 62–76. https://doi.org/10.1057/palgrave.ivs.9500116
4. Ren, D., Marusich, L. R., O'Donovan, J., Bakdash, J. Z., Schaffer, J., Cassenti, D. N., Kase, S. E., Roy, H., Lin, W. S., & Höllerer, T. (2019). Understanding node-link and matrix visualizations of networks: A large-scale online experiment. *Network Science*. https://donghaoren.org/publications/netsci19-nodelink-matrix.pdf
5. Lee, B., Plaisant, C., Parr, C. S., Fekete, J.-D., & Henry, N. (2006). Task taxonomy for graph visualization. *Proceedings of the 2006 AVI Workshop on BEyond Time and Errors: Novel Evaluation Methods for Information Visualization (BELIV)*. https://doi.org/10.1145/1168149.1168168
6. Behrisch, M., Bach, B., Riche, N. H., Schreck, T., & Fekete, J.-D. (2016). Matrix reordering methods for table and network visualization. *Computer Graphics Forum*, 35(3), 693–716. https://doi.org/10.1111/cgf.12935
7. Henry, N., Fekete, J.-D., & McGuffin, M. J. (2007). NodeTrix: A hybrid visualization of social networks. *IEEE Transactions on Visualization and Computer Graphics*, 13(6), 1302–1309. https://doi.org/10.1109/TVCG.2007.70582
8. Saket, B., Simonetto, P., Kobourov, S., & Börner, K. (2014). Node, node-link, and node-link-group diagrams: An evaluation. *IEEE Transactions on Visualization and Computer Graphics*, 20(12), 2231–2240. https://arxiv.org/abs/1404.1911
9. Bar-Joseph, Z., Gifford, D. K., & Jaakkola, T. S. (2001). Fast optimal leaf ordering for hierarchical clustering. *Bioinformatics*, 17(Suppl. 1), S22–S29. https://doi.org/10.1093/bioinformatics/17.suppl_1.S22
10. Lancichinetti, A., Fortunato, S., & Radicchi, F. (2008). Benchmark graphs for testing community detection algorithms. *Physical Review E*, 78, 046110. https://arxiv.org/abs/0805.4770
11. Combe, C. W., Fischer, L., & Rappsilber, J. (2015). xiNET: Cross-link network maps with residue resolution. *Molecular & Cellular Proteomics*, 14(4), 1137–1147. https://doi.org/10.1074/mcp.O114.042259
12. xiVIEW: Visualisation of crosslinking mass spectrometry data (2024). *Journal of Molecular Biology*. https://www.sciencedirect.com/science/article/pii/S0022283624002511
13. Hart, S. G., & Staveland, L. E. (1988). Development of NASA-TLX (Task Load Index): Results of empirical and theoretical research. In P. A. Hancock & N. Meshkati (Eds.), *Human Mental Workload* (pp. 139–183). North-Holland.

<!-- Slide number: 1 -->

TX-OMICS REVISION
Seven gene sets.
One auditable curation.

7
295
stocker CSVs
unique FBgns

4
3
A concise record of the publication rules, pathway selection, identity checks, and overlap that produced the stocker-ready inputs.
publication sets
pathway sets
Rosensweig–Shah 2026  ·  curation summary

### Notes:

<!-- Slide number: 2 -->
02
LIBRARY MAP
The final library at a glance
Four evidence-defined sets and three pathway-hit sets; 295 unique genes after FBgn deduplication.

01
02
03
04
6
20
97
4
genes
genes
genes
genes
Mechanistic
Homeostatic H×R
CSW 4+
HLH regulators
Publication core + screen ranks 1–6
15 History+/Rebound− + 5 opposite
Wake frequency ≥4, threshold-scoped
Four bHLH factors named in paper

05
06
07
99
167
5
genes
genes
genes
Ribosomal
Mitochondrial
Immune
DAVID term hits; approved keywords
DAVID mito/metabolism term hits
DAVID immune/Toll/Imd term hits
Rosensweig–Shah 2026  ·  flybase_gene_id is the membership key

### Notes:

<!-- Slide number: 3 -->
03
CURATION PIPELINE
Two curation tracks converge on one controlled folder
The shared gates are identity resolution, FBgn deduplication, and a clean seven-file handoff.

Publication categories
TRACK A

Paper + breakdown tables
Apply set-specific rules
Resolve FlyBase identity
Publish sets 01–04

1
2
3
4

Pathway hits
TRACK B

Four CSW master tables
DAVID · FDR ≤10%
Review term keywords
Publish sets 05–07

1
2
3
4

902
791
≥2
7
source observations
unique source FBgns
minimum pathway frequency
experiments
Thresholds remain separate: FC0.5 and FC1 are never collapsed into a max-frequency field.

### Notes:

<!-- Slide number: 4 -->
04
SETS 01–02
Publication sets use explicit, reproducible rules
Sets 01–02 preserve publication evidence and resolve one biologically important identity conflict.

SET 01  ·  6 GENES
SET 02  ·  20 GENES
Mechanistic
Homeostatic history × rebound
Publication core joined by FBgn to mechanistic screen ranks 1–6.
15
5
History+ / Rebound−
History− / Rebound+
Carry forward mechanism, literature status, Tx-Omics evidence, and fl.ai confidence.
Evidence labels: consistent, correlated, or curated.
No rank or annotation is inferred outside the source tables.

Identity correction
Source row trh / FBgn0262139 is resolved from publication context to Trhn / FBgn0035187 (tryptophan hydroxylase). The wrong FBgn is never retained.
unc79  ·  SIFa  ·  rumpel  ·  AstA-R2  ·  Trhn  ·  RpL23
Rosensweig–Shah 2026  ·  flybase_gene_id is the membership key

### Notes:

<!-- Slide number: 5 -->
05
SETS 03–04
Screen frequency and direct-regulator evidence stay separate
Set 03 is threshold-scoped screen evidence; Set 04 comes directly from the paper.

SET 03  ·  97 GENES
SET 04  ·  4 GENES
CSW 4+ wake-frequency union
HLH regulators
Named in Results / Table S1 after HOMER identified a CAGCTG E-box upstream of wake-induced genes.
97
7
0
FC0.5 wake frequency 4–6
also qualify at FC1
sleep genes at frequency 4
bigmax
FBgn0039509
Union rule
HLH3B
FBgn0011276
Take FC0.5 wake genes at frequency 4, 5, or 6.
Add FC1 wake frequency-4 genes; all seven are already in the FC0.5 set.
Keep FC0.5 and FC1 metadata in separate columns.
E(spl)m3-HLH
FBgn0002609
E(spl)mbeta-HLH
FBgn0002733
Rosensweig–Shah 2026  ·  flybase_gene_id is the membership key

### Notes:

<!-- Slide number: 6 -->
06
SETS 05–07
Pathway sets contain DAVID hits—not every annotated gene
Only hit genes from approved FDR-passing terms enter sets 05–07.

Submit four CSW masters
Keep FDR ≤10% terms
Review term keywords
Take hit FBgns
1

2

3

4

![Venn diagram of pathway gene-set overlap](Image0.jpg)
Why the overlap is intentional
167
Mitochondrial / metabolism

Mitochondrial translation and mitochondrial ribosomal-subunit terms match both keyword groups. Their hit genes are retained in every matched bucket.
99
Ribosomal / translation

5
Immune

31
ribo ∩ mito
All three lists are wake-only in this run.
Rosensweig–Shah 2026  ·  flybase_gene_id is the membership key

### Notes:

<!-- Slide number: 7 -->
07
FLYBASE AUDIT
Identity resolution is conservative by design
Current symbols can pass directly; ambiguous synonyms require biological context or a case-supported source match.

Trhn
Current symbol · tryptophan hydroxylase neuronal · publication-resolved in Homeostatic.
FBgn0035187
Approval gates

PASS
Exact current symbol
Publication-supported correction
Unique case-insensitive current-symbol match
Otherwise: stop for review

trh
Current symbol · trachealess · must not enter any generated identity field for this publication row.
FBgn0262139

REJECT

Trh
Synonym of Hn, Trhn, and trh · never auto-approved without biological context.
ambiguous

REVIEW
Membership key
flybase_gene_id

Symbols are display labels. All overlap, deduplication, and set membership use FBgn only.
Rosensweig–Shah 2026  ·  flybase_gene_id is the membership key

### Notes:

<!-- Slide number: 8 -->
08
OVERLAP  ·  295 UNIQUE FBGNS
Most genes are set-specific; overlap concentrates in CSW and pathways
Exact intersections expose the reusable biological core without inflating membership.

Largest exact intersection blocks
What matters
102
Mitochondrial only

200
genes in exactly one set
40
Ribosomal only

87
genes in exactly two sets
34
CSW 4+ ∩ mitochondrial

8
genes in exactly three sets
28
CSW 4+ only

24
Ribosomal ∩ mitochondrial only

Largest pairwise overlaps

41  CSW 4+ ∩ mito
35  CSW 4+ ∩ ribo
31  ribo ∩ mito
The immune and HLH sets remain fully distinct from the pathway core.
Rosensweig–Shah 2026  ·  flybase_gene_id is the membership key

### Notes:

<!-- Slide number: 9 -->
STOCKER HANDOFF
The contract is simple—and strict.

Seven inputs
01  Mechanistic

Exactly seven CSVs in the discovery folder
Required keys: ext_gene and flybase_gene_id
No overlap, GO, or audit CSVs beside the inputs
7
2
0
02  Homeostatic H×R

03  CSW 4+

04  HLH regulators

05  Ribosomal

Run
06  Mitochondrial

07  Immune

python -m fl_ai_reagent_stocker run ./data/gene_sets/Tx-Omics_Revision/ \
  --config ./data/config/stock_split_config_priority_example.json
Discovery is recursive: every extra *.csv changes the run. Keep analysis outputs elsewhere.
Stability rule
Rosensweig–Shah 2026  ·  Tx-Omics Revision

### Notes:

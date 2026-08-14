<!-- Slide number: 1 -->

TX-OMICS REVISION
Seven gene sets.
One auditable curation.

7
2
stocker CSVs
analysis families

4
3
A concise record of the publication rules, pathway selection, identity checks, and overlap that produced the stocker-ready inputs.
primary sets
GO term sets
Rosensweig–Shah 2026  ·  curation summary

### Notes:

<!-- Slide number: 2 -->
02
LIBRARY MAP
Two curation families—analyzed separately
Primary evidence sets and CSW-derived GO term–grouped sets use separate overlap analyses.
FAMILY 1  ·  PRIMARY CURATION SETS

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
CSW 4+ (FC0.5)
HLH regulators
Publication core + screen ranks 1–6
15 History+/Rebound− + 5 opposite
FC0.5 wake frequency ≥4
Four bHLH factors named in paper
FAMILY 2  ·  CSW-DERIVED GO TERM–GROUPED SETS (FC0.5 / FC1)

05
06
07
RULE
99
167
5
Keep them separate
genes
genes
genes
Ribosomal
Mitochondrial
Immune
Report overlap within each family—never as one seven-set comparison.
DAVID term hits; approved keywords
DAVID mito/metabolism term hits
DAVID immune / Toll / Imd term hits
Rosensweig–Shah 2026  ·  flybase_gene_id is the membership key

### Notes:

<!-- Slide number: 3 -->
03
CURATION PIPELINE
Two curation tracks converge on one controlled folder
They share identity gates and file standards—but their overlap analyses remain separate.

Primary curation sets
FAMILY 1

Paper + screen evidence
Apply set-specific rules
Resolve FlyBase identity
Publish sets 01–04

1
2
3
4

GO term–grouped sets
FAMILY 2

CSW FC0.5 + FC1 tables
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
FAMILY 2  ·  THRESHOLD PROVENANCE
CSW-derived GO term–grouped genes retain threshold provenance
Each gene records its source as FC0.5, FC1, or both; only hits from approved FDR-passing terms enter sets 05–07.

Submit four CSW masters
Keep FDR ≤10% terms
Review term keywords
Take hit FBgns
1

2

3

4

GO term–grouped outputs
Threshold provenance · 240 unique FBgns
167
Mitochondrial / metabolism
203
FC0.5 only
Ribosomal and mitochondrial dominate

99
Ribosomal / translation
5
FC1 only
BomBc2 · BomS1 · BomS2 · IM14 · Thor

5
Immune
32
both thresholds
FC1 ribosomal/mitochondrial genes are also at FC0.5

All three outputs are wake-direction genes.
Only the five immune genes are FC1-exclusive.
Key takeaway
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
FAMILY 1  ·  PRIMARY-SET OVERLAP
Primary-set overlap is limited to three genes
Mechanistic, Homeostatic, CSW 4+ (FC0.5), and HLH are compared only with one another.

Within Family 1
6
20
97
4
Mechanistic
Homeostatic
CSW 4+ · FC0.5
HLH
124
unique FBgns
Observed shared genes
121
in exactly one set
2

Mechanistic  ×  Homeostatic
Trhn  ·  AstA-R2
3
in exactly two sets
1

Mechanistic  ×  CSW 4+ (FC0.5)
RpL23
0
in three or four sets
0

All other pairings
No shared FBgns

HLH is fully disjoint from the other primary sets.
Rosensweig–Shah 2026  ·  flybase_gene_id is the membership key

### Notes:

<!-- Slide number: 9 -->
09
FAMILY 2  ·  GO OVERLAP
GO term–grouped overlap is confined to ribosomal × mitochondrial
The three pathway buckets are compared only within Family 2; threshold provenance remains a separate attribute.

GO term–group membership
What the overlap means
Schematic · counts shown

31
genes occur in both ribosomal and mitochondrial buckets

5
68
31
136
5
immune genes are fully disjoint
Immune
5
0
immune overlaps with either other GO bucket

Ribosomal
99
Mitochondrial
167
The 31-gene overlap is intentional: mitochondrial translation terms match both keyword groups.
240 unique FBgns in Family 2
No cross-family overlap is shown: Family 1 and Family 2 answer different curation questions.

### Notes:

<!-- Slide number: 10 -->
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

RUN

06  Mitochondrial

python -m fl_ai_reagent_stocker run ./data/gene_sets/Tx-Omics_Revision/ \
  --config ./data/config/stock_split_config_priority_example.json
07  Immune

Discovery is recursive: every extra *.csv changes the run. Keep analysis outputs elsewhere.
Stability rule
Rosensweig–Shah 2026  ·  Tx-Omics Revision

### Notes:

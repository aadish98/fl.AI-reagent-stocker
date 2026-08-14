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
Homeostatic genes
CSW 4+ genes
HLH genes
Consistent/correlated; published effect fits homeostatic feedback
Opposing history ↔ rebound correlations: 15 + 5
Consistent in ≥4 of 7 datasets
Potential upstream regulators
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

Paper-defined evidence
Apply four distinct rules
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
FAMILY 1  ·  SETS 01–02
Mechanistic and Homeostatic are different evidence categories
The six-gene functional set is distinct from the 20-gene opposing history–rebound correlation set.

SET 01  ·  6 GENES
SET 02  ·  20 GENES
Mechanistic
Homeostatic genes
Consistent or correlated sleep/wake genes with published effects in a direction consistent with potential homeostatic function.
15
5
History+ / Rebound−
History− / Rebound+
Expression is consistent across conditions or correlated with sleep/wake history.
Published perturbation affects sleep/wake in a direction compatible with negative feedback.
The paper highlights six such genes.

Genes correlate with sleep/wake history before sampling and with rebound in the opposite direction afterward.
Identity correction
Source row trh / FBgn0262139 is resolved from publication context to Trhn / FBgn0035187 (tryptophan hydroxylase). The wrong FBgn is never retained.
unc79  ·  SIFa  ·  rumpel  ·  AstA-R2  ·  Trhn  ·  RpL23
Rosensweig–Shah 2026  ·  flybase_gene_id is the membership key

### Notes:

<!-- Slide number: 5 -->
05
FAMILY 1  ·  SETS 03–04
CSW consistency and HLH regulation are separate categories
Set 03 measures recurrence across seven datasets; Set 04 identifies potential upstream regulators.

SET 03  ·  97 GENES
SET 04  ·  4 GENES
CSW 4+ genes
HLH genes
Named in Results / Table S1 after HOMER identified a CAGCTG E-box upstream of wake-induced genes.
97
7
0
Current output · FC0.5 wake frequency 4–6
also qualify at FC1
sleep genes at frequency 4
bigmax
FBgn0039509
Category rule
HLH3B
FBgn0011276
Consistent sleep/wake genes evident across four or more of seven datasets.
The current output contains 97 FC0.5 wake genes; seven also qualify at FC1.
No sleep-direction gene reaches frequency 4 in this run.
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
Mechanistic, Homeostatic genes, CSW 4+ genes, and HLH genes are compared only with one another.

Within Family 1
6
20
97
4
Mechanistic
Homeostatic genes
CSW 4+ genes
HLH genes
124
unique FBgns
Observed shared genes
121
in exactly one set
2

Mechanistic  ×  Homeostatic genes
Trhn  ·  AstA-R2
3
in exactly two sets
1

Mechanistic  ×  CSW 4+ genes
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
02  Homeostatic genes

03  CSW 4+ genes

04  HLH genes

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

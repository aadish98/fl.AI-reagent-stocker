<!-- Slide number: 1 -->

TX-OMICS REVISION
How each gene set was curated
Seven stocker-ready CSVs from Rosensweig–Shah 2026, with publication rules, DAVID pathway hits, FlyBase identity audit, and FBgn-only overlap.

4
3
295
7
publication categories
pathway-hit lists
unique FBgns
stocker CSVs
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
1  /  14

### Notes:

<!-- Slide number: 2 -->

OVERVIEW
The seven stocker inputs

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
Homeostatic HxR
CSW 4+
HLH regulators
Publication genes with in-vivo homeostatic support
15 Hist+/Reb−  +  5 Hist−/Reb+
Wake genes at frequency ≥4 (FC0.5)
Named bHLH factors from Results / Table S1

05
06
07
Σ
99
167
5
295
genes
genes
genes
genes
Ribosomal
Mitochondrial
Immune
Seven-set union
CSW genes hit for ribosome/translation GO terms
CSW genes for mito/metabolism GO terms
CSW genes hit for immune GO terms
Unique FBgns
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
2  /  14

### Notes:

<!-- Slide number: 3 -->

SCHEMA
Curation pipeline — two tracks, one stocker folder
TRACK A  ·  PUBLICATION CATEGORIES

A1
A2
A3
A4
Publication +
breakdown tables
Four category
CSVs
FlyBase identity
audit
Promote four
category files

TRACK B  ·  PATHWAY HITS

B1
B2
B3
B4
Four CSW master
tables (freq ≥ 2)
DAVID via
GO_Analysis
Keyword term
review
Three pathway-hit
CSVs

Promote only after both audits are clean
Copy exactly seven CSVs — and no analysis files — into data/gene_sets/Tx-Omics_Revision/. Stocker recursively discovers every *.csv, so overlap tables, GO workbooks, and this deck must stay out of that discovery set (the .pptx is safe; extra CSVs are not).
Membership key for overlap is flybase_gene_id only. Symbols are display labels.
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
3  /  14

### Notes:

<!-- Slide number: 4 -->

SOURCE UNIVERSE
CSW observations before any category cut

647
98
141
16
Wake FC0.5
Wake FC1
Sleep FC0.5
Sleep FC1

902
791
≥2
7
observation rows
unique FBgns
minimum frequency
experiments
keyed by FBgn × threshold × direction
never collapsed across thresholds
every master-table row
Baseline through THIP
Masters live in private/lab-materials/Tx-Omics_FollowUp/Data/Transcriptomics_Input/. Frequency, correlation, and cycling stay per threshold — never a max_frequency field.
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
4  /  14

### Notes:

<!-- Slide number: 5 -->

SET 01  ·  n = 6
Mechanistic: publication core joined to screen ranks 1–6

Curation rule
CSV schema
Start from Tx-Omics_homeostatic_genes_n=6genes.csv. Join by FBgn to ranks 1–6 of Tx-Omics_Mechanistic_Screen_Candidates_nsyb_elav.csv so mechanism, literature, and fl.ai fields survive.
| ext\_gene | required stocker symbol |
| --- | --- |
| flybase\_gene\_id | required stocker FBgn |
| mechanistic\_evidence\_type | consistent / correlated / curated |
| Mechanistic Category | from screen ranks 1–6 |
| Literature Status | published vs transcriptomic |
| Tx-Omics Evidence | Rosensweig–Shah 2026 |
| Proposed Mechanism | screen text |
| fl.ai Category / Confidence | joined, not invented |
| identity\_status | source |

unc79
NALCN auxiliary
consistent

SIFa
neuropeptide
consistent

rumpel
glial transport
correlated

AstA-R2
AstA receptor
correlated

Trhn
tryptophan hydroxylase
curated

RpL23
ribosome / wake translation
consistent
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
5  /  14

### Notes:

<!-- Slide number: 6 -->

SET 02  ·  n = 20
Homeostatic history × rebound, with one publication ID fix

15
5
Trhn, not trachealess
Source row trh / FBgn0262139 is publication-resolved to Trhn / FBgn0035187 (tryptophan hydroxylase). Results, Discussion, and Figure 14. The generated CSV never keeps the wrong FBgn.
History+  /  Rebound−
History−  /  Rebound+
Sleep-history positive and rebound negative — the paper’s homeostatic signature.
The opposite quadrant. Union is 20 unique FBgns.

Schema  ·  identity columns only; literature annotations follow. No history/rebound coefficients exist in-repo and none were fabricated.
ext_gene  ·  flybase_gene_id  ·  history_direction  ·  rebound_direction  ·  source_file  ·  identity_status  ·  identity_resolution_basis=RosensweigShah_2026_results_discussion_figure14
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
6  /  14

### Notes:

<!-- Slide number: 7 -->

SET 03  ·  n = 97
CSW 4+: wake frequency union, never a collapsed max

Union rule
97
1    FC0.5 wake frequency 4 + 5 + 6  →  97 unique genes
2    FC1 wake frequency 4  →  7 genes, all already inside the FC0.5 4+ set
wake genes
7 also qualify at FC1
0 sleep genes at freq 4
3    Union remains 97. Sleep never reaches frequency 4.
4    Keep frequency_FC0.5_wake and frequency_FC1_wake as separate fields.

Schema  ·  threshold-scoped metadata
ext_gene, flybase_gene_id, frequency_FC0.5_wake, frequency_FC1_wake, qualifies_FC0.5_4plus, qualifies_FC1_4plus, qualifying_thresholds, source_files, then per-threshold cycling, wake_corr_exps, and the seven experiment columns. No unsuffixed frequency or wake_corr_exps.
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
7  /  14

### Notes:

<!-- Slide number: 8 -->

SET 04  ·  n = 4
HLH regulators: named in the paper, not taken from CSW tables

bigmax
HLH3B
E(spl)m3-HLH
E(spl)mbeta-HLH
FBgn0039509
FBgn0011276
FBgn0002609
FBgn0002733

Membership evidence
Schema  ·  publication only
| ext\_gene / flybase\_gene\_id | stocker keys |
| --- | --- |
| publication\_section | Results |
| publication\_table | Table S1 |
| motif | CAGCTG E-box |
| publication\_evidence | direct-regulator statement |
| identity\_status | publication\_named |
HOMER found the CAGCTG E-box upstream of wake-induced genes. Results and Table S1 name these four bHLH factors as potential direct upstream regulators. FlyBase IDs are resolved from those symbols on the 2026_01 synonym table.
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
8  /  14

### Notes:

<!-- Slide number: 9 -->

SETS 05–07
Pathway genes are DAVID hits, not all annotated genes

1

2

3

4
Submit four CSW masters
Keep FDR ≤ 10% terms
Keyword-match term names
Take hit FBgns in those terms
GO_Analysis fetch_go_report
species 7227  ·  EASE 0.1
Sleep FC1: no passing terms
(explicit outcome)
ribosom / mitochondri /
immune, Toll, Imd…
One row per unique gene
per approved bucket

Conflict policy
Shared pathway schema
Mitochondrial translation and the mitochondrial large/small ribosomal subunits match both ribosomal and mitochondrial keywords. They were assigned to every matched bucket (include_all_matched_buckets). That is why 31 genes sit in both pathway lists.
ext_gene, source_flybase_gene_id, flybase_gene_id, pathway_bucket, source_tables, directions, thresholds, go_term_ids / names / FDR min, per-threshold frequency and corr, in_CSW_4plus, in_Mechanistic_6, in_Homeostatic_20, in_HLH_4.
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
9  /  14

### Notes:

<!-- Slide number: 10 -->

PATHWAY SETS
All three pathway lists are wake-only in this run

### Chart

| Category | Unique genes |
|---|---|
| Ribosomal /
translation | 99.0 |
| Mitochondrial /
metabolism | 167.0 |
| Immune | 5.0 |
![/Users/aadishms/Desktop/Projects.nosync/fl.AI-reagent-stocker/audit_outputs/Tx-Omics_Revision/Figures/pathway_overlap_venn.png](Image0.jpg)
Ribosomal 99  ·  Mitochondrial 167  ·  Immune 5. Intersection is mitochondrial ribosomal proteins from the three conflict terms. Immune (BomBc2, BomS1, BomS2, IM14, Thor) is disjoint.
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
10  /  14

### Notes:

<!-- Slide number: 11 -->

FLYBASE AUDIT
Identity rules that were allowed to pass

Trhn
FBgn0035187
Tryptophan hydroxylase neuronal. Exact current symbol. Publication-resolved in Homeostatic.

trh
FBgn0262139
trachealess. Exact current symbol. Must not appear in any generated identity field.

Trh
ambiguous
Synonym of Hn, Trhn, and trh. Never auto-approved without biological context.

RpL37a / arg
kept source IDs
Synonym-ambiguous, but exactly one candidate matches the current symbol case-insensitively (RpL37A, Arg).
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
11  /  14

### Notes:

<!-- Slide number: 12 -->

OVERLAP  ·  295 UNIQUE FBGNS
Pairwise overlapping-gene counts
|  | Mech | Homeo | CSW 4+ | HLH | Ribo | Mito | Immune |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Mech | 6 | 2 | 1 | 0 | 1 | 0 | 0 |
| Homeo | 2 | 20 | 0 | 0 | 0 | 0 | 0 |
| CSW 4+ | 1 | 0 | 97 | 0 | 35 | 41 | 0 |
| HLH | 0 | 0 | 0 | 4 | 0 | 0 | 0 |
| Ribo | 1 | 0 | 35 | 0 | 99 | 31 | 0 |
| Mito | 0 | 0 | 41 | 0 | 31 | 167 | 0 |
| Immune | 0 | 0 | 0 | 0 | 0 | 0 | 5 |

How to read it
Diagonal is set size. Off-diagonal is shared FBgns. Dark terracotta = largest overlaps.
200
in exactly one set
87
in exactly two sets
8
in exactly three sets
41
CSW 4+ ∩ mito
35
CSW 4+ ∩ ribo
31
ribo ∩ mito
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
12  /  14

### Notes:

<!-- Slide number: 13 -->

EXACT INTERSECTIONS
UpSet of the seven named sets

![/Users/aadishms/Desktop/Projects.nosync/fl.AI-reagent-stocker/audit_outputs/Tx-Omics_Revision/Figures/seven_set_upset.png](Image0.jpg)

Largest exclusive blocks
102
mitochondrial only
40
ribosomal only
34
CSW 4+ ∩ mito
28
CSW 4+ only
24
ribo ∩ mito only
5
immune only
4
HLH only
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
13  /  14

### Notes:

<!-- Slide number: 14 -->

WHAT CROSSES SETS
High-priority genes and the stocker contract
| Gene | FBgn | Sets |
| --- | --- | --- |
| Trhn | FBgn0035187 | Mechanistic · Homeostatic |
| AstA-R2 | FBgn0039595 | Mechanistic · Homeostatic |
| RpL23 | FBgn0010078 | Mechanistic · CSW 4+ · Ribosomal |
| unc79, SIFa, rumpel | FBgn0038693, FBgn0053527, FBgn0029950 | Mechanistic only |
| bigmax, HLH3B, E(spl)m3-HLH, E(spl)mbeta-HLH | four publication FBgns | HLH only |
| BomBc2, BomS1, BomS2, IM14, Thor | five immune FBgns | Immune only |

Stocker contract
Exactly seven CSVs. Every file must have ext_gene and flybase_gene_id. No analysis CSVs in this folder.
python -m fl_ai_reagent_stocker run ./data/gene_sets/Tx-Omics_Revision/ \
  --config ./data/config/stock_split_config_priority_example.json
Tx-Omics Revision  ·  gene-set curation  ·  Rosensweig–Shah 2026
14  /  14

### Notes:

<!-- Slide number: 1 -->

TX-OMICS REVISION
Seven gene sets.
Seven clear definitions.

4
3
primary sets
GO-derived sets

2
7
Family 1 captures four primary evidence categories. Family 2 contains three CSW-derived GO term groups, with FC0.5 and FC1 provenance shown separately.
thresholds
stocker CSVs

### Notes:

<!-- Slide number: 2 -->
02
CURATION MAP
Two families, seven non-interchangeable sets
Each set now receives one slide; overlap is reported separately within each family.

FAMILY 1  ·  PRIMARY EVIDENCE
FAMILY 2  ·  CSW-DERIVED GO
6
Mechanistic
Consistent/correlated + published homeostatic-consistent effect
99
01
Ribosomal
FC0.5 99  ·  FC1 12
05
20
Homeostatic genes
Opposing history versus rebound correlations
02
167
Mitochondrial
FC0.5 167  ·  FC1 23
06
97
CSW 4+ genes
Consistent in four or more of seven datasets
03
5
Immune
FC0.5 0  ·  FC1 5
07
4
HLH genes
Potential upstream regulators
04
FC0.5 and FC1 counts are source-provenance counts. A gene may appear at both thresholds.

### Notes:

<!-- Slide number: 3 -->
03
FAMILY 1  ·  SET 01
Mechanistic · 6 genes
Consistent or correlated sleep/wake genes with published effects compatible with homeostatic feedback.

unc79
SIFa
6
NALCN complex
neuropeptide
paper-highlighted genes
Inclusion requires transcriptomic evidence plus published sleep/wake function.

rumpel
AstA-R2
glial transporter
neuropeptide receptor

Consistent across sleep/wake conditions or correlated with sleep/wake history
Published perturbation changes sleep/wake
Effect direction fits a potential homeostatic loop

Trhn
RpL23
serotonin synthesis
ribosome / translation

These six are not the 20 opposing history–rebound correlation genes.

### Notes:

<!-- Slide number: 4 -->
04
FAMILY 1  ·  SET 02
Homeostatic genes · 20 genes
Genes correlated with pre-sampling sleep/wake history and post-sampling rebound in the opposite direction.

15  ·  HISTORY+ / REBOUND−
15
AstA-R2  ·  bumpel  ·  CG33080  ·  CG41378  ·  CG7460
CG7601  ·  Chrac-14  ·  Cyp28c1  ·  dpr21  ·  kumpel
Obp44a  ·  Sk1  ·  Ssk  ·  Trhn  ·  VGlut2
History+ / Rebound−
Positive correlation with prior sleep; negative correlation with rebound.
5

History− / Rebound+

5  ·  HISTORY− / REBOUND+
Negative correlation with prior sleep; positive correlation with rebound.
CG42323  ·  CG7079  ·  CG9313  ·  CG9377  ·  Gtpbp1

### Notes:

<!-- Slide number: 5 -->
05
FAMILY 1  ·  SET 03
CSW 4+ genes · 97 genes
Consistent sleep/wake genes evident across four or more of the seven datasets.

Recurrence across seven datasets
97
82
Frequency 4
FC0.5 wake genes

The category rule is recurrence across datasets—not GO enrichment.
14
Frequency 5

1
Frequency 6

7
also qualify at FC1

All seven are already contained within the 97-gene FC0.5 set.
No sleep-direction gene reaches frequency 4 in this run.

### Notes:

<!-- Slide number: 6 -->
06
FAMILY 1  ·  SET 04
HLH genes · 4 genes
Four bHLH factors identified as potential upstream regulators of wake-induced gene expression.

REGULATORY LOGIC
bigmax
HLH3B
CAGCTG
FBgn0039509
FBgn0011276
E-box enriched upstream of wake-induced genes

E(spl)m3-HLH
E(spl)mbeta-HLH
HOMER motif enrichment → candidate direct regulators

FBgn0002609
FBgn0002733
Results · Table S1

This set is regulatory evidence—not frequency or correlation evidence.

### Notes:

<!-- Slide number: 7 -->
07
FAMILY 2  ·  SET 05
Ribosomal / translation · 99 genes
GO term–grouped hit genes with FC0.5 and FC1 source provenance shown separately.

Most frequent approved terms
THRESHOLD PROVENANCE
99
93
Ribosome

92
Structural constituent of ribosome

unique genes in this set
99
77
FC0.5
Translation
genes supported

12
65
Cytoplasmic translation

FC1
genes supported
58
Cytosolic ribosome
12

BOTH
appear at both thresholds

87 FC0.5-only  ·  0 FC1-only  ·  12 at both thresholds
99 FC0.5 + 12 FC1 − 12 both = 99 unique

### Notes:

<!-- Slide number: 8 -->
08
FAMILY 2  ·  SET 06
Mitochondrial / metabolism · 167 genes
GO term–grouped hit genes with FC0.5 and FC1 source provenance shown separately.

Most frequent approved terms
THRESHOLD PROVENANCE
167
121
Mitochondrion

72
Mitochondrial inner membrane

unique genes in this set
167
46
FC0.5
Oxidative phosphorylation
genes supported

23
30
Mitochondrial translation

FC1
genes supported
24
Mitochondrial matrix
23

BOTH
appear at both thresholds

144 FC0.5-only  ·  0 FC1-only  ·  23 at both thresholds
167 FC0.5 + 23 FC1 − 23 both = 167 unique

### Notes:

<!-- Slide number: 9 -->
09
FAMILY 2  ·  SET 07
Immune · 5 genes
The immune GO term group is FC1-exclusive in this analysis.

Five FC1-supported genes
Approved GO terms
THRESHOLD PROVENANCE
5
BomBc2

Defense response
4

BomS1

unique genes in this set
Antibacterial humoral response
4

0
BomS2
FC0.5

genes supported
Response to bacterium
4

5
IM14

FC1
genes supported
Thor

0
BOTH
appear at both thresholds

0 FC0.5-only  ·  5 FC1-only  ·  0 at both

0 FC0.5 + 5 FC1 − 0 both = 5 unique

### Notes:

<!-- Slide number: 10 -->
10
FAMILY 1  ·  OVERLAP
Only three genes cross primary-set boundaries
Overlap is calculated among Mechanistic, Homeostatic, CSW 4+, and HLH only.

Observed shared genes
Family 1 union
124
2

unique FBgns
Mechanistic × Homeostatic
Trhn  ·  AstA-R2
121
in exactly one set
1

Mechanistic × CSW 4+
RpL23
3
in exactly two sets
0

All other pairings
No shared FBgns
0
in three or four sets

HLH genes are fully disjoint.

### Notes:

<!-- Slide number: 11 -->
11
FAMILY 2  ·  OVERLAP
The GO-derived union contains 240 unique genes
99 ribosomal + 167 mitochondrial + 5 immune − 31 shared ribosomal/mitochondrial genes = 240 unique FBgns.

GO term-group membership · 240 unique genes
Threshold provenance
across the 240-gene union

203
FC0.5 only
68

Ribosomal only

32
both thresholds
Category overlap
31
Ribosomal ∩ mitochondrial

genes retained in both ribosomal and mitochondrial buckets
31
136
5
Mitochondrial only

FC1 only
5
Immune only

Only the five immune genes are FC1-exclusive.

### Notes:

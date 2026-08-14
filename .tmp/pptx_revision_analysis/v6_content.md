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
Consistent/correlated + published homeostatic-consistent effect
Mechanistic
99
01
Ribosomal
FC0.5 99  ·  FC1 12
05
20
Opposing history versus rebound correlations
Homeostatic genes
02
167
Mitochondrial
FC0.5 167  ·  FC1 23
06
97
Consistent in four or more of seven datasets
CSW 4+ genes
03
5
Immune
FC0.5 0  ·  FC1 5
07
4
Potential upstream regulators
HLH genes
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
No sleep-upregulated gene met the ≥4-of-7 recurrence criterion.

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
GO-report assignment is shown here; experiment-level FC0.5 and FC1 regulation is shown on slide 12.

Most frequent approved terms
GO ASSIGNMENT SOURCE
99
93
Ribosome

92
Structural constituent of ribosome

unique genes in this set
99
77
FC0.5 GO
Translation
genes assigned

12
65
Cytoplasmic translation

FC1 GO
genes assigned
58
Cytosolic ribosome
12

BOTH
assigned by both reports

87 FC0.5-only  ·  0 FC1-only  ·  12 at both thresholds
99 FC0.5 + 12 FC1 − 12 both = 99 unique

### Notes:

<!-- Slide number: 8 -->
08
FAMILY 2  ·  SET 06
Mitochondrial / metabolism · 167 genes
GO-report assignment is shown here; experiment-level FC0.5 and FC1 regulation is shown on slide 13.

Most frequent approved terms
GO ASSIGNMENT SOURCE
167
121
Mitochondrion

72
Mitochondrial inner membrane

unique genes in this set
167
46
FC0.5 GO
Oxidative phosphorylation
genes assigned

23
30
Mitochondrial translation

FC1 GO
genes assigned
24
Mitochondrial matrix
23

BOTH
assigned by both reports

144 FC0.5-only  ·  0 FC1-only  ·  23 at both thresholds
167 FC0.5 + 23 FC1 − 23 both = 167 unique

### Notes:

<!-- Slide number: 9 -->
09
FAMILY 2  ·  SET 07
Immune · 5 genes
The set was assigned from FC1 GO enrichment; all five genes also meet FC0.5 regulation criteria.

Five FC1-supported genes
Approved GO terms · genes per term
GO ASSIGNMENT SOURCE
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
FC0.5 GO

genes assigned
Response to bacterium
4

5
IM14

FC1 GO
genes assigned
Genes may occur in more than one approved term.
Thor

0
BOTH
assigned by both reports

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

![Venn diagram showing 68 ribosomal-only, 31 shared ribosomal and mitochondrial, 136 mitochondrial-only, and 5 immune-only genes](Image0.jpg)
GO-assignment provenance
across the 240-gene union
203
FC0.5 only
32
both thresholds
5
FC1 only

240 = union of Ribosomal, Mitochondrial and Immune genes; each FBgn counted once.
Only the five Immune genes are FC1-exclusive.

### Notes:

<!-- Slide number: 12 -->
12
FAMILY 2  ·  SET 05  ·  EXPERIMENTS
Ribosomal / translation · experiment-level regulation
Bars show the share of category genes with logFC > threshold and q < 0.05 in each experiment; genes may appear in multiple bars.

FC0.5  ·  N=99
FC1  ·  N=12
40/99 · 40%
10/12 · 83%
Baseline
Baseline

68/99 · 69%
5/12 · 42%
MechSD 3h
MechSD 3h

97/99 · 98%
12/12 · 100%
MechSD 6h
MechSD 6h

48/99 · 48%
5/12 · 42%
MechSD 12h
MechSD 12h

0/99 · 0%
0/12 · 0%
R85C10>TrpA1
R85C10>TrpA1

25/99 · 25%
0/12 · 0%
R23E10>ChR
R23E10>ChR

37/99 · 37%
0/12 · 0%
THIP
THIP

n = category genes meeting this threshold in at least one experiment
n = category genes meeting this threshold in at least one experiment

MechSD6 regulates nearly all ribosomal genes at both thresholds; FC1 signal is concentrated in Baseline and mechanical deprivation.

### Notes:

<!-- Slide number: 13 -->
13
FAMILY 2  ·  SET 06  ·  EXPERIMENTS
Mitochondrial / metabolism · experiment-level regulation
Bars show the share of category genes with logFC > threshold and q < 0.05 in each experiment; genes may appear in multiple bars.

FC0.5  ·  N=167
FC1  ·  N=24
88/167 · 53%
20/24 · 83%
Baseline
Baseline

126/167 · 75%
10/24 · 42%
MechSD 3h
MechSD 3h

163/167 · 98%
24/24 · 100%
MechSD 6h
MechSD 6h

86/167 · 51%
6/24 · 25%
MechSD 12h
MechSD 12h

1/167 · 1%
0/24 · 0%
R85C10>TrpA1
R85C10>TrpA1

9/167 · 5%
0/24 · 0%
R23E10>ChR
R23E10>ChR

2/167 · 1%
0/24 · 0%
THIP
THIP

n = category genes meeting this threshold in at least one experiment
n = category genes meeting this threshold in at least one experiment

MechSD6 is the dominant mitochondrial response at both thresholds; FC1 hits are absent from R85C10, ChR and THIP.

### Notes:

<!-- Slide number: 14 -->
14
FAMILY 2  ·  SET 07  ·  EXPERIMENTS
Immune · experiment-level regulation
Bars show the share of category genes with logFC > threshold and q < 0.05 in each experiment; genes may appear in multiple bars.

FC0.5  ·  N=5
FC1  ·  N=5
1/5 · 20%
1/5 · 20%
Baseline
Baseline

4/5 · 80%
4/5 · 80%
MechSD 3h
MechSD 3h

3/5 · 60%
3/5 · 60%
MechSD 6h
MechSD 6h

4/5 · 80%
4/5 · 80%
MechSD 12h
MechSD 12h

0/5 · 0%
0/5 · 0%
R85C10>TrpA1
R85C10>TrpA1

0/5 · 0%
0/5 · 0%
R23E10>ChR
R23E10>ChR

0/5 · 0%
0/5 · 0%
THIP
THIP

n = category genes meeting this threshold in at least one experiment
n = category genes meeting this threshold in at least one experiment

All five immune genes meet both regulation thresholds; MechSD3 and MechSD12 regulate four of five genes.

### Notes:

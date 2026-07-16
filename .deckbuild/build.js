const fs = require("fs");
const pptxgen = require("pptxgenjs");
const intro = require("./slides_intro.js");
const ds = require("./slides_datasets.js");
const sum = require("./slides_summary.js");

const data = JSON.parse(fs.readFileSync(__dirname + "/data.json", "utf8"));
const TOTAL = 16;

function shortTier(name) {
  let base = name.split("_n=")[0];
  base = base.split(" Table")[0];
  base = base.replace(/_/g, " ").replace("frequency", "freq").replace("Tx-Omics ", "");
  return base;
}
function tierRows(keys) {
  const out = [];
  keys.forEach(k => (data.tx_tiers[k] || []).forEach(it => out.push([shortTier(it[0]), String(it[1])])));
  return out;
}

const pres = new pptxgen();
pres.defineLayout({ name: "WIDE", width: 13.333, height: 7.5 });
pres.layout = "WIDE";
pres.author = "Allada Lab";
pres.title = "FlyBase Reagent Stocker — Generated Stock Tables";

// 1-4 intro
intro.title(pres, data, TOTAL);
intro.pipeline(pres, data, TOTAL);
intro.overview(pres, data, TOTAL);
intro.anatomy(pres, data, TOTAL);

// 5 vGAT
ds.datasetSlide(pres, {
  kicker: "Run 1 · vGAT screen", title: "GABA gene list (n = 140)",
  n_genes: 140, total: data.targeted["vGAT — GABA gene list"].total, buckets: data.targeted["vGAT — GABA gene list"].buckets,
  schemaNote: "Curated symbol list with FlyBase IDs and a source-dataset tag.",
  cols: [["Gene", "submitted gene symbol"], ["primary_symbol", "current FlyBase symbol"], ["flybase_gene_id", "FBgn key used for matching"], ["corrected_Gene", "manual symbol fix (if any)"], ["Data Set", "source list (e.g. bottom-10% vGAT hits)"]],
  foot: "Source: Aldeb bulk-seq vGAT hits.",
}, 5, TOTAL);

// 6 BPD
ds.datasetSlide(pres, {
  kicker: "Run 2 · Priority v2", title: "BPD-GWAS DIOPT orthologs (n = 168)",
  n_genes: 168, total: data.targeted["BPD-GWAS orthologs"].total, buckets: data.targeted["BPD-GWAS orthologs"].buckets,
  schemaNote: "Human bipolar-disorder GWAS genes mapped to fly orthologs via DIOPT.",
  cols: [["Human_Symbol", "human GWAS gene (Nature 2025 S31)"], ["ext_gene", "fly ortholog symbol"], ["flybase_gene_id", "FBgn key used for matching"], ["DIOPT_Score / Weighted", "ortholog confidence"], ["Rank", "high / moderate / low"], ["Prediction Derived From", "supporting ortholog tools"]],
}, 6, TOTAL);

// 7 DIOPT AD
ds.datasetSlide(pres, {
  kicker: "Run 2 · Priority v2", title: "DIOPT Bellenguez/Lambert — Alzheimer's (n = 227)",
  n_genes: 227, total: data.targeted["DIOPT Bellenguez/Lambert (AD)"].total, buckets: data.targeted["DIOPT Bellenguez/Lambert (AD)"].buckets,
  schemaNote: "Alzheimer's GWAS orthologs (DIOPT) plus pre-pulled BDSC RNAi hints.",
  cols: [["Human.Symbol", "human AD GWAS gene"], ["primary_symbol", "fly ortholog symbol"], ["flybase_gene_id", "FBgn key used for matching"], ["DIOPT.Score / Rank", "ortholog confidence"], ["total_stock_count / BDSC_RNAi_stock_ids", "pre-pulled stock hints"], ["genotype_list / symbol_list", "candidate RNAi genotypes"]],
}, 7, TOTAL);

// 8 Jordan
ds.datasetSlide(pres, {
  kicker: "Run 2 · Priority v2", title: "Jordan VCM candidates (n = 35)",
  n_genes: 35, total: data.targeted["Jordan VCM"].total, buckets: data.targeted["Jordan VCM"].buckets,
  schemaNote: "Simple curated symbol-to-FBgn list (shared by the curated short-lists).",
  cols: [["ext_gene", "submitted gene symbol"], ["primary_symbol", "current FlyBase symbol"], ["flybase_gene_id", "FBgn key used for matching"], ["corrected_ext_gene", "manual symbol fix (if any)"]],
}, 8, TOTAL);

// 9 Slumber NP
ds.datasetSlide(pres, {
  kicker: "Run 2 · Priority v2", title: "Slumber neuropeptides + transmitters (n = 33)",
  n_genes: 33, total: data.targeted["Slumber Neuropeptides+transmitters"].total, buckets: data.targeted["Slumber Neuropeptides+transmitters"].buckets,
  schemaNote: "Simple curated symbol-to-FBgn list of neuropeptide / transmitter genes.",
  cols: [["ext_gene", "submitted gene symbol"], ["primary_symbol", "current FlyBase symbol"], ["flybase_gene_id", "FBgn key used for matching"], ["corrected_ext_gene", "manual symbol fix (if any)"]],
}, 9, TOTAL);

// 10 Slumber RNAi
ds.datasetSlide(pres, {
  kicker: "Run 2 · Priority v2", title: "Slumber RNAi candidates (n = 16)",
  n_genes: 16, total: data.targeted["Slumber RNAi candidates"].total, buckets: data.targeted["Slumber RNAi candidates"].buckets,
  schemaNote: "Simple curated symbol-to-FBgn list of RNAi candidate genes.",
  cols: [["ext_gene", "submitted gene symbol"], ["primary_symbol", "current FlyBase symbol"], ["flybase_gene_id", "FBgn key used for matching"], ["corrected_ext_gene", "manual symbol fix (if any)"]],
}, 10, TOTAL);

// 11 tx cycling
ds.txSlide(pres, {
  kicker: "Run 3 · Transcriptomics", title: "CSW sleep/wake cycling sets (FC0.5 & FC1)",
  schemaNote: "Per-gene cycling + sleep/wake correlation metrics, with AI literature annotations (flai_*).",
  cols: [["ext_gene / gene_id", "fly gene + FBgn"], ["frequency / is_cycling", "cycling evidence + recurrence"], ["sleep_corr_exps / wake_corr_exps", "correlated experiments"], ["Baseline, MechSD3/6/12, THIP", "per-condition statistics"], ["flai_Category / flai_Confidence", "AI sleep-wake call + confidence"], ["flai_Supporting_Refs / flai_Ref_*", "supporting publications"]],
  rightTitle: "Frequency tiers (11)", tiers: tierRows(["CSW FC0.5", "CSW FC1"]),
  foot: "Tier = recurrence frequency across experiments; higher tier = stronger / more reproducible signal.",
}, 11, TOTAL);

// 12 tx literature
ds.txSlide(pres, {
  kicker: "Run 3 · Transcriptomics", title: "SleepHistory correlates, overlaps & homeostatic",
  schemaNote: "Literature-derived correlate sets; AI category / confidence / rationale per gene (per-axis for overlaps).",
  cols: [["ext_gene / flybase_gene_id", "fly gene + FBgn"], ["Category / Confidence", "AI sleep-relevance call"], ["Rationale", "short justification"], ["Supporting_Refs / Ref_*", "publications + metadata"], ["History_* / Rebound_*", "per-axis calls (overlap sets)"], ["homeostatic set", "simple symbol-to-FBgn list"]],
  rightTitle: "Gene-set tables (5)", tiers: tierRows(["SleepHistory Correlates", "Overlapping SleepHistory<>SleepRebound Correlates", "Homeostatic Genes"]),
  foot: "Five tables across three categories; counts are gene-set sizes as submitted.",
}, 12, TOTAL);

// 13-16
sum.summaryTargeted(pres, data, 13, TOTAL);
sum.summaryTx(pres, data, 14, TOTAL);
sum.examples(pres, data, 15, TOTAL);
sum.closing(pres, data, 16, TOTAL);

const out = __dirname + "/Allada_Stocker_Tables_Overview.pptx";
pres.writeFile({ fileName: out }).then(() => console.log("WROTE " + out));

const { P, BC, HEAD, BODY, W, H, shadow, bg, pageFooter, header } = require("./lib.js");
const { card, stat, chip, dataTable, shortBucket } = require("./comp.js");

const BUCKETS = [
  ["Bloomington", "UAS / RNAi", "Stock-center, overexpression / knockdown"],
  ["Vienna (VDRC)", "UAS / RNAi", "VDRC knockdown collection"],
  ["Bloomington", "Allele / Insertion", "Stock-center classical alleles"],
  ["Stock Center", "Allele / Insertion", "Other orderable centers (Kyoto, etc.)"],
  ["Non-Stock-Center", "UAS / RNAi", "Custom phenotype reagents (no FBst)"],
  ["Non-Stock-Center", "Allele / Insertion", "Custom phenotype reagents (no FBst)"],
];

function title(pres, data, TOTAL) {
  const s = pres.addSlide();
  s.background = { color: P.INK };
  s.addShape("rect", { x: 0, y: 0, w: 0.28, h: H, fill: { color: P.TEAL } });
  s.addShape("rect", { x: 0.28, y: 0, w: 0.08, h: H, fill: { color: P.AMBER } });
  s.addText("FLYBASE REAGENT STOCKER", { x: 1.0, y: 1.6, w: 11, h: 0.4, fontFace: BODY, fontSize: 14, color: P.TEALL, charSpacing: 4, bold: true, margin: 0 });
  s.addText("Generated Stock Tables", { x: 0.95, y: 2.05, w: 11.6, h: 1.0, fontFace: HEAD, fontSize: 46, bold: true, color: P.WHITE, margin: 0 });
  s.addText("Prioritized Drosophila reagent tables across three input runs", { x: 1.0, y: 3.15, w: 11.4, h: 0.5, fontFace: BODY, fontSize: 18, color: "C7D2DE", margin: 0 });
  // quick stat ribbon
  const items = [["22", "input gene-set tables"], ["3", "screening runs"], ["6", "priority buckets / set"], ["8,300+", "stocks organized"]];
  let x = 1.0;
  items.forEach((it, i) => {
    s.addText(it[0], { x, y: 4.5, w: 2.7, h: 0.6, fontFace: HEAD, fontSize: 30, bold: true, color: P.AMBER, margin: 0 });
    s.addText(it[1], { x, y: 5.1, w: 2.7, h: 0.5, fontFace: BODY, fontSize: 12, color: "C7D2DE", margin: 0 });
    x += 2.85;
  });
  s.addText("Prepared for Dr. Ravi Allada", { x: 1.0, y: 6.5, w: 8, h: 0.4, fontFace: BODY, fontSize: 13, italic: true, color: P.TEALL, margin: 0 });
}

function pipeline(pres, data, TOTAL) {
  const s = pres.addSlide(); bg(s);
  header(s, "How the tables are built", "From a gene list to six prioritized stock tables");
  // process flow
  const steps = [
    ["Input gene set", "Curated / GWAS / transcriptomics list of fly genes (FBgn)"],
    ["FlyBase match", "Pull all stocks + alleles, constructs, phenotypes, references"],
    ["Filter + tier", "Collection, reagent type, phenotype & sleep-term relevance"],
    ["6 priority tables", "Each stock lands in the first bucket it qualifies for"],
  ];
  let x = 0.55; const y = 1.55, cw = 2.95, gap = 0.27;
  steps.forEach((st, i) => {
    card(s, x, y, cw, 1.5);
    s.addShape("oval", { x: x + 0.18, y: y + 0.2, w: 0.5, h: 0.5, fill: { color: P.TEAL } });
    s.addText(String(i + 1), { x: x + 0.18, y: y + 0.2, w: 0.5, h: 0.5, fontFace: HEAD, fontSize: 18, bold: true, color: P.WHITE, align: "center", valign: "middle", margin: 0 });
    s.addText(st[0], { x: x + 0.8, y: y + 0.22, w: cw - 0.95, h: 0.5, fontFace: BODY, fontSize: 13, bold: true, color: P.INK, valign: "middle", margin: 0 });
    s.addText(st[1], { x: x + 0.2, y: y + 0.82, w: cw - 0.4, h: 0.6, fontFace: BODY, fontSize: 10, color: P.MUT, margin: 0 });
    if (i < 3) s.addText("›", { x: x + cw - 0.02, y: y + 0.35, w: gap + 0.05, h: 0.7, fontFace: HEAD, fontSize: 26, bold: true, color: P.AMBER, align: "center", valign: "middle", margin: 0 });
    x += cw + gap;
  });
  // bucket legend
  s.addText("The six priority buckets (applied top-to-bottom; a stock appears in only one)", { x: 0.55, y: 3.35, w: 12, h: 0.34, fontFace: BODY, fontSize: 13, bold: true, color: P.INK, margin: 0 });
  let bx = 0.55, by = 3.78; const bw = 4.0, bh = 0.78;
  BUCKETS.forEach((b, i) => {
    const col = i % 3, rowi = Math.floor(i / 3);
    const X = 0.55 + col * (bw + 0.18), Y = 3.78 + rowi * (bh + 0.16);
    s.addShape("rect", { x: X, y: Y, w: 0.13, h: bh, fill: { color: BC[i] } });
    s.addShape("rect", { x: X + 0.13, y: Y, w: bw - 0.13, h: bh, fill: { color: P.WHITE }, line: { color: P.LINE, width: 1 } });
    s.addText(b[0] + "  ·  " + b[1], { x: X + 0.28, y: Y + 0.08, w: bw - 0.45, h: 0.32, fontFace: BODY, fontSize: 11.5, bold: true, color: P.TX, margin: 0 });
    s.addText(b[2], { x: X + 0.28, y: Y + 0.4, w: bw - 0.45, h: 0.3, fontFace: BODY, fontSize: 9, color: P.MUT, margin: 0 });
  });
  s.addText([{ text: "Phenotype gate: ", options: { bold: true } }, { text: "every table keeps only stocks with an associated phenotype. Sleep / circadian / rhythm / locomotor terms are used to flag literature-supported reagents.", options: {} }],
    { x: 0.55, y: 6.5, w: 12.2, h: 0.5, fontFace: BODY, fontSize: 10.5, color: P.MUT, margin: 0 });
  pageFooter(s, 2, TOTAL);
}

function overview(pres, data, TOTAL) {
  const s = pres.addSlide(); bg(s);
  header(s, "Inputs", "The input datasets at a glance");
  const rows = [
    [{ text: "vGAT screen", options: { bold: true, color: P.TEAL } }, "GABA gene list", "140", "Bulk-seq vGAT hits (Aldeb)", "Curated list (5 col)"],
    [{ text: "Priority v2", options: { bold: true, color: P.TEAL } }, "BPD-GWAS orthologs", "168", "Bipolar GWAS, Nature 2025 S31", "DIOPT ortholog (19 col)"],
    ["", "DIOPT Bellenguez/Lambert", "227", "Alzheimer's GWAS orthologs", "DIOPT + stocks (36 col)"],
    ["", "Jordan VCM", "35", "Curated candidate genes", "Simple list (4 col)"],
    ["", "Slumber neuropeptides", "33", "Curated neuropeptides/transmitters", "Simple list (4 col)"],
    ["", "Slumber RNAi candidates", "16", "Curated RNAi candidates", "Simple list (4 col)"],
    [{ text: "Transcriptomics", options: { bold: true, color: P.TEAL } }, "CSW FC0.5 (7 tiers)", "788", "Sleep/wake cycling + correlation", "Cycling + AI refs (31 col)"],
    ["", "CSW FC1 (4 tiers)", "114", "Sleep/wake cycling (stricter FC)", "Cycling + AI refs (31 col)"],
    ["", "SleepHistory correlates (2)", "215", "History +/- literature correlates", "Lit. category (18 col)"],
    ["", "History x Rebound overlaps (2)", "20", "Overlapping history/rebound sets", "Lit. category (14 col)"],
    ["", "Homeostatic genes (1)", "6", "Sleep homeostasis set", "Simple list (4 col)"],
  ];
  dataTable(s, 0.55, 1.5, 12.25, ["Run", "Dataset", "Genes", "Source / meaning", "Input schema"], rows, [1.6, 3.3, 0.9, 4.0, 2.45], { rowH: 0.42, bf: 10.5, hf: 11.5 });
  s.addText("22 input tables in total. Gene counts shown as submitted; the transcriptomics run groups 16 frequency tiers under five dataset categories.", { x: 0.55, y: 6.55, w: 12.2, h: 0.4, fontFace: BODY, fontSize: 10, italic: true, color: P.MUT, margin: 0 });
  pageFooter(s, 3, TOTAL);
}

function anatomy(pres, data, TOTAL) {
  const s = pres.addSlide(); bg(s);
  header(s, "Output", "Anatomy of every generated workbook");
  // left: sheets
  card(s, 0.55, 1.5, 5.7, 5.0);
  s.addText("What each workbook contains", { x: 0.78, y: 1.66, w: 5.2, h: 0.4, fontFace: BODY, fontSize: 14, bold: true, color: P.INK, margin: 0 });
  const sheets = [
    ["Sheet 1 – 6", "One per priority bucket (color-coded above)"],
    ["References", "Supporting publications per stock / allele"],
    ["All Phenotypic Stocks", "Full union of phenotype-backed stocks"],
    ["Contents", "Prioritization rules + per-sheet count breakdown"],
  ];
  let y = 2.2;
  sheets.forEach(sh => {
    s.addShape("rect", { x: 0.78, y: y + 0.05, w: 0.1, h: 0.62, fill: { color: P.TEAL } });
    s.addText(sh[0], { x: 0.98, y: y, w: 5.1, h: 0.34, fontFace: BODY, fontSize: 12.5, bold: true, color: P.TX, margin: 0 });
    s.addText(sh[1], { x: 0.98, y: y + 0.32, w: 5.1, h: 0.34, fontFace: BODY, fontSize: 10, color: P.MUT, margin: 0 });
    y += 0.82;
  });
  s.addText([{ text: "Priority rule:  ", options: { bold: true, color: P.AMBER } }, { text: "a stock is placed in the first bucket it qualifies for, so the six tables never double-count a stock.", options: { color: P.TX } }],
    { x: 0.78, y: 5.7, w: 5.25, h: 0.7, fontFace: BODY, fontSize: 10.5, margin: 0 });
  // right: key columns grouped
  card(s, 6.45, 1.5, 6.35, 5.0);
  s.addText("Key columns (≈59 total), grouped", { x: 6.68, y: 1.66, w: 5.9, h: 0.4, fontFace: BODY, fontSize: 14, bold: true, color: P.INK, margin: 0 });
  const groups = [
    ["Identity", "FBst, stock_number, collection, genotype"],
    ["Reagent type", "RNAi, UAS, GAL4, mutant, sgRNA, custom_stock"],
    ["Gene / allele", "gene_symbols, flybase_gene_ids, AlleleSymbol, class"],
    ["Construct parts", "FBti / FBtp / FBal symbols + IDs, Balancers"],
    ["Literature", "PMID, total_refs, sleep/circadian ref counts"],
    ["Phenotype + similarity", "PHENOTYPE_RELEVANCE_SCORE, cosine similarities"],
  ];
  y = 2.2;
  groups.forEach((g, i) => {
    s.addShape("oval", { x: 6.68, y: y + 0.02, w: 0.16, h: 0.16, fill: { color: BC[i] } });
    s.addText(g[0], { x: 6.95, y: y - 0.04, w: 5.6, h: 0.3, fontFace: BODY, fontSize: 12, bold: true, color: P.TX, margin: 0 });
    s.addText(g[1], { x: 6.95, y: y + 0.26, w: 5.7, h: 0.3, fontFace: BODY, fontSize: 9.5, color: P.MUT, margin: 0 });
    y += 0.72;
  });
  pageFooter(s, 4, TOTAL);
}

module.exports = { title, pipeline, overview, anatomy };

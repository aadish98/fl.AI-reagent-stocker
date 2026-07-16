#!/usr/bin/env node
/** Build Tx-Omics Follow-Up v3 gene-set overview deck for PI presentation. */
const fs = require("fs");
const path = require("path");
const { execFileSync } = require("child_process");
const pptxgen = require("../Tx-Omics_FollowUp/node_modules/pptxgenjs");
const { P, BC, HEAD, BODY, W, H, bg, pageFooter, header } = require("./lib.js");
const { card, stat, dataTable, hbars } = require("./comp.js");

const BREAKDOWN = path.resolve(__dirname, "../data/gene_sets/Tx-Omics-FollowUp_v3/Breakdown");
const OUT = path.join(BREAKDOWN, "Tx-Omics_GeneSets_Overview.pptx");
const TOTAL = 8;

function parseCSVLine(line) {
  const out = [];
  let cur = "";
  let inQ = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (c === '"') { inQ = !inQ; continue; }
    if (c === "," && !inQ) { out.push(cur); cur = ""; continue; }
    cur += c;
  }
  out.push(cur);
  return out;
}

function parseGeneSets() {
  const comboPath = path.join(BREAKDOWN, "combination_counts_summary.csv");
  const combo = {};
  fs.readFileSync(comboPath, "utf8").trim().split("\n").slice(1).forEach(line => {
    const p = parseCSVLine(line);
    const gs = p[1];
    if (!combo[gs]) combo[gs] = { genes: 0, stocks: 0, alleles: 0 };
    combo[gs].genes += parseInt(p[3], 10) || 0;
    combo[gs].stocks += parseInt(p[4], 10) || 0;
    combo[gs].alleles += parseInt(p[5], 10) || 0;
  });

  const sets = fs.readdirSync(BREAKDOWN)
    .filter(f => f.endsWith(".csv") && f !== "combination_counts_summary.csv")
    .map(f => {
      const m = f.match(/n=(\d+)genes/);
      const genes = m ? parseInt(m[1], 10) : fs.readFileSync(path.join(BREAKDOWN, f), "utf8").trim().split("\n").length - 1;
      const name = f.replace(".csv", "");
      let cat = "Other";
      if (name.startsWith("FC0.5_Sleep")) cat = "FC0.5 · Sleep rebound";
      else if (name.startsWith("FC0.5_Wake")) cat = "FC0.5 · Wake rebound";
      else if (name.startsWith("FC1_Sleep")) cat = "FC1 · Sleep rebound";
      else if (name.startsWith("FC1_Wake")) cat = "FC1 · Wake rebound";
      else if (name.includes("History") && name.includes("overlap")) cat = "History × Rebound overlap";
      else if (name.includes("History")) cat = "Sleep history correlates";
      else if (name.includes("homeostatic")) cat = "Curated · Homeostatic core";
      else if (name.includes("Mechanistic")) cat = "Curated · Mechanistic screen";
      return { file: f, name, genes, cat, stocks: (combo[name] || {}).stocks || 0, alleles: (combo[name] || {}).alleles || 0 };
    })
    .sort((a, b) => b.genes - a.genes);

  const mechPath = path.join(BREAKDOWN, "Tx-Omics_Mechanistic_Screen_Candidates_nsyb_elav.csv");
  const mech = fs.readFileSync(mechPath, "utf8").trim().split("\n").slice(1).map(line => {
    const cols = parseCSVLine(line);
    return {
      rank: cols[0],
      gene: cols[1],
      category: (cols[3] || "").split(" / ")[0],
      priority: cols[8] || "",
    };
  });

  return { sets, mech, combo };
}

function shortName(name) {
  return name
    .replace(/_n=\d+genes/, "")
    .replace(/ Table/, "")
    .replace(/frequency_/, "freq ")
    .replace(/_/g, " ")
    .replace("Tx-Omics ", "");
}

function titleSlide(pres, data) {
  const s = pres.addSlide();
  s.background = { color: P.INK };
  s.addShape("rect", { x: 0, y: 0, w: 0.28, h: H, fill: { color: P.TEAL } });
  s.addShape("rect", { x: 0.28, y: 0, w: 0.08, h: H, fill: { color: P.AMBER } });
  s.addText("TX-OMICS FOLLOW-UP v3", { x: 1.0, y: 1.35, w: 11, h: 0.38, fontFace: BODY, fontSize: 13, color: P.TEALL, charSpacing: 4, bold: true, margin: 0 });
  s.addText("Gene Set Breakdown", { x: 0.95, y: 1.85, w: 11.6, h: 1.05, fontFace: HEAD, fontSize: 44, bold: true, color: P.WHITE, margin: 0 });
  s.addText("Shared transcriptomic responses to sleep–wake manipulations", { x: 1.0, y: 2.95, w: 11.2, h: 0.45, fontFace: BODY, fontSize: 17, color: "C7D2DE", margin: 0 });
  s.addText("Rosensweig, Shah et al. · bioRxiv 2026.02.28.708752", { x: 1.0, y: 3.45, w: 11, h: 0.35, fontFace: BODY, fontSize: 12, italic: true, color: P.TEALL, margin: 0 });

  const stats = [
    [String(data.sets.length), "gene-set tiers"],
    [String(data.sets.reduce((a, x) => a + x.genes, 0)), "genes (tier sum)"],
    [String(data.sets.reduce((a, x) => a + x.stocks, 0).toLocaleString()), "FlyBase stocks mapped"],
    ["6", "homeostatic core + 19 screen candidates"],
  ];
  let x = 1.0;
  stats.forEach(([v, l]) => {
    s.addText(v, { x, y: 4.55, w: 2.65, h: 0.58, fontFace: HEAD, fontSize: 28, bold: true, color: P.AMBER, margin: 0 });
    s.addText(l, { x, y: 5.12, w: 2.65, h: 0.45, fontFace: BODY, fontSize: 11, color: "C7D2DE", margin: 0 });
    x += 2.85;
  });
  s.addText("Prepared for Dr. Ravi Allada · July 2026", { x: 1.0, y: 6.55, w: 8, h: 0.35, fontFace: BODY, fontSize: 12, italic: true, color: P.TEALL, margin: 0 });
}

function contextSlide(pres) {
  const s = pres.addSlide(); bg(s);
  header(s, "Biological framing", "Why these gene sets exist");
  const blocks = [
    ["Question", "Is there a homeostatic “sleeper” gene battery—analogous to per for circadian timing—that tracks sleep–wake history across diverse manipulations?"],
    ["Approach", "Brain transcriptomics across 7 datasets: mechanical SD, thermo/optogenetic activation, THIP, and baseline cycling. Genes ranked by fold-change consistency (FC ≥ 0.5 or ≥ 1) and cross-dataset frequency."],
    ["This deck", "17 tiered CSV gene lists in Breakdown/, each with per-set FlyBase stock workbooks. Sets move from broad rebound tiers → history correlates → curated mechanistic priorities."],
  ];
  blocks.forEach((b, i) => {
    const y = 1.55 + i * 1.65;
    card(s, 0.55, y, 12.2, 1.45);
    s.addShape("rect", { x: 0.55, y, w: 0.12, h: 1.45, fill: { color: BC[i] } });
    s.addText(b[0], { x: 0.85, y: y + 0.12, w: 2.2, h: 0.35, fontFace: BODY, fontSize: 13, bold: true, color: BC[i], margin: 0 });
    s.addText(b[1], { x: 0.85, y: y + 0.48, w: 11.6, h: 0.85, fontFace: BODY, fontSize: 12, color: P.TX, margin: 0 });
  });
  pageFooter(s, 2, TOTAL);
}

function taxonomySlide(pres, data) {
  const s = pres.addSlide(); bg(s);
  header(s, "Organization", "Six logical categories · 17 gene-set files");

  const cats = [
    ["FC0.5 · Sleep rebound", "Genes consistently up/down with sleep rebound (log₂FC ≥ 0.5) across manipulations", "2 tiers · 141 genes"],
    ["FC0.5 · Wake rebound", "Wake-enriched rebound responders; largest tiered collection", "5 tiers · 647 genes"],
    ["FC1 · Sleep / Wake", "Stricter FC ≥ 1 filter — higher-confidence, smaller lists", "4 tiers · 114 genes"],
    ["Sleep history correlates", "Expression tracks prior sleep/wake history (± literature via fl.ai)", "2 sets · 215 genes"],
    ["History × Rebound overlap", "Dual homeostatic signature: history correlate ∩ rebound responder", "2 sets · 20 genes"],
    ["Curated priorities", "6-gene homeostatic core + 19 mechanistic screen candidates (nsyb/elav)", "2 sets · 25 genes"],
  ];

  let y = 1.48;
  cats.forEach((c, i) => {
    const col = i % 2, row = Math.floor(i / 2);
    const X = 0.55 + col * 6.35, Y = 1.48 + row * 1.72;
    card(s, X, Y, 6.1, 1.55);
    s.addShape("oval", { x: X + 0.18, y: Y + 0.22, w: 0.42, h: 0.42, fill: { color: BC[i] } });
    s.addText(String(i + 1), { x: X + 0.18, y: Y + 0.22, w: 0.42, h: 0.42, fontFace: HEAD, fontSize: 14, bold: true, color: P.WHITE, align: "center", valign: "middle", margin: 0 });
    s.addText(c[0], { x: X + 0.72, y: Y + 0.18, w: 5.2, h: 0.38, fontFace: BODY, fontSize: 12.5, bold: true, color: P.INK, margin: 0 });
    s.addText(c[2], { x: X + 0.72, y: Y + 0.52, w: 5.2, h: 0.28, fontFace: BODY, fontSize: 10, bold: true, color: P.TEAL, margin: 0 });
    s.addText(c[1], { x: X + 0.22, y: Y + 0.88, w: 5.65, h: 0.55, fontFace: BODY, fontSize: 10, color: P.MUT, margin: 0 });
  });

  s.addText("Frequency tiers (freq 2–6): number of datasets in which a gene passes the FC threshold — higher freq = more robust.", { x: 0.55, y: 6.55, w: 12.2, h: 0.32, fontFace: BODY, fontSize: 10, italic: true, color: P.MUT, margin: 0 });
  pageFooter(s, 3, TOTAL);
}

function reboundSlide(pres, data) {
  const s = pres.addSlide(); bg(s);
  header(s, "Rebound tiers", "Consistent sleep/wake responders by FC threshold & frequency");

  const byFreq = (a, b) => {
    const fa = parseInt(a.name.match(/frequency_(\d+)/)[1], 10);
    const fb = parseInt(b.name.match(/frequency_(\d+)/)[1], 10);
    return fa - fb;
  };
  const sleepFc05 = data.sets.filter(x => x.cat === "FC0.5 · Sleep rebound").sort(byFreq);
  const wakeFc05 = data.sets.filter(x => x.cat === "FC0.5 · Wake rebound").sort(byFreq);
  const fc1 = data.sets.filter(x => x.cat.startsWith("FC1")).sort(byFreq);

  card(s, 0.55, 1.48, 6.0, 4.85);
  s.addText("FC0.5 · genes per frequency tier", { x: 0.78, y: 1.62, w: 5.5, h: 0.32, fontFace: BODY, fontSize: 12.5, bold: true, color: P.INK, margin: 0 });
  const barSleep = sleepFc05.map(x => ["Sleep freq " + x.name.match(/frequency_(\d+)/)[1], x.genes]);
  const barWake = wakeFc05.map(x => ["Wake freq " + x.name.match(/frequency_(\d+)/)[1], x.genes]);
  hbars(s, 0.78, 2.05, 5.5, barSleep, P.TEAL, { rh: 0.42, labelW: 1.35, valW: 0.55, lf: 9 });
  hbars(s, 0.78, 3.55, 5.5, barWake, "1C6FB0", { rh: 0.42, labelW: 1.35, valW: 0.55, lf: 9 });

  card(s, 6.75, 1.48, 6.0, 4.85);
  s.addText("FC1 (stricter) · all tiers", { x: 6.98, y: 1.62, w: 5.5, h: 0.32, fontFace: BODY, fontSize: 12.5, bold: true, color: P.INK, margin: 0 });
  const rows = fc1.map(x => {
    const label = x.name.includes("Sleep") ? "Sleep" : "Wake";
    const freq = x.name.match(/frequency_(\d+)/)[1];
    return [label + " · freq " + freq, String(x.genes), String(x.stocks)];
  });
  dataTable(s, 6.98, 2.05, 5.5, ["Tier", "Genes", "Stocks"], rows, [2.1, 1.0, 1.2], { rowH: 0.38, bf: 10, hf: 10.5 });

  s.addText([
    { text: "Key pattern: ", options: { bold: true, color: P.INK } },
    { text: "Wake rebound sets dominate at FC0.5 (353 genes at freq 2 alone). Stricter FC1 collapses to 114 genes — useful for high-confidence follow-up.", options: { color: P.MUT } },
  ], { x: 0.55, y: 6.55, w: 12.2, h: 0.32, fontFace: BODY, fontSize: 10, margin: 0 });
  pageFooter(s, 4, TOTAL);
}

function historySlide(pres, data) {
  const s = pres.addSlide(); bg(s);
  header(s, "History & overlap", "Homeostatic signatures beyond raw rebound");

  const hist = data.sets.filter(x => x.cat === "Sleep history correlates");
  const ov = data.sets.filter(x => x.cat === "History × Rebound overlap");

  card(s, 0.55, 1.5, 5.85, 2.35);
  stat(s, 0.55, 1.65, 5.85, hist.reduce((a, x) => a + x.genes, 0), "history-correlate genes", P.TEAL);
  s.addText("Pos correlates (n=149): genes whose expression rises with prior wake.\nNeg correlates (n=66): genes tracking prior sleep.", { x: 0.78, y: 2.85, w: 5.4, h: 0.75, fontFace: BODY, fontSize: 10.5, color: P.MUT, margin: 0 });

  card(s, 6.65, 1.5, 6.1, 2.35);
  stat(s, 6.65, 1.65, 6.1, ov.reduce((a, x) => a + x.genes, 0), "overlap genes (dual signature)", P.AMBER);
  s.addText("History(−) × Rebound(+) → 5 genes (e.g. CG9377)\nHistory(+) × Rebound(−) → 15 genes (e.g. AstA-R2, Mip pathway)", { x: 6.88, y: 2.85, w: 5.6, h: 0.75, fontFace: BODY, fontSize: 10.5, color: P.MUT, margin: 0 });

  card(s, 0.55, 4.05, 12.2, 2.15);
  s.addText("Interpretation for mechanistic follow-up", { x: 0.78, y: 4.18, w: 11.5, h: 0.32, fontFace: BODY, fontSize: 12.5, bold: true, color: P.INK, margin: 0 });
  const bullets = [
    "Overlap sets enrich for genes with opposing history vs. rebound dynamics — prime homeostatic regulators.",
    "fl.ai literature categories (sleep / circadian / both) annotate each gene; high-confidence overlaps feed the mechanistic screen list.",
    "Per-gene CSVs include cycling status, manipulation correlations (MechSD, TrpA1, THIP), and PMC-backed rationales.",
  ];
  s.addText(bullets.map((t, i) => ({
    text: t, options: { bullet: true, breakLine: i < bullets.length - 1, color: P.TX, fontSize: 11, fontFace: BODY },
  })), { x: 0.78, y: 4.55, w: 11.5, h: 1.45, margin: 0, paraSpaceAfter: 6 });
  pageFooter(s, 5, TOTAL);
}

function prioritySlide(pres, data) {
  const s = pres.addSlide(); bg(s);
  header(s, "Curated priorities", "Homeostatic core + mechanistic screen candidates");

  card(s, 0.55, 1.48, 4.0, 4.9);
  s.addText("6-gene homeostatic core", { x: 0.78, y: 1.62, w: 3.5, h: 0.32, fontFace: BODY, fontSize: 12, bold: true, color: P.INK, margin: 0 });
  s.addText("unc79 · SIFa · rumpel · AstA-R2 · Trhn · RpL23", { x: 0.78, y: 2.05, w: 3.5, h: 0.55, fontFace: BODY, fontSize: 10.5, color: P.TEAL, bold: true, margin: 0 });
  s.addText("NALCN channel auxiliaries, neuropeptide signaling, glial transport, serotonin biosynthesis, ribosome biogenesis — convergent homeostatic biology.", { x: 0.78, y: 2.65, w: 3.5, h: 1.1, fontFace: BODY, fontSize: 10, color: P.MUT, margin: 0 });
  stat(s, 0.78, 3.85, 3.5, "25", "FlyBase stocks", P.AMBER);

  card(s, 4.75, 1.48, 8.0, 4.9);
  s.addText("Top mechanistic screen candidates (nsyb-GAL4 / elav-GAL4)", { x: 4.98, y: 1.62, w: 7.5, h: 0.32, fontFace: BODY, fontSize: 12, bold: true, color: P.INK, margin: 0 });
  const top = data.mech.slice(0, 10).map(m => [
    m.rank, m.gene,
    { text: m.category, options: { fontSize: 9, color: P.MUT } },
    { text: m.priority, options: { bold: m.priority.startsWith("High"), color: m.priority.startsWith("High") ? P.TEAL : P.TX, fontSize: 9.5 } },
  ]);
  dataTable(s, 4.98, 2.0, 7.5, ["#", "Gene", "Mechanism", "Priority"], top, [0.35, 0.85, 3.5, 0.95], { rowH: 0.33, bf: 9.5, hf: 10, margin: 2 });
  s.addText("Full ranked list: 19 genes · Tx-Omics_Mechanistic_Screen_Candidates_nsyb_elav.csv", { x: 4.98, y: 5.55, w: 7.5, h: 0.3, fontFace: BODY, fontSize: 9, italic: true, color: P.MUT, margin: 0 });
  pageFooter(s, 6, TOTAL);
}

function stockSlide(pres, data) {
  const s = pres.addSlide(); bg(s);
  header(s, "Reagent coverage", "FlyBase stocks per gene-set tier");

  const topSets = [...data.sets].sort((a, b) => b.stocks - a.stocks).slice(0, 8);
  card(s, 0.55, 1.48, 7.0, 4.9);
  s.addText("Largest tiers by stock count", { x: 0.78, y: 1.62, w: 6.5, h: 0.32, fontFace: BODY, fontSize: 12.5, bold: true, color: P.INK, margin: 0 });
  hbars(s, 0.78, 2.05, 6.5, topSets.map(x => [shortName(x.name).slice(0, 28), x.stocks]), P.TEAL, { rh: 0.48, labelW: 2.8, valW: 0.7, lf: 9 });

  card(s, 7.75, 1.48, 5.0, 4.9);
  const totalStocks = data.sets.reduce((a, x) => a + x.stocks, 0);
  const totalAlleles = data.sets.reduce((a, x) => a + x.alleles, 0);
  s.addText(String(totalStocks.toLocaleString()), { x: 7.95, y: 1.85, w: 4.6, h: 0.55, fontFace: HEAD, fontSize: 28, bold: true, color: P.AMBER, align: "center", margin: 0 });
  s.addText("stocks across all tiers", { x: 7.95, y: 2.38, w: 4.6, h: 0.28, fontFace: BODY, fontSize: 10, color: P.MUT, align: "center", margin: 0 });
  s.addText(String(totalAlleles.toLocaleString()), { x: 7.95, y: 2.95, w: 4.6, h: 0.55, fontFace: HEAD, fontSize: 28, bold: true, color: P.TEAL, align: "center", margin: 0 });
  s.addText("alleles mapped", { x: 7.95, y: 3.48, w: 4.6, h: 0.28, fontFace: BODY, fontSize: 10, color: P.MUT, align: "center", margin: 0 });
  s.addText("Each tier → 6 priority stock tables\n(Bloomington / VDRC RNAi, alleles, non-stock-center)", { x: 7.95, y: 4.2, w: 4.6, h: 0.85, fontFace: BODY, fontSize: 10, color: P.MUT, margin: 0 });

  s.addText("Workbooks live in Breakdown/Per Gene Set Runs/ · summary in combination_counts_summary.csv", { x: 0.55, y: 6.55, w: 12.2, h: 0.32, fontFace: BODY, fontSize: 10, italic: true, color: P.MUT, margin: 0 });
  pageFooter(s, 7, TOTAL);
}

function closingSlide(pres) {
  const s = pres.addSlide();
  s.background = { color: P.INK };
  s.addShape("rect", { x: 0, y: H - 0.06, w: W, h: 0.06, fill: { color: P.AMBER } });
  s.addText("Takeaways", { x: 0.75, y: 1.2, w: 11, h: 0.55, fontFace: HEAD, fontSize: 32, bold: true, color: P.WHITE, margin: 0 });
  const pts = [
    "17 tiered gene lists span rebound consistency (FC0.5/FC1), sleep-history correlates, and dual-signature overlaps.",
    "Frequency bins prioritize robust responders across 7 manipulation datasets — wake tiers are largest; FC1 narrows to high-confidence targets.",
    "6-gene homeostatic core + 19 ranked mechanistic candidates are ready for nsyb/elav RNAi screening with mapped FlyBase reagents.",
    "Next: prioritize overlap + high-freq FC1 wake genes for behavioral validation; stock tables are generated per tier.",
  ];
  s.addText(pts.map((t, i) => ({
    text: t, options: { bullet: true, breakLine: i < pts.length - 1, color: "C7D2DE", fontSize: 13, fontFace: BODY },
  })), { x: 0.75, y: 2.0, w: 11.5, h: 3.5, margin: 0, paraSpaceAfter: 10 });
  s.addText("data/gene_sets/Tx-Omics-FollowUp_v3/Breakdown/", { x: 0.75, y: 6.2, w: 11, h: 0.35, fontFace: BODY, fontSize: 11, color: P.TEALL, margin: 0 });
  pageFooter(s, 8, TOTAL);
}

async function main() {
  execFileSync("python3", [path.join(__dirname, "build_txomics_breakdown_pptx.py")], { stdio: "inherit" });
}

main().catch(err => { console.error(err); process.exit(1); });

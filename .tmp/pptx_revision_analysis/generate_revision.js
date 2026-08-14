const pptxgen = require("../pptx_tool/node_modules/pptxgenjs");
const path = require("path");

const pptx = new pptxgen();
pptx.layout = "LAYOUT_WIDE";
pptx.author = "fl.ai";
pptx.subject = "Tx-Omics Revision gene-set curation";
pptx.title = "Tx-Omics Revision — Gene-set curation";
pptx.company = "fl.ai";
pptx.lang = "en-US";
const FONT = "Avenir Next";
pptx.theme = {
  headFontFace: FONT,
  bodyFontFace: FONT,
  lang: "en-US",
};
pptx.defineLayout({ name: "CUSTOM_WIDE", width: 13.333, height: 7.5 });
pptx.layout = "CUSTOM_WIDE";

const C = {
  ink: "17212B",
  navy: "122634",
  teal: "1D6F69",
  teal2: "5DA69D",
  mint: "DCEDE8",
  coral: "D96C4C",
  gold: "D8A238",
  cream: "F7F3EA",
  white: "FFFFFF",
  paper: "FFFCF6",
  muted: "66727A",
  line: "D7D8D2",
  pale: "ECEBE4",
  redPale: "F5DED6",
};

const OUT = path.resolve(__dirname, "Tx-Omics_Revision_gene_set_curation_revised.pptx");
const FIG = path.resolve(__dirname, "../../audit_outputs/Tx-Omics_Revision/Figures");

function addPage(slide, n, section = "GENE-SET CURATION") {
  slide.addText(section, {
    x: 0.68, y: 0.26, w: 3.4, h: 0.24, margin: 0,
    fontFace: FONT, fontSize: 10.5, bold: true, charSpacing: 1.8, color: C.teal,
  });
  slide.addText(String(n).padStart(2, "0"), {
    x: 12.1, y: 0.24, w: 0.55, h: 0.28, margin: 0,
    fontFace: FONT, fontSize: 11, bold: true, align: "right", color: C.muted,
  });
}

function addTitle(slide, title, subtitle, n, section) {
  addPage(slide, n, section);
  slide.addText(title, {
    x: 0.68, y: 0.66, w: 11.8, h: 0.58, margin: 0,
    fontFace: FONT, fontSize: 29, bold: true, color: C.ink,
    breakLine: false, fit: "shrink",
  });
  if (subtitle) {
    slide.addText(subtitle, {
      x: 0.68, y: 1.27, w: 11.8, h: 0.36, margin: 0,
      fontFace: FONT, fontSize: 15.5, color: C.muted, fit: "shrink",
    });
  }
}

function addFooter(slide, text = "Rosensweig–Shah 2026  ·  flybase_gene_id is the membership key") {
  slide.addText(text, {
    x: 0.68, y: 7.12, w: 11.2, h: 0.18, margin: 0,
    fontFace: FONT, fontSize: 10.5, color: C.muted,
  });
}

function card(slide, x, y, w, h, opts = {}) {
  slide.addShape(pptx.ShapeType.roundRect, {
    x, y, w, h, rectRadius: 0.08,
    fill: { color: opts.fill || C.white },
    line: { color: opts.line || C.line, width: opts.lineWidth || 0.8 },
    shadow: opts.shadow ? { type: "outer", color: "000000", blur: 1.5, offset: 1, angle: 135, opacity: 0.09 } : undefined,
  });
}

function label(slide, text, x, y, w, color = C.teal, fill = C.mint) {
  slide.addShape(pptx.ShapeType.roundRect, {
    x, y, w, h: 0.32, rectRadius: 0.08,
    fill: { color: fill }, line: { color: fill },
  });
  slide.addText(text, {
    x: x + 0.08, y: y + 0.055, w: w - 0.16, h: 0.18, margin: 0,
    fontSize: 9.5, bold: true, charSpacing: 1.1, align: "center", color,
  });
}

function metric(slide, value, caption, x, y, w, color = C.teal) {
  slide.addText(String(value), {
    x, y, w, h: 0.48, margin: 0, fontSize: 31, bold: true,
    fontFace: FONT, color,
  });
  slide.addText(caption, {
    x, y: y + 0.51, w, h: 0.35, margin: 0, fontSize: 12.5,
    color: C.muted, fit: "shrink",
  });
}

function bulletList(slide, items, x, y, w, h, fontSize = 15, color = C.ink) {
  const runs = [];
  items.forEach((item, i) => runs.push({
    text: item,
    options: { bullet: { indent: fontSize * 1.2 }, breakLine: i !== items.length - 1, paraSpaceAfterPt: 10 },
  }));
  slide.addText(runs, {
    x, y, w, h, margin: 0.04, fontFace: FONT, fontSize,
    color, breakLine: false, valign: "top", fit: "shrink",
  });
}

// 1 — Title
{
  const s = pptx.addSlide();
  s.background = { color: C.navy };
  s.addShape(pptx.ShapeType.rect, { x: 0, y: 0, w: 0.18, h: 7.5, fill: { color: C.coral }, line: { color: C.coral } });
  s.addText("TX-OMICS REVISION", {
    x: 0.82, y: 0.72, w: 4.4, h: 0.34, margin: 0,
    fontSize: 13, bold: true, charSpacing: 2.8, color: C.teal2,
  });
  s.addText("Seven gene sets.\nOne auditable curation.", {
    x: 0.82, y: 1.45, w: 7.9, h: 1.65, margin: 0,
    fontFace: FONT, fontSize: 39, bold: true, color: C.white,
    breakLine: false, fit: "shrink",
  });
  s.addText("A concise record of the publication rules, pathway selection, identity checks, and overlap that produced the stocker-ready inputs.", {
    x: 0.84, y: 3.35, w: 7.05, h: 0.92, margin: 0,
    fontSize: 18, color: "D6E0E3", breakLine: false, fit: "shrink",
  });
  const stats = [["7", "stocker CSVs"], ["2", "analysis families"], ["4", "primary sets"], ["3", "GO term sets"]];
  stats.forEach((d, i) => {
    const x = 8.65 + (i % 2) * 2.0;
    const y = 1.48 + Math.floor(i / 2) * 1.65;
    s.addShape(pptx.ShapeType.roundRect, {
      x, y, w: 1.7, h: 1.35, rectRadius: 0.08,
      fill: { color: i === 0 ? C.coral : "1B3443" }, line: { color: i === 0 ? C.coral : "35505E", width: 0.8 },
    });
    s.addText(d[0], { x: x + 0.18, y: y + 0.18, w: 1.35, h: 0.52, margin: 0, fontSize: 30, bold: true, color: C.white });
    s.addText(d[1], { x: x + 0.18, y: y + 0.79, w: 1.34, h: 0.32, margin: 0, fontSize: 11.5, color: "D6E0E3", fit: "shrink" });
  });
  s.addText("Rosensweig–Shah 2026  ·  curation summary", {
    x: 0.84, y: 6.75, w: 5, h: 0.25, margin: 0, fontSize: 10.5, color: "9CB0B8",
  });
}

// 2 — Two analytical families
{
  const s = pptx.addSlide();
  s.background = { color: C.paper };
  addTitle(s, "Two curation families—analyzed separately", "Primary evidence sets and CSW-derived GO term–grouped sets use separate overlap analyses.", 2, "LIBRARY MAP");
  const sets = [
    ["01", "6", "Mechanistic", "Consistent/correlated; published effect fits homeostatic feedback", C.coral],
    ["02", "20", "Homeostatic genes", "Opposing history ↔ rebound correlations: 15 + 5", C.gold],
    ["03", "97", "CSW 4+ genes", "Consistent in ≥4 of 7 datasets", C.teal],
    ["04", "4", "HLH genes", "Potential upstream regulators", "665C86"],
    ["05", "99", "Ribosomal", "DAVID term hits; approved keywords", C.coral],
    ["06", "167", "Mitochondrial", "DAVID mito/metabolism term hits", C.teal],
    ["07", "5", "Immune", "DAVID immune / Toll / Imd term hits", C.gold],
  ];
  const COLS = [0.68, 3.71, 6.74, 9.77];
  const ROWS = [2.05, 4.52];
  const CW = 2.70, CH = 1.90;
  s.addText("FAMILY 1  ·  PRIMARY CURATION SETS", {
    x: 0.68, y: 1.74, w: 4.2, h: 0.20, margin: 0,
    fontSize: 10.5, bold: true, charSpacing: 1.3, color: C.coral,
  });
  s.addText("FAMILY 2  ·  CSW-DERIVED GO TERM–GROUPED SETS (FC0.5 / FC1)", {
    x: 0.68, y: 4.21, w: 6.2, h: 0.20, margin: 0,
    fontSize: 10.5, bold: true, charSpacing: 1.3, color: C.teal,
  });
  sets.forEach((d, i) => {
    const x = COLS[i % 4];
    const y = ROWS[Math.floor(i / 4)];
    card(s, x, y, CW, CH, { fill: C.white, shadow: true, line: "E1DED4" });
    s.addShape(pptx.ShapeType.rect, { x, y, w: 0.11, h: CH, fill: { color: d[4] }, line: { color: d[4] } });
    s.addText(d[0], { x: x + 0.30, y: y + 0.22, w: 0.9, h: 0.22, margin: 0, fontSize: 10, bold: true, charSpacing: 1.2, color: d[4] });
    s.addText(d[1], { x: x + 0.30, y: y + 0.46, w: 1.1, h: 0.48, margin: 0, fontSize: 29, bold: true, color: C.ink });
    s.addText("genes", { x: x + 1.18, y: y + 0.65, w: 0.65, h: 0.20, margin: 0, fontSize: 10.5, color: C.muted });
    s.addText(d[2], { x: x + 0.30, y: y + 1.05, w: CW - 0.56, h: 0.26, margin: 0, fontSize: 14.5, bold: true, color: C.ink, fit: "shrink" });
    s.addText(d[3], { x: x + 0.30, y: y + 1.37, w: CW - 0.56, h: 0.38, margin: 0, fontSize: 11.5, color: C.muted, valign: "top", fit: "shrink" });
  });
  const ux = COLS[3], uy = ROWS[1];
  card(s, ux, uy, CW, CH, { fill: C.navy, line: C.navy, shadow: true });
  s.addText("RULE", { x: ux + 0.30, y: uy + 0.22, w: 0.9, h: 0.22, margin: 0, fontSize: 10, bold: true, charSpacing: 1.2, color: C.teal2 });
  s.addText("Keep them separate", { x: ux + 0.30, y: uy + 0.55, w: 2.05, h: 0.34, margin: 0, fontSize: 19, bold: true, color: C.white, fit: "shrink" });
  s.addText("Report overlap within each family—never as one seven-set comparison.", {
    x: ux + 0.30, y: uy + 1.10, w: CW - 0.56, h: 0.54, margin: 0,
    fontSize: 12.5, color: "C7D4D8", valign: "top", fit: "shrink",
  });
  addFooter(s);
}

// 3 — Pipeline
{
  const s = pptx.addSlide();
  s.background = { color: C.paper };
  addTitle(s, "Two curation tracks converge on one controlled folder", "They share identity gates and file standards—but their overlap analyses remain separate.", 3, "CURATION PIPELINE");
  const tracks = [
    { y: 1.92, tag: "FAMILY 1", title: "Primary curation sets", color: C.coral, steps: ["Paper-defined evidence", "Apply four distinct rules", "Resolve FlyBase identity", "Publish sets 01–04"] },
    { y: 3.72, tag: "FAMILY 2", title: "GO term–grouped sets", color: C.teal, steps: ["CSW FC0.5 + FC1 tables", "DAVID · FDR ≤10%", "Review term keywords", "Publish sets 05–07"] },
  ];
  tracks.forEach((t) => {
    label(s, t.tag, 0.7, t.y, 1.0, t.color, t.color === C.coral ? C.redPale : C.mint);
    s.addText(t.title, { x: 1.88, y: t.y + 0.02, w: 2.0, h: 0.28, margin: 0, fontSize: 14.5, bold: true, color: C.ink });
    t.steps.forEach((step, i) => {
      const x = 0.72 + i * 3.08;
      card(s, x, t.y + 0.48, 2.64, 0.92, { fill: C.white, line: "DFDDD4" });
      s.addShape(pptx.ShapeType.ellipse, { x: x + 0.16, y: t.y + 0.70, w: 0.42, h: 0.42, fill: { color: t.color }, line: { color: t.color } });
      s.addText(String(i + 1), { x: x + 0.16, y: t.y + 0.80, w: 0.42, h: 0.14, margin: 0, fontSize: 9.5, bold: true, align: "center", color: C.white });
      s.addText(step, { x: x + 0.72, y: t.y + 0.62, w: 1.68, h: 0.56, margin: 0, fontSize: 13.5, bold: i === 3, color: C.ink, fit: "shrink", valign: "mid" });
      if (i < 3) s.addShape(pptx.ShapeType.chevron, { x: x + 2.73, y: t.y + 0.79, w: 0.22, h: 0.24, fill: { color: C.line }, line: { color: C.line } });
    });
  });
  card(s, 0.7, 5.75, 11.94, 0.92, { fill: C.navy, line: C.navy });
  const bottom = [["902", "source observations"], ["791", "unique source FBgns"], ["≥2", "minimum pathway frequency"], ["7", "experiments"]];
  bottom.forEach((d, i) => {
    const x = 0.98 + i * 2.9;
    s.addText(d[0], { x, y: 5.94, w: 0.9, h: 0.34, margin: 0, fontSize: 22, bold: true, color: i === 2 ? C.gold : C.white });
    s.addText(d[1], { x: x + 0.82, y: 5.98, w: 1.78, h: 0.24, margin: 0, fontSize: 11, color: "C7D4D8", fit: "shrink" });
  });
  addFooter(s, "Thresholds remain separate: FC0.5 and FC1 are never collapsed into a max-frequency field.");
}

// 4 — Mechanistic + homeostatic
{
  const s = pptx.addSlide();
  s.background = { color: C.paper };
  addTitle(s, "Mechanistic and Homeostatic are different evidence categories", "The six-gene functional set is distinct from the 20-gene opposing history–rebound correlation set.", 4, "FAMILY 1  ·  SETS 01–02");
  card(s, 0.68, 1.85, 5.82, 4.85, { fill: C.white, shadow: true });
  label(s, "SET 01  ·  6 GENES", 0.98, 2.12, 1.72, C.coral, C.redPale);
  s.addText("Mechanistic", { x: 0.98, y: 2.62, w: 2.6, h: 0.42, margin: 0, fontSize: 24, bold: true, color: C.ink });
  s.addText("Consistent or correlated sleep/wake genes with published effects in a direction consistent with potential homeostatic function.", {
    x: 0.98, y: 3.12, w: 4.85, h: 0.66, margin: 0, fontSize: 15, color: C.muted, fit: "shrink",
  });
  bulletList(s, [
    "Expression is consistent across conditions or correlated with sleep/wake history.",
    "Published perturbation affects sleep/wake in a direction compatible with negative feedback.",
    "The paper highlights six such genes.",
  ], 0.98, 3.92, 4.85, 1.64, 14);
  s.addText("unc79  ·  SIFa  ·  rumpel  ·  AstA-R2  ·  Trhn  ·  RpL23", { x: 0.98, y: 5.90, w: 4.85, h: 0.36, margin: 0, fontSize: 13.5, bold: true, color: C.coral, fit: "shrink" });

  card(s, 6.82, 1.85, 5.82, 4.85, { fill: C.white, shadow: true });
  label(s, "SET 02  ·  20 GENES", 7.12, 2.12, 1.78, C.gold, "F4E8C9");
  s.addText("Homeostatic genes", { x: 7.12, y: 2.62, w: 4.55, h: 0.42, margin: 0, fontSize: 24, bold: true, color: C.ink, fit: "shrink" });
  metric(s, "15", "History+ / Rebound−", 7.12, 3.28, 2.15, C.gold);
  metric(s, "5", "History− / Rebound+", 9.72, 3.28, 2.15, C.coral);
  s.addShape(pptx.ShapeType.rect, { x: 7.12, y: 4.45, w: 4.95, h: 0.06, fill: { color: C.pale }, line: { color: C.pale } });
  s.addText("Genes correlate with sleep/wake history before sampling and with rebound in the opposite direction afterward.", {
    x: 7.12, y: 4.72, w: 4.85, h: 0.52, margin: 0, fontSize: 14.5, color: C.muted, fit: "shrink",
  });
  s.addText("Identity correction", { x: 7.12, y: 5.38, w: 2.1, h: 0.26, margin: 0, fontSize: 14, bold: true, color: C.ink });
  s.addText("Source row trh / FBgn0262139 is resolved from publication context to Trhn / FBgn0035187 (tryptophan hydroxylase). The wrong FBgn is never retained.", {
    x: 7.12, y: 5.72, w: 4.85, h: 0.66, margin: 0, fontSize: 13.5, color: C.muted, fit: "shrink",
  });
  addFooter(s);
}

// 5 — CSW + HLH
{
  const s = pptx.addSlide();
  s.background = { color: C.paper };
  addTitle(s, "CSW consistency and HLH regulation are separate categories", "Set 03 measures recurrence across seven datasets; Set 04 identifies potential upstream regulators.", 5, "FAMILY 1  ·  SETS 03–04");
  card(s, 0.68, 1.86, 7.45, 4.84, { fill: C.white, shadow: true });
  label(s, "SET 03  ·  97 GENES", 0.98, 2.13, 1.78, C.teal, C.mint);
  s.addText("CSW 4+ genes", { x: 0.98, y: 2.62, w: 4.9, h: 0.42, margin: 0, fontSize: 24, bold: true, color: C.ink });
  metric(s, "97", "Current output · FC0.5 wake frequency 4–6", 1.0, 3.33, 2.6, C.teal);
  metric(s, "7", "also qualify at FC1", 3.82, 3.33, 2.15, C.coral);
  metric(s, "0", "sleep genes at frequency 4", 6.1, 3.33, 1.5, C.muted);
  s.addText("Category rule", { x: 0.98, y: 4.78, w: 1.4, h: 0.26, margin: 0, fontSize: 14, bold: true, color: C.ink });
  bulletList(s, [
    "Consistent sleep/wake genes evident across four or more of seven datasets.",
    "The current output contains 97 FC0.5 wake genes; seven also qualify at FC1.",
    "No sleep-direction gene reaches frequency 4 in this run.",
  ], 0.98, 5.18, 6.42, 1.22, 13.5);

  card(s, 8.43, 1.86, 4.21, 4.84, { fill: C.navy, line: C.navy, shadow: true });
  label(s, "SET 04  ·  4 GENES", 8.76, 2.13, 1.68, C.white, C.coral);
  s.addText("HLH genes", { x: 8.76, y: 2.64, w: 3.1, h: 0.42, margin: 0, fontSize: 24, bold: true, color: C.white });
  s.addText("Named in Results / Table S1 after HOMER identified a CAGCTG E-box upstream of wake-induced genes.", {
    x: 8.76, y: 3.20, w: 3.25, h: 0.9, margin: 0, fontSize: 14.5, color: "CFDADD", fit: "shrink",
  });
  const genes = [["bigmax", "FBgn0039509"], ["HLH3B", "FBgn0011276"], ["E(spl)m3-HLH", "FBgn0002609"], ["E(spl)mbeta-HLH", "FBgn0002733"]];
  genes.forEach((g, i) => {
    s.addText(g[0], { x: 8.76, y: 4.40 + i * 0.43, w: 1.7, h: 0.22, margin: 0, fontSize: 13.5, bold: true, color: C.white, fit: "shrink" });
    s.addText(g[1], { x: 10.55, y: 4.40 + i * 0.43, w: 1.35, h: 0.22, margin: 0, fontSize: 11.5, color: "BFE1DC" });
  });
  addFooter(s);
}

// 6 — Pathways
{
  const s = pptx.addSlide();
  s.background = { color: C.paper };
  addTitle(s, "CSW-derived GO term–grouped genes retain threshold provenance", "Each gene records its source as FC0.5, FC1, or both; only hits from approved FDR-passing terms enter sets 05–07.", 6, "FAMILY 2  ·  THRESHOLD PROVENANCE");
  const stages = [["1", "Submit four CSW masters"], ["2", "Keep FDR ≤10% terms"], ["3", "Review term keywords"], ["4", "Take hit FBgns"]];
  stages.forEach((d, i) => {
    const x = 0.72 + i * 2.95;
    s.addShape(pptx.ShapeType.ellipse, { x, y: 1.84, w: 0.48, h: 0.48, fill: { color: C.teal }, line: { color: C.teal } });
    s.addText(d[0], { x, y: 1.96, w: 0.48, h: 0.14, margin: 0, align: "center", fontSize: 10, bold: true, color: C.white });
    s.addText(d[1], { x: x + 0.62, y: 1.90, w: 1.92, h: 0.32, margin: 0, fontSize: 13.5, bold: true, color: C.ink, fit: "shrink" });
    if (i < 3) s.addShape(pptx.ShapeType.chevron, { x: x + 2.60, y: 1.96, w: 0.18, h: 0.20, fill: { color: C.line }, line: { color: C.line } });
  });
  card(s, 0.68, 2.72, 5.22, 3.62, { fill: C.white, shadow: true });
  s.addText("GO term–grouped outputs", {
    x: 1.0, y: 3.02, w: 3.0, h: 0.28, margin: 0,
    fontSize: 17, bold: true, color: C.ink,
  });
  const bars = [
    ["Mitochondrial / metabolism", 167, C.teal],
    ["Ribosomal / translation", 99, C.coral],
    ["Immune", 5, C.gold],
  ];
  bars.forEach((b, i) => {
    const y = 3.58 + i * 0.78;
    s.addText(b[0], { x: 1.0, y, w: 2.02, h: 0.28, margin: 0, fontSize: 12.5, bold: true, color: C.ink, fit: "shrink" });
    s.addShape(pptx.ShapeType.roundRect, { x: 3.25, y: y + 0.03, w: 1.93, h: 0.28, rectRadius: 0.06, fill: { color: C.pale }, line: { color: C.pale } });
    s.addShape(pptx.ShapeType.roundRect, { x: 3.25, y: y + 0.03, w: Math.max(0.15, (b[1] / 167) * 1.93), h: 0.28, rectRadius: 0.06, fill: { color: b[2] }, line: { color: b[2] } });
    s.addText(String(b[1]), { x: 5.27, y: y - 0.01, w: 0.42, h: 0.28, margin: 0, fontSize: 14, bold: true, align: "right", color: b[2] });
  });
  s.addText("All three outputs are wake-direction genes.", { x: 1.0, y: 5.98, w: 4.2, h: 0.25, margin: 0, fontSize: 13.5, bold: true, color: C.teal });

  card(s, 6.18, 2.72, 6.46, 3.62, { fill: C.navy, line: C.navy, shadow: true });
  s.addText("Threshold provenance · 240 unique FBgns", {
    x: 6.52, y: 3.02, w: 4.9, h: 0.30, margin: 0,
    fontSize: 17, bold: true, color: C.white,
  });
  const thresholdStats = [
    ["203", "FC0.5 only", "Ribosomal and mitochondrial dominate", C.teal2],
    ["5", "FC1 only", "BomBc2 · BomS1 · BomS2 · IM14 · Thor", C.coral],
    ["32", "both thresholds", "FC1 ribosomal/mitochondrial genes are also at FC0.5", C.gold],
  ];
  thresholdStats.forEach((d, i) => {
    const y = 3.58 + i * 0.78;
    s.addText(d[0], { x: 6.52, y, w: 0.95, h: 0.40, margin: 0, fontSize: 25, bold: true, color: d[3] });
    s.addText(d[1], { x: 7.48, y: y + 0.02, w: 1.55, h: 0.25, margin: 0, fontSize: 13.5, bold: true, color: C.white });
    s.addText(d[2], { x: 9.14, y: y + 0.02, w: 2.92, h: 0.34, margin: 0, fontSize: 11.5, color: "C7D4D8", fit: "shrink" });
  });
  s.addShape(pptx.ShapeType.rect, { x: 6.52, y: 5.93, w: 5.52, h: 0.04, fill: { color: "35505E" }, line: { color: "35505E" } });
  s.addText("Key takeaway", { x: 6.52, y: 6.06, w: 1.25, h: 0.20, margin: 0, fontSize: 11, bold: true, color: C.teal2 });
  s.addText("Only the five immune genes are FC1-exclusive.", {
    x: 7.75, y: 6.04, w: 4.28, h: 0.23, margin: 0,
    fontSize: 13.5, bold: true, color: C.white, fit: "shrink",
  });
  addFooter(s);
}

// 7 — Identity QA
{
  const s = pptx.addSlide();
  s.background = { color: C.paper };
  addTitle(s, "Identity resolution is conservative by design", "Current symbols can pass directly; ambiguous synonyms require biological context or a case-supported source match.", 7, "FLYBASE AUDIT");
  const rows = [
    { y: 1.94, symbol: "Trhn", id: "FBgn0035187", status: "PASS", color: C.teal, fill: C.mint, text: "Current symbol · tryptophan hydroxylase neuronal · publication-resolved in Homeostatic." },
    { y: 3.26, symbol: "trh", id: "FBgn0262139", status: "REJECT", color: C.coral, fill: C.redPale, text: "Current symbol · trachealess · must not enter any generated identity field for this publication row." },
    { y: 4.58, symbol: "Trh", id: "ambiguous", status: "REVIEW", color: C.gold, fill: "F4E8C9", text: "Synonym of Hn, Trhn, and trh · never auto-approved without biological context." },
  ];
  rows.forEach((r) => {
    card(s, 0.68, r.y, 8.15, 1.08, { fill: C.white, line: "DFDDD4", shadow: true });
    label(s, r.status, 0.96, r.y + 0.36, 0.88, r.color, r.fill);
    s.addText(r.symbol, { x: 2.10, y: r.y + 0.22, w: 1.15, h: 0.34, margin: 0, fontSize: 21, bold: true, color: C.ink });
    s.addText(r.id, { x: 3.35, y: r.y + 0.27, w: 1.45, h: 0.26, margin: 0, fontSize: 12.5, bold: true, color: r.color });
    s.addText(r.text, { x: 4.98, y: r.y + 0.24, w: 3.42, h: 0.56, margin: 0, fontSize: 13.5, color: C.muted, fit: "shrink" });
  });
  card(s, 9.18, 1.94, 3.46, 3.72, { fill: C.navy, line: C.navy, shadow: true });
  s.addText("Approval gates", { x: 9.50, y: 2.27, w: 2.6, h: 0.38, margin: 0, fontSize: 22, bold: true, color: C.white });
  bulletList(s, [
    "Exact current symbol",
    "Publication-supported correction",
    "Unique case-insensitive current-symbol match",
    "Otherwise: stop for review",
  ], 9.50, 2.95, 2.60, 1.80, 14.5, C.white);
  s.addText("Membership key", { x: 9.50, y: 5.02, w: 1.60, h: 0.22, margin: 0, fontSize: 11.5, color: C.teal2 });
  s.addText("flybase_gene_id", { x: 9.50, y: 5.31, w: 2.20, h: 0.28, margin: 0, fontSize: 14.5, bold: true, color: C.white });
  card(s, 0.68, 6.03, 11.96, 0.62, { fill: "EEF2EF", line: "D9E1DC" });
  s.addText("Symbols are display labels. All overlap, deduplication, and set membership use FBgn only.", {
    x: 0.98, y: 6.22, w: 11.3, h: 0.22, margin: 0, fontSize: 14.5, bold: true, align: "center", color: C.teal,
  });
  addFooter(s);
}

// 8 — Family 1 overlap
{
  const s = pptx.addSlide();
  s.background = { color: C.paper };
  addTitle(s, "Primary-set overlap is limited to three genes", "Mechanistic, Homeostatic genes, CSW 4+ genes, and HLH genes are compared only with one another.", 8, "FAMILY 1  ·  PRIMARY-SET OVERLAP");
  card(s, 0.68, 1.86, 7.72, 4.86, { fill: C.white, shadow: true });
  const primarySets = [
    ["Mechanistic", "6", C.coral],
    ["Homeostatic genes", "20", C.gold],
    ["CSW 4+ genes", "97", C.teal],
    ["HLH genes", "4", "665C86"],
  ];
  primarySets.forEach((d, i) => {
    const x = 0.98 + i * 1.78;
    card(s, x, 2.18, 1.52, 0.88, { fill: "F8F7F2", line: "E1DED4" });
    s.addText(d[1], { x: x + 0.16, y: 2.35, w: 0.55, h: 0.34, margin: 0, fontSize: 22, bold: true, color: d[2] });
    s.addText(d[0], { x: x + 0.16, y: 2.70, w: 1.18, h: 0.25, margin: 0, fontSize: 11, bold: true, color: C.ink, fit: "shrink" });
  });
  s.addText("Observed shared genes", { x: 0.98, y: 3.46, w: 2.6, h: 0.28, margin: 0, fontSize: 17, bold: true, color: C.ink });
  const shared = [
    ["Mechanistic  ×  Homeostatic genes", "2", "Trhn  ·  AstA-R2", C.coral],
    ["Mechanistic  ×  CSW 4+ genes", "1", "RpL23", C.teal],
    ["All other pairings", "0", "No shared FBgns", C.muted],
  ];
  shared.forEach((d, i) => {
    const y = 3.98 + i * 0.72;
    s.addShape(pptx.ShapeType.rect, { x: 0.98, y, w: 0.08, h: 0.50, fill: { color: d[3] }, line: { color: d[3] } });
    s.addText(d[0], { x: 1.25, y: y + 0.02, w: 3.25, h: 0.24, margin: 0, fontSize: 13.5, bold: true, color: C.ink, fit: "shrink" });
    s.addText(d[1], { x: 4.72, y: y - 0.01, w: 0.55, h: 0.36, margin: 0, fontSize: 22, bold: true, color: d[3] });
    s.addText(d[2], { x: 5.40, y: y + 0.04, w: 2.35, h: 0.25, margin: 0, fontSize: 13.5, bold: i < 2, color: i < 2 ? C.ink : C.muted, fit: "shrink" });
  });
  card(s, 8.72, 1.86, 3.92, 4.86, { fill: C.navy, line: C.navy, shadow: true });
  s.addText("Within Family 1", { x: 9.08, y: 2.18, w: 2.9, h: 0.36, margin: 0, fontSize: 22, bold: true, color: C.white });
  const family1Stats = [["124", "unique FBgns"], ["121", "in exactly one set"], ["3", "in exactly two sets"], ["0", "in three or four sets"]];
  family1Stats.forEach((d, i) => {
    const y = 2.92 + i * 0.72;
    s.addText(d[0], { x: 9.08, y, w: 0.95, h: 0.40, margin: 0, fontSize: 25, bold: true, color: i === 0 ? C.teal2 : C.coral });
    s.addText(d[1], { x: 10.14, y: y + 0.06, w: 1.92, h: 0.26, margin: 0, fontSize: 13, color: "D4DFE2", fit: "shrink" });
  });
  s.addShape(pptx.ShapeType.rect, { x: 9.08, y: 5.92, w: 2.92, h: 0.04, fill: { color: "35505E" }, line: { color: "35505E" } });
  s.addText("HLH is fully disjoint from the other primary sets.", {
    x: 9.08, y: 6.12, w: 2.92, h: 0.38, margin: 0,
    fontSize: 13.5, bold: true, color: C.white, fit: "shrink",
  });
  addFooter(s);
}

// 9 — Family 2 overlap
{
  const s = pptx.addSlide();
  s.background = { color: C.paper };
  addTitle(s, "GO term–grouped overlap is confined to ribosomal × mitochondrial", "The three pathway buckets are compared only within Family 2; threshold provenance remains a separate attribute.", 9, "FAMILY 2  ·  GO OVERLAP");
  card(s, 0.68, 1.86, 8.00, 4.86, { fill: C.white, shadow: true });
  s.addText("GO term–group membership", { x: 1.02, y: 2.18, w: 3.2, h: 0.30, margin: 0, fontSize: 17, bold: true, color: C.ink });
  s.addText("Schematic · counts shown", { x: 6.40, y: 2.22, w: 1.82, h: 0.20, margin: 0, fontSize: 10.5, align: "right", color: C.muted });
  s.addShape(pptx.ShapeType.ellipse, {
    x: 1.22, y: 2.74, w: 3.18, h: 3.18,
    fill: { color: C.coral, transparency: 24 }, line: { color: C.coral, width: 1.4 },
  });
  s.addShape(pptx.ShapeType.ellipse, {
    x: 3.48, y: 2.74, w: 3.18, h: 3.18,
    fill: { color: C.teal, transparency: 24 }, line: { color: C.teal, width: 1.4 },
  });
  s.addShape(pptx.ShapeType.ellipse, {
    x: 6.92, y: 3.48, w: 1.10, h: 1.10,
    fill: { color: C.gold, transparency: 12 }, line: { color: C.gold, width: 1.2 },
  });
  s.addText("68", { x: 1.84, y: 3.89, w: 0.70, h: 0.46, margin: 0, fontSize: 28, bold: true, align: "center", color: C.ink });
  s.addText("31", { x: 3.57, y: 3.89, w: 0.72, h: 0.46, margin: 0, fontSize: 28, bold: true, align: "center", color: C.ink });
  s.addText("136", { x: 5.27, y: 3.89, w: 0.78, h: 0.46, margin: 0, fontSize: 28, bold: true, align: "center", color: C.ink });
  s.addText("5", { x: 7.16, y: 3.82, w: 0.62, h: 0.36, margin: 0, fontSize: 22, bold: true, align: "center", color: C.ink });
  s.addText("Ribosomal\n99", { x: 1.64, y: 5.85, w: 1.35, h: 0.50, margin: 0, fontSize: 12.5, bold: true, align: "center", color: C.coral });
  s.addText("Mitochondrial\n167", { x: 4.92, y: 5.85, w: 1.50, h: 0.50, margin: 0, fontSize: 12.5, bold: true, align: "center", color: C.teal });
  s.addText("Immune\n5", { x: 6.95, y: 4.72, w: 1.05, h: 0.48, margin: 0, fontSize: 11.5, bold: true, align: "center", color: C.gold });
  s.addText("240 unique FBgns in Family 2", {
    x: 1.02, y: 6.35, w: 7.25, h: 0.24, margin: 0,
    fontSize: 13.5, bold: true, align: "center", color: C.teal,
  });
  card(s, 9.00, 1.86, 3.64, 4.86, { fill: C.navy, line: C.navy, shadow: true });
  s.addText("What the overlap means", { x: 9.34, y: 2.18, w: 2.66, h: 0.58, margin: 0, fontSize: 21, bold: true, color: C.white, fit: "shrink" });
  const goTakeaways = [
    ["31", "genes occur in both ribosomal and mitochondrial buckets"],
    ["5", "immune genes are fully disjoint"],
    ["0", "immune overlaps with either other GO bucket"],
  ];
  goTakeaways.forEach((d, i) => {
    const y = 3.02 + i * 0.88;
    s.addText(d[0], { x: 9.34, y, w: 0.70, h: 0.40, margin: 0, fontSize: 25, bold: true, color: i === 0 ? C.coral : C.teal2 });
    s.addText(d[1], { x: 10.12, y: y + 0.02, w: 1.92, h: 0.46, margin: 0, fontSize: 12.5, color: "D4DFE2", fit: "shrink" });
  });
  s.addShape(pptx.ShapeType.rect, { x: 9.34, y: 5.78, w: 2.66, h: 0.04, fill: { color: "35505E" }, line: { color: "35505E" } });
  s.addText("The 31-gene overlap is intentional: mitochondrial translation terms match both keyword groups.", {
    x: 9.34, y: 5.98, w: 2.66, h: 0.54, margin: 0,
    fontSize: 12.5, bold: true, color: C.white, fit: "shrink",
  });
  addFooter(s, "No cross-family overlap is shown: Family 1 and Family 2 answer different curation questions.");
}

// 10 — Handoff
{
  const s = pptx.addSlide();
  s.background = { color: C.navy };
  s.addText("STOCKER HANDOFF", {
    x: 0.72, y: 0.34, w: 3.2, h: 0.24, margin: 0,
    fontSize: 11, bold: true, charSpacing: 2.0, color: C.teal2,
  });
  s.addText("The contract is simple—and strict.", {
    x: 0.72, y: 0.88, w: 7.4, h: 0.62, margin: 0,
    fontFace: FONT, fontSize: 32, bold: true, color: C.white,
  });
  const rules = [
    ["7", "Exactly seven CSVs in the discovery folder"],
    ["2", "Required keys: ext_gene and flybase_gene_id"],
    ["0", "No overlap, GO, or audit CSVs beside the inputs"],
  ];
  rules.forEach((r, i) => {
    const x = 0.72 + i * 3.03;
    card(s, x, 1.82, 2.72, 1.38, { fill: i === 0 ? C.coral : "1A3442", line: i === 0 ? C.coral : "35505E" });
    s.addText(r[0], { x: x + 0.24, y: 2.07, w: 0.62, h: 0.52, margin: 0, fontSize: 30, bold: true, color: C.white });
    s.addText(r[1], { x: x + 0.88, y: 2.04, w: 1.52, h: 0.68, margin: 0, fontSize: 12.5, color: "D6E0E3", fit: "shrink" });
  });
  card(s, 10.10, 0.88, 2.52, 4.62, { fill: C.cream, line: C.cream, shadow: true });
  s.addText("Seven inputs", { x: 10.42, y: 1.20, w: 1.85, h: 0.32, margin: 0, fontSize: 20, bold: true, color: C.ink });
  const names = ["01  Mechanistic", "02  Homeostatic genes", "03  CSW 4+ genes", "04  HLH genes", "05  Ribosomal", "06  Mitochondrial", "07  Immune"];
  names.forEach((n, i) => {
    s.addShape(pptx.ShapeType.ellipse, { x: 10.42, y: 1.80 + i * 0.46, w: 0.18, h: 0.18, fill: { color: i < 4 ? C.coral : C.teal }, line: { color: i < 4 ? C.coral : C.teal } });
    s.addText(n, { x: 10.74, y: 1.77 + i * 0.46, w: 1.44, h: 0.24, margin: 0, fontSize: 12.5, bold: true, color: C.ink, fit: "shrink" });
  });
  s.addText("RUN", { x: 0.72, y: 3.66, w: 1.0, h: 0.24, margin: 0, fontSize: 11, bold: true, charSpacing: 1.6, color: C.teal2 });
  s.addShape(pptx.ShapeType.roundRect, { x: 0.72, y: 4.02, w: 8.93, h: 1.18, rectRadius: 0.06, fill: { color: "0A1821" }, line: { color: "35505E", width: 0.8 } });
  s.addText("python -m fl_ai_reagent_stocker run ./data/gene_sets/Tx-Omics_Revision/ \\\n  --config ./data/config/stock_split_config_priority_example.json", {
    x: 1.02, y: 4.35, w: 8.30, h: 0.55, margin: 0,
    fontFace: FONT, fontSize: 12.5, color: "D7E3E6", breakLine: false, fit: "shrink",
  });
  s.addText("Stability rule", { x: 0.72, y: 5.72, w: 1.55, h: 0.24, margin: 0, fontSize: 12, bold: true, color: C.coral });
  s.addText("Discovery is recursive: every extra *.csv changes the run. Keep analysis outputs elsewhere.", {
    x: 2.25, y: 5.70, w: 7.4, h: 0.30, margin: 0, fontSize: 14, color: "D6E0E3", fit: "shrink",
  });
  s.addText("Rosensweig–Shah 2026  ·  Tx-Omics Revision", { x: 0.72, y: 7.06, w: 4.2, h: 0.18, margin: 0, fontSize: 9.5, color: "8FA6AF" });
}

pptx.writeFile({ fileName: OUT });

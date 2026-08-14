const pptxgen = require("../pptx_tool/node_modules/pptxgenjs");
const path = require("path");

const pptx = new pptxgen();
pptx.defineLayout({ name: "CUSTOM_WIDE", width: 13.333, height: 7.5 });
pptx.layout = "CUSTOM_WIDE";
pptx.author = "fl.ai";
pptx.title = "Tx-Omics Revision — Gene-set curation";
pptx.subject = "Rosensweig–Shah 2026 gene-set curation";
pptx.lang = "en-US";
const FONT = "Avenir Next";
pptx.theme = { headFontFace: FONT, bodyFontFace: FONT, lang: "en-US" };

const OUT = path.resolve(__dirname, "Tx-Omics_Revision_gene_set_curation_revised.pptx");
const FIG = path.resolve(__dirname, "../../audit_outputs/Tx-Omics_Revision/Figures");
const C = {
  ink: "17212B",
  navy: "122634",
  teal: "1D6F69",
  tealLight: "67AFA6",
  coral: "DC6B4D",
  gold: "D8A238",
  purple: "665C86",
  red: "A44A3F",
  white: "FFFFFF",
  paper: "FFFCF6",
  cream: "F7F3EA",
  muted: "5E6971",
  line: "D8D9D3",
  pale: "ECEBE4",
  mint: "DCEDE8",
  coralPale: "F5DED6",
  goldPale: "F4E8C9",
};

function card(slide, x, y, w, h, fill = C.white, line = C.line) {
  slide.addShape(pptx.ShapeType.roundRect, {
    x, y, w, h, rectRadius: 0.07,
    fill: { color: fill },
    line: { color: line, width: 0.8 },
  });
}

function header(slide, section, title, subtitle, page) {
  slide.background = { color: C.paper };
  slide.addText(section, {
    x: 0.72, y: 0.30, w: 5.2, h: 0.22, margin: 0,
    fontFace: FONT, fontSize: 10.5, bold: true, charSpacing: 1.6, color: C.teal,
  });
  slide.addText(String(page).padStart(2, "0"), {
    x: 12.0, y: 0.29, w: 0.58, h: 0.22, margin: 0,
    fontFace: FONT, fontSize: 10.5, bold: true, align: "right", color: C.muted,
  });
  slide.addText(title, {
    x: 0.72, y: 0.72, w: 11.85, h: 0.55, margin: 0,
    fontFace: FONT, fontSize: 28, bold: true, color: C.ink, fit: "shrink",
  });
  slide.addText(subtitle, {
    x: 0.72, y: 1.34, w: 11.85, h: 0.38, margin: 0,
    fontFace: FONT, fontSize: 15.5, color: C.muted, fit: "shrink",
  });
}

function pill(slide, text, x, y, w, color, fill) {
  slide.addShape(pptx.ShapeType.roundRect, {
    x, y, w, h: 0.34, rectRadius: 0.08,
    fill: { color: fill }, line: { color: fill },
  });
  slide.addText(text, {
    x: x + 0.10, y: y + 0.08, w: w - 0.20, h: 0.16, margin: 0,
    fontFace: FONT, fontSize: 9.5, bold: true, charSpacing: 1.0,
    align: "center", color,
  });
}

function bigStat(slide, value, label, x, y, w, color = C.teal, note = "") {
  slide.addText(String(value), {
    x, y, w, h: 0.62, margin: 0,
    fontFace: FONT, fontSize: 38, bold: true, color,
  });
  slide.addText(label, {
    x, y: y + 0.66, w, h: 0.30, margin: 0,
    fontFace: FONT, fontSize: 14, bold: true, color: C.white, fit: "shrink",
  });
  if (note) {
    slide.addText(note, {
      x, y: y + 1.00, w, h: 0.45, margin: 0,
      fontFace: FONT, fontSize: 13.5, color: "D5E0E3", fit: "shrink",
    });
  }
}

function bulletList(slide, items, x, y, w, h, fontSize = 15, color = C.ink) {
  const runs = items.map((item, i) => ({
    text: item,
    options: { bullet: { indent: fontSize * 1.2 }, breakLine: i < items.length - 1, paraSpaceAfterPt: 9 },
  }));
  slide.addText(runs, {
    x, y, w, h, margin: 0.03,
    fontFace: FONT, fontSize, color, valign: "top", fit: "shrink",
  });
}

function termBar(slide, name, value, max, x, y, w, color) {
  const labelW = 3.12;
  slide.addText(name, {
    x, y, w: labelW, h: 0.28, margin: 0,
    fontFace: FONT, fontSize: 14, bold: true, color: C.ink, fit: "shrink",
  });
  slide.addShape(pptx.ShapeType.roundRect, {
    x: x + labelW + 0.20, y: y + 0.02, w: w - labelW - 0.72, h: 0.27,
    rectRadius: 0.05, fill: { color: C.pale }, line: { color: C.pale },
  });
  slide.addShape(pptx.ShapeType.roundRect, {
    x: x + labelW + 0.20, y: y + 0.02,
    w: Math.max(0.10, (value / max) * (w - labelW - 0.72)), h: 0.27,
    rectRadius: 0.05, fill: { color }, line: { color },
  });
  slide.addText(String(value), {
    x: x + w - 0.40, y: y - 0.01, w: 0.40, h: 0.28, margin: 0,
    fontFace: FONT, fontSize: 13.5, bold: true, align: "right", color,
  });
}

function addThresholdPanel(slide, total, fc05, fc1, both, color, equation) {
  card(slide, 0.72, 1.98, 4.20, 4.78, C.navy, C.navy);
  pill(slide, "GO ASSIGNMENT SOURCE", 1.02, 2.28, 2.16, C.white, color);
  slide.addText(String(total), {
    x: 1.02, y: 2.86, w: 1.65, h: 0.78, margin: 0,
    fontFace: FONT, fontSize: 48, bold: true, color: C.white,
  });
  slide.addText("unique genes in this set", {
    x: 1.02, y: 3.67, w: 2.5, h: 0.28, margin: 0,
    fontFace: FONT, fontSize: 14, color: "D5E0E3",
  });
  const rows = [
    ["FC0.5 GO", fc05, "genes assigned"],
    ["FC1 GO", fc1, "genes assigned"],
    ["BOTH", both, "assigned by both reports"],
  ];
  rows.forEach((r, i) => {
    const y = 4.18 + i * 0.62;
    slide.addText(r[0], {
      x: 1.02, y, w: 0.92, h: 0.23, margin: 0,
      fontFace: FONT, fontSize: 10.5, bold: true, color: C.tealLight, fit: "shrink",
    });
    slide.addText(String(r[1]), {
      x: 2.02, y: y - 0.07, w: 0.65, h: 0.38, margin: 0,
      fontFace: FONT, fontSize: 23, bold: true, color: i === 1 ? C.coral : C.white,
    });
    slide.addText(r[2], {
      x: 2.78, y: y + 0.01, w: 1.46, h: 0.25, margin: 0,
      fontFace: FONT, fontSize: 12.5, color: "D5E0E3", fit: "shrink",
    });
  });
  slide.addShape(pptx.ShapeType.rect, {
    x: 1.02, y: 6.12, w: 3.55, h: 0.04,
    fill: { color: "35505E" }, line: { color: "35505E" },
  });
  slide.addText(equation, {
    x: 1.02, y: 6.31, w: 3.55, h: 0.30, margin: 0,
    fontFace: FONT, fontSize: 13, bold: true, color: C.white, fit: "shrink",
  });
}

function addExperimentDistributionSlide(
  page,
  setNumber,
  title,
  n05,
  counts05,
  n1,
  counts1,
  takeaway,
) {
  const s = pptx.addSlide();
  header(
    s,
    `FAMILY 2  ·  SET ${setNumber}  ·  EXPERIMENTS`,
    `${title} · experiment-level regulation`,
    "Within each threshold-responsive subset, bars show genes with logFC > threshold and q < 0.05 per experiment; genes may appear in multiple bars.",
    page,
  );
  const experiments = [
    "Baseline",
    "MechSD 3h",
    "MechSD 6h",
    "MechSD 12h",
    "R85C10>TrpA1",
    "R23E10>ChR",
    "THIP",
  ];
  const panels = [
    { x: 0.72, title: `FC0.5  ·  N=${n05}`, n: n05, counts: counts05, color: C.teal },
    { x: 6.82, title: `FC1  ·  N=${n1}`, n: n1, counts: counts1, color: C.coral },
  ];
  panels.forEach((panel) => {
    card(s, panel.x, 1.98, 5.80, 4.42, C.white, "E2DED4");
    pill(s, panel.title.toUpperCase(), panel.x + 0.30, 2.26, 3.20, C.white, panel.color);
    experiments.forEach((experiment, i) => {
      const y = 2.90 + i * 0.47;
      const count = panel.counts[i];
      const pct = panel.n ? Math.round((count / panel.n) * 100) : 0;
      s.addText(experiment, {
        x: panel.x + 0.30, y, w: 1.36, h: 0.24, margin: 0,
        fontFace: FONT, fontSize: 12.5, bold: true, color: C.ink, fit: "shrink",
      });
      s.addShape(pptx.ShapeType.roundRect, {
        x: panel.x + 1.74, y: y + 0.02, w: 2.54, h: 0.22,
        rectRadius: 0.04, fill: { color: C.pale }, line: { color: C.pale },
      });
      if (pct > 0) {
        s.addShape(pptx.ShapeType.roundRect, {
          x: panel.x + 1.74, y: y + 0.02, w: (pct / 100) * 2.54, h: 0.22,
          rectRadius: 0.04, fill: { color: panel.color }, line: { color: panel.color },
        });
      }
      s.addText(`${count} · ${pct}%`, {
        x: panel.x + 4.42, y: y - 0.01, w: 1.00, h: 0.26, margin: 0,
        fontFace: FONT, fontSize: 12, bold: true, align: "right", color: panel.color,
      });
    });
    s.addText("n = category genes meeting threshold in ≥1 experiment; independent of GO-assignment count", {
      x: panel.x + 0.30, y: 6.10, w: 5.18, h: 0.20, margin: 0,
      fontFace: FONT, fontSize: 11.5, color: C.muted, fit: "shrink",
    });
  });
  card(s, 0.72, 6.58, 11.90, 0.50, C.navy, C.navy);
  s.addText(takeaway, {
    x: 1.02, y: 6.72, w: 11.30, h: 0.21, margin: 0,
    fontFace: FONT, fontSize: 14, bold: true, align: "center", color: C.white, fit: "shrink",
  });
}

// 1 — Title
{
  const s = pptx.addSlide();
  s.background = { color: C.navy };
  s.addShape(pptx.ShapeType.rect, {
    x: 0, y: 0, w: 0.18, h: 7.5,
    fill: { color: C.coral }, line: { color: C.coral },
  });
  s.addText("TX-OMICS REVISION", {
    x: 0.82, y: 0.72, w: 4.4, h: 0.30, margin: 0,
    fontFace: FONT, fontSize: 13, bold: true, charSpacing: 2.8, color: C.tealLight,
  });
  s.addText("Seven gene sets.\nSeven clear definitions.", {
    x: 0.82, y: 1.48, w: 7.2, h: 1.55, margin: 0,
    fontFace: FONT, fontSize: 39, bold: true, color: C.white, fit: "shrink",
  });
  s.addText("Family 1 captures four primary evidence categories. Family 2 contains three CSW-derived GO term groups, with FC0.5 and FC1 provenance shown separately.", {
    x: 0.84, y: 3.38, w: 7.10, h: 1.10, margin: 0,
    fontFace: FONT, fontSize: 18, color: "D6E0E3", fit: "shrink",
  });
  const stats = [["4", "primary sets"], ["3", "GO-derived sets"], ["2", "thresholds"], ["7", "stocker CSVs"]];
  stats.forEach((d, i) => {
    const x = 8.65 + (i % 2) * 2.0;
    const y = 1.52 + Math.floor(i / 2) * 1.62;
    card(s, x, y, 1.72, 1.32, i === 0 ? C.coral : "1B3443", i === 0 ? C.coral : "35505E");
    s.addText(d[0], {
      x: x + 0.18, y: y + 0.18, w: 1.2, h: 0.50, margin: 0,
      fontFace: FONT, fontSize: 29, bold: true, color: C.white,
    });
    s.addText(d[1], {
      x: x + 0.18, y: y + 0.80, w: 1.30, h: 0.28, margin: 0,
      fontFace: FONT, fontSize: 13, color: "D6E0E3", fit: "shrink",
    });
  });
}

// 2 — Overview
{
  const s = pptx.addSlide();
  header(s, "CURATION MAP", "Two families, seven non-interchangeable sets", "Each set now receives one slide; overlap is reported separately within each family.", 2);
  card(s, 0.72, 1.98, 5.82, 4.72, C.white, "E2DED4");
  pill(s, "FAMILY 1  ·  PRIMARY EVIDENCE", 1.02, 2.28, 2.82, C.coral, C.coralPale);
  const primary = [
    ["01", "Mechanistic", "6", "Consistent/correlated + published homeostatic-consistent effect"],
    ["02", "Homeostatic genes", "20", "Opposing history versus rebound correlations"],
    ["03", "CSW 4+ genes", "97", "Consistent in four or more of seven datasets"],
    ["04", "HLH genes", "4", "Potential upstream regulators"],
  ];
  primary.forEach((d, i) => {
    const y = 2.93 + i * 0.78;
    s.addText(d[0], { x: 1.02, y, w: 0.40, h: 0.22, margin: 0, fontFace: FONT, fontSize: 10.5, bold: true, color: C.coral });
    s.addText(d[2], { x: 1.52, y: y - 0.08, w: 0.70, h: 0.38, margin: 0, fontFace: FONT, fontSize: 23, bold: true, color: C.ink });
    s.addText(d[1], { x: 2.28, y: y - 0.02, w: 1.72, h: 0.27, margin: 0, fontFace: FONT, fontSize: 14, bold: true, color: C.ink, fit: "shrink" });
    s.addText(d[3], { x: 4.02, y: y - 0.03, w: 2.10, h: 0.40, margin: 0, fontFace: FONT, fontSize: 12.5, color: C.muted, fit: "shrink" });
  });
  card(s, 6.82, 1.98, 5.80, 4.72, C.navy, C.navy);
  pill(s, "FAMILY 2  ·  CSW-DERIVED GO", 7.12, 2.28, 2.86, C.white, C.teal);
  const go = [
    ["05", "Ribosomal", "99", "FC0.5 99  ·  FC1 12"],
    ["06", "Mitochondrial", "167", "FC0.5 167  ·  FC1 23"],
    ["07", "Immune", "5", "FC0.5 0  ·  FC1 5"],
  ];
  go.forEach((d, i) => {
    const y = 3.02 + i * 1.00;
    s.addText(d[0], { x: 7.12, y, w: 0.40, h: 0.22, margin: 0, fontFace: FONT, fontSize: 10.5, bold: true, color: C.tealLight });
    s.addText(d[2], { x: 7.62, y: y - 0.10, w: 0.86, h: 0.43, margin: 0, fontFace: FONT, fontSize: 25, bold: true, color: C.white });
    s.addText(d[1], { x: 8.60, y: y - 0.03, w: 2.45, h: 0.28, margin: 0, fontFace: FONT, fontSize: 14, bold: true, color: C.white, fit: "shrink" });
    s.addText(d[3], { x: 10.96, y: y - 0.03, w: 1.40, h: 0.38, margin: 0, fontFace: FONT, fontSize: 12.5, color: "E2EAEC", fit: "shrink" });
  });
  s.addText("FC0.5 and FC1 counts are GO-report provenance; a gene may be assigned by both reports.", {
    x: 7.12, y: 6.04, w: 5.00, h: 0.40, margin: 0,
    fontFace: FONT, fontSize: 13, bold: true, color: C.tealLight, fit: "shrink",
  });
}

// 3 — Mechanistic
{
  const s = pptx.addSlide();
  header(s, "FAMILY 1  ·  SET 01", "Mechanistic · 6 genes", "Consistent or correlated sleep/wake genes with published effects compatible with homeostatic feedback.", 3);
  card(s, 0.72, 1.98, 4.10, 4.72, C.navy, C.navy);
  bigStat(s, "6", "paper-highlighted genes", 1.06, 2.38, 2.50, C.coral, "Inclusion requires transcriptomic evidence plus published sleep/wake function.");
  s.addShape(pptx.ShapeType.rect, { x: 1.06, y: 4.22, w: 3.38, h: 0.04, fill: { color: "35505E" }, line: { color: "35505E" } });
  bulletList(s, [
    "Consistent across sleep/wake conditions or correlated with sleep/wake history",
    "Published perturbation changes sleep/wake",
    "Effect direction fits a potential homeostatic loop",
  ], 1.06, 4.52, 3.36, 1.62, 14, C.white);
  const genes = [
    ["unc79", "NALCN complex"],
    ["SIFa", "neuropeptide"],
    ["rumpel", "glial transporter"],
    ["AstA-R2", "neuropeptide receptor"],
    ["Trhn", "serotonin synthesis"],
    ["RpL23", "ribosome / translation"],
  ];
  genes.forEach((d, i) => {
    const col = i % 2, row = Math.floor(i / 2);
    const x = 5.18 + col * 3.70, y = 2.12 + row * 1.35;
    card(s, x, y, 3.34, 1.06, C.white, "E2DED4");
    s.addText(d[0], { x: x + 0.25, y: y + 0.20, w: 2.75, h: 0.32, margin: 0, fontFace: FONT, fontSize: 20, bold: true, color: C.teal });
    s.addText(d[1], { x: x + 0.25, y: y + 0.62, w: 2.75, h: 0.23, margin: 0, fontFace: FONT, fontSize: 13.5, color: C.muted, fit: "shrink" });
  });
  card(s, 5.18, 6.20, 7.04, 0.50, C.coralPale, C.coralPale);
  s.addText("These six are not the 20 opposing history–rebound correlation genes.", {
    x: 5.48, y: 6.35, w: 6.44, h: 0.20, margin: 0,
    fontFace: FONT, fontSize: 14, bold: true, align: "center", color: C.coral,
  });
}

// 4 — Homeostatic
{
  const s = pptx.addSlide();
  header(s, "FAMILY 1  ·  SET 02", "Homeostatic genes · 20 genes", "Genes correlated with pre-sampling sleep/wake history and post-sampling rebound in the opposite direction.", 4);
  card(s, 0.72, 1.98, 3.35, 4.72, C.navy, C.navy);
  bigStat(s, "15", "History+ / Rebound−", 1.06, 2.38, 2.50, C.gold, "Positive correlation with prior sleep; negative correlation with rebound.");
  bigStat(s, "5", "History− / Rebound+", 1.06, 4.52, 2.50, C.coral, "Negative correlation with prior sleep; positive correlation with rebound.");
  card(s, 4.38, 1.98, 8.24, 2.62, C.white, "E2DED4");
  pill(s, "15  ·  HISTORY+ / REBOUND−", 4.70, 2.28, 2.52, C.gold, C.goldPale);
  s.addText("AstA-R2  ·  bumpel  ·  CG33080  ·  CG41378  ·  CG7460\nCG7601  ·  Chrac-14  ·  Cyp28c1  ·  dpr21  ·  kumpel\nObp44a  ·  Sk1  ·  Ssk  ·  Trhn  ·  VGlut2", {
    x: 4.70, y: 2.90, w: 7.56, h: 1.18, margin: 0,
    fontFace: FONT, fontSize: 15.5, bold: true, color: C.ink, breakLine: false, fit: "shrink",
  });
  card(s, 4.38, 4.92, 8.24, 1.78, C.white, "E2DED4");
  pill(s, "5  ·  HISTORY− / REBOUND+", 4.70, 5.22, 2.48, C.coral, C.coralPale);
  s.addText("CG42323  ·  CG7079  ·  CG9313  ·  CG9377  ·  Gtpbp1", {
    x: 4.70, y: 5.92, w: 7.48, h: 0.30, margin: 0,
    fontFace: FONT, fontSize: 15.5, bold: true, color: C.ink, fit: "shrink",
  });
}

// 5 — CSW 4+
{
  const s = pptx.addSlide();
  header(s, "FAMILY 1  ·  SET 03", "CSW 4+ genes · 97 genes", "Consistent sleep/wake genes evident across four or more of the seven datasets.", 5);
  card(s, 0.72, 1.98, 4.00, 4.72, C.navy, C.navy);
  bigStat(s, "97", "FC0.5 wake genes", 1.06, 2.38, 2.45, C.tealLight, "The category rule is recurrence across datasets—not GO enrichment.");
  s.addShape(pptx.ShapeType.rect, { x: 1.06, y: 4.24, w: 3.20, h: 0.04, fill: { color: "35505E" }, line: { color: "35505E" } });
  bigStat(s, "7", "also qualify at FC1", 1.06, 4.54, 2.45, C.coral, "All seven are already contained within the 97-gene FC0.5 set.");
  card(s, 5.04, 1.98, 7.58, 4.72, C.white, "E2DED4");
  s.addText("Recurrence across seven datasets", {
    x: 5.38, y: 2.28, w: 4.5, h: 0.32, margin: 0,
    fontFace: FONT, fontSize: 20, bold: true, color: C.ink,
  });
  const freq = [
    ["Frequency 4", 82, C.teal],
    ["Frequency 5", 14, C.gold],
    ["Frequency 6", 1, C.coral],
  ];
  freq.forEach((d, i) => termBar(s, d[0], d[1], 82, 5.38, 3.02 + i * 0.72, 6.62, d[2]));
  card(s, 5.38, 5.44, 6.62, 0.90, C.mint, C.mint);
  s.addText("No sleep-upregulated gene met the ≥4-of-7 recurrence criterion.", {
    x: 5.68, y: 5.72, w: 6.02, h: 0.26, margin: 0,
    fontFace: FONT, fontSize: 15, bold: true, align: "center", color: C.teal,
  });
}

// 6 — HLH
{
  const s = pptx.addSlide();
  header(s, "FAMILY 1  ·  SET 04", "HLH genes · 4 genes", "Four bHLH factors identified as potential upstream regulators of wake-induced gene expression.", 6);
  card(s, 0.72, 1.98, 4.18, 4.72, C.navy, C.navy);
  pill(s, "REGULATORY LOGIC", 1.04, 2.30, 1.78, C.white, C.purple);
  s.addText("CAGCTG", {
    x: 1.04, y: 3.00, w: 2.70, h: 0.64, margin: 0,
    fontFace: FONT, fontSize: 34, bold: true, color: C.white,
  });
  s.addText("E-box enriched upstream of wake-induced genes", {
    x: 1.04, y: 3.72, w: 3.15, h: 0.58, margin: 0,
    fontFace: FONT, fontSize: 15, color: "D5E0E3", fit: "shrink",
  });
  s.addShape(pptx.ShapeType.chevron, { x: 1.04, y: 4.62, w: 0.46, h: 0.40, fill: { color: C.tealLight }, line: { color: C.tealLight } });
  s.addText("HOMER motif enrichment → candidate direct regulators", {
    x: 1.72, y: 4.59, w: 2.75, h: 0.48, margin: 0,
    fontFace: FONT, fontSize: 14.5, bold: true, color: C.white, fit: "shrink",
  });
  s.addText("Results · Table S1", {
    x: 1.04, y: 5.74, w: 2.2, h: 0.24, margin: 0,
    fontFace: FONT, fontSize: 12.5, bold: true, color: C.tealLight,
  });
  const genes = [
    ["bigmax", "FBgn0039509"],
    ["HLH3B", "FBgn0011276"],
    ["E(spl)m3-HLH", "FBgn0002609"],
    ["E(spl)mbeta-HLH", "FBgn0002733"],
  ];
  genes.forEach((d, i) => {
    const col = i % 2, row = Math.floor(i / 2);
    const x = 5.26 + col * 3.68, y = 2.24 + row * 1.70;
    card(s, x, y, 3.30, 1.36, C.white, "E2DED4");
    s.addText(d[0], { x: x + 0.24, y: y + 0.28, w: 2.76, h: 0.34, margin: 0, fontFace: FONT, fontSize: 19, bold: true, color: C.teal, fit: "shrink" });
    s.addText(d[1], { x: x + 0.24, y: y + 0.82, w: 2.76, h: 0.25, margin: 0, fontFace: FONT, fontSize: 14, color: C.muted });
  });
  card(s, 5.26, 5.78, 6.98, 0.70, C.coralPale, C.coralPale);
  s.addText("This set is regulatory evidence—not frequency or correlation evidence.", {
    x: 5.58, y: 5.99, w: 6.34, h: 0.24, margin: 0,
    fontFace: FONT, fontSize: 14, bold: true, align: "center", color: C.coral,
  });
}

// 7 — Family 2 categorization workflow
{
  const s = pptx.addSlide();
  header(
    s,
    "FAMILY 2  ·  CATEGORIZATION WORKFLOW",
    "How CSW genes became three GO term–grouped sets",
    "The workflow selects genes from enriched DAVID terms; it does not assign every annotated gene to a category.",
    7,
  );
  const steps = [
    {
      title: "Prepare CSW lists",
      body: "Four master lists; every gene is consistent in ≥2 of 7 datasets.\n\nWake FC0.5 · 647\nWake FC1 · 98\nSleep FC0.5 · 141\nSleep FC1 · 16",
      color: C.coral,
    },
    {
      title: "Run DAVID",
      body: "Run each list separately.\n\nDrosophila · 7227\nDefault species background\nGO BP / CC / MF Direct + KEGG\nEASE cutoff · 0.1\nMin count · 1",
      color: C.teal,
    },
    {
      title: "Filter enriched terms",
      body: "Retain terms with:\n\nFDR ≤ 10%\n\nSleep FC1 produced no passing terms.",
      color: C.gold,
    },
    {
      title: "Assign term buckets",
      body: "Literal term-name patterns:\n\nribosom* / translat*\nmitochondri* / metabolism\nimmune / defense / Toll / Imd\n\nA term may match >1 bucket.",
      color: C.purple,
    },
    {
      title: "Extract gene hits",
      body: "Take only FBgn hits listed in matched terms.\n\nDedupe by FBgn.\nPreserve GO-report source.\nKeep every matched bucket.",
      color: C.teal,
    },
  ];
  steps.forEach((step, i) => {
    const x = 0.72 + i * 2.43;
    card(s, x, 2.05, 2.05, 3.30, C.white, "E2DED4");
    s.addShape(pptx.ShapeType.rect, {
      x, y: 2.05, w: 2.05, h: 0.10,
      fill: { color: step.color }, line: { color: step.color },
    });
    s.addShape(pptx.ShapeType.ellipse, {
      x: x + 0.22, y: 2.37, w: 0.48, h: 0.48,
      fill: { color: step.color }, line: { color: step.color },
    });
    s.addText(String(i + 1), {
      x: x + 0.22, y: 2.49, w: 0.48, h: 0.14, margin: 0,
      fontFace: FONT, fontSize: 10, bold: true, align: "center", color: C.white,
    });
    s.addText(step.title, {
      x: x + 0.82, y: 2.39, w: 1.00, h: 0.42, margin: 0,
      fontFace: FONT, fontSize: 14, bold: true, color: C.ink, fit: "shrink",
    });
    s.addText(step.body, {
      x: x + 0.22, y: 3.06, w: 1.62, h: 1.94, margin: 0,
      fontFace: FONT, fontSize: 12.5, color: C.muted, valign: "top", fit: "shrink",
    });
    if (i < steps.length - 1) {
      s.addShape(pptx.ShapeType.chevron, {
        x: x + 2.12, y: 3.42, w: 0.22, h: 0.34,
        fill: { color: C.line }, line: { color: C.line },
      });
    }
  });
  card(s, 0.72, 5.68, 11.90, 1.02, C.navy, C.navy);
  s.addText("OUTPUT", {
    x: 1.02, y: 5.96, w: 0.78, h: 0.22, margin: 0,
    fontFace: FONT, fontSize: 10.5, bold: true, charSpacing: 1.3, color: C.tealLight,
  });
  s.addText("99 Ribosomal  ·  167 Mitochondrial  ·  5 Immune", {
    x: 1.92, y: 5.88, w: 4.45, h: 0.34, margin: 0,
    fontFace: FONT, fontSize: 18, bold: true, color: C.white,
  });
  s.addText("All outputs are wake-direction. Mitochondrial translation terms match both keyword groups, so 31 genes are retained in both Ribosomal and Mitochondrial.", {
    x: 6.58, y: 5.84, w: 5.58, h: 0.48, margin: 0,
    fontFace: FONT, fontSize: 13.5, color: "D5E0E3", fit: "shrink",
  });
}

// 8 — Ribosomal GO
{
  const s = pptx.addSlide();
  header(s, "FAMILY 2  ·  SET 05", "Ribosomal / translation · 99 genes", "GO-report assignment is shown here; experiment-level FC0.5 and FC1 regulation is shown on slide 13.", 8);
  addThresholdPanel(s, 99, 99, 12, 12, C.coral, "99 FC0.5 + 12 FC1 − 12 both = 99 unique");
  card(s, 5.24, 1.98, 7.38, 4.78, C.white, "E2DED4");
  s.addText("Most frequent approved terms", {
    x: 5.58, y: 2.30, w: 5.60, h: 0.32, margin: 0,
    fontFace: FONT, fontSize: 20, bold: true, color: C.ink,
  });
  const terms = [
    ["Ribosome", 93],
    ["Structural constituent of ribosome", 92],
    ["Translation", 77],
    ["Cytoplasmic translation", 65],
    ["Cytosolic ribosome", 58],
  ];
  terms.forEach((d, i) => termBar(s, d[0], d[1], 93, 5.58, 3.02 + i * 0.58, 6.46, C.coral));
  card(s, 5.58, 6.08, 6.46, 0.44, C.coralPale, C.coralPale);
  s.addText("87 FC0.5-GO-only  ·  0 FC1-GO-only  ·  12 in both GO reports", {
    x: 5.84, y: 6.20, w: 5.94, h: 0.19, margin: 0,
    fontFace: FONT, fontSize: 13.5, bold: true, align: "center", color: C.coral,
  });
}

// 9 — Mitochondrial GO
{
  const s = pptx.addSlide();
  header(s, "FAMILY 2  ·  SET 06", "Mitochondrial / metabolism · 167 genes", "GO-report assignment is shown here; experiment-level FC0.5 and FC1 regulation is shown on slide 14.", 9);
  addThresholdPanel(s, 167, 167, 23, 23, C.teal, "167 FC0.5 + 23 FC1 − 23 both = 167 unique");
  card(s, 5.24, 1.98, 7.38, 4.78, C.white, "E2DED4");
  s.addText("Most frequent approved terms", {
    x: 5.58, y: 2.30, w: 5.60, h: 0.32, margin: 0,
    fontFace: FONT, fontSize: 20, bold: true, color: C.ink,
  });
  const terms = [
    ["Mitochondrion", 121],
    ["Mitochondrial inner membrane", 72],
    ["Oxidative phosphorylation", 46],
    ["Mitochondrial translation", 30],
    ["Mitochondrial matrix", 24],
  ];
  terms.forEach((d, i) => termBar(s, d[0], d[1], 121, 5.58, 3.02 + i * 0.58, 6.46, C.teal));
  card(s, 5.58, 6.08, 6.46, 0.44, C.mint, C.mint);
  s.addText("144 FC0.5-GO-only  ·  0 FC1-GO-only  ·  23 in both GO reports", {
    x: 5.84, y: 6.20, w: 5.94, h: 0.19, margin: 0,
    fontFace: FONT, fontSize: 13.5, bold: true, align: "center", color: C.teal,
  });
}

// 10 — Immune GO
{
  const s = pptx.addSlide();
  header(s, "FAMILY 2  ·  SET 07", "Immune · 5 genes", "Assigned from FC1 GO enrichment; all five genes also meet FC0.5 regulation criteria (experiment view on slide 15).", 10);
  addThresholdPanel(s, 5, 0, 5, 0, C.gold, "0 FC0.5 + 5 FC1 − 0 both = 5 unique");
  card(s, 5.24, 1.98, 7.38, 4.78, C.white, "E2DED4");
  s.addText("Five genes assigned by the FC1 GO report", {
    x: 5.58, y: 2.30, w: 2.35, h: 0.62, margin: 0,
    fontFace: FONT, fontSize: 17.5, bold: true, color: C.ink, fit: "shrink",
  });
  const genes = ["BomBc2", "BomS1", "BomS2", "IM14", "Thor"];
  genes.forEach((g, i) => {
    const y = 3.00 + i * 0.58;
    s.addShape(pptx.ShapeType.ellipse, { x: 5.58, y: y + 0.02, w: 0.28, h: 0.28, fill: { color: C.gold }, line: { color: C.gold } });
    s.addText(g, { x: 6.04, y, w: 1.60, h: 0.28, margin: 0, fontFace: FONT, fontSize: 15, bold: true, color: C.ink });
  });
  s.addText("Approved GO terms · genes per term", {
    x: 8.45, y: 2.30, w: 3.25, h: 0.32, margin: 0,
    fontFace: FONT, fontSize: 17.5, bold: true, color: C.ink, fit: "shrink",
  });
  const terms = [["Defense response", 4], ["Antibacterial humoral response", 4], ["Response to bacterium", 4]];
  terms.forEach((d, i) => {
    const y = 3.04 + i * 0.76;
    s.addShape(pptx.ShapeType.rect, { x: 8.45, y: y + 0.04, w: 0.12, h: 0.28, fill: { color: C.gold }, line: { color: C.gold } });
    s.addText(d[0], { x: 8.76, y, w: 2.32, h: 0.34, margin: 0, fontFace: FONT, fontSize: 14, bold: true, color: C.ink, fit: "shrink" });
    s.addText(String(d[1]), { x: 11.32, y, w: 0.42, h: 0.34, margin: 0, fontFace: FONT, fontSize: 17, bold: true, align: "right", color: C.gold });
  });
  s.addText("Genes may occur in more than one approved term.", {
    x: 8.45, y: 5.18, w: 3.57, h: 0.26, margin: 0,
    fontFace: FONT, fontSize: 12.5, color: C.muted, fit: "shrink",
  });
  card(s, 8.45, 5.66, 3.57, 0.70, C.goldPale, C.goldPale);
  s.addText("0 FC0.5-GO-only  ·  5 FC1-GO-only  ·  0 in both reports", {
    x: 8.66, y: 5.88, w: 3.15, h: 0.24, margin: 0,
    fontFace: FONT, fontSize: 13.5, bold: true, align: "center", color: "9A6E13", fit: "shrink",
  });
}

// 11 — Family 1 overlap
{
  const s = pptx.addSlide();
  header(s, "FAMILY 1  ·  OVERLAP", "Only three genes cross primary-set boundaries", "Overlap is calculated among Mechanistic, Homeostatic, CSW 4+, and HLH only.", 11);
  card(s, 0.72, 1.98, 7.54, 4.72, C.white, "E2DED4");
  s.addText("Observed shared genes", {
    x: 1.06, y: 2.30, w: 3.2, h: 0.32, margin: 0,
    fontFace: FONT, fontSize: 20, bold: true, color: C.ink,
  });
  const rows = [
    ["Mechanistic × Homeostatic", "2", "Trhn  ·  AstA-R2", C.coral],
    ["Mechanistic × CSW 4+", "1", "RpL23", C.teal],
    ["All other pairings", "0", "No shared FBgns", C.muted],
  ];
  rows.forEach((d, i) => {
    const y = 3.04 + i * 0.90;
    s.addShape(pptx.ShapeType.rect, { x: 1.06, y, w: 0.08, h: 0.58, fill: { color: d[3] }, line: { color: d[3] } });
    s.addText(d[0], { x: 1.36, y: y + 0.06, w: 3.05, h: 0.28, margin: 0, fontFace: FONT, fontSize: 14.5, bold: true, color: C.ink, fit: "shrink" });
    s.addText(d[1], { x: 4.62, y: y - 0.02, w: 0.60, h: 0.44, margin: 0, fontFace: FONT, fontSize: 27, bold: true, color: d[3] });
    s.addText(d[2], { x: 5.40, y: y + 0.06, w: 2.20, h: 0.28, margin: 0, fontFace: FONT, fontSize: 14.5, bold: i < 2, color: i < 2 ? C.ink : C.muted, fit: "shrink" });
  });
  card(s, 8.58, 1.98, 4.04, 4.72, C.navy, C.navy);
  s.addText("Family 1 union", { x: 8.94, y: 2.30, w: 2.8, h: 0.34, margin: 0, fontFace: FONT, fontSize: 21, bold: true, color: C.white });
  const stats = [["124", "unique FBgns"], ["121", "in exactly one set"], ["3", "in exactly two sets"], ["0", "in three or four sets"]];
  stats.forEach((d, i) => {
    const y = 3.00 + i * 0.70;
    s.addText(d[0], { x: 8.94, y, w: 0.95, h: 0.40, margin: 0, fontFace: FONT, fontSize: 25, bold: true, color: i === 0 ? C.tealLight : C.coral });
    s.addText(d[1], { x: 10.02, y: y + 0.06, w: 1.85, h: 0.27, margin: 0, fontFace: FONT, fontSize: 13.5, color: "D5E0E3", fit: "shrink" });
  });
  card(s, 8.94, 6.02, 3.12, 0.44, "1B3443", "35505E");
  s.addText("HLH genes are fully disjoint.", {
    x: 9.16, y: 6.14, w: 2.68, h: 0.20, margin: 0,
    fontFace: FONT, fontSize: 13.5, bold: true, align: "center", color: C.white,
  });
}

// 12 — Family 2 overlap
{
  const s = pptx.addSlide();
  header(s, "FAMILY 2  ·  OVERLAP", "The GO-derived union contains 240 unique genes", "99 ribosomal + 167 mitochondrial + 5 immune − 31 shared ribosomal/mitochondrial genes = 240 unique FBgns.", 12);
  card(s, 0.72, 1.98, 7.50, 4.72, C.white, "E2DED4");
  s.addImage({
    path: path.join(FIG, "pathway_overlap_venn.png"),
    x: 0.98, y: 2.20, w: 6.98, h: 3.82,
    sizing: { type: "contain", w: 6.98, h: 3.82 },
    altText: "Venn diagram showing 68 ribosomal-only, 31 shared ribosomal and mitochondrial, 136 mitochondrial-only, and 5 immune-only genes",
  });
  card(s, 1.02, 6.08, 6.90, 0.42, C.mint, C.mint);
  s.addText("240 = union of Ribosomal, Mitochondrial and Immune genes; each FBgn counted once.", {
    x: 1.24, y: 6.19, w: 6.46, h: 0.19, margin: 0,
    fontFace: FONT, fontSize: 13, bold: true, align: "center", color: C.teal, fit: "shrink",
  });
  card(s, 8.54, 1.98, 4.08, 4.72, C.navy, C.navy);
  s.addText("GO-assignment provenance\nacross the 240-gene union", {
    x: 8.90, y: 2.30, w: 3.08, h: 0.68, margin: 0,
    fontFace: FONT, fontSize: 20, bold: true, color: C.white, fit: "shrink",
  });
  const threshold = [
    ["203", "FC0.5 GO only", C.tealLight],
    ["32", "both GO reports", C.gold],
    ["5", "FC1 GO only", C.coral],
  ];
  threshold.forEach((d, i) => {
    const y = 3.34 + i * 0.90;
    s.addText(d[0], { x: 8.90, y, w: 0.96, h: 0.48, margin: 0, fontFace: FONT, fontSize: 29, bold: true, color: d[2] });
    s.addText(d[1], { x: 9.98, y: y + 0.09, w: 1.84, h: 0.28, margin: 0, fontFace: FONT, fontSize: 14, bold: true, color: C.white });
  });
  s.addShape(pptx.ShapeType.rect, { x: 8.90, y: 6.10, w: 3.02, h: 0.04, fill: { color: "35505E" }, line: { color: "35505E" } });
  s.addText("Only the five Immune genes are assigned by FC1 GO alone.", {
    x: 8.90, y: 6.28, w: 3.02, h: 0.28, margin: 0,
    fontFace: FONT, fontSize: 14, bold: true, color: C.white, fit: "shrink",
  });
}

// 13–15 — Experiment-level regulation within GO categories
addExperimentDistributionSlide(
  13,
  "05",
  "Ribosomal / translation",
  99,
  [40, 68, 97, 48, 0, 25, 37],
  12,
  [10, 5, 12, 5, 0, 0, 0],
  "MechSD6 regulates 97/99 FC0.5-responsive genes and all 12 FC1-responsive genes; FC1 signal is mechanical-deprivation focused.",
);

addExperimentDistributionSlide(
  14,
  "06",
  "Mitochondrial / metabolism",
  167,
  [88, 126, 163, 86, 1, 9, 2],
  24,
  [20, 10, 24, 6, 0, 0, 0],
  "MechSD6 regulates 163/167 FC0.5-responsive genes and all 24 FC1-responsive genes; FC1 has no R85C10, ChR or THIP hits.",
);

addExperimentDistributionSlide(
  15,
  "07",
  "Immune",
  5,
  [1, 4, 3, 4, 0, 0, 0],
  5,
  [1, 4, 3, 4, 0, 0, 0],
  "All five immune genes meet both thresholds in at least one experiment; MechSD3 and MechSD12 each regulate four of five.",
);

pptx.writeFile({ fileName: OUT });

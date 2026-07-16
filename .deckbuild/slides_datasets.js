const { P, BC, HEAD, BODY, W, H, bg, pageFooter, header } = require("./lib.js");
const { card, stat, dataTable, shortBucket } = require("./comp.js");

function bucketBreakdown(s, x, y, w, buckets) {
  const head = [
    { text: " ", options: { fill: { color: P.INK }, color: P.WHITE } },
    { text: "Priority bucket", options: { fill: { color: P.INK }, color: P.WHITE, bold: true, fontSize: 10.5, fontFace: BODY } },
    { text: "Stocks", options: { fill: { color: P.INK }, color: P.WHITE, bold: true, fontSize: 10.5, fontFace: BODY, align: "center" } },
    { text: "Alleles", options: { fill: { color: P.INK }, color: P.WHITE, bold: true, fontSize: 10.5, fontFace: BODY, align: "center" } },
    { text: "Genes", options: { fill: { color: P.INK }, color: P.WHITE, bold: true, fontSize: 10.5, fontFace: BODY, align: "center" } },
  ];
  const body = buckets.map((b, i) => {
    const f = i % 2 ? "EEF2F5" : P.WHITE;
    return [
      { text: " ", options: { fill: { color: BC[i] } } },
      { text: shortBucket(b[0]).replace("  ·  Phenotype", ""), options: { fill: { color: f }, color: P.TX, fontSize: 9.5, fontFace: BODY, valign: "middle" } },
      { text: b[1], options: { fill: { color: f }, color: P.TX, fontSize: 10, fontFace: BODY, align: "center", bold: true, valign: "middle" } },
      { text: b[2], options: { fill: { color: f }, color: P.MUT, fontSize: 10, fontFace: BODY, align: "center", valign: "middle" } },
      { text: b[3], options: { fill: { color: f }, color: P.MUT, fontSize: 10, fontFace: BODY, align: "center", valign: "middle" } },
    ];
  });
  s.addTable([head, ...body], { x, y, w, colW: [0.16, 3.04, 0.9, 0.9, 0.85], border: { type: "solid", pt: 0.5, color: P.LINE }, rowH: 0.34, margin: 2 });
}

function datasetSlide(pres, cfg, n, TOTAL) {
  const s = pres.addSlide(); bg(s);
  header(s, cfg.kicker, cfg.title);
  // stat callouts
  const t = cfg.total;
  const stats = [["Input genes", cfg.n_genes, P.AMBER], ["Unique stocks", t[1], P.TEAL], ["Unique alleles", t[2], P.TEAL], ["Genes w/ stocks", t[3], P.TEAL]];
  let x = 0.55;
  stats.forEach(st => {
    card(s, x, 1.5, 2.92, 1.05);
    stat(s, x, 1.62, 2.92, st[1], st[0], st[2]);
    x += 3.05;
  });
  // left: schema
  card(s, 0.55, 2.78, 5.55, 3.7);
  s.addText("Input schema", { x: 0.78, y: 2.92, w: 5.1, h: 0.34, fontFace: BODY, fontSize: 13.5, bold: true, color: P.INK, margin: 0 });
  s.addText(cfg.schemaNote, { x: 0.78, y: 3.24, w: 5.1, h: 0.5, fontFace: BODY, fontSize: 9.5, italic: true, color: P.MUT, margin: 0 });
  let y = 3.78;
  cfg.cols.forEach(c => {
    s.addShape("oval", { x: 0.8, y: y + 0.04, w: 0.13, h: 0.13, fill: { color: P.TEAL } });
    s.addText([{ text: c[0] + "  ", options: { bold: true, color: P.TX } }, { text: "— " + c[1], options: { color: P.MUT } }],
      { x: 1.02, y: y - 0.04, w: 5.0, h: 0.32, fontFace: BODY, fontSize: 9.8, margin: 0 });
    y += 0.42;
  });
  // right: bucket breakdown
  card(s, 6.3, 2.78, 6.5, 3.7);
  s.addText("Stocks per priority bucket", { x: 6.53, y: 2.92, w: 6.0, h: 0.34, fontFace: BODY, fontSize: 13.5, bold: true, color: P.INK, margin: 0 });
  bucketBreakdown(s, 6.53, 3.34, 5.85, cfg.buckets);
  if (cfg.foot) s.addText(cfg.foot, { x: 6.53, y: 6.12, w: 6.0, h: 0.32, fontFace: BODY, fontSize: 9, italic: true, color: P.MUT, margin: 0 });
  pageFooter(s, n, TOTAL);
}

// transcriptomics category slide: shared schema + tier list
function txSlide(pres, cfg, n, TOTAL) {
  const s = pres.addSlide(); bg(s);
  header(s, cfg.kicker, cfg.title);
  // left: schema
  card(s, 0.55, 1.5, 5.55, 4.98);
  s.addText("Input schema", { x: 0.78, y: 1.64, w: 5.1, h: 0.34, fontFace: BODY, fontSize: 13.5, bold: true, color: P.INK, margin: 0 });
  s.addText(cfg.schemaNote, { x: 0.78, y: 1.97, w: 5.15, h: 0.62, fontFace: BODY, fontSize: 9.5, italic: true, color: P.MUT, margin: 0 });
  let y = 2.66;
  cfg.cols.forEach(c => {
    s.addShape("oval", { x: 0.8, y: y + 0.04, w: 0.13, h: 0.13, fill: { color: P.TEAL } });
    s.addText([{ text: c[0] + "  ", options: { bold: true, color: P.TX } }, { text: "— " + c[1], options: { color: P.MUT } }],
      { x: 1.02, y: y - 0.04, w: 5.0, h: 0.34, fontFace: BODY, fontSize: 9.8, margin: 0 });
    y += 0.46;
  });
  // right: tiers table
  card(s, 6.3, 1.5, 6.5, 4.98);
  s.addText(cfg.rightTitle, { x: 6.53, y: 1.64, w: 6.0, h: 0.34, fontFace: BODY, fontSize: 13.5, bold: true, color: P.INK, margin: 0 });
  dataTable(s, 6.53, 2.04, 5.85, ["Gene-set table (tier)", "Genes"], cfg.tiers, [4.85, 1.0], { rowH: 0.3, bf: 9.8, hf: 10.5 });
  if (cfg.foot) s.addText(cfg.foot, { x: 6.53, y: 6.05, w: 6.0, h: 0.36, fontFace: BODY, fontSize: 9, italic: true, color: P.MUT, margin: 0 });
  pageFooter(s, n, TOTAL);
}

module.exports = { datasetSlide, txSlide };

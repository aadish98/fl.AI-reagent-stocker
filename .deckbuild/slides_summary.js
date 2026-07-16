const { P, BC, HEAD, BODY, W, H, bg, pageFooter, header } = require("./lib.js");
const { card, stat, dataTable, shortBucket, hbars } = require("./comp.js");

function summaryTargeted(pres, data, n, TOTAL) {
  const s = pres.addSlide(); bg(s);
  header(s, "Summary counts", "Targeted gene-list runs — totals");
  const order = ["vGAT — GABA gene list", "BPD-GWAS orthologs", "DIOPT Bellenguez/Lambert (AD)", "Jordan VCM", "Slumber Neuropeptides+transmitters", "Slumber RNAi candidates"];
  const rows = order.map(name => {
    const d = data.targeted[name]; const t = d.total;
    return [name.replace("Slumber Neuropeptides+transmitters", "Slumber neuropeptides"), String(d.n_genes), t[1], t[2], t[3]];
  });
  // total
  let st = 0, sa = 0;
  order.forEach(name => { st += parseInt(data.targeted[name].total[1]); sa += parseInt(data.targeted[name].total[2]); });
  rows.push([{ text: "All targeted runs", options: { bold: true } }, "619", { text: String(st), options: { bold: true } }, { text: String(sa), options: { bold: true } }, "—"]);
  // left table
  card(s, 0.55, 1.5, 6.55, 4.55);
  s.addText("Unique stocks / alleles per dataset", { x: 0.78, y: 1.64, w: 6.1, h: 0.34, fontFace: BODY, fontSize: 13, bold: true, color: P.INK, margin: 0 });
  dataTable(s, 0.78, 2.04, 6.1, ["Dataset", "Genes", "Stocks", "Alleles", "w/ stocks"], rows, [2.65, 0.78, 0.95, 0.95, 0.92], { rowH: 0.42, bf: 10, hf: 10.5 });
  // right chart
  card(s, 7.3, 1.5, 5.5, 4.55);
  s.addText("Stocks organized per dataset", { x: 7.5, y: 1.64, w: 5.1, h: 0.34, fontFace: BODY, fontSize: 13, bold: true, color: P.INK, margin: 0 });
  const labels = ["GABA", "BPD-GWAS", "DIOPT AD", "Jordan", "Slumber NP", "Slumber RNAi"];
  const barRows = order.map((name, i) => [labels[i], parseInt(data.targeted[name].total[1])]);
  hbars(s, 7.55, 2.2, 5.05, barRows, P.TEAL, { rh: 0.6, labelW: 1.55, valW: 0.85, lf: 9.5 });
  s.addText("Counts are de-duplicated within each run (a stock is counted once across the six buckets).", { x: 0.55, y: 6.2, w: 12.2, h: 0.34, fontFace: BODY, fontSize: 9.5, italic: true, color: P.MUT, margin: 0 });
  pageFooter(s, n, TOTAL);
}

function summaryTx(pres, data, n, TOTAL) {
  const s = pres.addSlide(); bg(s);
  header(s, "Summary counts", "Transcriptomics run — totals");
  // stat ribbon
  const ribbon = [["16", "gene-set tiers", P.AMBER], ["5", "dataset categories", P.TEAL], [String(data.tx_total_stocks), "stocks organized", P.TEAL]];
  let x = 0.55;
  ribbon.forEach(r => { card(s, x, 1.5, 2.6, 1.0); stat(s, x, 1.62, 2.6, r[0], r[1], r[2]); x += 2.75; });
  // left: per dataset
  card(s, 0.55, 2.75, 6.05, 3.7);
  s.addText("By dataset category", { x: 0.78, y: 2.89, w: 5.6, h: 0.32, fontFace: BODY, fontSize: 12.5, bold: true, color: P.INK, margin: 0 });
  const nm = { "CSW FC0.5": "CSW FC0.5", "CSW FC1": "CSW FC1", "Homeostatic Genes": "Homeostatic genes", "Overlapping SleepHistory<>SleepRebound Correlates": "History x Rebound overlaps", "SleepHistory Correlates": "SleepHistory correlates" };
  const rows = data.tx_per_dataset.map(r => [nm[r[0]] || r[0], String(r[1]), String(r[2])]);
  dataTable(s, 0.78, 3.26, 5.6, ["Dataset category", "Tiers", "Stocks"], rows, [3.5, 0.95, 1.15], { rowH: 0.42, bf: 10, hf: 10.5 });
  // right: by bucket
  card(s, 6.8, 2.75, 6.0, 3.7);
  s.addText("By priority bucket (all 16 tiers)", { x: 7.0, y: 2.89, w: 5.6, h: 0.32, fontFace: BODY, fontSize: 12.5, bold: true, color: P.INK, margin: 0 });
  const bucketShort = ["Bloom · UAS", "Vienna · UAS", "Bloom · Allele", "Stock Ctr · Allele", "Non-SC · UAS", "Non-SC · Allele"];
  const bRows = data.tx_by_bucket.map((b, i) => [bucketShort[i] || shortBucket(b[0]), b[1]]);
  hbars(s, 7.0, 3.35, 5.55, bRows, P.AMBER, { rh: 0.48, labelW: 1.75, valW: 0.8, lf: 9 });
  s.addText("Transcriptomics totals are summed across tiers; genes/alleles recur across tiers, so stocks are the headline metric.", { x: 0.55, y: 6.55, w: 12.2, h: 0.34, fontFace: BODY, fontSize: 9.5, italic: true, color: P.MUT, margin: 0 });
  pageFooter(s, n, TOTAL);
}

function examples(pres, data, n, TOTAL) {
  const s = pres.addSlide(); bg(s);
  header(s, "Worked examples", "One example reagent from each category");
  s.addText("Representative rows from the GABA gene-list workbook (n=140) — every dataset follows the same six-table structure.", { x: 0.55, y: 1.34, w: 12.2, h: 0.32, fontFace: BODY, fontSize: 10.5, color: P.MUT, margin: 0 });
  const ex = data.examples;
  const cw = 6.05, ch = 1.70, gx = 0.18, gy = 0.10, startY = 1.74;
  ex.forEach((e, i) => {
    const col = i % 2, rowi = Math.floor(i / 2);
    const X = 0.55 + col * (cw + gx), Y = startY + rowi * (ch + gy);
    s.addShape("rect", { x: X, y: Y, w: 0.14, h: ch, fill: { color: BC[i] } });
    s.addShape("roundRect", { x: X + 0.14, y: Y, w: cw - 0.14, h: ch, rectRadius: 0.06, fill: { color: P.WHITE }, line: { color: P.LINE, width: 1 } });
    s.addText(e.label, { x: X + 0.32, y: Y + 0.1, w: cw - 0.5, h: 0.3, fontFace: BODY, fontSize: 11, bold: true, color: BC[i], margin: 0 });
    s.addText(e.rows + " stocks in this bucket", { x: X + 0.32, y: Y + 0.4, w: cw - 0.5, h: 0.24, fontFace: BODY, fontSize: 8.5, italic: true, color: P.MUT, margin: 0 });
    const line = (lab, val, yy) => {
      s.addText([{ text: lab + "  ", options: { bold: true, color: P.INK } }, { text: val, options: { color: P.TX } }],
        { x: X + 0.32, y: yy, w: cw - 0.55, h: 0.24, fontFace: BODY, fontSize: 9, margin: 0 });
    };
    line("ID:", e.id, Y + 0.64);
    line("Collection:", e.collection, Y + 0.88);
    line("Gene:", e.gene, Y + 1.12);
    line("Genotype:", e.genotype, Y + 1.36);
  });
  pageFooter(s, n, TOTAL);
}

function closing(pres, data, n, TOTAL) {
  const s = pres.addSlide();
  s.background = { color: P.INK };
  s.addShape("rect", { x: 0, y: 0, w: 0.28, h: H, fill: { color: P.TEAL } });
  s.addShape("rect", { x: 0.28, y: 0, w: 0.08, h: H, fill: { color: P.AMBER } });
  s.addText("TAKEAWAYS", { x: 1.0, y: 1.0, w: 11, h: 0.4, fontFace: BODY, fontSize: 13, color: P.TEALL, charSpacing: 4, bold: true, margin: 0 });
  s.addText("What the tables give you", { x: 0.95, y: 1.4, w: 11.6, h: 0.8, fontFace: HEAD, fontSize: 34, bold: true, color: P.WHITE, margin: 0 });
  const pts = [
    ["22 gene sets, one consistent format", "Every input list becomes the same six prioritized, phenotype-filtered stock tables."],
    ["Ordering by orderability + reagent type", "Bloomington / Vienna / other centers first; UAS vs classical alleles separated; custom reagents last."],
    ["Built-in literature relevance", "Sleep / circadian / rhythm / locomotor references and phenotype-similarity scores travel with each stock."],
    ["8,300+ phenotype-backed stocks organized", "Ready to filter down to an orderable shortlist per gene."],
  ];
  let y = 2.7;
  pts.forEach((p, i) => {
    s.addShape("oval", { x: 1.0, y: y, w: 0.42, h: 0.42, fill: { color: P.TEAL } });
    s.addText(String(i + 1), { x: 1.0, y: y, w: 0.42, h: 0.42, fontFace: HEAD, fontSize: 15, bold: true, color: P.WHITE, align: "center", valign: "middle", margin: 0 });
    s.addText(p[0], { x: 1.6, y: y - 0.04, w: 11, h: 0.36, fontFace: BODY, fontSize: 15, bold: true, color: P.WHITE, margin: 0 });
    s.addText(p[1], { x: 1.6, y: y + 0.32, w: 11, h: 0.4, fontFace: BODY, fontSize: 11.5, color: "C7D2DE", margin: 0 });
    y += 0.95;
  });
  pageFooter(s, n, TOTAL);
}

module.exports = { summaryTargeted, summaryTx, examples, closing };

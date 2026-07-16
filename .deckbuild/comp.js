const { P, BC, HEAD, BODY, shadow } = require("./lib.js");

// short label for a bucket combination string
function shortBucket(combo) {
  return combo.replace(/ >> /g, "  ·  ")
    .replace("AlleleOrInsertion", "Allele/Insertion")
    .replace("Non Stock Center", "Non-Stock-Center")
    .replace("Vienna", "Vienna (VDRC)");
}

// rounded card with soft shadow
function card(slide, x, y, w, h, fill) {
  slide.addShape("roundRect", { x, y, w, h, rectRadius: 0.09, fill: { color: fill || P.CARD }, line: { color: P.LINE, width: 1 }, shadow: shadow() });
}

// big stat callout
function stat(slide, x, y, w, value, label, color) {
  slide.addText(String(value), { x, y, w, h: 0.62, fontFace: HEAD, fontSize: 30, bold: true, color: color || P.TEAL, align: "center", margin: 0 });
  slide.addText(label, { x, y: y + 0.62, w, h: 0.34, fontFace: BODY, fontSize: 10.5, color: P.MUT, align: "center", margin: 0 });
}

// a colored chip with text (for bucket legend)
function chip(slide, x, y, w, h, color, title, sub) {
  slide.addShape("rect", { x, y, w: 0.12, h, fill: { color } });
  slide.addShape("rect", { x: x + 0.12, y, w: w - 0.12, h, fill: { color: P.WHITE }, line: { color: P.LINE, width: 1 } });
  slide.addText(title, { x: x + 0.26, y: y + 0.04, w: w - 0.4, h: h - 0.06, fontFace: BODY, fontSize: 10.5, bold: true, color: P.TX, valign: "middle", margin: 0 });
  if (sub) slide.addText(sub, { x: x + 0.26, y: y + h - 0.3, w: w - 0.4, h: 0.26, fontFace: BODY, fontSize: 8.5, color: P.MUT, margin: 0 });
}

// styled table from header + rows
function dataTable(slide, x, y, w, header, rows, colW, opts) {
  opts = opts || {};
  const head = header.map(t => ({ text: t, options: { fill: { color: P.INK }, color: P.WHITE, bold: true, fontSize: opts.hf || 11, align: "left", fontFace: BODY, valign: "middle" } }));
  const body = rows.map((r, i) => r.map(c => {
    const cell = (typeof c === "object") ? c : { text: String(c), options: {} };
    cell.options = Object.assign({ fill: { color: i % 2 ? "EEF2F5" : P.WHITE }, color: P.TX, fontSize: opts.bf || 10.5, fontFace: BODY, valign: "middle" }, cell.options);
    return cell;
  }));
  slide.addTable([head, ...body], { x, y, w, colW, border: { type: "solid", pt: 0.5, color: P.LINE }, rowH: opts.rowH || 0.3, autoPage: false, margin: opts.margin == null ? 3 : opts.margin });
}

// hand-built horizontal bars (avoids pptxgenjs chart XML entirely)
function hbars(slide, x, y, w, rows, color, opts) {
  opts = opts || {};
  const rh = opts.rh || 0.5;
  const labelW = opts.labelW || 1.4;
  const valW = opts.valW || 0.85;
  const barMax = w - labelW - valW - 0.12;
  const maxV = Math.max.apply(null, rows.map(r => r[1]));
  rows.forEach((r, i) => {
    const ry = y + i * rh;
    slide.addText(r[0], { x, y: ry, w: labelW - 0.1, h: rh, fontFace: BODY, fontSize: opts.lf || 9.5, color: P.MUT, align: "right", valign: "middle", margin: 0 });
    const bw = Math.max(0.05, barMax * (r[1] / maxV));
    const by = ry + rh * 0.2;
    slide.addShape("rect", { x: x + labelW, y: by, w: bw, h: rh * 0.6, fill: { color } });
    slide.addText(r[1].toLocaleString(), { x: x + labelW + bw + 0.06, y: ry, w: valW, h: rh, fontFace: BODY, fontSize: opts.lf || 9.5, bold: true, color: P.TX, valign: "middle", margin: 0 });
  });
}

module.exports = { shortBucket, card, stat, chip, dataTable, hbars };

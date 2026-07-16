// Shared style + helpers for the Allada stocker deck
const P = {
  INK: "15233B", INK2: "203250", TEAL: "0E7C7B", TEALD: "0A5C5B",
  TEALL: "8FCFCD", AMBER: "E0A13A", BG: "F4F6F8", WHITE: "FFFFFF",
  TX: "1F2A37", MUT: "667385", LINE: "E2E8F0", CARD: "FFFFFF"
};
// 6 priority-bucket colors (consistent across deck)
const BC = ["0E7C7B", "1C6FB0", "6B4FA0", "C2792B", "B2415A", "4F7A2E"];
const HEAD = "Georgia";
const BODY = "Calibri";
const W = 13.333, H = 7.5;

function shadow() { return { type: "outer", color: "9AA6B2", blur: 7, offset: 3, angle: 135, opacity: 0.28 }; }

function bg(slide) { slide.background = { color: P.BG }; }

function pageFooter(slide, n, total) {
  slide.addShape("rect", { x: 0, y: H - 0.32, w: W, h: 0.32, fill: { color: P.INK } });
  slide.addText("Allada Lab  ·  FlyBase Reagent Stocker", { x: 0.5, y: H - 0.34, w: 8, h: 0.32, fontFace: BODY, fontSize: 9, color: P.TEALL, valign: "middle", margin: 0 });
  slide.addText(n + " / " + total, { x: W - 1.5, y: H - 0.34, w: 1.0, h: 0.32, fontFace: BODY, fontSize: 9, color: P.TEALL, align: "right", valign: "middle", margin: 0 });
}

function header(slide, kicker, title) {
  slide.addShape("rect", { x: 0, y: 0, w: W, h: 1.18, fill: { color: P.INK } });
  slide.addShape("rect", { x: 0, y: 1.18, w: W, h: 0.06, fill: { color: P.AMBER } });
  if (kicker) slide.addText(kicker.toUpperCase(), { x: 0.55, y: 0.18, w: 12, h: 0.3, fontFace: BODY, fontSize: 11, color: P.TEALL, charSpacing: 3, bold: true, margin: 0 });
  slide.addText(title, { x: 0.55, y: 0.44, w: 12.2, h: 0.62, fontFace: HEAD, fontSize: 25, color: P.WHITE, bold: true, margin: 0 });
}

module.exports = { P, BC, HEAD, BODY, W, H, shadow, bg, pageFooter, header };

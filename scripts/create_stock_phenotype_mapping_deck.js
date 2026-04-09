const path = require("path");
const PptxGenJS = require("pptxgenjs");

const pptx = new PptxGenJS();
pptx.layout = "LAYOUT_WIDE";
pptx.author = "OpenAI";
pptx.company = "OpenAI";
pptx.subject = "Stock Phenotype Sheet ID mapping flow";
pptx.title = "From input gene symbols to Stock Phenotype Sheet rows";
pptx.lang = "en-US";
pptx.theme = {
  headFontFace: "Georgia",
  bodyFontFace: "Calibri",
  lang: "en-US",
};

const OUT = path.resolve(
  process.cwd(),
  "stock_phenotype_sheet_mapping_flow.pptx",
);

const W = 13.333;
const H = 7.5;

const C = {
  navy: "16324F",
  teal: "0E7490",
  aqua: "14B8A6",
  gold: "F0B429",
  cream: "F7FBFC",
  ice: "E6F4F7",
  ink: "14212B",
  slate: "51606E",
  line: "C8D7DE",
  white: "FFFFFF",
  paleTeal: "DDF3F1",
  paleGold: "FEF3D0",
  paleBlue: "E8EFFA",
  rose: "FCEAEA",
};

function shadow() {
  return {
    type: "outer",
    color: "000000",
    blur: 3,
    offset: 1.5,
    angle: 45,
    opacity: 0.12,
  };
}

function addSlideFrame(slide, title, kicker, opts = {}) {
  slide.background = { color: opts.bg || C.cream };
  slide.addShape(pptx.ShapeType.rect, {
    x: 0,
    y: 0,
    w: W,
    h: 0.82,
    fill: { color: opts.header || C.navy },
    line: { color: opts.header || C.navy, transparency: 100 },
  });
  slide.addText(kicker, {
    x: 0.64,
    y: 0.18,
    w: 5.3,
    h: 0.18,
    fontFace: "Calibri",
    fontSize: 10,
    color: "CDE7EF",
    bold: true,
    charSpacing: 0.6,
    margin: 0,
  });
  slide.addText(title, {
    x: 0.62,
    y: 0.36,
    w: 9.4,
    h: 0.28,
    fontFace: "Georgia",
    fontSize: 24,
    color: C.white,
    bold: true,
    margin: 0,
  });
}

function addFooter(slide, text, pageNo, opts = {}) {
  slide.addText(text, {
    x: 0.68,
    y: 7.08,
    w: 11.5,
    h: 0.18,
    fontFace: "Calibri",
    fontSize: opts.fontSize || 9.5,
    color: opts.color || C.slate,
    italic: true,
    margin: 0,
  });
  slide.addText(String(pageNo), {
    x: 12.3,
    y: 7.04,
    w: 0.35,
    h: 0.22,
    fontFace: "Calibri",
    fontSize: 10,
    color: opts.color || C.slate,
    bold: true,
    align: "right",
    margin: 0,
  });
}

function addPill(slide, x, y, w, label, fill, color = C.ink, opts = {}) {
  const pillH = opts.h || 0.42;
  slide.addShape(pptx.ShapeType.roundRect, {
    x,
    y,
    w,
    h: pillH,
    rectRadius: 0.07,
    fill: { color: fill },
    line: { color: fill },
  });
  slide.addText(label, {
    x: x + 0.12,
    y: y + (pillH - 0.18) / 2,
    w: w - 0.24,
    h: 0.18,
    fontFace: "Calibri",
    fontSize: opts.fontSize || 11,
    bold: true,
    color,
    align: "center",
    margin: 0,
    fit: "shrink",
  });
}

function addCard(slide, cfg) {
  const accentColor = cfg.accent || C.teal;
  slide.addShape(pptx.ShapeType.rect, {
    x: cfg.x,
    y: cfg.y,
    w: cfg.w,
    h: cfg.h,
    fill: { color: cfg.fill || C.white },
    line: { color: cfg.line || C.line, pt: 1 },
    shadow: shadow(),
  });
  slide.addShape(pptx.ShapeType.rect, {
    x: cfg.x,
    y: cfg.y,
    w: 0.1,
    h: cfg.h,
    fill: { color: accentColor },
    line: { color: accentColor, transparency: 100 },
  });
  slide.addText(cfg.title, {
    x: cfg.x + 0.2,
    y: cfg.y + 0.14,
    w: cfg.w - 0.32,
    h: 0.28,
    fontFace: "Georgia",
    fontSize: cfg.titleSize || 16,
    color: C.ink,
    bold: true,
    margin: 0,
  });
  slide.addText(cfg.body, {
    x: cfg.x + 0.2,
    y: cfg.y + 0.52,
    w: cfg.w - 0.32,
    h: cfg.h - 0.64,
    fontFace: cfg.bodyFont || "Calibri",
    fontSize: cfg.bodySize || 12,
    color: cfg.bodyColor || C.slate,
    breakLine: false,
    valign: "top",
    margin: 0,
    fit: "shrink",
  });
}

function addConnector(slide, x, y, w, h = 0.08, color = C.aqua) {
  slide.addShape(pptx.ShapeType.rect, {
    x,
    y,
    w,
    h,
    fill: { color },
    line: { color, transparency: 100 },
  });
}

function addColumnHeader(slide, x, y, w, label) {
  slide.addShape(pptx.ShapeType.rect, {
    x,
    y,
    w,
    h: 0.46,
    fill: { color: C.navy },
    line: { color: C.navy, pt: 1 },
  });
  slide.addText(label, {
    x: x + 0.08,
    y: y + 0.12,
    w: w - 0.16,
    h: 0.16,
    fontFace: "Calibri",
    fontSize: 11,
    bold: true,
    color: C.white,
    margin: 0,
    align: "center",
  });
}

function addCell(slide, x, y, w, h, text, fill, opts = {}) {
  slide.addShape(pptx.ShapeType.rect, {
    x,
    y,
    w,
    h,
    fill: { color: fill },
    line: { color: C.line, pt: 0.75 },
  });
  slide.addText(text, {
    x: x + 0.08,
    y: y + 0.08,
    w: w - 0.16,
    h: h - 0.16,
    fontFace: opts.fontFace || "Calibri",
    fontSize: opts.fontSize || 10.5,
    color: opts.color || C.ink,
    bold: opts.bold || false,
    margin: 0,
    fit: "shrink",
    valign: "mid",
  });
}

function addFlowBox(slide, cfg) {
  const pillH = cfg.pillHeight || 0.34;
  const pillY = cfg.source ? cfg.y + cfg.h - pillH - 0.1 : null;
  const bodyBottomInset = cfg.sourceInline
    ? 0.28
    : cfg.source
      ? pillH + 0.16
      : 0.16;
  slide.addShape(pptx.ShapeType.rect, {
    x: cfg.x,
    y: cfg.y,
    w: cfg.w,
    h: cfg.h,
    fill: { color: cfg.fill || C.white },
    line: { color: cfg.line || C.line, pt: 1 },
    shadow: shadow(),
  });
  slide.addShape(pptx.ShapeType.ellipse, {
    x: cfg.x + 0.15,
    y: cfg.y + 0.14,
    w: 0.42,
    h: 0.42,
    fill: { color: cfg.badgeFill || C.teal },
    line: { color: cfg.badgeFill || C.teal, transparency: 100 },
  });
  slide.addText(String(cfg.step), {
    x: cfg.x + 0.15,
    y: cfg.y + 0.24,
    w: 0.42,
    h: 0.12,
    align: "center",
    fontFace: "Calibri",
    fontSize: 11,
    bold: true,
    color: C.white,
    margin: 0,
  });
  slide.addText(cfg.title, {
    x: cfg.x + 0.68,
    y: cfg.y + 0.14,
    w: cfg.w - 0.82,
    h: 0.22,
    fontFace: "Georgia",
    fontSize: 15,
    color: C.ink,
    bold: true,
    margin: 0,
    fit: "shrink",
  });
  slide.addText(cfg.body, {
    x: cfg.x + (cfg.bodyIndentUnderBadge ? 0.68 : 0.18),
    y: cfg.y + 0.38,
    w: cfg.w - (cfg.bodyIndentUnderBadge ? 0.82 : 0.34),
    h: cfg.h - 0.38 - bodyBottomInset,
    fontFace: cfg.bodyFont || "Calibri",
    fontSize: cfg.bodySize || 11.5,
    color: C.slate,
    margin: 0,
    fit: "shrink",
  });
  if (cfg.sourceInline && cfg.source) {
    slide.addText(cfg.source, {
      x: cfg.x + 0.18,
      y: cfg.y + cfg.h - 0.2,
      w: cfg.w - 0.36,
      h: 0.14,
      fontFace: "Calibri",
      fontSize: cfg.sourceFontSize || 9.2,
      bold: true,
      color: C.teal,
      margin: 0,
      fit: "shrink",
    });
  } else if (cfg.source) {
    addPill(
      slide,
      cfg.x + 0.18,
      pillY,
      cfg.w - 0.36,
      cfg.source,
      cfg.pillFill || C.paleBlue,
      cfg.pillColor || C.navy,
      { h: pillH, fontSize: cfg.sourceFontSize || 9.5 },
    );
  }
}

function buildTitleSlide() {
  const slide = pptx.addSlide();
  slide.background = { color: C.navy };
  slide.addShape(pptx.ShapeType.rect, {
    x: 0,
    y: 0,
    w: W,
    h: H,
    fill: { color: C.navy },
    line: { color: C.navy, transparency: 100 },
  });
  slide.addShape(pptx.ShapeType.rect, {
    x: 8.55,
    y: 0,
    w: 4.78,
    h: 7.5,
    fill: { color: C.teal, transparency: 18 },
    line: { color: C.teal, transparency: 100 },
  });
  slide.addShape(pptx.ShapeType.rect, {
    x: 8.95,
    y: 0.48,
    w: 3.95,
    h: 6.55,
    fill: { color: C.aqua, transparency: 72 },
    line: { color: C.aqua, transparency: 100 },
  });

  slide.addText("fl_ai_reagent_stocker", {
    x: 0.7,
    y: 0.62,
    w: 3.1,
    h: 0.22,
    fontFace: "Calibri",
    fontSize: 12,
    bold: true,
    color: "CDE7EF",
    charSpacing: 1.2,
    margin: 0,
  });
  slide.addText("From input gene symbols to Stock Phenotype Sheet rows", {
    x: 0.7,
    y: 1.18,
    w: 6.25,
    h: 1.25,
    fontFace: "Georgia",
    fontSize: 28,
    bold: true,
    color: C.white,
    margin: 0,
    fit: "shrink",
  });
  slide.addText(
    "Exact FlyBase ID mapping, source tables, and implementation comments for the soft-run phenotype path.",
    {
      x: 0.72,
      y: 2.56,
      w: 5.85,
      h: 0.8,
      fontFace: "Calibri",
      fontSize: 17,
      color: "DCEBF0",
      margin: 0,
      fit: "shrink",
    },
  );

  addPill(slide, 0.72, 3.72, 1.7, "7 mapping edges", C.paleGold, C.ink);
  addPill(slide, 2.58, 3.72, 1.7, "6 source tables", C.paleTeal, C.ink);
  addPill(slide, 4.44, 3.72, 2.0, "4 caveats to state", C.paleBlue, C.ink);

  slide.addShape(pptx.ShapeType.rect, {
    x: 8.12,
    y: 1.0,
    w: 4.38,
    h: 5.18,
    fill: { color: "F8FCFD" },
    line: { color: "9EC3D0", pt: 1.2 },
    shadow: shadow(),
  });
  slide.addText("Presenter shorthand", {
    x: 8.38,
    y: 1.26,
    w: 2.6,
    h: 0.24,
    fontFace: "Calibri",
    fontSize: 12,
    bold: true,
    color: C.teal,
    margin: 0,
  });
  slide.addText(
    "symbol -> FBgn -> FBal\n-> FBtp/FBti -> FBst\n-> relevant_component_ids\n-> genotype_FBids overlap\n-> phenotype-sheet row",
    {
      x: 8.38,
      y: 1.72,
      w: 3.66,
      h: 2.2,
      fontFace: "Consolas",
      fontSize: 16,
      color: C.ink,
      bold: true,
      margin: 0,
      fit: "shrink",
      valign: "mid",
    },
  );
  slide.addText(
    "Downstream inclusion is driven by the unioned Stage 1 component IDs, not by the original user symbol column.",
    {
      x: 8.38,
      y: 4.42,
      w: 3.55,
      h: 0.9,
      fontFace: "Calibri",
      fontSize: 14,
      color: C.ink,
      margin: 0,
      fit: "shrink",
    },
  );
  addFooter(
    slide,
    "Sources: scripts/fetch_fbgn_ids.py, stock_finding.py, stock_splitting.py, FlyBase Genes/alleles_and_stocks/transgenic tables",
    1,
    { color: "DCEBF0", fontSize: 9.5 },
  );
}

function buildFlowSlide() {
  const slide = pptx.addSlide();
  addSlideFrame(
    slide,
    "End-to-End Flowchart",
    "Operational path used by the soft-run Stock Phenotype Sheet",
  );

  const leftX = 0.72;
  const rightX = 7.02;
  const boxW = 5.58;
  const boxH = 1.04;
  const ys = [1.18, 2.34, 3.5, 4.66];

  const leftSteps = [
    {
      step: 1,
      title: "Input symbols",
      body: "Read user CSV values from ext_gene or the chosen input-gene column.",
      source: "Source: CSV input",
      badgeFill: C.navy,
      bodySize: 11.6,
    },
    {
      step: 2,
      title: "Convert to FBgn",
      body: "Map current symbols first, then expanded synonyms, keeping only Dmel primary FBgn IDs.",
      source: "Source: fb_synonym",
      badgeFill: C.teal,
      bodySize: 11.6,
    },
    {
      step: 3,
      title: "Seed source alleles",
      body: "Use direct GeneID -> AlleleID joins and construct-linked allele seeding through regulatory-region or encoded-product matches.",
      source: "Source: fbal_to_fbgn + construct tables",
      badgeFill: C.aqua,
      bodySize: 11.3,
    },
    {
      step: 4,
      title: "Expand FBtp / FBti",
      body: "Collect construct and insertion IDs so a later stock match can happen through more than one component namespace.",
      source: "Source: construct + insertion expansion",
      badgeFill: C.gold,
      bodySize: 11.4,
    },
  ];

  const rightSteps = [
    {
      step: 5,
      title: "Match to FBst",
      body: "Any FBal, FBtp, or FBti component can hit a stock via derived_stock_component.",
      source: "Source: derived stock components",
      badgeFill: C.navy,
      bodySize: 11.4,
    },
    {
      step: 6,
      title: "Build join key",
      body: "Store relevant_component_ids = FBal union FBtp union FBti for the matched source-allele family.",
      source: "Source: Stage 1 stock row",
      badgeFill: C.teal,
      bodySize: 11.2,
    },
    {
      step: 7,
      title: "Filter phenotype rows",
      body: "Keep genotype_phenotype_data rows whose genotype_FBids overlap the relevant component set.",
      source: "Source: genotype_phenotype_data",
      badgeFill: C.aqua,
      bodySize: 11.4,
    },
    {
      step: 8,
      title: "Write phenotype rows",
      body: "Resolve stock label, gene, reagent symbol, phenotype, qualifier, and reference metadata into the final sheet.",
      source: "Source: Stock Phenotype Sheet",
      badgeFill: C.gold,
      bodySize: 11.3,
    },
  ];

  leftSteps.forEach((cfg, idx) => {
    addFlowBox(slide, {
      ...cfg,
      x: leftX,
      y: ys[idx],
      w: boxW,
      h: boxH,
      fill: C.white,
      bodyIndentUnderBadge: true,
      sourceInline: true,
      sourceFontSize: 9.2,
    });
  });
  rightSteps.forEach((cfg, idx) => {
    addFlowBox(slide, {
      ...cfg,
      x: rightX,
      y: ys[idx],
      w: boxW,
      h: boxH,
      fill: C.white,
      bodyIndentUnderBadge: true,
      sourceInline: true,
      sourceFontSize: 9.2,
    });
  });

  addConnector(slide, 3.46, 2.24, 0.08, 0.1);
  addConnector(slide, 3.46, 3.4, 0.08, 0.1);
  addConnector(slide, 3.46, 4.56, 0.08, 0.1);
  addConnector(slide, 9.76, 2.24, 0.08, 0.1);
  addConnector(slide, 9.76, 3.4, 0.08, 0.1);
  addConnector(slide, 9.76, 4.56, 0.08, 0.1);

  slide.addShape(pptx.ShapeType.rect, {
    x: 6.46,
    y: 2.62,
    w: 0.84,
    h: 1.5,
    fill: { color: C.ice },
    line: { color: C.line, pt: 1 },
    shadow: shadow(),
  });
  slide.addText("Then", {
    x: 6.6,
    y: 2.84,
    w: 0.56,
    h: 0.18,
    fontFace: "Georgia",
    fontSize: 13,
    bold: true,
    color: C.ink,
    margin: 0,
    align: "center",
  });
  slide.addText("switch to stock + phenotype recovery", {
    x: 6.56,
    y: 3.18,
    w: 0.66,
    h: 0.56,
    fontFace: "Calibri",
    fontSize: 9.2,
    color: C.slate,
    margin: 0,
    fit: "shrink",
    align: "center",
    valign: "mid",
  });

  addCard(slide, {
    x: 0.72,
    y: 5.96,
    w: 11.92,
    h: 0.74,
    title: "What actually drives inclusion",
    body: "Inclusion depends on the chain FBgn -> source allele family -> relevant_component_ids -> genotype_FBids overlap.",
    fill: C.ice,
    accent: C.navy,
    titleSize: 14,
    bodySize: 12.6,
  });

  addFooter(
    slide,
    "Implementation anchors: fetch_fbgn_ids.py, stock_finding.py gene-component helpers, stock_splitting.py phenotype-sheet builder",
    2,
    { fontSize: 9.4 },
  );
}

function buildExactMapSlide() {
  const slide = pptx.addSlide();
  addSlideFrame(
    slide,
    "Exact ID Map and Data Sources",
    "Each edge the phenotype-sheet path can traverse",
    { bg: C.cream, header: C.teal },
  );

  const x0 = 0.62;
  const y0 = 1.18;
  const rowH = 0.67;
  const colW = [1.45, 2.18, 3.32, 5.1];
  const headers = ["From", "To / key field", "Data source", "Comment"];

  let x = x0;
  headers.forEach((header, i) => {
    addColumnHeader(slide, x, y0, colW[i], header);
    x += colW[i];
  });

  const rows = [
    [
      "gene symbol",
      "FBgn via flybase_gene_id",
      "Genes/fb_synonym*.tsv(.gz)",
      "Maps current symbols first, then expanded fullname/symbol synonyms; restricted to Dmel and primary FBgn IDs.",
    ],
    [
      "FBgn",
      "FBal via GeneID -> AlleleID",
      "fbal_to_fbgn*.tsv(.gz)",
      "Main direct gene-to-allele edge used to seed source alleles.",
    ],
    [
      "FBgn or FBsf helper",
      "extra FBal seeds",
      "transgenic_construct_descriptions*.tsv + fbsf_to_fbgn.csv",
      "Construct rows can seed component alleles when the input gene matches a regulatory region or encoded product/tool.",
    ],
    [
      "FBal",
      "FBtp",
      "transgenic_construct_descriptions*.tsv(.gz)",
      "Expands a source allele into construct reagent IDs plus symbols, class terms, and supporting references.",
    ],
    [
      "FBal or FBtp",
      "FBti",
      "insertion_allele_descriptions*.tsv + fbtp_to_fbti.csv",
      "Insertion IDs are recovered directly from allele descriptions or indirectly from construct-to-insertion expansion.",
    ],
    [
      "FBal / FBtp / FBti",
      "FBst via\nderived_stock_component",
      "fbst_to_derived_stock_component.csv",
      "A stock can match through any of the three component namespaces, not only through alleles.",
    ],
    [
      "matched stock row",
      "relevant_component_ids",
      "Stage 1 workbook columns",
      "Stage 1 stores the union of gene-relevant FBal, FBtp, and FBti IDs; this union is the downstream phenotype-sheet join key.",
    ],
  ];

  rows.forEach((row, idx) => {
    const fill = idx % 2 === 0 ? C.white : C.ice;
    let cx = x0;
    row.forEach((cell, i) => {
      addCell(
        slide,
        cx,
        y0 + 0.46 + idx * rowH,
        colW[i],
        rowH,
        cell,
        fill,
        {
          fontFace: i < 2 ? "Consolas" : "Calibri",
          fontSize: i < 2 ? 10 : 10.2,
          color: i === 2 ? C.navy : C.ink,
          bold: i === 0,
        },
      );
      cx += colW[i];
    });
  });

  slide.addShape(pptx.ShapeType.rect, {
    x: 0.62,
    y: 6.28,
    w: 12.01,
    h: 0.54,
    fill: { color: C.paleGold },
    line: { color: C.gold, pt: 1 },
    shadow: shadow(),
  });
  slide.addText("Key presenter line: the phenotype sheet is keyed by component overlap, not by a direct gene-to-phenotype lookup.", {
    x: 0.82,
    y: 6.45,
    w: 11.58,
    h: 0.18,
    fontFace: "Georgia",
    fontSize: 13.5,
    bold: true,
    color: C.ink,
    margin: 0,
    fit: "shrink",
  });

  addFooter(
    slide,
    "Files referenced here come from data/flybase/Genes, data/flybase/alleles_and_stocks, data/flybase/transgenic_constructs, and data/flybase/transgenic_insertions",
    3,
  );
}

function buildJoinKeySlide() {
  const slide = pptx.addSlide();
  addSlideFrame(
    slide,
    "How Stage 1 Builds the Join Key",
    "Why the union is broader than the direct stock match",
  );

  addCard(slide, {
    x: 0.72,
    y: 1.22,
    w: 5.25,
    h: 4.95,
    title: "Stored formula",
    body:
      "source_allele_ids\n\n" +
      "-> relevant_fbal_ids\n" +
      "-> relevant_fbtp_ids gathered from constructs_by_allele\n" +
      "-> relevant_fbti_ids gathered from insertions_by_allele\n\n" +
      "relevant_component_ids = FBal union FBtp union FBti\n\n" +
      "This value is written into each Stage 1 stock row and later reused by the phenotype-sheet builder.",
    fill: C.white,
    accent: C.teal,
    titleSize: 18,
    bodySize: 15,
    bodyFont: "Consolas",
    bodyColor: C.ink,
  });

  addPill(slide, 1.0, 5.56, 1.18, "FBal IDs", C.paleBlue, C.navy);
  addPill(slide, 2.32, 5.56, 1.18, "FBtp IDs", C.paleTeal, C.navy);
  addPill(slide, 3.64, 5.56, 1.18, "FBti IDs", C.paleGold, C.navy);

  addCard(slide, {
    x: 6.34,
    y: 1.22,
    w: 6.28,
    h: 1.5,
    title: "Match provenance can come from more than one route",
    body: "direct_allele, construct_regulatory_region, and construct_encoded_product all feed the same source-allele family before the union is built.",
    fill: C.white,
    accent: C.navy,
    titleSize: 16,
    bodySize: 13.5,
  });
  addCard(slide, {
    x: 6.34,
    y: 2.98,
    w: 6.28,
    h: 1.56,
    title: "Why the reverse index uses both Chado and Stage 1 IDs",
    body: "A phenotype row may mention an FBal ID even when the visible Chado stock components are only FBtp or FBti. The code therefore indexes both derived_stock_component and relevant_component_ids.",
    fill: C.white,
    accent: C.aqua,
    titleSize: 16,
    bodySize: 13.5,
  });
  addCard(slide, {
    x: 6.34,
    y: 4.82,
    w: 6.28,
    h: 1.34,
    title: "Custom phenotype reagent path",
    body: "If an input-linked allele has phenotype evidence but never lands on an FBst stock, Stage 1 creates a synthetic Custom phenotype reagent row keyed by the allele symbol.",
    fill: C.white,
    accent: C.gold,
    titleSize: 16,
    bodySize: 13.5,
  });

  addFooter(
    slide,
    "Code anchors: stock_finding.py:_build_gene_component_tables, _build_stock_mapping, _build_custom_phenotype_rows",
    4,
    { fontSize: 9.4, color: "465766" },
  );
}

function buildPhenotypeJoinSlide() {
  const slide = pptx.addSlide();
  addSlideFrame(
    slide,
    "How a phenotype row gets into the sheet",
    "Stage 2 soft-run filtering, reverse resolution, and row construction",
    { bg: C.cream, header: C.teal },
  );

  addFlowBox(slide, {
    step: 1,
    x: 0.72,
    y: 1.28,
    w: 3.22,
    h: 1.42,
    title: "Extract row IDs",
    body: "Parse FlyBase IDs out of genotype_FBids with a regex-based extractor.",
    source: "genotype_FBids field",
    fill: C.white,
    badgeFill: C.navy,
    pillFill: C.paleBlue,
    sourceFontSize: 9.2,
    pillHeight: 0.3,
  });
  addFlowBox(slide, {
    step: 2,
    x: 0.72,
    y: 3.0,
    w: 3.22,
    h: 1.42,
    title: "Apply inclusion filter",
    body: "Keep the phenotype row only if its extracted IDs intersect all_relevant_ids from stock-backed and custom rows.",
    source: "intersection test",
    fill: C.white,
    badgeFill: C.teal,
    pillFill: C.paleTeal,
    sourceFontSize: 9.2,
    pillHeight: 0.3,
  });
  addFlowBox(slide, {
    step: 3,
    x: 0.72,
    y: 4.72,
    w: 3.22,
    h: 1.42,
    title: "Resolve stock keys",
    body: "Map matched component IDs back to FBst rows or custom stock_number keys, then expand FBrf metadata into reference fields.",
    source: "reverse index + references",
    fill: C.white,
    badgeFill: C.gold,
    pillFill: C.paleGold,
    sourceFontSize: 9.2,
    pillHeight: 0.3,
  });
  addConnector(slide, 2.27, 2.76, 0.08, 0.18);
  addConnector(slide, 2.27, 4.48, 0.08, 0.18);

  addCard(slide, {
    x: 4.38,
    y: 1.28,
    w: 4.1,
    h: 3.24,
    title: "Output fields emitted per resolved phenotype row",
    body:
      "Gene\n" +
      "Reagent Type or Allele Symbol\n" +
      "Source/ Stock #\n" +
      "Genotype\n" +
      "Phenotype\n" +
      "Qualifier\n" +
      "PMID / PMCID\n" +
      "Reference\n" +
      "Authors / Journal / Year of Publication",
    fill: C.white,
    accent: C.navy,
    titleSize: 17,
    bodySize: 14,
    bodyFont: "Calibri",
    bodyColor: C.ink,
  });

  addCard(slide, {
    x: 8.78,
    y: 1.28,
    w: 3.86,
    h: 1.55,
    title: "Important nuance",
    body: "The sheet builder is currently passed the full Stage 1 stocks_df, even though the contents-sheet prose says phenotype rows are limited to stocks present in output sheets.",
    fill: C.white,
    accent: C.gold,
    titleSize: 15,
    bodySize: 12.8,
  });
  addCard(slide, {
    x: 8.78,
    y: 3.12,
    w: 3.86,
    h: 1.72,
    title: "Actual uniqueness key",
    body: "Rows are grouped by Source/ Stock #, Genotype, Phenotype, Qualifier, and reference metadata fields. That is more specific than source/stock + genotype + reference alone.",
    fill: C.white,
    accent: C.aqua,
    titleSize: 15,
    bodySize: 12.8,
  });
  addCard(slide, {
    x: 4.38,
    y: 5.18,
    w: 8.26,
    h: 0.88,
    title: "Presenter takeaway",
    body: "A phenotype row survives because one of its genotype_FBids overlaps the precomputed relevant component set and can be resolved back to a stock key.",
    fill: C.ice,
    accent: C.teal,
    titleSize: 14,
    bodySize: 12,
  });

  addFooter(
    slide,
    "Code anchors: stock_splitting.py:_extract_flybase_ids, _build_stock_phenotype_sheet, grouping + row-order logic",
    5,
    { fontSize: 9.4, color: "465766" },
  );
}

function buildCaveatSlide() {
  const slide = pptx.addSlide();
  slide.background = { color: C.navy };
  slide.addShape(pptx.ShapeType.rect, {
    x: 0,
    y: 0,
    w: W,
    h: 0.86,
    fill: { color: C.teal },
    line: { color: C.teal, transparency: 100 },
  });
  slide.addText("Comments worth saying out loud", {
    x: 0.68,
    y: 0.26,
    w: 6.4,
    h: 0.3,
    fontFace: "Georgia",
    fontSize: 24,
    bold: true,
    color: C.white,
    margin: 0,
  });

  const cards = [
    {
      x: 0.82,
      y: 1.34,
      w: 5.9,
      h: 1.54,
      accent: C.gold,
      title: "Workbook prose vs code path",
      body: "The contents-sheet note says the phenotype sheet covers stocks present in output sheets, but the code calls the builder with the full Stage 1 stocks_df.",
      fill: "F9FCFD",
    },
    {
      x: 6.92,
      y: 1.34,
      w: 5.6,
      h: 1.54,
      accent: C.aqua,
      title: "Grouping is more specific than the prose says",
      body: "Implementation groups by phenotype and qualifier in addition to stock, genotype, and reference metadata, so distinct phenotypes remain on separate rows.",
      fill: "F9FCFD",
    },
    {
      x: 0.82,
      y: 3.22,
      w: 5.9,
      h: 1.54,
      accent: C.teal,
      title: "Regex caution",
      body: "_extract_flybase_ids() uses the loose pattern FB[a-zA-Z]{2,4}[0-9]+. The audit notes this can overmatch non-FlyBase strings such as FBXO11.",
      fill: "F9FCFD",
    },
    {
      x: 6.92,
      y: 3.22,
      w: 5.6,
      h: 1.54,
      accent: C.gold,
      title: "Gene labels are reconstructed",
      body: "The phenotype-sheet Gene column is rebuilt from matched component metadata, so it is not a direct copy of the user's original input symbol column.",
      fill: "F9FCFD",
    },
  ];

  cards.forEach((card) => addCard(slide, card));

  slide.addShape(pptx.ShapeType.rect, {
    x: 0.82,
    y: 5.22,
    w: 11.7,
    h: 1.08,
    fill: { color: C.aqua, transparency: 12 },
    line: { color: C.aqua, pt: 1.2 },
  });
  slide.addText("Bottom line", {
    x: 1.08,
    y: 5.62,
    w: 1.8,
    h: 0.2,
    fontFace: "Calibri",
    fontSize: 12,
    color: C.white,
    bold: true,
    margin: 0,
  });
  slide.addText(
    "If you describe this workflow as 'component-driven phenotype row recovery using precomputed gene-relevant FBal/FBtp/FBti unions,' you will be faithful to the implementation.",
    {
      x: 2.3,
      y: 5.44,
      w: 9.3,
      h: 0.42,
      fontFace: "Georgia",
      fontSize: 16,
      bold: true,
      color: C.white,
      margin: 0,
      fit: "shrink",
    },
  );

  addFooter(
    slide,
    "Caveats sourced from stock_splitting.py plus scripts/audit_stock_phenotype_pipeline.py",
    6,
    { color: "DCEBF0", fontSize: 9.5 },
  );
}

buildTitleSlide();
buildFlowSlide();
buildExactMapSlide();
buildJoinKeySlide();
buildPhenotypeJoinSlide();
buildCaveatSlide();

pptx.writeFile({ fileName: OUT });

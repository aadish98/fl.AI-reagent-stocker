#!/usr/bin/env python3
"""Build the shareable Stocker SOP PDF (docs/SOP_Stocker_Cursor_CLI.pdf)."""

from __future__ import annotations

from pathlib import Path

from reportlab.lib.colors import HexColor, white
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Flowable,
    KeepTogether,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_PATH = REPO_ROOT / "docs" / "SOP_Stocker_Cursor_CLI.pdf"

INK = HexColor("#1C1915")
MUTED = HexColor("#5C564C")
RULE = HexColor("#8B1E1E")
RULE_BG = HexColor("#F6EDE6")
NAVY = HexColor("#1F3A3A")
NAVY_SOFT = HexColor("#E7EEEE")
CREAM = HexColor("#F7F3EC")
CODE_BG = HexColor("#EFEBE3")
LINE = HexColor("#D4CFC4")
PAPER = HexColor("#FBF8F2")

FLYBASE_URL = "https://flybase.org/convert/id"
FLYBASE_REPORT_EXAMPLE = "https://flybase.org/reports/FBgn0023076"
GITHUB_URL = "https://github.com/Allada-Lab/fl-AI-reagent-stocker"
CURSOR_URL = "https://cursor.com"
PYTHON_URL = "https://www.python.org/downloads/"
DEFAULT_CONFIG = "data/config/stock_split_config_priority_example.json"
SOP_VERSION = "v3"
SOP_DATE = "27 August 2026"

# Body is 11 pt. Other sizes are the previous scale multiplied by 11 / 9.1.
BODY = 11.0
SCALE = BODY / 9.1


def _pt(previous: float) -> float:
    return round(previous * SCALE, 1)


def _register_fonts() -> None:
    avenir = "/System/Library/Fonts/Avenir Next.ttc"
    pdfmetrics.registerFont(TTFont("AvenirNextB", avenir, subfontIndex=0))
    pdfmetrics.registerFont(TTFont("AvenirNextBI", avenir, subfontIndex=1))
    pdfmetrics.registerFont(TTFont("AvenirNextD", avenir, subfontIndex=2))
    pdfmetrics.registerFont(TTFont("AvenirNextDI", avenir, subfontIndex=3))
    pdfmetrics.registerFont(TTFont("AvenirNextI", avenir, subfontIndex=4))
    pdfmetrics.registerFont(TTFont("AvenirNextM", avenir, subfontIndex=5))
    pdfmetrics.registerFont(TTFont("AvenirNext", avenir, subfontIndex=7))
    pdfmetrics.registerFontFamily(
        "AvenirNext",
        normal="AvenirNext",
        bold="AvenirNextD",
        italic="AvenirNextI",
        boldItalic="AvenirNextDI",
    )
    pdfmetrics.registerFont(TTFont("Menlo", "/System/Library/Fonts/Menlo.ttc", subfontIndex=0))
    pdfmetrics.registerFont(TTFont("MenloB", "/System/Library/Fonts/Menlo.ttc", subfontIndex=1))


class Hairline(Flowable):
    def __init__(self, color=LINE, thickness=0.6, space_before=2, space_after=6):
        super().__init__()
        self.color = color
        self.thickness = thickness
        self.spaceBefore = space_before
        self.spaceAfter = space_after
        self.height = thickness

    def wrap(self, availWidth, availHeight):
        self.width = availWidth
        return availWidth, self.height

    def draw(self):
        self.canv.setStrokeColor(self.color)
        self.canv.setLineWidth(self.thickness)
        self.canv.line(0, 0, self.width, 0)


class Callout(Flowable):
    def __init__(self, title, body_flowables, border=RULE, fill=RULE_BG, width=None):
        super().__init__()
        self.title = title
        self.body = body_flowables
        self.border = border
        self.fill = fill
        self._width = width
        self._inner = None

    def wrap(self, availWidth, availHeight):
        self.width = self._width or availWidth
        inner_w = self.width - 18
        story = [
            Paragraph(
                self.title,
                ParagraphStyle(
                    "calloutTitle",
                    fontName="AvenirNextD",
                    fontSize=_pt(9.2),
                    leading=_pt(12),
                    textColor=self.border,
                    spaceAfter=5,
                ),
            )
        ]
        story.extend(self.body)
        y = 0
        wrapped = []
        for item in story:
            w, h = item.wrap(inner_w, availHeight)
            wrapped.append((item, w, h))
            y += h + getattr(item, "spaceAfter", 0) + getattr(item, "spaceBefore", 0)
        self._inner = wrapped
        self.height = y + 18
        return self.width, self.height

    def draw(self):
        canv = self.canv
        canv.setFillColor(self.fill)
        canv.roundRect(0, 0, self.width, self.height, 3, fill=1, stroke=0)
        canv.setFillColor(self.border)
        canv.rect(0, 0, 4.5, self.height, fill=1, stroke=0)
        y = self.height - 9
        for item, _w, h in self._inner:
            y -= getattr(item, "spaceBefore", 0)
            y -= h
            item.drawOn(canv, 11, y)
            y -= getattr(item, "spaceAfter", 0)


class CodeBlock(Flowable):
    def __init__(self, lines):
        super().__init__()
        self.raw_lines = lines
        self._line_h = _pt(11.2)
        self._font_size = _pt(7.6)
        self.lines = lines

    def wrap(self, availWidth, availHeight):
        self.width = availWidth
        max_chars = max(8, int((availWidth - 18) / (self._font_size * 0.62)))
        wrapped = []
        for line in self.raw_lines:
            if not line:
                wrapped.append("")
                continue
            while len(line) > max_chars:
                wrapped.append(line[:max_chars])
                line = line[max_chars:]
            wrapped.append(line)
        self.lines = wrapped
        self.height = 12 + self._line_h * len(self.lines)
        return self.width, self.height

    def draw(self):
        canv = self.canv
        canv.setFillColor(CODE_BG)
        canv.roundRect(0, 0, self.width, self.height, 2, fill=1, stroke=0)
        canv.setFillColor(INK)
        canv.setFont("Menlo", self._font_size)
        y = self.height - 14
        for line in self.lines:
            canv.drawString(8, y, line)
            y -= self._line_h


def _styles():
    base = getSampleStyleSheet()
    styles = {
        "kicker": ParagraphStyle(
            "kicker",
            parent=base["Normal"],
            fontName="AvenirNextD",
            fontSize=_pt(8),
            leading=_pt(10),
            textColor=NAVY,
            spaceAfter=3,
        ),
        "title": ParagraphStyle(
            "title",
            parent=base["Normal"],
            fontName="AvenirNextD",
            fontSize=_pt(16.5),
            leading=_pt(19),
            textColor=INK,
            spaceAfter=4,
        ),
        "subtitle": ParagraphStyle(
            "subtitle",
            parent=base["Normal"],
            fontName="AvenirNextI",
            fontSize=_pt(9.4),
            leading=_pt(12.2),
            textColor=MUTED,
            spaceAfter=10,
        ),
        "h": ParagraphStyle(
            "h",
            parent=base["Normal"],
            fontName="AvenirNextD",
            fontSize=_pt(10),
            leading=_pt(13),
            textColor=NAVY,
            spaceBefore=9,
            spaceAfter=4,
        ),
        "body": ParagraphStyle(
            "body",
            parent=base["Normal"],
            fontName="AvenirNext",
            fontSize=BODY,
            leading=_pt(12.1),
            textColor=INK,
            alignment=TA_JUSTIFY,
            spaceAfter=5,
        ),
        "bodyL": ParagraphStyle(
            "bodyL",
            parent=base["Normal"],
            fontName="AvenirNext",
            fontSize=BODY,
            leading=_pt(12.1),
            textColor=INK,
            alignment=TA_LEFT,
            spaceAfter=4,
        ),
        "bullet": ParagraphStyle(
            "bullet",
            parent=base["Normal"],
            fontName="AvenirNext",
            fontSize=_pt(9.0),
            leading=_pt(12.0),
            textColor=INK,
            leftIndent=14,
            bulletIndent=0,
            spaceAfter=2,
        ),
        "step": ParagraphStyle(
            "step",
            parent=base["Normal"],
            fontName="AvenirNext",
            fontSize=_pt(9.0),
            leading=_pt(12.0),
            textColor=INK,
            leftIndent=18,
            firstLineIndent=-14,
            spaceAfter=3,
        ),
        "meta": ParagraphStyle(
            "meta",
            parent=base["Normal"],
            fontName="AvenirNextM",
            fontSize=_pt(7.6),
            leading=_pt(10),
            textColor=MUTED,
        ),
        "ruleBody": ParagraphStyle(
            "ruleBody",
            parent=base["Normal"],
            fontName="AvenirNext",
            fontSize=_pt(8.7),
            leading=_pt(11.5),
            textColor=INK,
            spaceAfter=3,
        ),
        "ruleBodyB": ParagraphStyle(
            "ruleBodyB",
            parent=base["Normal"],
            fontName="AvenirNextD",
            fontSize=_pt(8.7),
            leading=_pt(11.5),
            textColor=INK,
            spaceAfter=4,
        ),
        "tiny": ParagraphStyle(
            "tiny",
            parent=base["Normal"],
            fontName="AvenirNext",
            fontSize=_pt(8.0),
            leading=_pt(10.5),
            textColor=MUTED,
        ),
        "th": ParagraphStyle(
            "th",
            parent=base["Normal"],
            fontName="AvenirNextD",
            fontSize=_pt(8.0),
            leading=_pt(10.5),
            textColor=INK,
        ),
        "td": ParagraphStyle(
            "td",
            parent=base["Normal"],
            fontName="AvenirNext",
            fontSize=_pt(8.0),
            leading=_pt(10.5),
            textColor=INK,
        ),
        "cmd": ParagraphStyle(
            "cmd",
            parent=base["Normal"],
            fontName="Menlo",
            fontSize=_pt(7.5),
            leading=_pt(10),
            textColor=INK,
        ),
    }
    return styles


def _header_footer(canvas, doc):
    canvas.saveState()
    width, height = letter
    canvas.setFillColor(PAPER)
    canvas.rect(0, 0, width, height, fill=1, stroke=0)
    canvas.setFillColor(NAVY)
    canvas.rect(0, height - 32, width, 32, fill=1, stroke=0)
    canvas.setFillColor(white)
    canvas.setFont("AvenirNextM", _pt(8))
    canvas.drawString(0.65 * inch, height - 20, "ALLADA LAB  ·  fl.AI-reagent-stocker")
    canvas.drawRightString(width - 0.65 * inch, height - 20, f"SOP  ·  {SOP_VERSION}  ·  {SOP_DATE}")
    canvas.setFillColor(NAVY)
    canvas.rect(0, 0, width, 30, fill=1, stroke=0)
    canvas.setFillColor(white)
    canvas.setFont("AvenirNext", _pt(7.5))
    canvas.drawString(0.65 * inch, 12, "Never edit the user gene list. Convert with scripts/fetch_fbgn_ids.py only.")
    canvas.drawRightString(width - 0.65 * inch, 12, f"Page {doc.page}")
    canvas.restoreState()


def _meta_table(styles):
    rows = [
        [
            Paragraph("<b>Audience</b><br/>Anyone installing or running Stocker for the first time", styles["meta"]),
            Paragraph("<b>Default config</b><br/>priority_example.json", styles["meta"]),
            Paragraph(
                "<b>Validate IDs</b><br/><link href='https://flybase.org/convert/id' color='#1F3A3A'>"
                "flybase.org/convert/id</link>",
                styles["meta"],
            ),
        ]
    ]
    t = Table(rows, colWidths=[2.5 * inch, 2.4 * inch, 2.2 * inch])
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), NAVY_SOFT),
                ("BOX", (0, 0), (-1, -1), 0, NAVY_SOFT),
                ("LEFTPADDING", (0, 0), (-1, -1), 8),
                ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                ("TOPPADDING", (0, 0), (-1, -1), 7),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ]
        )
    )
    return t


def _kv_table(rows, styles, col_widths=None):
    data = [
        [Paragraph(a, styles["th" if i == 0 else "td"]), Paragraph(b, styles["th" if i == 0 else "td"])]
        for i, (a, b) in enumerate(rows)
    ]
    t = Table(data, colWidths=col_widths or [2.45 * inch, 4.65 * inch])
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), NAVY_SOFT),
                ("BACKGROUND", (0, 1), (-1, -1), CREAM),
                ("BOX", (0, 0), (-1, -1), 0.4, LINE),
                ("INNERGRID", (0, 0), (-1, -1), 0.3, LINE),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ]
        )
    )
    return t


def _code(text: str) -> str:
    return f"<font name='Menlo' size='{_pt(7.6)}'>{text}</font>"


def build():
    _register_fonts()
    styles = _styles()
    doc = SimpleDocTemplate(
        str(OUT_PATH),
        pagesize=letter,
        leftMargin=0.65 * inch,
        rightMargin=0.65 * inch,
        topMargin=0.58 * inch,
        bottomMargin=0.48 * inch,
        title="SOP: Gene list to Stocker stock sheets",
        author="Allada Lab",
        subject="Standalone SOP for fl.AI-reagent-stocker",
    )

    story = []
    story.append(Paragraph("STANDARD OPERATING PROCEDURE", styles["kicker"]))
    story.append(Paragraph("Gene list to Stocker stock sheets", styles["title"]))
    story.append(
        Paragraph(
            "This PDF is enough to install Stocker from GitHub, understand what it does, "
            "and run the lab default. The output is only as good as the gene list, the "
            "config you leave in place, and the stock limits that config applies.",
            styles["subtitle"],
        )
    )
    story.append(_meta_table(styles))
    story.append(Spacer(1, 8))

    story.append(Paragraph("What this document covers", styles["h"]))
    story.append(Hairline(NAVY, 0.8, 0, 5))
    story.append(
        Paragraph(
            "Three decisions determine whether the workbook is usable. Everything else "
            "(install, commands, troubleshooting) exists so those three decisions are done correctly.",
            styles["body"],
        )
    )
    for text in (
        "1. <b>The gene list is solid.</b> Every symbol maps to the intended current "
        "<i>D. melanogaster</i> FlyBase gene. Wrong FBgn IDs silently query the wrong genes.",
        "2. <b>The default config is understood.</b> The JSON file is the sorting protocol. "
        "It does not add genes or stocks. Changing a keyword or a sheet order changes the science.",
        "3. <b>Stock limits are understood.</b> The default keeps at most 5 stocks per gene "
        "and 3 per allele, counted across the whole workbook. Higher-priority sheets spend the quota first.",
    ):
        story.append(Paragraph(text, styles["step"]))

    story.append(Paragraph("Copy-paste URLs", styles["h"]))
    story.append(Hairline(NAVY, 0.8, 0, 4))
    story.append(
        CodeBlock(
            [
                f"Project (GitHub):     {GITHUB_URL}",
                f"ID validator:         {FLYBASE_URL}",
                f"Gene report example:  {FLYBASE_REPORT_EXAMPLE}",
                f"Python:               {PYTHON_URL}",
                f"Cursor (optional):    {CURSOR_URL}",
            ]
        )
    )
    story.append(Spacer(1, 4))
    story.append(
        Paragraph(
            f"Every FBgn in the review workbook uses this report pattern: "
            f"{_code('https://flybase.org/reports/FBgn########')} "
            f"(example Clock: {_code(FLYBASE_REPORT_EXAMPLE)}).",
            styles["body"],
        )
    )

    # --- How it works ---
    story.append(Paragraph("1. How Stocker works", styles["h"]))
    story.append(Hairline(NAVY, 0.8, 0, 5))
    story.append(
        Paragraph(
            "Stocker is a <i>Drosophila</i> reagent-lookup tool. You give it gene symbols. "
            "It returns the FlyBase stocks that perturb those genes, together with the papers "
            "already attached to those stocks, sorted into experimental categories.",
            styles["body"],
        )
    )
    story.append(
        Paragraph(
            "It does not invent stocks, call DEGs, or decide that a knockdown “worked.” "
            "It organizes records that already exist in a local FlyBase snapshot and in PubMed "
            "titles and abstracts. Availability and stock health still have to be checked at "
            "BDSC, VDRC, or the relevant center before you order.",
            styles["body"],
        )
    )
    story.append(
        Paragraph(
            "The path, in biological order:",
            styles["bodyL"],
        )
    )
    for text in (
        "1. <b>Identity.</b> Each symbol is matched to one current fly FBgn (for example Clock → FBgn0023076).",
        "2. <b>Reagents.</b> Alleles, constructs, and insertions of that FBgn are collected.",
        "3. <b>Stocks.</b> Orderable (and some custom) genotypes at BDSC, VDRC, and other centers.",
        "4. <b>Literature and phenotypes.</b> FlyBase citations are attached; PubMed supplies title and abstract. "
        "Curated FlyBase phenotypes are scored against the config’s phenotype words.",
        "5. <b>Sorting.</b> Each stock is tested against the config’s category list, in order, and kept in the "
        "<b>first</b> category it qualifies for.",
        "6. <b>Limits.</b> Across the whole workbook, at most 5 stocks per gene and 3 per allele are kept. "
        "Leftovers that hit the cap are dropped, not spilled to later sheets.",
    ):
        story.append(Paragraph(text, styles["step"]))
    story.append(Spacer(1, 3))
    story.append(
        Paragraph(
            "One gene-list CSV becomes one workbook. Two lists (for example two DEG contrasts) "
            "become two independent workbooks. You provide the list and the config. Stocker "
            "consults FlyBase and PubMed. You receive sorted sheets plus a counts table.",
            styles["body"],
        )
    )

    # --- Install ---
    story.append(Paragraph("2. One-time install from GitHub", styles["h"]))
    story.append(Hairline(NAVY, 0.8, 0, 5))
    story.append(
        Paragraph(
            "You need Python 3 and git. Every later command runs from the repo root "
            f"({_code('fl-AI-reagent-stocker')}), not from a Downloads folder or a PDF-only Cursor project. "
            f"If the GitHub repo is private, ask the maintainer for access first. "
            f"Source: {_code(GITHUB_URL)}.",
            styles["body"],
        )
    )
    story.append(
        Paragraph(
            "<b>A. Clone and install Python packages.</b> "
            f"{_code('python')} means Python 3; use {_code('python3')} if {_code('python')} is missing.",
            styles["bodyL"],
        )
    )
    story.append(
        CodeBlock(
            [
                f"git clone {GITHUB_URL}.git",
                "cd fl-AI-reagent-stocker",
                "python3 --version",
                "python3 -m pip install -r requirements.txt",
            ]
        )
    )
    story.append(Spacer(1, 4))
    story.append(
        Paragraph(
            "<b>B. Download the FlyBase citation table.</b> A GitHub clone includes most stock, "
            "allele, phenotype, and synonym tables, plus the helper files "
            f"{_code('fbst_to_derived_stock_component.csv')} and {_code('fbtp_to_fbti.csv')}. "
            "The large citation file "
            f"{_code('data/flybase/references/entity_publication_fb_*.tsv.gz')} is not in git. "
            "Without it, Stocker cannot attach papers. Run this once after clone (and again when "
            "you want a newer FlyBase release):",
            styles["body"],
        )
    )
    story.append(
        CodeBlock(
            [
                "python3 scripts/refresh_flybase_data.py",
            ]
        )
    )
    story.append(Spacer(1, 4))
    story.append(
        Paragraph(
            "That command pulls the current FlyBase TSV families into "
            f"{_code('data/flybase/')}. You do <b>not</b> need "
            f"{_code('--with-xml')} for a first run. The XML flag is only for rebuilding the helper CSVs, "
            "which already ship in the repo.",
            styles["body"],
        )
    )
    story.append(
        Paragraph(
            "<b>C. Optional keys in a repo-root</b> "
            f"{_code('.env')}. The default config turns <b>embeddings on</b> and "
            "<b>GPT validation off</b>. Stock sheets still build without a key. Cosine-similarity "
            "columns need "
            f"{_code('OPENAI_API_KEY')}. {_code('NCBI_API_KEY')} only raises PubMed rate limits.",
            styles["body"],
        )
    )
    story.append(
        CodeBlock(
            [
                "OPENAI_API_KEY=sk-...",
                "NCBI_API_KEY=...",
            ]
        )
    )
    story.append(Spacer(1, 4))
    story.append(
        Paragraph(
            f"Confirm {_code('data')}, {_code('scripts')}, {_code('fl_ai_reagent_stocker')}, and "
            f"{_code(DEFAULT_CONFIG)} are visible. If you were handed a complete folder on a shared drive, "
            "skip clone; still run the pip install and confirm "
            f"{_code('entity_publication')} exists under {_code('data/flybase/references/')}.",
            styles["body"],
        )
    )
    story.append(
        Paragraph(
            "Optional: open that same repo folder in Cursor and paste the commands below into "
            f"Agent chat or <b>View &gt; Terminal</b>. Cursor itself is at {_code(CURSOR_URL)}. "
            "The agent must run every Python command from the repo root.",
            styles["body"],
        )
    )

    # --- Gene list ---
    gene_intro = [
        Paragraph("3. Make the gene list solid", styles["h"]),
        Hairline(RULE, 0.8, 0, 5),
        Paragraph(
            "This is the hard gate. Do not run Stocker until a human has reviewed the FBgn IDs. "
            "A clean install and a correct config cannot rescue a list that points at the wrong genes.",
            styles["body"],
        ),
    ]
    story.append(KeepTogether(gene_intro))
    rule_body = [
        Paragraph(
            "Wrong FBgn IDs silently send Stocker at the wrong genes. Conversion must "
            f"<b>always</b> go through {_code('scripts/fetch_fbgn_ids.py')}.",
            styles["ruleBodyB"],
        ),
        Paragraph(
            "Never edit the user-generated gene list. The script writes new files "
            f"({_code('validated_&lt;original&gt;.csv')}, {_code('needs-review.csv')}, and "
            f"{_code('validated_&lt;original&gt;.xlsx')}). It does not overwrite the original.",
            styles["ruleBody"],
        ),
        Paragraph(
            "Never invent, recall, autocomplete, or type "
            f"{_code('flybase_gene_id')} / FBgn values. Do not copy IDs from memory, chat, "
            "FlyBase pages, or an old spreadsheet. Do not “fix” a dash by filling an ID. "
            f"Edit {_code('validated_*.csv')} or FBgn columns <b>only if the user explicitly asks</b>.",
            styles["ruleBody"],
        ),
        Paragraph(
            f"After conversion, review the xlsx (each FBgn is a hyperlink). If anything is unmatched, "
            f"paste those symbols at {_code(FLYBASE_URL)}. Stocker runs <b>only</b> on "
            f"{_code('validated_*.csv')}, and <b>only after that review</b>.",
            styles["ruleBody"],
        ),
    ]
    story.append(Callout("NON-NEGOTIABLE  ·  HUMANS AND AGENTS", rule_body))
    story.append(Spacer(1, 6))

    story.append(
        Paragraph(
            "<b>A. Copy the user’s list unchanged.</b> One new folder per request. "
            "Do not edit the user’s CSV. If you started from Excel, export a CSV and keep the workbook.",
            styles["bodyL"],
        )
    )
    story.append(
        CodeBlock(
            [
                "python3 -c \"from pathlib import Path; Path(",
                "    'data/gene_sets/YourName-Project-MMDDYYYY'",
                ").mkdir(parents=True, exist_ok=True)\"",
            ]
        )
    )
    story.append(Spacer(1, 4))
    story.append(
        Paragraph(
            f"One gene set = one CSV. First row must be exactly {_code('ext_gene')}. "
            "One Drosophila gene symbol per data row. No empty header, no extra title rows, "
            "no human gene symbols mixed in unless you intend to leave those unmatched.",
            styles["body"],
        )
    )

    story.append(
        Paragraph(
            "<b>B. Convert symbols to FBgn IDs.</b> Agents run this script. They do not fill IDs.",
            styles["bodyL"],
        )
    )
    story.append(
        CodeBlock(
            [
                'python3 scripts/fetch_fbgn_ids.py "data/gene_sets/YourName-Project-MMDDYYYY"',
            ]
        )
    )
    story.append(Spacer(1, 4))
    story.append(
        Paragraph(
            f"Reads source {_code('*.csv')} files (skips its own sidecars). "
            "Leaves the original file bytes unchanged. Writes:",
            styles["body"],
        )
    )
    story.append(
        CodeBlock(
            [
                "validated_<original>.csv     mapped rows only (FBgn IDs present)",
                "needs-review.csv             unmatched rows from every input + source_file",
                "validated_<original>.xlsx    review workbook: Validated + Needs review sheets",
            ]
        )
    )
    story.append(Spacer(1, 4))
    story.append(
        Paragraph(
            f"In the xlsx, each {_code('flybase_gene_id')} is a hyperlink to a "
            f"{_code(FLYBASE_REPORT_EXAMPLE)}-style report. "
            "<b>Tell the requester this xlsx exists and wait for them to review it for precision.</b> "
            f"Optional second argument is the symbol column (default {_code('ext_gene')}).",
            styles["body"],
        )
    )

    review_gate = [
        Paragraph(
            "<b>C. Human review gate. Stop here until the list is confirmed.</b>",
            styles["bodyL"],
        ),
        Paragraph(
            f"Open {_code('validated_&lt;original&gt;.xlsx')}. Click every FBgn. "
            "Confirm it is the intended current <i>D. melanogaster</i> gene, not a paralog, "
            "a human gene, or an old secondary ID that now points elsewhere.",
            styles["body"],
        ),
        Paragraph(
            f"If {_code('needs-review.csv')} has any data rows, paste those symbols at "
            f"{_code(FLYBASE_URL)} before anyone treats the list as ready. "
            "Enable “Match synonyms / secondary IDs” if a current symbol fails.",
            styles["body"],
        ),
        Paragraph(
            "Typical unmatched causes (fix on a <b>copy</b> of the user’s list, then re-run B; "
            "never edit the original in place):",
            styles["bodyL"],
        ),
    ]
    for text in (
        f"• RNA-seq uniquifiers ({_code('sky1')}, {_code('RpL231')}) — put the parent symbol "
        f"on a new working CSV, then re-run the script. Do not strip a digit that is part of "
        f"the real name ({_code('Argk1')}, {_code('miple1')}).",
        f"• Not a fly gene (for example {_code('EGFP')}) — leave it in needs-review; do not invent an ID.",
        f"• Old CG numbers / synonyms — the script usually fills {_code('primary_symbol')}. That is expected.",
        "• Mixed species or outdated symbols — resolve on FlyBase, then put the current fly symbol on the working CSV.",
    ):
        review_gate.append(Paragraph(text, styles["bullet"]))
    review_gate.append(
        Paragraph(
            "Do not start Stocker while a real fly gene remains unmatched. "
            "If FlyBase disagrees with the script, keep the script output, flag the row, and ask the maintainer.",
            styles["body"],
        )
    )
    story.append(KeepTogether(review_gate))

    # --- Default config ---
    story.append(Paragraph("4. Default config — what the knobs mean", styles["h"]))
    story.append(Hairline(NAVY, 0.8, 0, 5))
    story.append(
        Paragraph(
            f"Unless the requester names another JSON file, use {_code(DEFAULT_CONFIG)}. "
            "This file is the sorting protocol, not a data source. It does not add genes or stocks. "
            "It only tells Stocker which columns to read, what “topic-relevant” means, which "
            "phenotype buckets to fill, and how many stocks to keep.",
            styles["body"],
        )
    )
    story.append(
        Paragraph(
            "Do not swap this file for "
            f"{_code('stock_split_config_example.json')} (publication-evidence Ref matrix, effectively "
            "no stock caps) or "
            f"{_code('stock_split_config_priority_all_phenotypes.json')} (exhaustive coverage) "
            "unless that is what was asked for. Those files change both the sheets and the limits.",
            styles["body"],
        )
    )
    story.append(
        Paragraph(
            "Settings in the default file, in plain English:",
            styles["bodyL"],
        )
    )
    story.append(
        _kv_table(
            [
                ("Setting", "Value and meaning"),
                (
                    f"{_code('input.geneCol')}<br/>{_code('input.inputGeneCol')}",
                    f"{_code('flybase_gene_id')} and {_code('ext_gene')}. Stocker looks up reagents by FBgn "
                    "and keeps your original symbol for display.",
                ),
                (
                    f"{_code('input.skipFbgnidConversion')}",
                    f"<b>true</b>. Stage 1 trusts the validated CSV IDs. That is why section 3 exists: "
                    "conversion already happened outside Stocker. Do not set this false and also skip review.",
                ),
                (
                    f"{_code('relevantSearchTerms')}",
                    "<b>sleep, circadian, rhythm, locomotor</b> (case-insensitive). A paper is topic-relevant "
                    "if its title or abstract contains one of these words. That is a text match, not a claim "
                    "that the paper demonstrated the phenotype.",
                ),
                (
                    f"{_code('phenotypeSimilarityTargets')}",
                    "The same four words, each with an embedding phrase. Used for curated-phenotype tiers "
                    "(Phenotype++ if a FlyBase phenotype name matches a target; Phenotype+ if any curated "
                    "phenotype exists) and for optional cosine columns when embeddings run.",
                ),
                (
                    f"{_code('maxStocksPerGene')}<br/>{_code('maxStocksPerAllele')}",
                    "<b>5</b> and <b>3</b>. See section 5. These are the shortlist caps.",
                ),
                (
                    f"{_code('pubmed.batchSize')}",
                    "<b>50</b>. How many PubMed records are requested per batch. Does not change which papers count.",
                ),
                (
                    f"{_code('embeddings.enabled')}",
                    "<b>true</b> in this file. Adds cosine-similarity columns against the phenotype targets. "
                    f"Needs {_code('OPENAI_API_KEY')} for those columns. Leave it as shipped unless asked to change it.",
                ),
                (
                    f"{_code('validation.enabled')}",
                    "<b>false</b>. GPT does not read full text to decide whether a knockdown “worked.” "
                    "Do not turn this on unless the requester asks.",
                ),
                (
                    f"{_code('output.preserveUnsplit')}<br/>{_code('Workbook')}",
                    "<b>false</b>. The primary run keeps the organized workbook only, not the intermediate Stage 1 file.",
                ),
            ],
            styles,
        )
    )
    story.append(Spacer(1, 8))
    story.append(
        Paragraph(
            "Publication-evidence tiers (computed as "
            f"{_code('ALLELE_PAPER_RELEVANCE_SCORE')}, inherited by every stock of that allele):",
            styles["bodyL"],
        )
    )
    story.append(
        _kv_table(
            [
                ("Tier", "Meaning"),
                ("<b>Ref++</b> (score 2)", "At least one linked paper mentions a configured keyword in the title or abstract."),
                ("<b>Ref+</b> (score 1)", "The allele has papers, but none of those titles or abstracts hit the keywords."),
                ("<b>Ref-</b> (score 0)", "FlyBase lists no publications for the allele."),
            ],
            styles,
            col_widths=[1.6 * inch, 5.5 * inch],
        )
    )
    story.append(Spacer(1, 8))
    story.append(
        Paragraph(
            "This default does <b>not</b> split sheets by Ref++ / Ref+ / Ref-. It keeps stocks that "
            "have <b>any curated FlyBase phenotype</b> (Phenotype+ or Phenotype++), then bins them "
            "into six named buckets. Combinations are tested top to bottom. Each "
            f"({_code('stock_id')}, collection) reagent lands in only the first matching sheet.",
            styles["body"],
        )
    )
    story.append(
        _kv_table(
            [
                ("Sheet (tab name)", "Who gets in"),
                ("1  BDSC·RNAi·Pheno", "Bloomington RNAi reagent with a curated phenotype."),
                ("2  VDRC·RNAi·Pheno", "Vienna / VDRC RNAi reagent with a curated phenotype."),
                ("3  BDSC·Allele·Pheno", "Bloomington allele or insertion (not the broad RNAi/UAS proxy) with a phenotype."),
                ("4  OtherSC·Allele·Pheno", "Any other orderable stock-center allele/insertion with a phenotype."),
                ("5  Custom·RNAi·Pheno", "Phenotype-backed RNAi with no orderable FBst record."),
                ("6  Custom·Allele·Pheno", "Phenotype-backed allele/insertion with no orderable FBst record."),
            ],
            styles,
            col_widths=[2.2 * inch, 4.9 * inch],
        )
    )
    story.append(Spacer(1, 6))
    story.append(
        Paragraph(
            "Empty categories do not appear as tabs. A Bloomington RNAi stock that also matches "
            "the custom RNAi rule appears only on sheet 1. Changing sheet order changes which "
            "stocks spend the per-gene quota.",
            styles["body"],
        )
    )
    story.append(
        Paragraph(
            "Useful filter meanings in this file: "
            "<b>RNAi_reagent</b> is FlyBase "
            f"{_code('transgenic_product_class_terms')} containing "
            f"{_code('RNAi_reagent')}. <b>AlleleOrInsertion</b> is the complement of the repo’s "
            "broad UAS/RNAi proxy (genotype has UAS or RNAi, or a VDRC GD/KK/VSH knockdown). "
            "<b>Stock Center</b> vs <b>Non Stock Center</b> is whether an orderable FBst record exists. "
            "<b>Phenotype</b> is any non-zero curated-phenotype score.",
            styles["body"],
        )
    )

    # --- Limits ---
    story.append(Paragraph("5. Stock limits and how sheets are filled", styles["h"]))
    story.append(Hairline(NAVY, 0.8, 0, 5))
    story.append(
        Paragraph(
            "After identity, reagent expansion, and first-match sorting, the default config "
            "keeps a shortlist. These two numbers are the ones people miss:",
            styles["body"],
        )
    )
    story.append(
        _kv_table(
            [
                ("Limit", "What it does in the default file"),
                (
                    f"{_code('maxStocksPerGene')} = 5",
                    "Once a gene already has 5 accepted stocks anywhere in the workbook, later matching stocks of that gene are dropped.",
                ),
                (
                    f"{_code('maxStocksPerAllele')} = 3",
                    "Once an allele already has 3 accepted stocks, further stocks of that same allele are dropped even if the gene still has room.",
                ),
            ],
            styles,
        )
    )
    story.append(Spacer(1, 8))
    story.append(
        Paragraph(
            "How the quota is spent:",
            styles["bodyL"],
        )
    )
    for text in (
        "1. Stocks are sorted so stronger publication evidence is considered first "
        "(Ref++ before Ref+, then more keyword-hitting papers, then more total papers).",
        "2. Sheet 1 is filled first, then sheet 2, and so on. A stock accepted on an earlier "
        "sheet never appears again.",
        "3. Gene and allele counters are <b>shared across all six sheets</b>. Five Bloomington "
        "RNAi stocks for <i>Clk</i> on sheet 1 means later sheets get no further <i>Clk</i> stocks.",
        "4. A stock that matched sheet 1 but was rejected because the gene or allele was already "
        "at the cap is not spilled onto sheet 2. Leftovers are omitted.",
        "5. The Contents sheet and "
        f"{_code('combination_counts_summary.xlsx')} count what survived these limits, not every "
        "FlyBase stock that theoretically matched.",
    ):
        story.append(Paragraph(text, styles["step"]))
    story.append(Spacer(1, 3))
    story.append(
        Paragraph(
            "Experimentally: you receive a prioritized shortlist, not an inventory. If the "
            "requester wants every phenotype-backed reagent, they must name a different config "
            "(the exhaustive "
            f"{_code('stock_split_config_priority_all_phenotypes.json')} file lifts these caps) "
            f"or raise {_code('maxStocksPerGene')} / {_code('maxStocksPerAllele')} on a copy of "
            "the default JSON. Do not silently lift the caps.",
            styles["body"],
        )
    )
    story.append(
        Paragraph(
            "Limits do not check whether BDSC or VDRC currently ships the stock. They only cap "
            "how many FlyBase records are written into the workbook.",
            styles["body"],
        )
    )

    # --- Run ---
    story.append(Paragraph("6. Run Stocker on the validated list", styles["h"]))
    story.append(Hairline(NAVY, 0.8, 0, 5))
    story.append(
        Paragraph(
            f"Run only after section 3 is done. Stay with {_code(DEFAULT_CONFIG)} unless another "
            f"JSON file was named. Discovery skips the original CSV and {_code('needs-review.csv')} "
            f"when {_code('validated_*.csv')} is present.",
            styles["body"],
        )
    )
    story.append(
        CodeBlock(
            [
                "python3 -m fl_ai_reagent_stocker run \\",
                '  "data/gene_sets/YourName-Project-MMDDYYYY" \\',
                f"  --config {DEFAULT_CONFIG}",
            ]
        )
    )
    story.append(Spacer(1, 4))
    story.append(
        Paragraph(
            "Optional recount after a partial rerun:",
            styles["bodyL"],
        )
    )
    story.append(
        CodeBlock(
            [
                "python3 scripts/summarize_combination_counts.py \\",
                '  "data/gene_sets/YourName-Project-MMDDYYYY" \\',
                f"  --config {DEFAULT_CONFIG}",
            ]
        )
    )

    story.append(Paragraph("7. What you should have when it finishes", styles["h"]))
    story.append(Hairline(NAVY, 0.8, 0, 5))
    story.append(
        Paragraph(
            "Keep the original CSV, the validated CSV, the review xlsx, and needs-review.csv. "
            "They are the record of what was queried and what was checked.",
            styles["body"],
        )
    )
    story.append(
        CodeBlock(
            [
                "data/gene_sets/<run>/<original>.csv                         (untouched user list)",
                "data/gene_sets/<run>/validated_<original>.csv               (stocker input)",
                "data/gene_sets/<run>/validated_<original>.xlsx              (human ID review)",
                "data/gene_sets/<run>/needs-review.csv",
                "data/gene_sets/<run>/Per Gene Set Runs/<validated name>/",
                "    Stocks/<validated name>_aggregated.xlsx",
                "data/gene_sets/<run>/combination_counts_summary.xlsx",
            ]
        )
    )
    story.append(Spacer(1, 4))
    story.append(
        Paragraph(
            "Share the aggregated workbook and the counts table. Inside the workbook: "
            "<b>Contents</b>, one tab per populated category from section 4, <b>References</b> "
            "(papers cited by stocks that made it through the limits), and "
            "<b>Stock Sheet by Gene</b>. Cosine columns appear when embeddings ran with a key. "
            "Check BDSC, VDRC, or the relevant center before ordering.",
            styles["body"],
        )
    )

    story.append(Paragraph("8. If something fails", styles["h"]))
    story.append(Hairline(NAVY, 0.8, 0, 5))
    fail_rows = [
        [
            Paragraph("<b>Symptom</b>", styles["th"]),
            Paragraph("<b>What to do</b>", styles["th"]),
        ],
        [
            Paragraph("No project on disk", styles["td"]),
            Paragraph(
                f"Install from GitHub (section 2). Source: {GITHUB_URL}. A code-only zip still needs "
                f"{_code('refresh_flybase_data.py')} so {_code('entity_publication')} is present.",
                styles["td"],
            ),
        ],
        [
            Paragraph("git clone fails (private repo)", styles["td"]),
            Paragraph("Ask the maintainer for GitHub access, or take a complete folder from a shared drive and skip clone.", styles["td"]),
        ],
        [
            Paragraph("Run dies looking up papers / missing entity_publication", styles["td"]),
            Paragraph(
                f"The citation table is not in git. From the repo root: {_code('python3 scripts/refresh_flybase_data.py')}.",
                styles["td"],
            ),
        ],
        [
            Paragraph("“No gene-list CSV files found”", styles["td"]),
            Paragraph("Wrong path, or files are not named .csv. Run from the repo root.", styles["td"]),
        ],
        [
            Paragraph("“Missing required columns”", styles["td"]),
            Paragraph(
                f"Stocker needs {_code('ext_gene')} and {_code('flybase_gene_id')} on the validated CSV. Re-run conversion; do not type IDs.",
                styles["td"],
            ),
        ],
        [
            Paragraph("needs-review.csv has rows", styles["td"]),
            Paragraph(
                f"Do not run Stocker. Review the xlsx and {FLYBASE_URL}. Fix symbols on a new working CSV, then re-run the script.",
                styles["td"],
            ),
        ],
        [
            Paragraph("Workbook has far fewer stocks than FlyBase shows", styles["td"]),
            Paragraph(
                "Expected under the default caps (5 / gene, 3 / allele) and first-match sheets. "
                "See section 5. Do not lift limits unless asked.",
                styles["td"],
            ),
        ],
        [
            Paragraph("python: command not found", styles["td"]),
            Paragraph(
                f"Use {_code('python3')} or {_code('python3 -m fl_ai_reagent_stocker')} so the same interpreter is used.",
                styles["td"],
            ),
        ],
        [
            Paragraph("Embedding / OpenAI errors", styles["td"]),
            Paragraph(
                f"Stock sheets can still be used. Cosine columns need {_code('OPENAI_API_KEY')} in {_code('.env')}. "
                "Do not turn on validation to “fix” this.",
                styles["td"],
            ),
        ],
    ]
    ft = Table(fail_rows, colWidths=[2.35 * inch, 4.75 * inch])
    ft.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), NAVY_SOFT),
                ("BACKGROUND", (0, 1), (-1, -1), CREAM),
                ("BOX", (0, 0), (-1, -1), 0.4, LINE),
                ("INNERGRID", (0, 0), (-1, -1), 0.3, LINE),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ]
        )
    )
    story.append(ft)

    agent_block = [
        Paragraph("9. Agent command policy", styles["h"]),
        Hairline(NAVY, 0.8, 0, 5),
    ]
    for text in (
        "• Run only the commands in this SOP, from the Stocker repo root.",
        f"• Never edit the user-generated gene list. Never write FBgn IDs by hand. Always call {_code('scripts/fetch_fbgn_ids.py')}.",
        f"• Tell the user that {_code('validated_*.xlsx')} exists and wait for them to review it. "
        f"If {_code('needs-review.csv')} is non-empty, send them to {_code(FLYBASE_URL)} first.",
        f"• Run Stocker only on {_code('validated_*.csv')}, and only after that review.",
        f"• Do not edit {_code('validated_*.csv')} or FBgn columns unless the user explicitly asks.",
        f"• Use {_code(DEFAULT_CONFIG)} unless another JSON file is named. Do not lift stock limits, "
        "reorder sheets, or turn on GPT validation unless asked. Leave embeddings as shipped.",
        "• Do not invent stocks or paper conclusions. If a command fails, show the error. "
        "Do not invent a workaround that writes gene IDs.",
    ):
        agent_block.append(Paragraph(text, styles["bullet"]))
    agent_block.append(Spacer(1, 6))
    agent_block.append(
        Paragraph(
            "Optional further reading in the repo: "
            f"{_code('docs/how_stocker_works.md')}, "
            f"{_code('docs/pipeline_usage.md')}, "
            f"{_code('docs/stock_split_config_priority_example.md')}, "
            f"{_code('docs/AGENTS.md')}.",
            styles["tiny"],
        )
    )
    story.append(KeepTogether(agent_block))

    doc.build(story, onFirstPage=_header_footer, onLaterPages=_header_footer)
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    build()

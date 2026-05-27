"""Render docs/stock_split_config_example.md to a clean printable PDF.

Output: docs/stock_split_config_example.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import markdown
from weasyprint import HTML, CSS


REPO_ROOT = Path(__file__).resolve().parents[1]
MD_PATH = REPO_ROOT / "docs" / "stock_split_config_example.md"
PDF_PATH = REPO_ROOT / "docs" / "stock_split_config_example.pdf"


CSS_STYLES = """
@page {
    size: Letter;
    margin: 1in 0.9in 1in 0.9in;
    @bottom-right {
        content: counter(page) " / " counter(pages);
        font-family: "Helvetica", "Arial", sans-serif;
        font-size: 9pt;
        color: #666;
    }
    @bottom-left {
        content: "Stock-splitting config walkthrough";
        font-family: "Helvetica", "Arial", sans-serif;
        font-size: 9pt;
        color: #666;
    }
}

html, body {
    font-family: "Helvetica", "Arial", sans-serif;
    font-size: 10.5pt;
    line-height: 1.45;
    color: #1a1a1a;
}

h1 {
    font-size: 20pt;
    margin-top: 0;
    margin-bottom: 0.4em;
    color: #111;
    border-bottom: 2px solid #888;
    padding-bottom: 0.2em;
}

h2 {
    font-size: 14.5pt;
    margin-top: 1.4em;
    margin-bottom: 0.3em;
    color: #111;
    page-break-after: avoid;
}

h3 {
    font-size: 12pt;
    margin-top: 1em;
    margin-bottom: 0.25em;
    color: #222;
    page-break-after: avoid;
}

p, ul, ol {
    margin-top: 0.3em;
    margin-bottom: 0.6em;
}

ul, ol {
    padding-left: 1.4em;
}

li {
    margin-bottom: 0.15em;
}

code {
    font-family: "Menlo", "Consolas", monospace;
    font-size: 9.5pt;
    background: #f3f3f3;
    padding: 0.05em 0.3em;
    border-radius: 3px;
    color: #b03030;
    word-break: break-word;
}

table {
    border-collapse: collapse;
    width: 100%;
    margin: 0.6em 0 1em 0;
    font-size: 9.5pt;
    page-break-inside: auto;
}

thead {
    display: table-header-group;
}

tr {
    page-break-inside: avoid;
    page-break-after: auto;
}

th, td {
    border: 1px solid #cccccc;
    padding: 5px 8px;
    text-align: left;
    vertical-align: top;
}

th {
    background: #f0f0f0;
    font-weight: 600;
}

tbody tr:nth-child(even) {
    background: #fafafa;
}

blockquote {
    border-left: 3px solid #888;
    margin: 0.6em 0 0.6em 0;
    padding: 0.1em 0.8em;
    color: #333;
    background: #f7f7f7;
}

hr {
    border: none;
    border-top: 1px solid #ddd;
    margin: 1.2em 0;
}

strong {
    color: #111;
}

.subtitle {
    font-size: 12pt;
    color: #555;
    margin-top: -0.4em;
    margin-bottom: 1.5em;
}
"""


def build_html(md_text: str) -> str:
    md = markdown.Markdown(
        extensions=[
            "extra",
            "tables",
            "sane_lists",
            "meta",
            "toc",
        ],
        output_format="html5",
    )
    body_html = md.convert(md_text)

    meta = getattr(md, "Meta", {}) or {}
    title = " ".join(meta.get("title", ["Stock-splitting config walkthrough"]))
    subtitle = " ".join(meta.get("subtitle", []))
    subtitle_block = f'<p class="subtitle">{subtitle}</p>' if subtitle else ""

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
</head>
<body>
{subtitle_block}
{body_html}
</body>
</html>
"""


def main() -> int:
    if not MD_PATH.exists():
        print(f"Source markdown not found: {MD_PATH}", file=sys.stderr)
        return 1

    md_text = MD_PATH.read_text(encoding="utf-8")
    html_text = build_html(md_text)

    HTML(string=html_text, base_url=str(MD_PATH.parent)).write_pdf(
        target=str(PDF_PATH),
        stylesheets=[CSS(string=CSS_STYLES)],
    )

    size_kb = PDF_PATH.stat().st_size / 1024
    print(f"Wrote {PDF_PATH.relative_to(REPO_ROOT)} ({size_kb:.1f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

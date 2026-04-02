from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "docs" / "images" / "pipeline-data-flowchart.png"

FONT_REGULAR_PATH = Path("/System/Library/Fonts/Supplemental/Arial.ttf")
FONT_BOLD_PATH = Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf")

CANVAS_WIDTH = 2600
CANVAS_HEIGHT = 1550

BG = "#f7f9fc"
TEXT = "#172033"
MUTED = "#4d5b75"
LINE = "#58708f"
SOFT_BORDER = "#b9c8dc"

BLUE = "#2d6fb7"
BLUE_FILL = "#dfeeff"
TEAL = "#1f8f98"
TEAL_FILL = "#ddf6f7"
PURPLE = "#7a4eb2"
PURPLE_FILL = "#efe5fb"
AMBER = "#d39229"
AMBER_FILL = "#fff0cf"
GREEN = "#3b9556"
GREEN_FILL = "#e1f4e6"
SLATE_FILL = "#edf2f8"


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    font_path = FONT_BOLD_PATH if bold else FONT_REGULAR_PATH
    if font_path.exists():
        return ImageFont.truetype(str(font_path), size=size)
    return ImageFont.load_default()


TITLE_FONT = load_font(54, bold=True)
SUBTITLE_FONT = load_font(28)
BANNER_FONT = load_font(30, bold=True)
CHIP_FONT = load_font(24, bold=True)
CARD_TITLE_FONT = load_font(28, bold=True)
CARD_BODY_FONT = load_font(22)
CARD_SMALL_FONT = load_font(19)
FOOTNOTE_FONT = load_font(20)


def wrap_text(draw: ImageDraw.ImageDraw, text: str, font, max_width: int) -> list[str]:
    words = text.split()
    if not words:
        return [""]
    lines: list[str] = []
    current = words[0]
    for word in words[1:]:
        candidate = f"{current} {word}"
        bbox = draw.textbbox((0, 0), candidate, font=font)
        if bbox[2] - bbox[0] <= max_width:
            current = candidate
        else:
            lines.append(current)
            current = word
    lines.append(current)
    return lines


def draw_centered_text(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], text: str, font, fill: str) -> None:
    lines = text.split("\n")
    metrics = [draw.textbbox((0, 0), line, font=font) for line in lines]
    line_heights = [bbox[3] - bbox[1] for bbox in metrics]
    total_height = sum(line_heights) + (len(lines) - 1) * 6
    x0, y0, x1, y1 = box
    cursor_y = y0 + (y1 - y0 - total_height) / 2
    for line, bbox, line_height in zip(lines, metrics, line_heights):
        width = bbox[2] - bbox[0]
        cursor_x = x0 + (x1 - x0 - width) / 2
        draw.text((cursor_x, cursor_y), line, font=font, fill=fill)
        cursor_y += line_height + 6


def draw_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    body_lines: list[str],
    outline: str,
    fill: str,
    header_fill: str,
    bullet: bool = True,
    footer_lines: list[str] | None = None,
) -> None:
    x0, y0, x1, y1 = box
    header_height = 58
    radius = 26
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=4)
    draw.rounded_rectangle((x0, y0, x1, y0 + header_height), radius=radius, fill=header_fill, outline=outline, width=4)
    draw.rectangle((x0, y0 + header_height - 10, x1, y0 + header_height), fill=header_fill, outline=header_fill)

    draw_centered_text(draw, (x0 + 12, y0 + 6, x1 - 12, y0 + header_height), title, CARD_TITLE_FONT, "white")

    cursor_y = y0 + header_height + 24
    body_width = x1 - x0 - 52
    bullet_indent = 28

    for line in body_lines:
        wrapped = wrap_text(draw, line, CARD_BODY_FONT, body_width - (bullet_indent if bullet else 0))
        for idx, segment in enumerate(wrapped):
            if bullet and idx == 0:
                draw.ellipse((x0 + 22, cursor_y + 9, x0 + 34, cursor_y + 21), fill=outline)
                draw.text((x0 + 44, cursor_y), segment, font=CARD_BODY_FONT, fill=TEXT)
            else:
                offset = 44 if bullet else 24
                draw.text((x0 + offset, cursor_y), segment, font=CARD_BODY_FONT, fill=TEXT)
            cursor_y += 34
        cursor_y += 8

    if footer_lines:
        cursor_y += 2
        for line in footer_lines:
            wrapped = wrap_text(draw, line, CARD_SMALL_FONT, body_width)
            for segment in wrapped:
                draw.text((x0 + 24, cursor_y), segment, font=CARD_SMALL_FONT, fill=MUTED)
                cursor_y += 28
            cursor_y += 4


def draw_file_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    lines: list[str],
    outline: str,
    fill: str,
) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=26, fill=fill, outline=outline, width=4)
    draw.text((x0 + 24, y0 + 20), title, font=CARD_TITLE_FONT, fill=TEXT)
    cursor_y = y0 + 76
    width = x1 - x0 - 44
    for line in lines:
        wrapped = wrap_text(draw, line, CARD_BODY_FONT, width)
        for segment in wrapped:
            draw.text((x0 + 24, cursor_y), segment, font=CARD_BODY_FONT, fill=TEXT)
            cursor_y += 32
        cursor_y += 6


def draw_chip(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], label: str, fill: str, outline: str) -> None:
    draw.rounded_rectangle(box, radius=22, fill=fill, outline=outline, width=3)
    draw_centered_text(draw, box, label, CHIP_FONT, outline)


def draw_poly_arrow(draw: ImageDraw.ImageDraw, points: list[tuple[int, int]], color: str, width: int = 6) -> None:
    draw.line(points, fill=color, width=width)
    if len(points) < 2:
        return
    (x0, y0), (x1, y1) = points[-2], points[-1]
    head = 16
    if abs(x1 - x0) >= abs(y1 - y0):
        if x1 >= x0:
            arrow = [(x1, y1), (x1 - head, y1 - 8), (x1 - head, y1 + 8)]
        else:
            arrow = [(x1, y1), (x1 + head, y1 - 8), (x1 + head, y1 + 8)]
    else:
        if y1 >= y0:
            arrow = [(x1, y1), (x1 - 8, y1 - head), (x1 + 8, y1 - head)]
        else:
            arrow = [(x1, y1), (x1 - 8, y1 + head), (x1 + 8, y1 + head)]
    draw.polygon(arrow, fill=color)


def main(output_path: Path = DEFAULT_OUTPUT) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    image = Image.new("RGB", (CANVAS_WIDTH, CANVAS_HEIGHT), BG)
    draw = ImageDraw.Draw(image)

    draw.text((88, 52), "CLI Data Flow: fl_ai_reagent_stocker", font=TITLE_FONT, fill=TEXT)
    draw.text(
        (92, 122),
        "How data moves through find-stocks, split-stocks, validate-stocks, and the run-full-pipeline wrapper",
        font=SUBTITLE_FONT,
        fill=MUTED,
    )

    banner = (640, 168, 1970, 248)
    draw.rounded_rectangle(banner, radius=28, fill=SLATE_FILL, outline=SOFT_BORDER, width=3)
    draw_centered_text(
        draw,
        banner,
        "python -m fl_ai_reagent_stocker <command> [options]",
        BANNER_FONT,
        TEXT,
    )

    chip_y0 = 274
    chip_y1 = 330
    draw_chip(draw, (690, chip_y0, 930, chip_y1), "find-stocks", BLUE_FILL, BLUE)
    draw_chip(draw, (960, chip_y0, 1200, chip_y1), "split-stocks", TEAL_FILL, TEAL)
    draw_chip(draw, (1230, chip_y0, 1510, chip_y1), "validate-stocks", PURPLE_FILL, PURPLE)
    draw_chip(draw, (1540, chip_y0, 1930, chip_y1), "run-full-pipeline = 1 -> 2 -> 3", AMBER_FILL, AMBER)

    gene_box = (70, 410, 320, 590)
    flybase_box = (70, 700, 320, 980)
    config_box = (1090, 360, 1520, 470)
    stage1_box = (380, 390, 1010, 980)
    workbook_box = (1090, 560, 1520, 860)
    stage2_box = (1590, 390, 2130, 780)
    stage3_box = (1590, 840, 2130, 1330)
    services_box = (1110, 1160, 1500, 1430)
    outputs_box = (2190, 560, 2550, 1080)
    notes_box = (70, 1110, 1010, 1450)

    draw_file_card(
        draw,
        gene_box,
        "Input gene list CSVs",
        [
            "./gene_lists/*.csv",
            "One or more CSVs with genes of interest.",
            "find-stocks reads these first.",
        ],
        BLUE,
        BLUE_FILL,
    )

    draw_file_card(
        draw,
        flybase_box,
        "Reference data tables",
        [
            "data/flybase/alleles_and_stocks/",
            "data/flybase/references/",
            "data/flybase/transgenic_*/",
            "PubMed metadata cache / API",
        ],
        GREEN,
        GREEN_FILL,
    )

    draw_file_card(
        draw,
        config_box,
        "Config JSON",
        [
            "data/config/*.json",
            "Stage 1 keywords plus Stage 2/3 split rules.",
        ],
        AMBER,
        AMBER_FILL,
    )

    draw_card(
        draw,
        stage1_box,
        "Stage 1: find-stocks",
        [
            "Read CSV gene lists from the input folder.",
            "Optionally convert gene symbols to FBgn IDs.",
            "Map FBgn -> FBal / FBtp / FBti -> FBst using FlyBase tables.",
            "Link matched stocks to publications and PMIDs.",
            "Score keyword relevance from title and abstract metadata.",
        ],
        BLUE,
        BLUE_FILL,
        BLUE,
    )

    draw_file_card(
        draw,
        workbook_box,
        "Stage 1 outputs",
        [
            "./gene_lists/Stocks/",
            "aggregated_stock_refs.xlsx",
            "Sheets: Stocks, References",
            "Sidecar report:",
            "references_without_pmid_fbrf.txt",
        ],
        BLUE,
        "#ebf3ff",
    )

    draw_card(
        draw,
        stage2_box,
        "Stage 2: split-stocks",
        [
            "Load the Stage 1 workbook from ./gene_lists/Stocks.",
            "Compute derived columns used by stock ranking and grouping.",
            "Apply JSON filters, combinations, and per-gene / per-allele limits.",
            "Keep only references cited by stocks that survive into output sheets.",
            "Write organized workbooks without GPT validation.",
        ],
        TEAL,
        TEAL_FILL,
        TEAL,
    )

    draw_card(
        draw,
        stage3_box,
        "Stage 3: validate-stocks",
        [
            "Reload the same Stage 1 workbook and the same JSON config.",
            "Rebuild the split outputs with the same filter and limit logic.",
            "Select only Ref++ stocks that actually appear in output sheets.",
            "Retrieve paper text when available, then run GPT functional validation.",
            "Merge validation columns back into the organized workbook(s).",
        ],
        PURPLE,
        PURPLE_FILL,
        PURPLE,
    )

    draw_file_card(
        draw,
        services_box,
        "Stage 3 services",
        [
            "OPENAI_API_KEY / OPENAI_MODEL",
            "UNPAYWALL_TOKEN",
            "PMC / Crossref / PubMed full-text lookups",
        ],
        PURPLE,
        "#f5effd",
    )

    draw_file_card(
        draw,
        outputs_box,
        "Organized outputs",
        [
            "./gene_lists/Stocks/",
            "Organized Stocks/",
            "*_aggregated.xlsx",
            "Output sheets from config combinations",
            "References sheet narrowed to cited PMIDs",
            "Validation columns appear only after Stage 3",
            "references_without_pmid_fbrf.txt is moved here",
        ],
        TEAL,
        "#ecfbfb",
    )

    draw_file_card(
        draw,
        notes_box,
        "Wrapper behavior and important CLI details",
        [
            "run-full-pipeline runs Stage 1 -> Stage 2 -> Stage 3 in order.",
            "split-stocks and validate-stocks both read ./gene_lists/Stocks as input.",
            "validate-stocks rebuilds split outputs instead of using Stage 2 as its source of truth.",
            "The split config controls both organization and Ref++ gating.",
        ],
        LINE,
        "#f1f5fb",
    )

    # Main arrows
    draw_poly_arrow(draw, [(320, 500), (380, 500)], BLUE)
    draw_poly_arrow(draw, [(320, 820), (350, 820), (350, 700), (380, 700)], GREEN)
    draw_poly_arrow(draw, [(1305, 470), (1305, 520), (1100, 520), (1100, 590)], AMBER)
    draw_poly_arrow(draw, [(1010, 670), (1090, 670)], BLUE)
    draw_poly_arrow(draw, [(1520, 710), (1590, 710)], TEAL)
    draw_poly_arrow(draw, [(1520, 760), (1550, 760), (1550, 980), (1590, 980)], PURPLE)
    draw_poly_arrow(draw, [(1500, 1290), (1590, 1290)], PURPLE)
    draw_poly_arrow(draw, [(2130, 670), (2190, 670)], TEAL)
    draw_poly_arrow(draw, [(2130, 1030), (2170, 1030), (2170, 790), (2190, 790)], PURPLE)

    # Config arrows into stages
    draw_poly_arrow(draw, [(1090, 415), (1040, 415), (1040, 450), (1010, 450)], AMBER)
    draw_poly_arrow(draw, [(1520, 415), (1555, 415), (1555, 470), (1590, 470)], AMBER)
    draw_poly_arrow(draw, [(1520, 440), (1565, 440), (1565, 920), (1590, 920)], AMBER)

    # Top wrapper line
    draw.line([(770, 352), (1930, 352)], fill=AMBER, width=4)
    draw_poly_arrow(draw, [(930, 352), (930, 390)], AMBER)
    draw_poly_arrow(draw, [(1080, 352), (1865, 352), (1865, 390)], AMBER)

    image.save(output_path, format="PNG")
    print(f"Saved flowchart to {output_path}")


if __name__ == "__main__":
    main()

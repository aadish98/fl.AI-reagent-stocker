#!/usr/bin/env python3
"""Fix python-pptx OOXML issues that trigger PowerPoint repair on macOS."""
import re
import zipfile
from pathlib import Path


def _patch_sldsz(pres_xml: str) -> str:
    m = re.search(r"<p:sldSz\s+([^/]+?)/>", pres_xml)
    if not m:
        return pres_xml
    attrs = m.group(1)
    cx = int(re.search(r'cx="(\d+)"', attrs).group(1))
    cy = int(re.search(r'cy="(\d+)"', attrs).group(1))
    if abs(cx - 12192000) <= 60000 and abs(cy - 6858000) <= 60000:
        new_el = '<p:sldSz cx="12192000" cy="6858000" type="screen16x9"/>'
    elif abs(cx - 9144000) <= 60000 and abs(cy - 6858000) <= 60000:
        new_el = '<p:sldSz cx="9144000" cy="6858000" type="screen4x3"/>'
    else:
        new_attrs = re.sub(r'\s+type="[^"]+"', "", attrs).strip()
        new_el = f"<p:sldSz {new_attrs}/>"
    return pres_xml.replace(m.group(0), new_el)


def _next_rid(rels_xml: str) -> str:
    ids = [int(m.group(1)) for m in re.finditer(r'Id="rId(\d+)"', rels_xml)]
    return f"rId{max(ids) + 1 if ids else 1}"


def _strip_notes_slides(entries: dict) -> dict:
    for name in [n for n in entries if n.startswith("ppt/notesSlides/")]:
        del entries[name]

    for name in list(entries.keys()):
        if name.startswith("ppt/slides/_rels/") and name.endswith(".rels"):
            rels = entries[name].decode("utf-8")
            new_rels = re.sub(r'<Relationship[^>]*notesSlide[^>]*/>\s*', "", rels)
            if new_rels != rels:
                entries[name] = new_rels.encode("utf-8")

    ct_path = "[Content_Types].xml"
    if ct_path in entries:
        ct = entries[ct_path].decode("utf-8")
        entries[ct_path] = re.sub(
            r'<Override PartName="/ppt/notesSlides/[^"]+"[^>]*/>\s*',
            "",
            ct,
        ).encode("utf-8")
    return entries


def _scrub_content_types(entries: dict) -> dict:
    """Drop [Content_Types] overrides for parts that are not in the package."""
    ct_path = "[Content_Types].xml"
    if ct_path not in entries:
        return entries
    ct = entries[ct_path].decode("utf-8")

    def _keep(m: re.Match) -> str:
        part = m.group(1).lstrip("/")
        return m.group(0) if part in entries else ""

    ct = re.sub(r'<Override PartName="([^"]+)"[^>]*/>\s*', _keep, ct)
    entries[ct_path] = ct.encode("utf-8")
    return entries


def _strip_notes_master(entries: dict) -> dict:
    for name in [n for n in entries if n.startswith("ppt/notesMasters/")]:
        del entries[name]

    rels_path = "ppt/_rels/presentation.xml.rels"
    if rels_path in entries:
        rels = entries[rels_path].decode("utf-8")
        entries[rels_path] = re.sub(
            r'<Relationship[^>]*notesMaster[^>]*/>', "", rels
        ).encode("utf-8")

    pres_path = "ppt/presentation.xml"
    if pres_path in entries:
        pres = entries[pres_path].decode("utf-8")
        entries[pres_path] = re.sub(
            r"<p:notesMasterIdLst>.*?</p:notesMasterIdLst>", "", pres
        ).encode("utf-8")

    ct_path = "[Content_Types].xml"
    if ct_path in entries:
        ct = entries[ct_path].decode("utf-8")
        entries[ct_path] = re.sub(
            r'<Override PartName="/ppt/notesMasters/notesMaster1\.xml"[^>]*/>',
            "",
            ct,
        ).encode("utf-8")
    return entries


def _strip_printer_settings(entries: dict) -> dict:
    to_drop = [n for n in entries if n.startswith("ppt/printerSettings/")]
    for name in to_drop:
        del entries[name]
    if not to_drop:
        return entries

    rels_path = "ppt/_rels/presentation.xml.rels"
    if rels_path in entries:
        rels = entries[rels_path].decode("utf-8")
        entries[rels_path] = re.sub(
            r'<Relationship[^>]*printerSettings[^>]*/>', "", rels
        ).encode("utf-8")

    ct_path = "[Content_Types].xml"
    if ct_path in entries:
        ct = entries[ct_path].decode("utf-8")
        entries[ct_path] = re.sub(
            r'<Default Extension="bin"[^>]*/>', "", ct
        ).encode("utf-8")
    return entries


def _fix_empty_table_paragraphs(xml: str) -> str:
    """pptxgenjs emits text-less table cells that trigger PowerPoint repair."""
    pat = re.compile(
        r"(<a:txBody>.*?<a:p><a:pPr[^>]*><a:buNone/></a:pPr>)<a:endParaRPr",
        re.S,
    )

    def _repl(m: re.Match) -> str:
        return m.group(1) + '<a:r><a:rPr lang="en-US" sz="100" dirty="0"/><a:t xml:space="preserve"> </a:t></a:r><a:endParaRPr'

    return pat.sub(_repl, xml)


def _fix_duplicate_shape_ids(xml: str) -> str:
    """pptxgenjs can reuse cNvPr ids when tables are added after text shapes."""
    seen: set[str] = set()
    max_id = 0

    def _fix_id(m: re.Match) -> str:
        nonlocal max_id
        sid = m.group(1)
        max_id = max(max_id, int(sid))
        if sid not in seen:
            seen.add(sid)
            return m.group(0)
        max_id += 1
        new_id = str(max_id)
        seen.add(new_id)
        return f'<p:cNvPr id="{new_id}"'

    return re.sub(r'<p:cNvPr id="(\d+)"', _fix_id, xml)


def _patch_slide_xml(entries: dict) -> dict:
    for name in list(entries.keys()):
        if not name.startswith("ppt/slides/slide") or not name.endswith(".xml"):
            continue
        xml = entries[name].decode("utf-8")
        xml = _fix_empty_table_paragraphs(xml)
        xml = _fix_duplicate_shape_ids(xml)
        entries[name] = xml.encode("utf-8")
    return entries


def _ordered_names(names: list[str]) -> list[str]:
    rest = sorted(n for n in names if n != "[Content_Types].xml")
    return ["[Content_Types].xml"] + rest if "[Content_Types].xml" in names else rest


def sanitize(path: str | Path) -> Path:
    path = Path(path)
    with zipfile.ZipFile(path, "r") as zin:
        entries = {n: zin.read(n) for n in zin.namelist()}

    entries = _strip_printer_settings(entries)
    entries = _strip_notes_master(entries)
    entries = _strip_notes_slides(entries)
    entries = _scrub_content_types(entries)
    entries = _patch_slide_xml(entries)

    pres_p, rels_p = "ppt/presentation.xml", "ppt/_rels/presentation.xml.rels"
    if pres_p in entries and rels_p in entries:
        pres = entries[pres_p].decode("utf-8")
        pres = _patch_sldsz(pres)
        entries[pres_p] = pres.encode("utf-8")

    tmp = path.with_suffix(".sanitized.pptx")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for name in _ordered_names(list(entries.keys())):
            zout.writestr(name, entries[name])
    tmp.replace(path)
    return path


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {sys.argv[0]} <file.pptx>")
    out = sanitize(sys.argv[1])
    print(f"SANITIZED {out}")

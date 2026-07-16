#!/usr/bin/env python3
import re
import shutil
import sys
import zipfile

TARGET = "ppt/presentation.xml"

def fix_xml(xml: str) -> str:
    m = re.search(r"\s*<p:notesMasterIdLst>.*?</p:notesMasterIdLst>", xml, re.S)
    if not m:
        return xml
    block = m.group(0).strip()
    xml = xml[: m.start()] + xml[m.end():]
    if re.search(r"</p:sldMasterIdLst>\s*<p:notesMasterIdLst>", xml):
        return xml
    return xml.replace("</p:sldMasterIdLst>", "</p:sldMasterIdLst>\n  " + block, 1)

def main(path: str) -> int:
    with zipfile.ZipFile(path) as zin:
        names = zin.namelist()
        data = {n: zin.read(n) for n in names}
        infos = {n: zin.getinfo(n) for n in names}
    xml = data[TARGET].decode("utf-8")
    fixed = fix_xml(xml)
    if fixed == xml:
        print("postfix: no change needed")
        return 0
    data[TARGET] = fixed.encode("utf-8")
    tmp = path + ".tmp"
    with zipfile.ZipFile(tmp, "w") as zout:
        for n in names:
            zi = zipfile.ZipInfo(n, date_time=infos[n].date_time)
            zi.compress_type = infos[n].compress_type
            zi.external_attr = infos[n].external_attr
            zout.writestr(zi, data[n])
    shutil.move(tmp, path)
    print("postfix: reordered notesMasterIdLst in " + TARGET)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))

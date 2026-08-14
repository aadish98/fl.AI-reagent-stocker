import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BACKUP = ROOT / ".tmp" / "family1_label_backup"
DEST = ROOT / "data" / "gene_sets" / "Tx-Omics_Revision"

DEFINITIONS = {
    "Mechanistic": (
        "Consistent or correlated sleep/wake genes with published effects on sleep/wake "
        "that are consistent with a potential homeostatic function."
    ),
    "Homeostatic genes": (
        "Genes correlated with sleep/wake history before sampling and rebound sleep/wake "
        "in the opposite direction after sampling (15 positive-history/negative-rebound "
        "and 5 negative-history/positive-rebound genes)."
    ),
    "CSW 4+ genes": (
        "Consistent sleep/wake genes evident across four or more of the seven datasets."
    ),
    "HLH genes": "HLH genes identified as potential upstream regulators.",
}


def transform(name: str, gene_set: str) -> None:
    source = BACKUP / name
    destination = DEST / name
    with source.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        original_fields = list(reader.fieldnames or [])

    fields = list(original_fields)
    if name.startswith("01_"):
        fields[fields.index("Mechanistic Category")] = "mechanistic_subcategory"
    insert_at = fields.index("flybase_gene_id") + 1
    fields[insert_at:insert_at] = ["gene_set", "gene_set_definition"]

    output_rows = []
    for row in rows:
        if name.startswith("01_"):
            category = row.pop("Mechanistic Category")
            row["mechanistic_subcategory"] = category.removeprefix("Homeostatic / ").strip()
        row["gene_set"] = gene_set
        row["gene_set_definition"] = DEFINITIONS[gene_set]
        output_rows.append(row)

    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"{name}: {len(output_rows)} rows, gene_set={gene_set}")


transform("01_Mechanistic_n=6genes.csv", "Mechanistic")
transform("02_Homeostatic_HistoryxRebound_n=20genes.csv", "Homeostatic genes")
transform("03_CSW_4plus_n=97genes.csv", "CSW 4+ genes")
transform("04_HLH_upstream_regulators_n=4genes.csv", "HLH genes")


screen_path = (
    ROOT
    / "data"
    / "gene_sets"
    / "Tx-Omics-FollowUp_v3"
    / "Breakdown"
    / "Tx-Omics_Mechanistic_Screen_Candidates_nsyb_elav.csv"
)
with screen_path.open(newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle)
    screen_rows = list(reader)
    screen_fields = list(reader.fieldnames or [])

old_header = "Mechanistic Category"
new_header = "Mechanistic Subcategory"
if old_header in screen_fields:
    screen_fields[screen_fields.index(old_header)] = new_header
for row in screen_rows:
    value = row.pop(old_header, row.get(new_header, ""))
    row[new_header] = value.removeprefix("Homeostatic / ").strip()

with screen_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=screen_fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(screen_rows)
print(f"{screen_path.name}: renamed category column to {new_header}")

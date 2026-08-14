"""Pinned paths, hashes, and expected invariants for Tx-Omics Revision."""

from __future__ import annotations

from pathlib import Path

from fl_ai_reagent_stocker.config import FLYBASE_DATA, REPO_ROOT

AUDIT_ROOT = REPO_ROOT / "audit_outputs" / "Tx-Omics_Revision"
STOCKER_DIR = REPO_ROOT / "data" / "gene_sets" / "Tx-Omics_Revision"
BREAKDOWN_DIR = (
    REPO_ROOT / "data" / "gene_sets" / "Tx-Omics-FollowUp_v3" / "Breakdown"
)
TX_INPUT_DIR = (
    REPO_ROOT
    / "private"
    / "lab-materials"
    / "Tx-Omics_FollowUp"
    / "Data"
    / "Transcriptomics_Input"
)
GO_ANALYSIS_DIR = REPO_ROOT / "private" / "lab-materials" / "Tx-Omics_FollowUp"
FLYBASE_SYNONYM_PATH = FLYBASE_DATA / "genes" / "fb_synonym_fb_2026_01.tsv.gz"
FETCH_FBGN_SCRIPT = REPO_ROOT / "scripts" / "fetch_fbgn_ids.py"

GO_SCRIPT_BLOBS = {
    "GenerateGOChartReport.py": "8430bf1fb0209ce2ff4c956750672d10981b3e34",
    "ProcessGOresults.py": "e8287d99e6163f111cdfa2dbfaf7832f520c9b80",
    "VisualizeGOResults.py": "56ea2c11817fa31aa589e76712ee27ed34065288",
}
GO_ANALYSIS_REMOTE_COMMIT = "95501af480473d33f906d0c72714ef7613ace722"

DAVID_WSDL = (
    "https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl"
)
DAVID_ENDPOINT = (
    "https://davidbioinformatics.nih.gov/webservice/services/"
    "DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/"
)
DAVID_EMAIL_DEFAULT = "aadishms@umich.edu"
DAVID_SPECIES = "7227"
DAVID_CATEGORIES = (
    "GOTERM_BP_DIRECT,GOTERM_CC_DIRECT,GOTERM_MF_DIRECT,KEGG_PATHWAY"
)
DAVID_EASE = 0.1
DAVID_MIN_COUNT = 1
DAVID_FDR_PERCENT_MAX = 10.0

MASTER_TABLES = {
    "Wake FC0.5": {
        "filename": "Wake gene list FC0.5 Table.csv",
        "threshold": "FC0.5",
        "direction": "wake",
        "expected_rows": 647,
    },
    "Wake FC1": {
        "filename": "Wake gene list FC1 Table.csv",
        "threshold": "FC1",
        "direction": "wake",
        "expected_rows": 98,
    },
    "Sleep FC0.5": {
        "filename": "Sleep gene list FC0.5 Table.csv",
        "threshold": "FC0.5",
        "direction": "sleep",
        "expected_rows": 141,
    },
    "Sleep FC1": {
        "filename": "Sleep gene list FC1 Table.csv",
        "threshold": "FC1",
        "direction": "sleep",
        "expected_rows": 16,
    },
}

EXPECTED_MASTER_ROWS = 902
EXPECTED_UNIQUE_FBGN = 791
MIN_FREQUENCY = 2

EXPERIMENT_COLUMNS = [
    "Baseline",
    "MechSD3",
    "MechSD6",
    "MechSD12",
    "R85C10>TrpA1",
    "R23E10>ChR",
    "THIP",
]

CSW_4PLUS_SOURCES = [
    "FC0.5_Wake_frequency_4_n=82genes.csv",
    "FC0.5_Wake_frequency_5_n=14genes.csv",
    "FC0.5_Wake_frequency_6_n=1genes.csv",
    "FC1_Wake_frequency_4_n=7genes.csv",
]

MECHANISTIC_IDENTITIES = [
    ("unc79", "FBgn0038693"),
    ("SIFa", "FBgn0053527"),
    ("rumpel", "FBgn0029950"),
    ("AstA-R2", "FBgn0039595"),
    ("Trhn", "FBgn0035187"),
    ("RpL23", "FBgn0010078"),
]

MECHANISTIC_DEFINITION = (
    "Consistent or correlated sleep/wake genes with published effects on sleep/wake "
    "that are consistent with a potential homeostatic function."
)
HOMEOSTATIC_DEFINITION = (
    "Genes correlated with sleep/wake history before sampling and rebound sleep/wake "
    "in the opposite direction after sampling (15 positive-history/negative-rebound "
    "and 5 negative-history/positive-rebound genes)."
)
CSW_4PLUS_DEFINITION = (
    "Consistent sleep/wake genes evident across four or more of the seven datasets."
)
HLH_DEFINITION = "HLH genes identified as potential upstream regulators."

HLH_PUBLICATION_GENES = [
    ("bigmax", "FBgn0039509"),
    ("HLH3B", "FBgn0011276"),
    ("E(spl)m3-HLH", "FBgn0002609"),
    ("E(spl)mbeta-HLH", "FBgn0002733"),
]

TRHN_SYMBOL = "Trhn"
TRHN_FBGN = "FBgn0035187"
TRACHEALESS_SYMBOL = "trh"
TRACHEALESS_FBGN = "FBgn0262139"
PUBLICATION_RESOLUTION_BASIS = "RosensweigShah_2026_results_discussion_figure14"

SET_NAMES = [
    "Mechanistic",
    "Homeostatic History/Rebound",
    "CSW 4+",
    "HLH Upstream Regulators",
    "CSW Ribosomal/Translation",
    "CSW Mitochondrial/Metabolism",
    "CSW Immune",
]

PATHWAY_BUCKETS = [
    "Ribosomal/translation",
    "Mitochondrial/metabolism",
    "Immune",
]

BUCKET_TO_SET_NAME = {
    "Ribosomal/translation": "CSW Ribosomal/Translation",
    "Mitochondrial/metabolism": "CSW Mitochondrial/Metabolism",
    "Immune": "CSW Immune",
}

BUCKET_TO_STEM = {
    "Ribosomal/translation": "05_CSW_Ribosomal_Translation",
    "Mitochondrial/metabolism": "06_CSW_Mitochondrial_Metabolism",
    "Immune": "07_CSW_Immune",
}

FBGN_RE = r"^FBgn\d{7}$"

KEYWORD_RULES = {
    "Ribosomal/translation": [
        ("ribosom", r"ribosom"),
        ("translat", r"translat"),
        ("rRNA", r"rrna"),
        ("polysome", r"polysome"),
        ("translational initiation", r"translational initiation"),
        ("translational elongation", r"translational elongation"),
        ("translational termination", r"translational termination"),
        ("eukaryotic initiation", r"eukaryotic initiation"),
    ],
    "Mitochondrial/metabolism": [
        ("mitochondri", r"mitochondri"),
        ("oxidative phosphorylation", r"oxidative phosphorylation"),
        ("respiratory chain", r"respiratory chain"),
        ("NADH dehydrogenase", r"nadh dehydrogenase"),
        ("ATP synthase", r"atp synthase"),
        ("metabolic process", r"metabolic process"),
        ("metabolism", r"metabolism"),
        ("TCA", r"\btca\b"),
        ("citrate cycle", r"citrate cycle"),
    ],
    "Immune": [
        ("immune", r"immune"),
        ("defense response", r"defense response"),
        ("antimicrobial", r"antimicrobial"),
        ("antibacterial", r"antibacterial"),
        ("response to bacterium", r"response to bacterium"),
        ("Toll", r"\btoll\b"),
        ("Imd", r"\bimd\b"),
        ("melanization", r"melanization"),
        ("peptidoglycan recognition", r"peptidoglycan recognition"),
        ("phagocytosis", r"phagocytosis"),
    ],
}


def audit_paths(root: Path | None = None) -> dict[str, Path]:
    base = Path(root) if root is not None else AUDIT_ROOT
    return {
        "root": base,
        "staging_categories": base / "Staging" / "Categories",
        "go_input": base / "GO" / "Input",
        "go_raw": base / "GO" / "Raw",
        "go_processed": base / "GO" / "Processed",
        "go_review": base / "GO" / "Review",
        "pathways_preaudit": base / "Pathways" / "PreAudit",
        "pathways_auditrun": base / "Pathways" / "AuditRun",
        "pathways_approved": base / "Pathways" / "Approved",
        "evidence": base / "Evidence",
        "audits": base / "Audits",
        "figures": base / "Figures",
        "figure_data": base / "FigureData",
        "overlap_report": base / "Overlap_Report.md",
        "readme": base / "README.md",
        "manifest": base / "run_manifest.json",
    }

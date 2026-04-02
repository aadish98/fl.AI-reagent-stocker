"""Generate a pipeline data-flow diagram as PNG using Graphviz."""

import graphviz

dot = graphviz.Digraph(
    "fl_ai_reagent_stocker_pipeline",
    format="png",
    engine="dot",
    graph_attr={
        "rankdir": "TB",
        "bgcolor": "#FFFFFF",
        "fontname": "Helvetica Neue",
        "pad": "0.4",
        "nodesep": "0.5",
        "ranksep": "0.7",
        "dpi": "200",
        "label": "<<B><FONT POINT-SIZE=\"26\">fl.AI Reagent Stocker — Data Flow</FONT></B>>",
        "labelloc": "t",
        "labeljust": "c",
        "fontcolor": "#1a1a2e",
        "splines": "polyline",
        "compound": "true",
    },
    node_attr={
        "fontname": "Helvetica Neue",
        "fontsize": "10",
        "style": "filled",
        "penwidth": "1.4",
        "margin": "0.15,0.08",
    },
    edge_attr={
        "fontname": "Helvetica Neue",
        "fontsize": "8",
        "color": "#666666",
        "arrowsize": "0.7",
        "penwidth": "1.0",
    },
)

INPUT_FILL = "#e3f2fd"
INPUT_BORDER = "#1976d2"
STAGE_FILL = "#fff8e1"
STAGE_BORDER = "#f9a825"
PROC_FILL = "#f5f5f5"
PROC_BORDER = "#757575"
OUT_FILL = "#e8f5e9"
OUT_BORDER = "#388e3c"
API_FILL = "#fce4ec"
API_BORDER = "#c62828"
CFG_FILL = "#f3e5f5"
CFG_BORDER = "#7b1fa2"

# ═══════════════════════════════════════════════════════════════
#  INPUT DATA
# ═══════════════════════════════════════════════════════════════
with dot.subgraph(name="cluster_inputs") as c:
    c.attr(
        label="<<B>Input Data</B>>",
        style="rounded,dashed",
        color="#1976d2",
        fontcolor="#1976d2",
        fontsize="13",
        margin="16",
    )
    c.node("gene_csv", "Gene List CSVs\n(*.csv with gene symbols\nor flybase_gene_id)",
           shape="folder", fillcolor=INPUT_FILL, color=INPUT_BORDER)
    c.node("fb_alleles", "FlyBase Alleles & Stocks\nfbal_to_fbgn · stocks ·\nallele descriptions",
           shape="folder", fillcolor=INPUT_FILL, color=INPUT_BORDER)
    c.node("fb_refs", "FlyBase References\nentity_publication ·\nfbrf_pmid_pmcid_doi",
           shape="folder", fillcolor=INPUT_FILL, color=INPUT_BORDER)
    c.node("fb_ti", "Transgenic Data\nconstruct descriptions ·\nfbtp_to_fbti",
           shape="folder", fillcolor=INPUT_FILL, color=INPUT_BORDER)
    c.node("fb_stcomp", "Derived Stock Map\nfbst_to_derived_\nstock_component.csv",
           shape="folder", fillcolor=INPUT_FILL, color=INPUT_BORDER)
    c.node("json_cfg", "Split Config JSON\nfilters · combinations ·\nrelevantSearchTerms",
           shape="hexagon", fillcolor=CFG_FILL, color=CFG_BORDER)

# ═══════════════════════════════════════════════════════════════
#  STAGE 1 — find-stocks
# ═══════════════════════════════════════════════════════════════
with dot.subgraph(name="cluster_s1") as c:
    c.attr(
        label="<<B>Stage 1 — find-stocks</B>>",
        style="filled,rounded",
        fillcolor="#fffde7",
        color="#f9a825",
        fontcolor="#e65100",
        fontsize="14",
        margin="16",
    )
    c.node("s1_fbgn", "Gene Symbol → FBgn\nconversion (optional)",
           shape="box", fillcolor=PROC_FILL, color=PROC_BORDER, style="filled,rounded")
    c.node("s1_map", "Build Stock Mapping\nFBgn → FBal → FBtp/FBti → FBst",
           shape="box", fillcolor=PROC_FILL, color=PROC_BORDER, style="filled,rounded")
    c.node("s1_ref", "Link References\nmap FBrf → PMID per\nstock & component",
           shape="box", fillcolor=PROC_FILL, color=PROC_BORDER, style="filled,rounded")
    c.node("s1_kw", "Keyword Scoring\ncount title/abstract\nkeyword hits",
           shape="box", fillcolor=PROC_FILL, color=PROC_BORDER, style="filled,rounded")
    c.node("pubmed", "NCBI PubMed\nEntrez API",
           shape="component", fillcolor=API_FILL, color=API_BORDER)

dot.edge("gene_csv", "s1_fbgn")
dot.edge("s1_fbgn", "s1_map")
dot.edge("fb_alleles", "s1_map")
dot.edge("fb_stcomp", "s1_map")
dot.edge("fb_ti", "s1_map")
dot.edge("s1_map", "s1_ref")
dot.edge("fb_refs", "s1_ref")
dot.edge("s1_ref", "pubmed", label="fetch metadata\nfor PMIDs")
dot.edge("pubmed", "s1_kw", label="title · abstract")
dot.edge("s1_ref", "s1_kw")
dot.edge("json_cfg", "s1_kw", label="keywords", style="dashed", color="#7b1fa2")

# Stage 1 outputs
dot.node("s1_xlsx", "Stocks/\naggregated_stock_refs.xlsx\n(Stocks + References sheets)",
         shape="note", fillcolor=OUT_FILL, color=OUT_BORDER)
dot.node("s1_nopmid", "references_without_\npmid_fbrf.txt",
         shape="note", fillcolor=OUT_FILL, color=OUT_BORDER)
dot.edge("s1_kw", "s1_xlsx")
dot.edge("s1_kw", "s1_nopmid")

# ═══════════════════════════════════════════════════════════════
#  STAGE 2 — split-stocks
# ═══════════════════════════════════════════════════════════════
with dot.subgraph(name="cluster_s2") as c:
    c.attr(
        label="<<B>Stage 2 — split-stocks</B>>",
        style="filled,rounded",
        fillcolor="#e8f5e9",
        color="#388e3c",
        fontcolor="#1b5e20",
        fontsize="14",
        margin="16",
    )
    c.node("s2_derive", "Compute Derived Columns\nBalancers · multiple_insertions ·\nALLELE_PAPER_RELEVANCE_SCORE",
           shape="box", fillcolor=PROC_FILL, color=PROC_BORDER, style="filled,rounded")
    c.node("s2_filter", "Apply JSON Filters\n& Combinations\n(per-sheet partitions)",
           shape="box", fillcolor=PROC_FILL, color=PROC_BORDER, style="filled,rounded")
    c.node("s2_limit", "Apply Stock Limits\nmaxStocksPerGene ·\nmaxStocksPerAllele",
           shape="box", fillcolor=PROC_FILL, color=PROC_BORDER, style="filled,rounded")

dot.edge("s1_xlsx", "s2_derive", style="bold", color="#388e3c", label="reads Stage 1\noutput")
dot.edge("json_cfg", "s2_filter", label="filters ·\ncombinations", style="dashed", color="#7b1fa2")
dot.edge("gene_csv", "s2_limit", label="input genes\n(for limits)", style="dashed", color="#1976d2")
dot.edge("s2_derive", "s2_filter")
dot.edge("s2_filter", "s2_limit")

dot.node("s2_xlsx", "Organized Stocks/\n*_aggregated.xlsx\n(Contents · per-combo sheets ·\nStock Sheet by Gene)",
         shape="note", fillcolor=OUT_FILL, color=OUT_BORDER)
dot.edge("s2_limit", "s2_xlsx")

# ═══════════════════════════════════════════════════════════════
#  STAGE 3 — validate-stocks
# ═══════════════════════════════════════════════════════════════
with dot.subgraph(name="cluster_s3") as c:
    c.attr(
        label="<<B>Stage 3 — validate-stocks</B>>",
        style="filled,rounded",
        fillcolor="#fce4ec",
        color="#c62828",
        fontcolor="#b71c1c",
        fontsize="14",
        margin="16",
    )
    c.node("s3_refpp", "Identify Ref++ Stocks\nfrom output sheets",
           shape="box", fillcolor=PROC_FILL, color=PROC_BORDER, style="filled,rounded")
    c.node("s3_fetch", "Fetch Full Text\nfor Ref++ PMIDs",
           shape="box", fillcolor=PROC_FILL, color=PROC_BORDER, style="filled,rounded")
    c.node("s3_gpt", "GPT Functional\nValidation\n(per allele-PMID pair)",
           shape="box", fillcolor=PROC_FILL, color=PROC_BORDER, style="filled,rounded")
    c.node("s3_merge", "Merge Validation\nColumns into Workbook",
           shape="box", fillcolor=PROC_FILL, color=PROC_BORDER, style="filled,rounded")
    c.node("fulltext_api", "Full-Text Sources\nPMC · Europe PMC ·\nUnpaywall",
           shape="component", fillcolor=API_FILL, color=API_BORDER)
    c.node("openai_api", "OpenAI API\n(GPT model)",
           shape="component", fillcolor=API_FILL, color=API_BORDER)

dot.edge("s1_xlsx", "s3_refpp", style="bold", color="#c62828",
         label="re-reads Stage 1\noutput (re-applies\nfilters + limits)")
dot.edge("json_cfg", "s3_refpp", label="Ref++ filter\nrules", style="dashed", color="#7b1fa2")
dot.edge("s3_refpp", "s3_fetch")
dot.edge("s3_fetch", "fulltext_api")
dot.edge("fulltext_api", "s3_gpt", label="paper text")
dot.edge("s3_gpt", "openai_api", label="prompt +\npaper text")
dot.edge("openai_api", "s3_gpt", label="validation\nresult", style="dashed")
dot.edge("s3_gpt", "s3_merge")

dot.node("s3_xlsx",
         "Organized Stocks/\n*_aggregated.xlsx\n(+ Functionally Valid Stock? ·\nvalidation_rationale ·\nphenotypes_observed)",
         shape="note", fillcolor=OUT_FILL, color=OUT_BORDER)
dot.edge("s3_merge", "s3_xlsx")

# ═══════════════════════════════════════════════════════════════
#  run-full-pipeline wrapper
# ═══════════════════════════════════════════════════════════════
dot.node("full_pipe", "run-full-pipeline\n(orchestrates Stages 1 → 2 → 3)",
         shape="doubleoctagon", fillcolor="#eceff1", color="#455a64",
         style="filled,bold", fontsize="10", penwidth="2")
dot.edge("full_pipe", "s1_fbgn", style="dashed", color="#455a64", label="①", penwidth="1.5")
dot.edge("full_pipe", "s2_derive", style="dashed", color="#455a64", label="②", penwidth="1.5")
dot.edge("full_pipe", "s3_refpp", style="dashed", color="#455a64", label="③", penwidth="1.5")

# ═══════════════════════════════════════════════════════════════
#  Legend
# ═══════════════════════════════════════════════════════════════
with dot.subgraph(name="cluster_legend") as c:
    c.attr(label="<<B>Legend</B>>", style="rounded", color="#bdbdbd",
           fontsize="11", fontcolor="#616161", margin="12")
    c.node("l1", "Input Data", shape="folder", fillcolor=INPUT_FILL,
           color=INPUT_BORDER, fontsize="9", width="1.2")
    c.node("l2", "Config", shape="hexagon", fillcolor=CFG_FILL,
           color=CFG_BORDER, fontsize="9", width="1.2")
    c.node("l3", "Processing Step", shape="box", fillcolor=PROC_FILL,
           color=PROC_BORDER, fontsize="9", style="filled,rounded", width="1.2")
    c.node("l4", "External API", shape="component", fillcolor=API_FILL,
           color=API_BORDER, fontsize="9", width="1.2")
    c.node("l5", "Output File", shape="note", fillcolor=OUT_FILL,
           color=OUT_BORDER, fontsize="9", width="1.2")
    c.edge("l1", "l2", style="invis")
    c.edge("l2", "l3", style="invis")
    c.edge("l3", "l4", style="invis")
    c.edge("l4", "l5", style="invis")

output_path = dot.render(
    filename="pipeline-data-flow",
    directory=".",
    cleanup=True,
)
print(f"Flowchart saved to: {output_path}")

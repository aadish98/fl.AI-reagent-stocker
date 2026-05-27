# Pipeline Flowcharts

These flowcharts show the project at two levels:

- the high-level command flow from input gene lists to final workbooks
- the low-level evidence flow that explains where stock references, phenotype references, keyword buckets, and validation decisions come from

## Data flowchart

```mermaid
flowchart TD
    subgraph DataPreprocessing["Data Preprocessing"]
        direction TB
        GeneSet(["Gene set"])
        FlyBaseTSV[("FlyBase TSV")]
        ChadoFBst[("FlyBase XML<br/>chado_FBst")]
        ChadoFBti[("FlyBase XML<br/>chado_FBti")]

        FBstLookup["Build FBst -> component lookup<br/>fbst_to_derived_stock_component.csv"]
        FBtiLookup["Build FBtp -> FBti lookup<br/>fbtp_to_fbti.csv"]

        Components["Find all associated<br/>components<br/>(FBal / FBtp / FBti)"]
        BranchInputs["Shared downstream inputs<br/>components + FBst lookup"]
    end

    subgraph DownstreamBranches["Downstream branches"]
        direction LR

        subgraph PhenotypicStocks["Phenotypic Stocks"]
            direction TB
            PhenoFBstRoute["FBst -> component lookup<br/>(full lookup)"]
            Phenotypes["Find all associated<br/>phenotypes for components"]
            KeywordSplit{"Split by keywords<br/>(true / false)"}
            EmbeddingSearch["Contextual search using<br/>OpenAI embedding model<br/>on phenotypes"]
            AttachPhenotypeStocks["Attach FBst stocks<br/>(consumes FBst -> component lookup)"]
            Phenotypic[["ALL PHENOTYPIC REAGENTS<br/>+ associated refs"]]
        end

        subgraph AvailableStocks["Available Stocks"]
            direction TB
            StockFBstRoute["component -> stock lookup<br/>(filtered lookup)"]
            Stocks["Find all associated<br/>stocks"]
            Categorize["Categorize by ★Ref++ <br/>(stock + reagent type)"]
            Validate["Validate Ref++ stocks<br/>using paper text + OpenAI"]
            Available[["ALL AVAILABLE REAGENTS<br/>+ associated refs"]]
        end
    end

    GeneSet --> Components
    FlyBaseTSV -. source .-> Components
    ChadoFBti -. source .-> FBtiLookup
    FBtiLookup -. expands FBtp components .-> Components
    ChadoFBst -. source .-> FBstLookup

    Components --> BranchInputs
    FBstLookup -.-> BranchInputs

    PhenoFBstRoute -.-> Phenotypes

    StockFBstRoute -.-> Stocks

    Stocks --> Categorize
    Categorize --> Validate
    Validate --> Available

    Phenotypes --> KeywordSplit
    Phenotypes --> EmbeddingSearch
    KeywordSplit --> AttachPhenotypeStocks
    EmbeddingSearch --> AttachPhenotypeStocks
    AttachPhenotypeStocks --> Phenotypic

    linkStyle 0 stroke:#dc2626,stroke-width:4px
    linkStyle 1 stroke:#2563eb,stroke-width:2px
    linkStyle 2 stroke:#db2777,stroke-width:2px
    linkStyle 3 stroke:#ea580c,stroke-width:2.5px
    linkStyle 4 stroke:#db2777,stroke-width:2px
    linkStyle 5 stroke:#334155,stroke-width:2.5px
    linkStyle 6 stroke:#ea580c,stroke-width:2.5px
    linkStyle 7 stroke:#f97316,stroke-width:2.5px
    linkStyle 8 stroke:#f97316,stroke-width:2.5px
    linkStyle 9 stroke:#0891b2,stroke-width:2.5px
    linkStyle 10 stroke:#0891b2,stroke-width:2.5px
    linkStyle 11 stroke:#059669,stroke-width:3px
    linkStyle 12 stroke:#0891b2,stroke-width:2.5px
    linkStyle 13 stroke:#4b5563,stroke-width:2.5px
    linkStyle 14 stroke:#0891b2,stroke-width:2.5px
    linkStyle 15 stroke:#4b5563,stroke-width:2.5px
    linkStyle 16 stroke:#059669,stroke-width:3px

    classDef gene fill:#0f172a,stroke:#334155,stroke-width:2px,color:#ffffff
    classDef tsv fill:#dbeafe,stroke:#2563eb,stroke-width:2px,color:#1e3a8a
    classDef xml fill:#fce7f3,stroke:#db2777,stroke-width:2px,color:#831843
    classDef lookup fill:#ffedd5,stroke:#ea580c,stroke-width:2px,color:#7c2d12
    classDef router fill:#fef9c3,stroke:#ca8a04,stroke-width:2px,color:#713f12
    classDef route fill:#fff7ed,stroke:#f97316,stroke-width:2px,color:#7c2d12
    classDef step fill:#ecfeff,stroke:#0891b2,stroke-width:1.5px,color:#164e63
    classDef ai fill:#f3f4f6,stroke:#4b5563,stroke-width:2px,color:#111827
    classDef output fill:#ecfdf5,stroke:#059669,stroke-width:3px,color:#064e3b

    class GeneSet gene
    class FlyBaseTSV tsv
    class ChadoFBst,ChadoFBti xml
    class FBstLookup,FBtiLookup lookup
    class BranchInputs router
    class PhenoFBstRoute,StockFBstRoute route
    class Components,Stocks,Categorize,Phenotypes,KeywordSplit,AttachPhenotypeStocks step
    class Validate,EmbeddingSearch ai
    class Available,Phenotypic output
```

- ★ `Ref++` = a publication that has a configured keyword in its title or abstract.
- Data Preprocessing turns the input gene set plus FlyBase TSV/XML-derived lookups into the shared component and stock-lookup context used downstream.
- ++Available Stocks++ is the all-stocks path: find associated stocks, categorize by `Ref++` and reagent type, then validate `Ref++` candidates with paper text + OpenAI.
- ++Phenotypic Stocks++ is the phenotype path: find phenotype-associated reagents, split/search those phenotypes, then attach FBst stocks using the full stock-component lookup.

---

## High level technical flowchart

The pipeline now has four explicit stages plus an end-to-end wrapper.
`split-stocks`, `phenotype-sheet`, and `validate-stocks` are sibling
downstream commands over Stage 1 output. `run-full-pipeline --mode normal`
runs Stage 1 + `split-stocks` + `validate-stocks`; `run-full-pipeline --mode
phenotype` runs Stage 1 + `phenotype-sheet` only and never invokes GPT
validation.

```mermaid
flowchart TD
    A[Input gene-list CSVs] --> S1[Stage 1<br/>find-stocks]
    S1 --> O1[(aggregated_stock_refs.xlsx<br/>Stocks + References + Gene Reagent Index)]

    O1 --> S2[split-stocks]
    S2 --> O2[(Organized Stocks workbook)]
    O2 --> S3[validate-stocks]
    S3 --> O3[(Validated workbook<br/>+ validation columns)]

    O1 --> PS[phenotype-sheet]
    PS --> P1[(Stock Phenotype Sheet workbook)]
    PS -. --similarity tiers .-> P2[(Cosine similarity tier workbook)]
    PS -. --similarity simple-buckets .-> P3[(Simple-bucket workbook)]
    PS -. --similarity keyword-buckets .-> P4[(Keyword Hits / No Keyword Hits workbook)]

    classDef stage fill:#e3f2fd,stroke:#1565c0,color:#0d47a1
    classDef artifact fill:#fff8e1,stroke:#f57f17,color:#e65100
    classDef pheno fill:#f3e5f5,stroke:#6a1b9a,color:#4a148c
    class S1,S2,S3 stage
    class PS pheno
    class O1,O2,O3,P1,P2,P3,P4 artifact
```

---

## Low level technical flowchart

Each stage is its own subgraph. Within a stage, every step flows top to
bottom. Cross-stage arrows only appear at clear hand-off points
(`Stage 1 → split-stocks`, `Stage 1 → phenotype-sheet`, `split-stocks → validate-stocks`).
FlyBase table names are written into the steps that use them so there are no
shared fan-out edges from a single "Inputs" node.

```mermaid
flowchart TD
    A[Input gene-list CSVs]

    subgraph Stage1["Stage 1: find-stocks"]
        direction TB
        A1[Map gene symbols → FBgn]
        A1 --> A2[gene → allele<br/>fbal_to_fbgn]
        A2 --> A3[allele → construct<br/>transgenic_construct_descriptions]
        A3 --> A4[construct → insertion<br/>fbtp_to_fbti]
        A4 --> A5[component → stock<br/>fbst_to_derived_stock_component]
        A5 --> A6[Stock rows with<br/>relevant FBal / FBtp / FBti]
        A6 --> A7[component → FBrf<br/>entity_publication]
        A7 --> A8[FBrf → PMID / PMCID / DOI<br/>fbrf_pmid_pmcid_doi]
        A8 --> A9[Score keywords against<br/>title / abstract metadata]
        A9 --> O1[(aggregated_stock_refs.xlsx)]
    end

    subgraph Stage2["split-stocks (organization only)"]
        direction TB
        B1[Load Stage 1 workbook]
        B1 --> B2[Compute derived columns<br/>Balancers, RNAi,<br/>ALLELE_PAPER_RELEVANCE_SCORE]
        B2 --> B3[Apply JSON filters and combinations]
        B3 --> B4[Apply maxStocksPerGene<br/>and maxStocksPerAllele]
        B4 --> O2[(Organized workbook sheets)]
    end

    subgraph PhenotypeSheet["phenotype-sheet (gene-first phenotype reagent flow)"]
        direction TB
        C1[Input-gene reagents from Gene Reagent Index]
        C1 --> C2[Match curated phenotype rows<br/>genotype_phenotype_data]
        CL[Global FBst lookup<br/>fbst_to_derived_stock_component.csv]
        CL --> C4
        C2 --> C3[Resolve phenotype FBrf → PMID / PMCID<br/>fbrf_pmid_pmcid_doi]
        C3 --> C4[(Stock Phenotype Sheet)]
        C4 --> C5[OpenAI embeddings vs<br/>phenotypeSimilarityTargets]
        C5 --> C6{"--similarity"}
        C6 -->|"none"| C10[(no sidecar)]
        C6 -->|"tiers"| C7[(Cosine similarity tiers)]
        C6 -->|"simple-buckets"| C8[(Simple buckets)]
        C6 -->|"keyword-buckets"| C9[(Keyword Hits /<br/>No Keyword Hits)]
    end

    subgraph Stage3["validate-stocks (GPT validation only)"]
        direction TB
        D1[Rebuild split outputs]
        D1 --> D2[Select Ref++ stocks<br/>in output sheets]
        D2 --> D3[Pair with keyword-hit PMIDs<br/>from Stage 1]
        D3 --> D4[Retrieve full text<br/>PubMed + PMC + Unpaywall]
        D4 --> D5[GPT functional validation]
        D5 --> D6[Validation status, phenotypes,<br/>confidence, rationale]
        D6 --> O3[(Workbook + validation columns)]
    end

    A --> A1
    O1 --> B1
    O1 --> C1
    O2 --> D1

    classDef artifact fill:#fff8e1,stroke:#f57f17,color:#e65100
    classDef external fill:#ede7f6,stroke:#5e35b1,color:#311b92
    class O1,O2,O3,C4,C7,C8,C9,C10 artifact
    class C5,D4,D5 external
```

---

## Important Distinction

The pipeline has two separate reference concepts:

- **Stock-level references** come from `entity_publication` and answer: "Which papers are associated with this stock's relevant allele, construct, or insertion components?"
- **Phenotype-sheet references** come from `genotype_phenotype_data.reference` and answer: "Which paper supports this specific curated phenotype row?"

A PMID can appear in Stage 1 but not in the Stock Phenotype Sheet when FlyBase associates the paper with the stock component but does not use that paper as the reference for a curated phenotype row.

---

## Stock Phenotype Sheet is Phenotype-Table-Driven

The Stock Phenotype Sheet is driven by every FBal / FBtp / FBti reagent associated with an input gene (the complete Stage 1 **Gene Reagent Index**), filtered against curated rows in `genotype_phenotype_data`. Stocks are then attached afterward via the full, unfiltered `fbst_to_derived_stock_component.csv`. This means:

- A reagent that has a phenotype row is included regardless of whether any FBst stock was found by Stage 1 stock matching or kept by Stage 2 JSON split limits.
- Every FBst FlyBase associates with a matched reagent surfaces as a Source/ Stock # row, not just the stocks that survived Stage 1 ranking.
- A phenotype-matched reagent with no FBst still produces a "no-stock phenotype reagent" row sourced from the reagent index.

# Citing Resources Used by `fl_ai_reagent_stocker`

If you publish work that depends on this pipeline, please cite the underlying
data sources, APIs, and models the tool relies on. The pipeline does not
generate the underlying biological knowledge: it organizes, links, ranks, and
text-mines records from the resources listed below.

This guide groups citations by category and provides:

- a one-line description of how the pipeline uses the resource
- the canonical reference to cite in a Methods section
- a ready-to-paste BibTeX entry where one is available
- the recommended URL for "Data availability" sections

A suggested Methods-section paragraph appears at the end.

> Always verify each resource's current "How to cite" page at the time of
> submission. Database papers, model versions, and access dates change.

---

## 1. FlyBase (primary genetic database)

**Used for:** all gene, allele (FBal), construct (FBtp), insertion (FBti),
stock (FBst), and reference (FBrf) identifiers; gene-symbol → FBgn conversion;
allele-to-construct, construct-to-insertion, and component-to-stock lookups;
genotype-phenotype curation; stock-to-publication links; and FBrf → PMID /
PMCID / DOI resolution. Specific report families consumed by the pipeline:

- `fbal_to_fbgn`
- `transgenic_construct_descriptions`
- `entity_publication`
- `fbrf_pmid_pmcid_doi`
- `genotype_phenotype_data`
- `chado_FBst*.xml(.gz)` and `chado_FBti*.xml(.gz)` (chado XML)

Files are pulled from
`https://s3ftp.flybase.org/releases/current/precomputed_files/` and
`https://s3ftp.flybase.org/releases/current/chado-xml/` via
`scripts/refresh_flybase_data.py`. **Record the FlyBase release version used
(e.g. `FB2025_03`)** in your Methods section; it is part of the on-disk
filenames.

**Cite:**

> Öztürk-Çolak A., Marygold S. J., Antonazzo G., Attrill H., Goutte-Gattat D.,
> Jenkins V. K., Matthews B. B., Millburn G., dos Santos G., Tabone C. J.,
> the FlyBase Consortium. FlyBase: updates to the *Drosophila* genes and
> genomes database. *Genetics* (current release citation).

```bibtex
@article{flybase2024,
  author  = {{\"O}zt{\"u}rk-{\c{C}}olak, Aya{\c{s}}e and Marygold, Steven J. and
             Antonazzo, Giulia and Attrill, Helen and Goutte-Gattat, Damien and
             Jenkins, Victoria K. and Matthews, Beverley B. and Millburn, Gillian
             and dos Santos, Gilberto and Tabone, Christopher J. and the FlyBase Consortium},
  title   = {{FlyBase: updates to the Drosophila genes and genomes database}},
  journal = {Genetics},
  year    = {2024},
  doi     = {10.1093/genetics/iyad211}
}
```

URL: <https://flybase.org/>

---

## 2. Stock centers behind the FBst records

The pipeline does not query stock-center catalogs directly, but every FBst
returned by the pipeline ultimately points back to a physical stock held by a
stock center. **You must cite the stock center for any stock you actually
order or use** in experiments. The most common ones are:

### Bloomington Drosophila Stock Center (BDSC)

> Cook K. R., Parks A. L., Jacobus L. M., Kaufman T. C., Matthews K. A. New
> research resources at the Bloomington *Drosophila* Stock Center.
> *Fly (Austin)* 4, 88–91 (2010).

NIH support: P40OD018537. Acknowledge the grant when ordering BDSC stocks.

URL: <https://bdsc.indiana.edu/>

### Vienna Drosophila Resource Center (VDRC)

> Dietzl G. *et al.* A genome-wide transgenic RNAi library for conditional
> gene inactivation in *Drosophila*. *Nature* 448, 151–156 (2007).

URL: <https://stockcenter.vdrc.at/>

### Transgenic RNAi Project (TRiP, Harvard Medical School)

> Perkins L. A. *et al.* The Transgenic RNAi Project at Harvard Medical
> School: resources and validation. *Genetics* 201, 843–852 (2015).

URL: <https://fgr.hms.harvard.edu/>

### Kyoto Stock Center (DGRC, Kyoto Institute of Technology)

If you use Kyoto-derived stocks, cite the Kyoto DGRC and credit the
respective depositor.

URL: <https://kyotofly.kit.jp/cgi-bin/stocks/index.cgi>

### NIG-Fly (National Institute of Genetics, Japan)

URL: <https://shigen.nig.ac.jp/fly/nigfly/>

---

## 3. Literature retrieval and metadata APIs

These services are queried by `fl_ai_reagent_stocker.integrations.*` to attach
PMIDs, abstracts, and full text to FlyBase references during keyword scoring
and validation.

### NCBI PubMed and PMC (E-utilities + ID Converter + OAI)

**Used for:** PMID/PMCID metadata, abstracts, OAI-PMH full-text records, and
the PMC ID converter at
`https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/`. Set `NCBI_API_KEY`
to raise the request rate limit.

> Sayers E. W. *et al.* Database resources of the National Center for
> Biotechnology Information. *Nucleic Acids Research* (current annual
> database issue).

URLs:

- PubMed: <https://pubmed.ncbi.nlm.nih.gov/>
- PMC: <https://pmc.ncbi.nlm.nih.gov/>
- E-utilities terms: <https://www.ncbi.nlm.nih.gov/books/NBK25497/>

### Europe PMC (search + full-text XML)

**Used for:** keyword-scoped search and `fullTextXML` retrieval at
`https://www.ebi.ac.uk/europepmc/webservices/rest/`.

> Europe PMC Consortium. Europe PMC: a full-text literature database for the
> life sciences and platform for innovation. *Nucleic Acids Research* 43,
> D1042–D1048 (2015).

URL: <https://europepmc.org/>

### Unpaywall

**Used for:** locating open-access PDFs by DOI when full text is paywalled.
Set `UNPAYWALL_TOKEN` (your contact email; default `aadishms@umich.edu`) per
Unpaywall's terms of service.

> Piwowar H. *et al.* The state of OA: a large-scale analysis of the
> prevalence and impact of Open Access articles. *PeerJ* 6, e4375 (2018).

URL: <https://unpaywall.org/>

### Crossref

**Used for:** DOI metadata fall-back at `https://api.crossref.org/works`.

> Hendricks G., Tkaczyk D., Lin J., Feeney P. Crossref: the sustainable
> source of community-owned scholarly metadata. *Quantitative Science Studies*
> 1, 414–427 (2020).

URL: <https://www.crossref.org/>

### OpenAlex

**Used for:** works-metadata enrichment at
`https://api.openalex.org/works`.

> Priem J., Piwowar H., Orr R. OpenAlex: a fully-open index of scholarly
> works, authors, venues, institutions, and concepts. *arXiv* 2205.01833
> (2022).

URL: <https://openalex.org/>

---

## 4. AI models (OpenAI)

**Used for:** (a) phenotype-text embeddings when the embedding analysis is
enabled with `settings.embeddings.enabled`, which produces cosine columns,
similarity tiers, and plots, and (b) functional
validation of `Ref++` stocks during `run` (or `validate-stocks`) when the
required JSON config block `settings.validation.enabled` is `true`. Model
defaults are configured in `fl_ai_reagent_stocker/config.py`:

- chat / validation model: `gpt-5-mini` (override with `OPENAI_MODEL`)
- embedding model: `text-embedding-3-large` (override with
  `OPENAI_EMBEDDING_MODEL`)

**Cite the specific model version you used at run time** (record it from
your run logs; the default may change with future releases).

> OpenAI. *GPT-5 / text-embedding-3-large*. Accessed via the OpenAI API at
> `https://api.openai.com/`, accessed YYYY-MM-DD.

```bibtex
@misc{openai_gpt5_2025,
  author       = {{OpenAI}},
  title        = {{GPT-5: model and API documentation}},
  year         = {2025},
  howpublished = {\url{https://platform.openai.com/docs/models}},
  note         = {Accessed via the OpenAI API}
}
```

URL: <https://platform.openai.com/docs/models>

---

## 5. The pipeline itself

If `fl_ai_reagent_stocker` materially shaped your analysis, please cite this
software repository in addition to the data sources above.

```bibtex
@software{fl_ai_reagent_stocker,
  author = {Allada Lab},
  title  = {{fl\_ai\_reagent\_stocker: modular Drosophila stock-processing pipeline}},
  url    = {https://github.com/Allada-Lab/fl-AI-reagent-stocker},
  note   = {Record the commit hash or release tag used in your analysis.}
}
```

When you cite the software, also report:

- the commit hash or release tag of `fl_ai_reagent_stocker`
- the FlyBase release version (e.g. `FB2025_03`)
- the OpenAI model versions resolved at run time
- the JSON config file used (for example
  `data/config/stock_split_config_example.json`)

These four items make the run reproducible.

---

## 6. Suggested Methods paragraph

A drop-in template you can adapt:

> Candidate *Drosophila* stocks were assembled with `fl_ai_reagent_stocker`
> (commit `<HASH>`), which links input gene symbols to FlyBase
> (release `<FBYYYY_NN>`; Öztürk-Çolak et al., 2024) records via FBgn,
> FBal, FBtp, and FBti identifiers and joins each component to FBst stocks
> through the FlyBase chado XML and precomputed report families. Stock
> publications were resolved through `entity_publication` and
> `fbrf_pmid_pmcid_doi`, then cross-referenced against PubMed (Sayers et al.,
> NAR), PMC, Europe PMC (Europe PMC Consortium, 2015), Crossref (Hendricks
> et al., 2020), and OpenAlex (Priem et al., 2022); open-access full text was
> located through Unpaywall (Piwowar et al., 2018). Where indicated,
> phenotype text was embedded with OpenAI `text-embedding-3-large` and
> functional validation of `Ref++` stocks was performed with OpenAI
> `<MODEL_VERSION>` against retrieved full text. Stocks were ordered from
> the Bloomington *Drosophila* Stock Center (NIH P40OD018537; Cook et al.,
> 2010), the Vienna *Drosophila* Resource Center (Dietzl et al., 2007),
> and/or the Transgenic RNAi Project at Harvard Medical School (Perkins et
> al., 2015) as appropriate.

---

## 7. Quick checklist before submission

- [ ] FlyBase paper cited and FlyBase release version reported.
- [ ] Stock-center paper(s) cited for every center an ordered stock came from.
- [ ] BDSC NIH grant `P40OD018537` acknowledged if BDSC stocks were used.
- [ ] PubMed / PMC, Europe PMC, Crossref, OpenAlex, and Unpaywall cited if
      their data appear in your reference table or full-text validation.
- [ ] OpenAI model versions and access date recorded if validation or
      embedding similarity were used.
- [ ] `fl_ai_reagent_stocker` commit hash and JSON config reported for
      reproducibility.

# Data handling: what leaves your machine

ClawBio is local-first. This page says exactly what that means, skill by skill,
so that an institution can decide which skills it may run on data it is
responsible for. It was written against `main` on 2026-09-04 by reading the
code, not the descriptions. `tests/test_data_handling_doc.py` scans every skill
for outbound-call code and fails if a networked skill is missing from this page,
so the list cannot silently fall behind the code.

## The short version

- The `clawbio` package (CLI, runner, shared helpers) makes no network calls.
  There is no telemetry, no update check, no analytics. The audit layer in
  `clawbio/common/audit.py` writes a local JSONL file and nothing else.
- Most skills never touch the network. They read your files, compute, and write
  to the output directory you name.
- The skills that do reach the network are all listed below, in five classes,
  from least to most sensitive. Nothing is sent unless you run the skill, and a
  `--demo` run is offline unless the table says otherwise.
- Credentials are read from the environment variables named below. ClawBio
  never stores them, and reproducibility bundles do not contain them.

## The classes

| Class | What is sent | Typical example |
|-------|--------------|-----------------|
| 1 | Nothing of yours. The skill downloads public reference data, summary statistics, package metadata or model weights. | `gwas-prs` fetching a PGS Catalog scoring file |
| 2 | Query terms you typed: gene symbols, rsIDs, regions, drug names, trial IDs, SQL, free text. | `clinpgx` looking up a gene and drug |
| 3 | Your variants: chromosome, position and alleles read from your VCF or genotype file, sometimes as HGVS strings or rsIDs. No sample identifier is sent, but a set of variants can identify a person on its own. | `vcf-annotator` sending each variant to Ensembl VEP |
| 4 | Whole sequences or files uploaded to a third-party service or a hosted model. | `gi-annotation` uploading a FASTA sequence |
| 5 | Your chat messages, photos and voice notes to a hosted language model. | The RoboTerri bot |

## Class 3: skills that send your variants

| Skill | Host(s) | What is sent | When | Offline |
|-------|---------|--------------|------|---------|
| `clinical-variant-reporter` | rest.ensembl.org, grch37.rest.ensembl.org (VEP) | Each variant as a region or HGVS string. No sample ID, no genotype for other loci. | Every non-demo run | `--demo` uses the bundled evidence cache and makes no call |
| `vcf-annotator` | rest.ensembl.org (VEP), gnomad.broadinstitute.org (GraphQL), eutils.ncbi.nlm.nih.gov (ClinVar) | Per variant: `chrom:g.posref>alt` to VEP, `chrom-pos-ref-alt` to gnomAD, the rsID to ClinVar. NCBI requests carry `NCBI_TOOL_EMAIL` if set. | Every non-demo run | `--demo` uses bundled fixtures |
| `variant-annotation` | rest.ensembl.org (VEP region endpoint) | A batch of variant strings in one POST body. | On a cache miss. Results are cached under `~/.clawbio/variant_annotation_cache` with a TTL. | Runs from the cache only for variants already seen |

There is currently no local Ensembl VEP backend for these three skills. An
air-gapped deployment must exclude them or supply its own annotation source.

## Class 4: skills that upload sequences or files

| Skill | Host(s) | What is sent | When | Credential |
|-------|---------|--------------|------|------------|
| `gi-annotation`, `gi-chromatin`, `gi-enhancer`, `gi-expression`, `gi-promoter`, `gi-splice` | api.genomicintelligence.ai (or `GI_BASE_URL`) | The entire nucleotide sequence from your FASTA, in the request body, through `clawbio/gi/gi_client.py`. | Every non-demo run. Fails without a key; there is no local fallback. | `GI_API_KEY` |
| `galaxy-bridge` | usegalaxy.org by default, or `GALAXY_URL` (usegalaxy.eu is the other documented server) | Only with `--run`: your input file is uploaded to a history on that Galaxy server and the tool runs there. Without `--run` the skill only fetches the public tool catalogue (class 1). | `--run` only | `GALAXY_API_KEY` |
| `flow-bio` | app.flow.bio (or `FLOW_URL`) | Queries against your own Flow account and, if you use the upload path, the sample files you name. | On run | `FLOW_TOKEN`, or `FLOW_USERNAME` and `FLOW_PASSWORD` |

Each `gi-*` SKILL.md carries the same warning the client enforces: do not submit
identifiable patient data to the hosted service without a data-use agreement.

## Class 2: skills that send what you asked about

| Skill | Host(s) | What is sent | Credential |
|-------|---------|--------------|------------|
| `gwas-lookup` | rest.ensembl.org, www.ebi.ac.uk (GWAS Catalog, eQTL Catalogue), gtexportal.org, r12.finngen.fi, pheweb.org, pheweb.jp, portaldev.sph.umich.edu, api.platform.opentargets.org, sashagusev.github.io | The rsID you pass with `--rsid`, and the gene and region derived from it. If that rsID is one of your own variants, it leaves the machine. Results cache locally; `--demo` uses pre-fetched data. | none |
| `fine-mapping` | rest.ensembl.org | Region and gene lookups for the locus being fine-mapped, used to annotate the report. Summary statistics are processed locally. | none |
| `locuscompare-region-render` | rest.ensembl.org (plus the class 1 downloads below) | Region and gene lookups for the plotted locus. | none |
| `pathway-enricher` | maayanlab.cloud (Enrichr) | Your gene list, in full. | none |
| `clinical-trial-finder` | clinicaltrials.gov, www.clinicaltrialsregister.eu, api.platform.opentargets.org, www.ebi.ac.uk, id.nlm.nih.gov and hl7.org terminology lookups | The gene, rsID, condition or free-text query you pass; a `--input` file of genes is sent term by term. | none |
| `clinpgx` | api.clinpgx.org | Gene, drug and star-allele names. | none |
| `omics-target-evidence-mapper` | api.platform.opentargets.org, rest.uniprot.org, clinicaltrials.gov, eutils.ncbi.nlm.nih.gov | Target and gene symbols. | none |
| `lit-synthesizer` | eutils.ncbi.nlm.nih.gov, api.biorxiv.org | Your free-text literature query. | none |
| `pubmed-summariser` | eutils.ncbi.nlm.nih.gov | Your free-text PubMed query. | none |
| `article-data-fetcher` | eutils, pmc, ftp and edata at ncbi.nlm.nih.gov, api.crossref.org, api.datacite.org, api.figshare.com, datadryad.org, www.ebi.ac.uk | Article identifiers (PMID, PMCID, DOI). | none |
| `protocols-io` | www.protocols.io | Search terms and protocol identifiers, against your account. | `PROTOCOLS_IO_ACCESS_TOKEN` |
| `labstep` | The Labstep service, through the `labstepPy` client | Reads and writes to your own Labstep workspace: experiment, protocol and inventory records you ask for or create. `--demo` is offline synthetic data. | `LABSTEP_API_KEY` |
| `illumina-bridge` | ica.illumina.com (or `ILLUMINA_ICA_BASE_URL`) | Metadata queries against your own Illumina Connected Analytics projects. | `ILLUMINA_ICA_API_KEY` |
| `xena-tcga-gene-query` | biotree.top:38123 over plain HTTP by default (or `UCSCXENA_API_BASE_URL`) | Gene symbols and cohort names. Note the default endpoint is unencrypted HTTP to a third-party host. | `UCSCXENA_API_KEY` |
| `bigquery-public` | Google BigQuery (www.googleapis.com) | The SQL you write against public datasets, under your Google credentials. | Google application default credentials |
| `ukb-navigator` | Voyage AI, only if `VOYAGE_API_KEY` is set | The text of your query, for embedding. Without the key the skill uses the local ChromaDB default embedding and makes no call. | `VOYAGE_API_KEY` (optional) |
| `bgpt-mcp` | bgpt.pro (remote MCP server at `/mcp/sse` and `/mcp/stream`) | Your free-text literature question. This skill has no code of its own; the agent connects to the remote server directly. Free for the first results, key thereafter. | optional bgpt.pro key |

## Class 1: skills that download reference data only

Nothing derived from your data is sent. These skills fetch public files and
cache them locally where noted.

| Skill | Host(s) | What is fetched | Credential |
|-------|---------|-----------------|------------|
| `gwas-prs` | www.pgscatalog.org, ftp.ebi.ac.uk | PGS Catalog score metadata and harmonised scoring files. Your genotypes are scored locally and never leave. | none |
| `just-prs-mcp` | PyPI, through `uvx`, on first run | The pinned `just-prs-mcp` package. The MCP server it starts is local stdio, not remote. | none |
| `eqtl-catalogue-region-fetch` | ftp.ebi.ac.uk, www.ebi.ac.uk | eQTL Catalogue summary statistics for a region. Cache: `EQTL_CATALOGUE_CACHE_DIR`. | none |
| `gwas-catalog-region-fetch` | ftp.ebi.ac.uk | GWAS Catalog summary statistics for a region. Cache: `GWAS_CATALOG_CACHE_DIR`. | none |
| `ld-1000g-region-compute` | ftp.1000genomes.ebi.ac.uk, www.cog-genomics.org | 1000 Genomes reference genotypes for the region, and the plink binary if `PLINK_BIN` is unset. LD is computed locally. | none |
| `ukb-ppp-region-fetch` | www.synapse.org, repo-prod.prod.sagebase.org | UKB-PPP summary statistics for a region. Cache: `UKB_PPP_CACHE_DIR`. | `SYNAPSE_AUTH_TOKEN` |
| `bioconductor-bridge` | bioconductor.org | Package metadata for recommendations. | none |
| `busco-assessor` | eutils.ncbi.nlm.nih.gov, ftp.ensemblgenomes.org | Taxonomy lookups and BUSCO lineage datasets. | `NCBI_API_KEY` (optional) |
| `deepspot-m` | huggingface.co | Model weights on first use. Inference is local. | none |
| `proteomics-clock` | raw.githubusercontent.com | Published coefficient tables, cached under `CLAWBIO_CACHE`. | none |
| `nfcore-rnaseq-wrapper`, `nfcore-sarek-wrapper`, `nfcore-scrnaseq-wrapper` | nf-co.re, github.com, and the container registries the pipeline declares | Nextflow pulls the pinned pipeline and its containers. Your samples stay in the local work directory. Set `NXF_OFFLINE=true` with pre-pulled assets to forbid all of it. | none (Sentieon licence variables for that sarek path only) |

## Class 5: the RoboTerri bot

`bot/` (Telegram, Discord and WhatsApp adapters) is an optional conversational
interface, not part of the skill library. Every message, photo and voice note
you send it goes to the language model at `LLM_BASE_URL` using `LLM_API_KEY`
(falling back to `OPENAI_API_KEY`); the default model is `gemini-2.0-flash`
through an OpenAI-compatible endpoint, and voice replies use the `tts-1` model
on the same key. The messaging platforms themselves see every message by
construction. `drug-photo` has no code of its own: the photo is interpreted by
the vision model of whichever agent is driving the skill.

## Viewer assets

Two skills reference public CDNs for display only. `struct-predictor` fetches
3Dmol.js from cdnjs.cloudflare.com when building its HTML viewer and falls back
to an inline copy; `analyze-fasta` HTML reports reference fonts.googleapis.com
when opened in a browser. Neither sends any data.

## Checked and found to make no outbound call

These were on the list of things to check because their names or descriptions
suggest a service, and they do not call one: `bioqc-mcp` (a local stdio MCP
server), `turingdb-graph` (localhost only), `affinity-proteomics`, `wgs-prs`,
`genome-compare`, `pharmgx-reporter` (the Ensembl URLs in its source are
comments and citation links), and `drug-photo` (no code).

## Running ClawBio with egress blocked

- Everything not in classes 2 to 5 works with outbound network blocked once
  class 1 assets are cached.
- Class 3 skills need Ensembl VEP, gnomAD or ClinVar. Until a local backend
  exists, exclude them or route their hosts through an internal mirror you
  control.
- Never put an API key in a skill's input file or output directory.
  Reproducibility bundles record commands and file paths, so review them
  before sharing.
- ClawBio is a research tool, not a medical device. This page is a factual
  description of network behaviour, not a compliance statement under GDPR,
  HIPAA or any other regime. Your institution's own review still applies.

## Keeping this page true

`tests/test_data_handling_doc.py` scans every skill for outbound-call code
(`requests`, `urllib`, `httpx`, the Genomic Intelligence client, Nextflow and
`uvx` launches, cloud SDKs) and for declared remote endpoints in SKILL.md
frontmatter, and fails when a skill it finds is not named on this page. When
you add a network call to a skill, add a row here in the same pull request and
say what is sent.

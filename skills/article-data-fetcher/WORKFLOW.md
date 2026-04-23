# article-data-fetcher — Workflow Diagram

```mermaid
flowchart TD
    A([User Input]) --> B{Identifier type?}

    B -->|doi.org URL\nor 10.xxx/...| C[resolve_doi\nCrossref API\n→ title, abstract, links, DOI]
    B -->|PMID or PMID:xxx| D[resolve_pmid\nNCBI esummary\n→ title, DOI from articleids]
    B -->|PMC123 or PMC URL| E[resolve_pmcid\nFetch PMC HTML\n→ title, DOI from meta tag,\nfull page text + all hrefs]
    B -->|Repository URL\ne.g. GEO / Zenodo| F[Pass URL text\nthrough pattern scan]
    B -->|--demo flag| DEMO([Demo mode\nGSE100866\nSkip to file listing])

    C & D & E & F --> G[discover_accessions\nScan text + links\nfor known patterns]

    G --> H{Accessions\nfound?}

    H -->|Yes| I[Fetch file listings]
    H -->|No| FB[Fallback: search\nexternal repositories]

    FB --> FB1[DataCite API\nrelatedIdentifier = DOI]
    FB --> FB2[Zenodo API\nrelated_identifier = DOI\nor title search]
    FB --> FB3[NCBI GEO eSearch\nby PMID / DOI / title\n→ GSE accessions]
    FB --> FB4[ENA portal search\nstudy_doi or pubmed_id]
    FB --> FB5[Dryad API\nq = DOI]

    FB1 & FB2 & FB3 & FB4 & FB5 --> FBM[Merge + deduplicate\nall fallback files]
    FBM --> H2{Any files\nfound?}
    H2 -->|No| STOP([⚠ Tell user to\nprovide accession manually])
    H2 -->|Yes| SHOW

    I --> I1{Repo type?}
    I1 -->|GSE...| GEO[fetch_geo_files\nParse GEO FTP\nApache listing\n→ filename + size]
    I1 -->|PRJNA / SRP / ERP| ENA[fetch_ena_files\nENA Portal API\n→ fastq / submitted FTP]
    I1 -->|10.5281/zenodo.NNN\nor zenodo.org/records/NNN| ZEN[fetch_zenodo_files\nZenodo REST API\n→ files + sizes + checksums]
    I1 -->|zenodo.org/communities/X| ZCOM[fetch_zenodo_community_files\nZenodo community API\n→ all records + files]
    I1 -->|E-XXXX-NNN| AE[fetch_ena_files\nArrayExpress via ENA]

    GEO & ENA & ZEN & ZCOM & AE --> SHOW

    SHOW[print_file_listing\n─────────────────────────\n# FILENAME SIZE DATA TYPE SOURCE\nwith infer_data_type tags]

    SHOW --> SEL{User selects\nfiles}
    SEL -->|number / range\ne.g. 1,3,5-7| NUM[Select by index]
    SEL -->|extension\ne.g. vcf,csv.gz| EXT[file_matches_types]
    SEL -->|all| ALL[All files]
    SEL -->|none / Enter| CANCEL([Download cancelled])

    NUM & EXT & ALL --> WARN{Total size\n> 10 GB?}
    WARN -->|Yes| CONF{User confirms?}
    CONF -->|No| CANCEL
    CONF -->|Yes| DL
    WARN -->|No| DL

    DL[download_file\nStream with tqdm progress bar\nVerify MD5 checksum if available]
    DL --> UZ{--no-unzip\nflag?}
    UZ -->|No = default| UNZIP[unzip_file\nDecompress .gz\nDelete original]
    UZ -->|Yes| KEEP[Keep compressed]

    UNZIP & KEEP --> OUT[Write outputs\nmanifest.json\nreport.md]
    OUT --> DONE([✅ Done])

    style DEMO fill:#e8f5e9,stroke:#388e3c
    style STOP fill:#ffebee,stroke:#c62828
    style CANCEL fill:#fff3e0,stroke:#e65100
    style DONE fill:#e8f5e9,stroke:#388e3c
```

## Summary of steps

| Step | What happens |
|---|---|
| **1. Parse identifier** | Normalise input → DOI, PMID, PMCID, or URL |
| **2. Resolve article** | Fetch title + DOI + full text from Crossref / NCBI / PMC |
| **3. Discover accessions** | Regex-scan text and hrefs for GEO, ENA, Zenodo, Figshare, Dryad, OSF patterns |
| **4a. Fetch file listings** | Hit repository APIs to get filename, size, URL, checksum |
| **4b. Fallback search** | If no accessions found, query DataCite, Zenodo, GEO, ENA, Dryad by DOI/title |
| **5. Tag data types** | Keyword rules infer `count matrix`, `metadata`, `variant calls`, etc. |
| **6. User selects files** | By number, range, extension, or `all` |
| **7. Download** | Stream with progress bar, verify MD5, auto-unzip `.gz` by default |
| **8. Write outputs** | `manifest.json` (paths + checksums) and `report.md` |

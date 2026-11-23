# PhageViz: Product Requirements Document

| | |
|---|---|
| **Product** | PhageViz, a bacteriophage genome annotation viewer |
| **Version** | 1.0 (MVP) |
| **Status** | Draft |
| **Owner** | _your name_ |
| **Last updated** | 2026-11-23 |

---

## 1. Problem statement

Phage researchers annotate genomes with command-line tools like Pharokka, or download them from NCBI. The result is a GenBank or GFF3 file, which is hard to read as text. To understand a phage genome's layout (where the lysis cassette sits, how the structural genes are arranged, how much of the genome is hypothetical), people use desktop tools such as SnapGene or Artemis, or write one-off plotting scripts.

Students and wet-lab biologists who don't use the command line have no quick way to answer basic questions about a phage:

- What does this genome look like, and how are its genes organised?
- Does it have the genes it needs to be a complete, viable phage?
- How many of its genes have unknown function?

**PhageViz** answers these in a browser. It needs no installation and no command line, and it loads a genome from a file or from NCBI in seconds.

## 2. Goals

1. Show any phage genome as an interactive, function-coloured map in under 10 seconds.
2. Accept the formats people already have: GenBank, GFF3, FASTA, or just an NCBI accession.
3. Give a quick completeness check through hallmark genes.
4. Let users take results away as files (CSV, FASTA, GenBank, figures).

### Non-goals (v1)

- Full functional annotation against large databases (that is Pharokka's job)
- User accounts, saved projects, or a database
- Comparing several genomes (planned for v2)
- Lifestyle prediction, host prediction, or prophage detection
- Bacterial genomes or metagenomic assemblies

## 3. Target users

| Persona | Description | Main need |
|---|---|---|
| **Graduate researcher** | Isolates and sequences phages and runs Pharokka | Quick visual check and publication-ready map figures |
| **Bioinformatics student** | Learning phage genomics | Exploring reference phages by accession without setting up tools |
| **Wet-lab biologist** | Works with phage therapy candidates, limited coding | Checking a candidate for lysogeny-related genes (integrase, repressor) |
| **Instructor** | Teaches microbiology or genomics | Live demos of lambda, T4 and T7 in class |

## 4. User stories

### Loading data

| ID | Story | Priority |
|---|---|---|
| US-01 | As a researcher, I want to upload my Pharokka GenBank file so I can see my annotated phage. | Must |
| US-02 | As a student, I want to type an NCBI accession like `NC_001416` so I can load a reference phage without downloading files. | Must |
| US-03 | As a student, I want to search NCBI by phage name so I can find a genome when I don't know its accession. | Should |
| US-04 | As a new visitor, I want a one-click demo genome so I can see what the app does before using my own data. | Must |
| US-05 | As a researcher, I want to upload a raw FASTA and have genes predicted so I can explore a genome I haven't annotated yet. | Should |
| US-06 | As a user, I want a clear error message when my file is invalid or the accession doesn't exist so I know what to fix. | Must |

### Exploring

| ID | Story | Priority |
|---|---|---|
| US-07 | As a user, I want genes coloured by function on a strand-aware map so I can see the genome's modular layout. | Must |
| US-08 | As a user, I want to zoom and pan along the map so I can inspect a specific region. | Must |
| US-09 | As a user, I want to click a gene and see its product, coordinates and protein sequence. | Must |
| US-10 | As a user, I want GC content and GC skew tracks under the map so I can spot unusual regions. | Should |
| US-11 | As a user, I want to filter the map and table by functional category. | Should |
| US-12 | As a user, I want to search the gene table by product name or locus tag. | Must |
| US-13 | As a researcher, I want summary statistics (length, GC%, coding density, % hypothetical) for my paper's methods section. | Must |
| US-14 | As a wet-lab biologist, I want to see whether hallmark genes (terminase, portal, capsid, endolysin, holin, integrase) are present. | Must |
| US-15 | As a researcher, I want to rotate the genome to start at the terminase large subunit, following convention. | Could |

### Exporting

| ID | Story | Priority |
|---|---|---|
| US-16 | As a user, I want to download the gene table as CSV. | Must |
| US-17 | As a user, I want to download all protein sequences as FASTA for downstream analysis (BLAST, HHpred). | Must |
| US-18 | As a researcher, I want to download the map as SVG or PNG for figures. | Should |
| US-19 | As a user, I want to download the genome as GenBank, including predicted genes and rotation. | Could |

## 5. Functional requirements

### FR-1: Input

- FR-1.1 Accept GenBank, GFF3 (with or without an embedded `##FASTA` section) and FASTA uploads up to 50 MB.
- FR-1.2 Detect the format from content, not just the file extension.
- FR-1.3 Handle files with several records and let the user pick one.
- FR-1.4 Provide at least three bundled demo genomes (lambda, T7, T4).
- FR-1.5 Reject genomes over 1,000,000 bp with an explanation.

### FR-2: NCBI integration

- FR-2.1 Fetch a nucleotide record in GenBank format from an accession, with or without a version number.
- FR-2.2 Validate the accession format before calling NCBI.
- FR-2.3 Before fetching, check the record's length and title through a summary call, and refuse records over the size limit.
- FR-2.4 Search NCBI Nucleotide by free text and return up to 20 results with accession, title and length.
- FR-2.5 Cache fetched records so repeated requests don't call NCBI again within 24 hours.
- FR-2.6 Respect NCBI usage rules: send an email address with every request, and stay within 3 requests per second (10 with an API key).
- FR-2.7 Show specific messages for: accession not found, network error, NCBI unavailable, record too large, missing configuration.
- FR-2.8 Read the email and optional API key from Streamlit secrets, never from code.

### FR-3: Processing

- FR-3.1 Turn every input into the same internal data model (see SPEC §4).
- FR-3.2 Assign each gene one of the 10 PHROG functional categories: from PHROG annotations when present, otherwise from keyword rules.
- FR-3.3 Call genes with `pyrodigal` for FASTA input and for GenBank or GFF3 files with no CDS features.
- FR-3.4 Compute GC content and GC skew in sliding windows (window 100–5,000 bp, default 1,000).
- FR-3.5 Compute summary statistics and the hallmark check.
- FR-3.6 Rotate the genome to a chosen gene, splitting features that cross the new origin.

### FR-4: Display

- FR-4.1 Interactive genome map: forward strand above the axis, reverse below, one colour per category, range slider.
- FR-4.2 Clicking a gene opens a detail panel.
- FR-4.3 GC and skew tracks aligned to the map's x-axis.
- FR-4.4 Searchable, sortable gene table.
- FR-4.5 Summary tab with metric cards, a category chart and a gene length histogram.
- FR-4.6 Hallmarks tab with present or absent status and the matching gene IDs.

### FR-5: Export

- FR-5.1 CSV of the gene table (columns defined in SPEC §8).
- FR-5.2 Protein FASTA of all CDS.
- FR-5.3 GenBank of the current genome state.
- FR-5.4 SVG and PNG of the genome map.

## 6. Non-functional requirements

| Category | Requirement |
|---|---|
| **Performance** | A 50 kb genome parses and renders in under 3 seconds; a 170 kb genome (T4) in under 8 seconds. An NCBI fetch finishes in under 10 seconds under normal conditions. |
| **Memory** | Peak use under 500 MB for a 1 Mb genome, so it fits free hosting tiers. |
| **Usability** | A first-time user loads a demo and clicks a gene within 1 minute with no instructions. |
| **Reliability** | Malformed input never crashes the app; it always shows a message. |
| **Privacy** | Uploads stay in memory for the session only and are never written to disk or logged. |
| **Security** | No execution of uploaded content. Secrets stored in Streamlit secrets only. |
| **Compatibility** | Latest Chrome, Firefox, Safari and Edge on desktop. Mobile is usable but not optimised. |
| **Maintainability** | Bioinformatics logic lives in `src/phageviz/` with no Streamlit imports, and has at least 80% test coverage. |

## 7. UI/UX flow

```
┌─────────────────────────────────────────────────────────────┐
│ SIDEBAR                     │ MAIN AREA                     │
│                             │                               │
│ Input source                │ [Map] [Table] [Summary]       │
│ ( ) Demo                    │ [Hallmarks]                   │
│ ( ) Upload                  │                               │
│ ( ) NCBI                    │ Genome title · length · GC%   │
│                             │                               │
│ ── depends on source ──     │ ┌───────────────────────────┐ │
│ Demo:   [lambda ▾]          │ │ genome map (+ strand)     │ │
│ Upload: [drop file]         │ │ ───────────────────────── │ │
│ NCBI:   [accession/name]    │ │ genome map (− strand)     │ │
│         [Search] [Fetch]    │ ├───────────────────────────┤ │
│         results list        │ │ GC content track          │ │
│                             │ ├───────────────────────────┤ │
│ Record [▾] (if several)     │ │ GC skew track             │ │
│                             │ └───────────────────────────┘ │
│ Display                     │                               │
│ Categories [multiselect]    │ Selected gene                 │
│ GC window  [slider]         │ A · terminase large subunit   │
│ Start at   [gene ▾]         │ 711–2,636 (+) · 641 aa        │
│                             │ MVSK... [copy]                │
└─────────────────────────────────────────────────────────────┘
```

### Main flows

**Flow A: first visit (demo)**
Open app → lambda demo loads by default → map shows → user clicks a gene → detail panel updates.

**Flow B: NCBI accession**
Choose NCBI → type `NC_001604` → Fetch → spinner ("Fetching from NCBI…") → summary check passes → record downloads → map shows. If something fails, an error box explains why and suggests a fix.

**Flow C: NCBI search**
Choose NCBI → type `Pseudomonas phage` → Search → table of up to 20 results (accession, title, length) → select one → Fetch → same as Flow B.

**Flow D: upload raw FASTA**
Upload `.fasta` → info banner ("No annotation found, predicting genes with pyrodigal") → map shows with keyword-based categories and a notice that categories are approximate.

### UX principles

- The app is never empty: a demo loads on first visit.
- Every slow step shows a spinner with a description of what it's doing.
- Errors say what happened and what to do next.
- Approximate results (keyword categories, predicted genes) are clearly labelled.

## 8. Success metrics

| Metric | Target | How measured |
|---|---|---|
| Time to first map (demo) | < 5 s | Manual timing on deployed app |
| Time to first map (NCBI fetch, 50 kb) | < 10 s | Manual timing |
| Parser correctness | CDS count matches source for all demo and test files | pytest |
| Gene-calling recall vs RefSeq | ≥ 90% on lambda, T7, T4 | Benchmark script |
| Crash rate on malformed input | 0 unhandled exceptions across the test corpus | pytest with bad-input fixtures |
| Test coverage of `src/` | ≥ 80% | `pytest --cov` |
| Usability | 3 of 3 test users load a phage by accession without help | Informal user testing |
| Adoption (portfolio) | Deployed with a public link and README screenshots | Checklist |

## 9. Release plan

| Milestone | Scope | Stories |
|---|---|---|
| **M1: Core viewer** | Data model, GenBank parser, demo genomes, map, gene table, click-to-inspect | US-01, 04, 07, 08, 09, 12 |
| **M2: NCBI** | Accession fetch, search, caching, error handling, secrets config | US-02, 03, 06 |
| **M3: Analysis** | GC/skew tracks, summary tab, hallmarks, category filter | US-10, 11, 13, 14 |
| **M4: Raw input and export** | FASTA path with pyrodigal, GFF3 parser, CSV, FASTA, figures | US-05, 16, 17, 18 |
| **M5: Polish and ship** | Rotation, GenBank export, tests, CI, Docker, deployment, README | US-15, 19 |

## 10. Risks and mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| NCBI downtime or rate limiting | NCBI feature fails | Caching, retries with backoff, clear message, demos and upload still work |
| Pharokka output format changes between versions | Categories lost | Parse defensively, fall back to keyword rules, test against several Pharokka versions |
| Large genomes exhaust free-tier memory | App crashes | Size checks before fetch and after upload, numpy-based metrics |
| Keyword categories mislead users | Wrong conclusions | Label them "approximate" and recommend Pharokka |
| User fetches a non-phage record | Confusing output | Show the record title and warn when "phage" or "virus" is absent from it, but still allow it |

## 11. Future work (v2+)

- Compare several genomes, with synteny links
- Lifestyle prediction (lytic vs temperate)
- Optional PHROG HMM search with `pyhmmer`
- Batch NCBI fetch for several accessions
- Circular map view

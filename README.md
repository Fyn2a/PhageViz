# PhageViz

**Interactive bacteriophage genome annotation viewer, built with Streamlit.**

Load a phage genome from a file or straight from NCBI. PhageViz shows it as an interactive map coloured by function, with GC content and GC skew tracks, a searchable gene table, genome statistics and a hallmark-gene check. You can download every result.

<!-- Replace with a real screenshot or GIF once the app runs -->
![PhageViz screenshot](docs/screenshot.png)

> **Live demo:** _add your Streamlit Community Cloud or Hugging Face Spaces link here_

---

## Features

- **Three ways to load a genome**
  - Upload a GenBank, GFF3 or FASTA file
  - Fetch a genome from **NCBI** by accession (e.g. `NC_001416`) or search by name
  - One-click demo genomes (phage lambda, T7, T4)
- **Interactive genome map:** forward and reverse strands, genes coloured by PHROG functional category, zoom, pan and click-to-inspect
- **GC content and GC skew tracks** with an adjustable window size
- **Gene table:** searchable and filterable, with protein sequences
- **Genome summary:** length, GC%, coding density, % hypothetical proteins, strand bias, category breakdown
- **Hallmark check:** whether terminase, portal, major capsid, endolysin and holin genes are present
- **Gene calling for raw FASTA** with `pyrodigal`, plus keyword-based functional categories
- **Genome rotation:** set any gene (usually the terminase large subunit) as position 1
- **Downloads:** CSV, protein FASTA, GenBank, and map figures as SVG or PNG

---

## Quick start

### 1. Clone and install

```bash
git clone https://github.com/<your-username>/phageviz.git
cd phageviz
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Requires **Python 3.10+**.

### 2. Configure NCBI access (needed for the NCBI feature only)

NCBI asks every Entrez user to give an email address. Create `.streamlit/secrets.toml`:

```toml
[ncbi]
email = "you@example.com"
api_key = ""   # optional: raises the rate limit from 3 to 10 requests per second
```

You can get a free API key from your NCBI account settings. **Never commit `secrets.toml`.** It is already listed in `.gitignore`.

Without this file, the app still runs, and the NCBI tab shows a configuration message.

### 3. Run

```bash
streamlit run app.py
```

The app opens at `http://localhost:8501`.

### Run with Docker

```bash
docker build -t phageviz .
docker run -p 8501:8501 -v $(pwd)/.streamlit:/app/.streamlit phageviz
```

---

## Usage

1. Choose an input source in the sidebar: **Demo**, **Upload** or **NCBI**.
2. For NCBI, type an accession (`NC_001416`) or a name (`Escherichia phage T7`), pick a result and click **Fetch**.
3. If the file holds more than one record, pick the one you want.
4. Look through the tabs: **Map**, **Table**, **Summary**, **Hallmarks**.
5. Click a gene on the map to see its details and protein sequence.
6. Download results from the **Table** and **Map** tabs.

---

## Inputs

| Input | Formats | Notes |
|---|---|---|
| Annotated genome | GenBank (`.gb`, `.gbk`, `.gbff`), GFF3 (`.gff`, `.gff3`) | Pharokka output is supported; PHROG categories are read when present |
| Raw genome | FASTA (`.fa`, `.fasta`, `.fna`) | Genes are called with `pyrodigal` |
| NCBI accession | RefSeq (`NC_001416`) or GenBank (`J02459`), version optional | Fetched as GenBank through NCBI E-utilities |
| NCBI search | Free-text name | Returns up to 20 matching complete genomes |
| Settings | Sidebar | GC window size, category filter, rotation start gene |

Limits: 50 MB upload, 1,000,000 bp per genome.

## Outputs

| Output | Where | Format |
|---|---|---|
| Genome map with GC and skew tracks | Map tab | Interactive; export as SVG or PNG |
| Gene detail panel | Map tab | On screen |
| Gene table | Table tab | Interactive; export as CSV |
| Protein sequences | Table tab | FASTA |
| Annotated genome | Table tab | GenBank |
| Genome statistics and category chart | Summary tab | On screen |
| Hallmark presence/absence | Hallmarks tab | On screen |

---

## Project structure

```
phageviz/
├── app.py                    # Streamlit UI
├── requirements.txt
├── Dockerfile
├── .streamlit/
│   ├── config.toml           # theme and upload size limit
│   └── secrets.toml          # NCBI email and key (not committed)
├── data/
│   ├── demo/                 # lambda, T7, T4 GenBank files
│   └── phrog_annot.tsv       # PHROG ID -> annotation and category
├── src/phageviz/
│   ├── models.py             # GenomeRecord data model
│   ├── parsers.py            # GenBank, GFF3, FASTA
│   ├── ncbi.py               # Entrez search and fetch
│   ├── gene_calling.py       # pyrodigal wrapper
│   ├── annotation.py         # functional categories, hallmarks
│   ├── metrics.py            # GC, skew, summary stats
│   ├── transforms.py         # genome rotation
│   ├── exporters.py          # CSV, FASTA, GenBank
│   └── plots.py              # Plotly figures
├── tests/
└── docs/
    ├── PRD.md
    └── SPEC.md
```

---

## Testing

```bash
pytest -q
```

The NCBI tests are marked `network` and skipped by default. Run them with `pytest -m network`.

## Validation

Gene calling compared with RefSeq annotations (a gene counts as a match when it has the same stop coordinate and strand):

| Phage | RefSeq CDS | Predicted | Precision | Recall |
|---|---|---|---|---|
| Lambda (NC_001416) | _tbd_ | _tbd_ | _tbd_ | _tbd_ |
| T7 (NC_001604) | _tbd_ | _tbd_ | _tbd_ | _tbd_ |
| T4 (NC_000866) | _tbd_ | _tbd_ | _tbd_ | _tbd_ |

---

## Limitations

- Keyword-based categories are rough. For publication-quality annotation, run [Pharokka](https://github.com/gbouras13/pharokka) and upload its GenBank output.
- `pyrodigal` is not phage-specific, so some gene boundaries will differ from PHANOTATE or RefSeq.
- Genetic code 11 is used by default.
- NCBI fetches depend on NCBI availability and rate limits.

## Data sources and credits

- **NCBI Nucleotide** via E-utilities and Biopython `Bio.Entrez`
- **PHROGs** functional categories (Terzian et al., 2021, *NAR Genomics and Bioinformatics*)
- **Pyrodigal** (Larralde, 2022, *JOSS*), Cython bindings to Prodigal (Hyatt et al., 2010)
- **Biopython** (Cock et al., 2009, *Bioinformatics*)

## License

MIT. See `LICENSE`.

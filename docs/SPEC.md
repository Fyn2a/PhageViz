# PhageViz: Technical Specification

| | |
|---|---|
| **Version** | 1.0 |
| **Related** | [PRD.md](PRD.md), [README.md](../README.md) |
| **Last updated** | 2026-11-23 |

---

## 1. Overview

PhageViz is a single-process Streamlit app. It has **no backend API and no database**. Genome data comes from three sources (bundled demo files, user uploads and NCBI E-utilities), is turned into one in-memory data model, processed with pandas and numpy, and rendered with Plotly. Results leave the app as downloaded files.

## 2. Architecture

```
                         ┌──────────────────────────────┐
                         │         Browser (UI)         │
                         └──────────────┬───────────────┘
                                        │ websocket (handled by Streamlit)
┌───────────────────────────────────────▼──────────────────────────────────────┐
│ app.py — Streamlit UI layer                                                  │
│  sidebar widgets · tabs · session_state · cache decorators · error display   │
└───────┬───────────────┬──────────────────┬─────────────────┬────────────────┘
        │               │                  │                 │
┌───────▼──────┐ ┌──────▼──────┐  ┌────────▼────────┐ ┌──────▼──────┐
│ data/demo/   │ │ file upload │  │ ncbi.py         │ │ secrets     │
│ (static)     │ │ (bytes)     │  │ Bio.Entrez      │─┤ email, key  │
└───────┬──────┘ └──────┬──────┘  └────────┬────────┘ └─────────────┘
        │               │                  │  HTTPS
        │               │          ┌───────▼────────┐
        │               │          │ NCBI E-utilities│
        │               │          └───────┬────────┘
        └───────────────┴─────────┬────────┘  GenBank text
                          ┌───────▼────────┐
                          │ parsers.py     │  → GenomeRecord
                          └───────┬────────┘
                ┌─────────────────┼──────────────────┐
        ┌───────▼───────┐ ┌───────▼───────┐ ┌────────▼───────┐
        │ gene_calling  │ │ annotation    │ │ metrics        │
        │ (pyrodigal)   │ │ categories,   │ │ GC, skew,      │
        │               │ │ hallmarks     │ │ summary        │
        └───────┬───────┘ └───────┬───────┘ └────────┬───────┘
                └─────────────────┼──────────────────┘
                          ┌───────▼────────┐
                          │ transforms     │  rotation
                          └───────┬────────┘
                    ┌─────────────┴─────────────┐
            ┌───────▼───────┐           ┌───────▼───────┐
            │ plots.py      │           │ exporters.py  │
            │ Plotly figs   │           │ CSV/FASTA/GBK │
            └───────────────┘           └───────────────┘
```

**Layering rule:** nothing in `src/phageviz/` imports Streamlit. Modules take and return plain Python, pandas or Plotly objects, so they can be tested without the UI. `app.py` does the caching, state and error display.

## 3. Tech stack

| Component | Package | Version |
|---|---|---|
| Runtime | Python | ≥ 3.10 |
| UI | `streamlit` | ≥ 1.35 (needed for `on_select` in `st.plotly_chart`) |
| Data | `pandas`, `numpy` | ≥ 2.0, ≥ 1.26 |
| Bioinformatics | `biopython` | ≥ 1.83 |
| Gene calling | `pyrodigal` | ≥ 3.0 |
| Plots | `plotly` | ≥ 5.20 |
| Figure export | `kaleido` | ≥ 0.2 |
| Testing | `pytest`, `pytest-cov` | latest |

Pin exact versions in `requirements.txt` before release.

## 4. Data model

### 4.1 `GenomeRecord`

```python
@dataclass
class GenomeRecord:
    id: str                    # accession or record ID, e.g. "NC_001416.1"
    description: str           # record title / definition line
    sequence: str              # uppercase DNA, A/C/G/T/N only
    features: pd.DataFrame     # schema in 4.2
    source: Literal["demo", "upload", "ncbi"]
    annotation_origin: Literal["file", "pharokka", "predicted"]
    topology: Literal["linear", "circular"] = "linear"
    rotation_offset: int = 0   # 0 = original coordinates

    @property
    def length(self) -> int: return len(self.sequence)
```

### 4.2 Features DataFrame schema

| Column | dtype | Nullable | Description |
|---|---|---|---|
| `id` | `string` | no | Locus tag, or `{record_id}_{n:04d}` when missing. Unique within the record. |
| `type` | `category` (fixed: `CDS`, `tRNA`, `tmRNA`) | no | Feature type |
| `start` | `int64` | no | **1-based, inclusive.** Always `1 ≤ start ≤ end`. |
| `end` | `int64` | no | 1-based, inclusive. `end ≤ L`, except origin-crossing features (§4.4). |
| `strand` | `int8` | no | `1` or `-1` |
| `product` | `string` | no | Default `"hypothetical protein"` |
| `category` | `category` (fixed: the 10 in §4.3, in that order) | no | Functional category |
| `category_source` | `category` (fixed: `phrog`, `keyword`) | no | How `category` was assigned |
| `phrog` | `string` | yes | PHROG ID, e.g. `phrog_1` |
| `aa` | `string` | yes | Protein sequence without the trailing `*`. `<NA>` for RNA genes. |
| `length_aa` | `Int64` (nullable integer) | yes | `len(aa)`; `<NA>` when `aa` is `<NA>` |
| `wraps_origin` | `bool` | no | True if the feature crosses position L→1 (circular only) |

**dtypes.** `string` is `pd.StringDtype()` and `Int64` is pandas' nullable integer; both behave the
same on pandas 2.x and 3.x. Missing values are `pd.NA`, never `""` or `NaN`. numpy `int64` can't
hold missing values, which is why `length_aa` is `Int64`. Categorical columns use fixed category
lists so frames from different parsers concatenate and compare cleanly. Parsers build frames with
`models.features_from_rows()`, which sets dtypes and defaults; `models.validate_features()` checks
columns, dtypes and the invariants above.

**RNA genes.** `tRNA` and `tmRNA` have no PHROG category. They are always `category="other"`,
`category_source="keyword"`, `aa=<NA>`, `length_aa=<NA>`.

### 4.3 Functional categories (PHROG)

```
head and packaging
connector
tail
lysis
integration and excision
DNA, RNA and nucleotide metabolism
transcription regulation
moron, auxiliary metabolic gene and host takeover
other
unknown function
```

Each category has a fixed colour in `plots.CATEGORY_COLORS`, using a colour-blind-safe palette.

### 4.4 Coordinate conventions

| Source | Native convention | Conversion to internal |
|---|---|---|
| Biopython `FeatureLocation` | 0-based, half-open | `start + 1`, `end` |
| GFF3 | 1-based, inclusive | none |
| pyrodigal `Gene` | 1-based, inclusive (`begin`, `end`) | none |
| Python string slicing | 0-based, half-open | `seq[start - 1 : end]` |

**Origin-crossing features** (circular genomes only, e.g. Biopython `join(48001..48502,1..300)` on
a 48,502 bp genome) are stored **unwrapped**: `start` is the first base before the origin and `end`
is the last base after the origin **plus L**, with `wraps_origin=True`.

| Feature | `start` | `end` | `wraps_origin` |
|---|---|---|---|
| `join(48001..48502,1..300)`, L = 48,502 | 48001 | 48802 | True |

Invariants for wrapped rows: `start ≤ L < end < start + L`. For all other rows `end ≤ L` and
`wraps_origin=False`. Linear records never have wrapped rows. This keeps `start ≤ end` everywhere,
so `length = end - start + 1` and sorting by `start` work without special cases.

- Splitting into real parts (for drawing two bars or writing GenBank): `models.feature_parts(start,
  end, L)` returns `[(start, L), (1, end - L)]` for wrapped rows, `[(start, end)]` otherwise.
- Slicing a wrapped feature: `seq[start - 1 :] + seq[: end - L]`.
- No module stores `start > end`, including `transforms.py`.

### 4.5 Validation

`GenomeRecord` is a frozen dataclass (`eq=False`, since `==` on a DataFrame field is ambiguous)
that validates itself in `__post_init__` and raises `errors.InvalidRecord` on any violation:

- `source`, `annotation_origin` and `topology` are one of their allowed values
- `sequence` is non-empty and matches `[ACGTN]+` (parsers upper-case and replace other IUPAC codes
  before building the record, so they can report how many were replaced)
- `0 ≤ rotation_offset < length`
- `validate_features(features, length)` passes (§4.2 dtypes and nulls, unique `id`, strand,
  `start ≤ end`, §4.4 bounds, RNA-gene rules, no trailing `*`, `length_aa` consistent with `aa`)
- linear records have no `wraps_origin` rows

Size limits (§6.5) are input policy and are enforced by `parsers.load`, not by the model.
Transforms return new records, so every record in the app has passed these checks.

## 5. Module interfaces

```python
# models.py
CATEGORIES: tuple[str, ...]                                      # §4.3, in order
FEATURE_COLUMNS: tuple[str, ...]                                 # §4.2, in order
def empty_features() -> pd.DataFrame
def features_from_rows(rows: list[dict]) -> pd.DataFrame          # fills defaults, sets dtypes
def validate_features(df: pd.DataFrame, length: int | None = None) -> None   # raises on violation
def feature_parts(start: int, end: int, length: int) -> list[tuple[int, int]]

# parsers.py
@dataclass(frozen=True)
class LoadResult:
    records: tuple[GenomeRecord, ...]   # at least one
    warnings: tuple[str, ...] = ()      # user-facing; app.py shows them with st.warning

def decode(data: bytes) -> str                                   # gzip, BOM, encoding, newlines
def detect_format(text: str) -> Literal["genbank", "gff3", "fasta"]
def parse_genbank(text: str, source: str) -> LoadResult
def parse_gff3(text: str, source: str) -> LoadResult             # M4
def parse_fasta(text: str, source: str) -> LoadResult            # M4; features empty
def load(data: bytes, source: str) -> LoadResult

# ncbi.py
def configure(email: str, api_key: str | None) -> None
def validate_accession(acc: str) -> str                          # normalised or raises
def summarize(accession: str) -> NcbiSummary
def search(term: str, retmax: int = 20) -> list[NcbiSummary]
def fetch_genbank(accession: str) -> str                         # raw GenBank text

# gene_calling.py
def call_genes(record: GenomeRecord, table: int = 11) -> pd.DataFrame

# annotation.py
def categorize(product: str) -> str
def assign_categories(types: pd.Series, products: pd.Series,
                      functions: pd.Series) -> pd.DataFrame  # category, category_source
def apply_categories(df: pd.DataFrame, phrog_map: pd.DataFrame) -> pd.DataFrame
def hallmark_check(df: pd.DataFrame) -> pd.DataFrame

# metrics.py
def gc_content(seq: str) -> float
def gc_windows(seq: str, window: int, step: int) -> pd.DataFrame  # pos, gc, skew
def summary(record: GenomeRecord) -> dict

# transforms.py
def rotate(record: GenomeRecord, new_start: int) -> GenomeRecord

# exporters.py
def to_csv(record: GenomeRecord) -> bytes
def to_protein_fasta(record: GenomeRecord) -> bytes
def to_genbank(record: GenomeRecord) -> bytes

# plots.py
def genome_map(record: GenomeRecord, categories: list[str] | None) -> go.Figure
def gc_tracks(windows: pd.DataFrame, length: int) -> go.Figure
def category_chart(df: pd.DataFrame) -> go.Figure
def length_histogram(df: pd.DataFrame) -> go.Figure
```

```python
@dataclass
class NcbiSummary:
    accession: str      # AccessionVersion, e.g. "NC_001416.1"
    title: str
    length: int
    uid: str
```

## 6. Input specifications

### 6.1 Format detection (`detect_format`)

| First non-blank line starts with | Format |
|---|---|
| `LOCUS` | GenBank |
| `##gff-version 3` | GFF3 |
| `>` | FASTA |
| anything else | error `UnsupportedFormat` |

The file extension is only a hint for the uploader's file picker.

### 6.2 GenBank

- Parse with `SeqIO.parse(StringIO(text), "genbank")`, with `BiopythonWarning` silenced
  (cosmetic issues such as odd LOCUS spacing). Any exception from Biopython becomes `ParseError`.
- Keep feature types `CDS`, `tRNA`, `tmRNA`. Skip `gene` (it duplicates `CDS`) and `source`.
- Qualifiers read: `locus_tag`, `product`, `translation`, `transl_table`, `codon_start`, `phrog`
  (Pharokka), `function` (Pharokka category). Only the first value of each is used.
- **Product:** RNA genes without a `product` use their type (`"tRNA"`); CDS default to
  `"hypothetical protein"`.
- **PHROG ID:** `"15"`, `"phrog15"` and `"phrog_15"` (any case) become `"phrog_15"`; anything else
  (e.g. `"No_PHROG"`) is missing.
- **Translation:** the `translation` qualifier (trailing `*` stripped) if present. Otherwise the
  feature is extracted (strand- and join-aware, so it works across the origin), offset by
  `codon_start`, trimmed to whole codons and translated with `table=transl_table or 11,
  cds=False, to_stop=True`. If the first codon is a start codon of that table (e.g. GTG, TTG),
  the first residue is set to `M`, because alternative start codons initiate with Met; this is
  skipped for 5'-partial genes (`<` on +, `>` on −) and when `codon_start != 1`, since those don't
  begin at their real start. An unknown table gives `aa=<NA>`. (Lambda gene `A` starts with GTG;
  plain `cds=False` translation would wrongly start it with `V`.)
- **Coordinates** (§4.4): simple locations use `start + 1`, `end`. Joins use outer bounds. A join
  whose parts, put in genome order (minus-strand parts are listed in reading order by Biopython,
  so they are reversed), jump back once to a lower position crosses the origin and is stored
  unwrapped (`end + L`). Biopython's own `location.start/end` can't be used for these: it reports
  `0..L`.
- **Skipped features** (left out, with one warning per record naming up to three): no location,
  no single strand, unknown position, parts that jump back more than once, or crossing the origin
  in a record whose topology is `linear` (topology is never inferred).
- **IDs:** missing `locus_tag` → `{record_id}_{n:04d}`, where `n` is the 1-based position among
  kept features; duplicates get `_2`, `_3`, … skipping any suffix already in use.
- Topology from `rec.annotations.get("topology", "linear")`; any other value → `linear`.
- **Pharokka detection:** any CDS has a `phrog` qualifier → `annotation_origin="pharokka"`.
- No CDS features → treat like FASTA and run gene calling. **Until M4** this raises
  `UnsupportedFormat` ("has no annotated genes …").
- **Record-level failures** (too short or long, no sequence, no CDS, `InvalidRecord`): the record
  is skipped with a warning when other records load. If no record loads, the error itself is
  raised for a single-record file, and a `ParseError` summarising the first problem otherwise.

### 6.3 GFF3

- Read lines up to `##FASTA`; the rest is the sequence (FASTA). If there's no embedded sequence, show "GFF3 needs an embedded ##FASTA section".
- Columns: `seqid, source, type, start, end, score, strand, phase, attributes`.
- Attributes are split on `;` then `=`, and values URL-decoded.
- `product` from `product=`, falling back to `Name=`.

### 6.4 FASTA

- Each record becomes one `GenomeRecord` with empty features, then `call_genes` runs.
- Sequence is upper-cased. Characters outside `ACGTN` (IUPAC ambiguity codes) become `N`, with a warning showing the count.

### 6.5 Validation rules

| Rule | Limit | Behaviour |
|---|---|---|
| Upload size | 50 MB (`server.maxUploadSize` in `config.toml`) | Blocked by Streamlit |
| Genome length | 1,000,000 bp | Reject record, show message |
| Minimum length | 1,000 bp | Reject: "too short to be a phage genome" |
| Encoding | UTF-8 (BOM stripped), falling back to latin-1 | latin-1 decodes any bytes, so decoding never fails; binary input is rejected by `detect_format` |
| Line endings | `\r\n` and `\r` | Normalised to `\n` on decode |
| Gzip | Detected by magic bytes `1f 8b`, not file name | Decompress, then parse; corrupt data → `ParseError` |
| Decompressed size | 100 MB | `ParseError` (guards against gzip bombs; a 1 Mb genome as GenBank is ~3 MB) |
| Records per file | 50 | Load the first 50 with a warning; the rest are not parsed |

Length limits are inclusive (1,000 and 1,000,000 bp are accepted) and are enforced in
`parse_genbank`, so uploads and NCBI fetches share them. `load` takes no file name: the format
comes from content and file names are never logged, so it has no use.

M1 status: `load` detects GFF3 and FASTA but raises `UnsupportedFormat` ("… aren't supported yet")
until M4.

### 6.6 Input errors

| Exception | Trigger |
|---|---|
| `UnsupportedFormat` | Empty file, unrecognised first line, GFF3/FASTA before M4, GenBank record without CDS before M4 |
| `ParseError` | Malformed GenBank, corrupt or oversized gzip, no records, all records unusable |
| `MissingSequence` | GenBank record without sequence (CON record) |
| `RecordTooLarge` / `RecordTooShort` | Genome outside 1,000–1,000,000 bp (also used by NCBI, §7.6) |
| `InvalidRecord` | Parsed data breaks the model (§4.5), e.g. a feature past the sequence end |

All subclass `PhageVizError`; messages name the record ID and say what to do next.

## 7. NCBI Entrez integration

### 7.1 Configuration

`.streamlit/secrets.toml`:

```toml
[ncbi]
email = "you@example.com"
api_key = ""
```

```python
from Bio import Entrez
Entrez.email = email
Entrez.api_key = api_key or None
Entrez.tool = "phageviz"
Entrez.max_tries = 3
Entrez.sleep_between_tries = 5
```

If `email` is missing, the NCBI option shows a setup message and makes no calls.

### 7.2 Calls

These are the external API calls the app makes. PhageViz exposes no endpoints of its own.

| Purpose | E-utility | Parameters | Returns |
|---|---|---|---|
| Search | `esearch` | `db="nuccore"`, `term=<user term> AND "complete genome"[Title]`, `retmax=20`, `idtype="acc"` | List of accessions |
| Summary | `esummary` | `db="nuccore"`, `id=<accession list>` | Title, length (`Length`), `AccessionVersion` |
| Fetch | `efetch` | `db="nuccore"`, `id=<accession>`, `rettype="gbwithparts"`, `retmode="text"` | GenBank text |

`gbwithparts` is used instead of `gb` so the sequence is always included, even for records assembled from parts.

### 7.3 Accession validation

```python
ACCESSION_RE = re.compile(
    r"^(?:[A-Z]{2}_\d{6,9}|[A-Z]{1,2}\d{5,6}|[A-Z]{4}\d{8,10})(?:\.\d+)?$"
)
```

Input is trimmed and upper-cased first. RefSeq (`NC_001416`), GenBank (`J02459`, `MN123456`) and WGS-style accessions pass. Anything else is treated as a search term.

### 7.4 Fetch sequence

```
user input
   │
   ├─ matches ACCESSION_RE? ── no ──▶ search(term) ──▶ results table ──▶ user picks one
   │                                                                          │
   yes ◀──────────────────────────────────────────────────────────────────────┘
   │
   ▼
summarize(acc) ── not found ──▶ AccessionNotFound
   │
   ├─ length > 1,000,000 ──▶ RecordTooLarge
   │
   ▼
fetch_genbank(acc) ──▶ parse_genbank(text, source="ncbi") ──▶ GenomeRecord
```

The summary check costs one small request and avoids downloading a whole bacterial chromosome by mistake.

### 7.5 Rate limiting and caching

- Biopython's `Entrez` module already spaces requests to 3/s without a key and 10/s with one. No extra throttle is needed for a single user.
- In `app.py`:
  ```python
  @st.cache_data(ttl=86_400, show_spinner="Fetching from NCBI…")
  def cached_fetch(accession: str) -> str: ...

  @st.cache_data(ttl=3_600)
  def cached_search(term: str) -> list[NcbiSummary]: ...
  ```
- The cache key is the normalised accession, so `nc_001416` and `NC_001416` share an entry. An unversioned accession is cached under the version NCBI returns.
- Raw GenBank text is cached, not the parsed record, so the parser can change without invalidating the cache.

### 7.6 Errors

| Exception | Trigger | Message shown |
|---|---|---|
| `NcbiNotConfigured` | No email in secrets | "NCBI access isn't set up. Add your email to `.streamlit/secrets.toml` (see README)." |
| `InvalidAccession` | Fails regex in fetch mode | "`XYZ` doesn't look like an accession. Try `NC_001416`, or search by name." |
| `AccessionNotFound` | esummary returns no record or an error field | "NCBI has no record for `XYZ`." |
| `RecordTooLarge` | Length > limit | "`XYZ` is 4,641,652 bp. PhageViz handles genomes up to 1 Mb." |
| `NcbiUnavailable` | HTTP 5xx, timeout, `URLError` after retries | "NCBI isn't responding. Try again in a minute, or upload a file instead." |
| `NcbiRateLimited` | HTTP 429 after retries | "NCBI is limiting requests. Wait a moment, or add an API key." |
| `NoResults` | esearch returns 0 IDs | "No complete genomes matched. Try a broader name." |

All are subclasses of `PhageVizError` and are caught in `app.py` and shown with `st.error`. Unexpected exceptions are caught at the top level, logged without user data, and shown as a generic message.

### 7.7 Non-phage records

If the record title doesn't contain `phage` or `virus` (case-insensitive), show `st.warning("This record may not be a phage genome.")` and continue.

## 8. Processing

### 8.1 Category assignment (in priority order)

1. Pharokka `function` qualifier, if it matches a PHROG category → `category_source="phrog"`
2. `phrog` ID looked up in `data/phrog_annot.tsv` → `"phrog"`
3. Keyword rules on `product` (first match wins) → `"keyword"`

**Where each step runs.** The `function` qualifier has no column in the features schema, so steps
1 and 3 run in the parser through `annotation.assign_categories(types, products, functions)`,
which is vectorised over all features. Step 2 runs afterwards in
`apply_categories(df, phrog_map)`, which only upgrades rows whose `category_source` is `keyword`
and whose `phrog` ID is in the map; rows set from `function` are never touched, so the priority
order holds. `categorize(product)` is the single-value form of step 3: a plain-Python loop over the same
`KEYWORD_RULES` table (a one-row pandas path costs ~25 ms per call); a test checks both paths
agree on every T4 product.

**`function` matching.** The qualifier matches a category when it equals one of the 10 names in
§4.3 after trimming, collapsing whitespace and ignoring case. Anything else (including `tRNA` or
a missing value) falls through to the keyword rules. A missing or blank `product` is treated as
`hypothetical protein`.

| Category | Regex (case-insensitive) |
|---|---|
| connector | `head.tail\|connector\|neck\|adaptor` |
| lysis | `holin\|endolysin\|lysozyme\|spanin\|lysis\|amidase` |
| head and packaging | `terminase\|portal\|capsid\|scaffold\|packag\|head` |
| tail | `tail\|fib(?:er\|re)\|baseplate\|sheath\|tape.measure` |
| integration and excision | `integrase\|excisionase\|recombinase\|transposase` |
| DNA, RNA and nucleotide metabolism | `polymerase\|helicase\|primase\|nuclease\|ligase\|reductase\|kinase\|topoisomerase\|methyltransferase\|replication\|ssDNA` |
| transcription regulation | `repressor\|regulator\|anti.?terminat\|sigma\|activator\|transcription` |
| moron, auxiliary metabolic gene and host takeover | `superinfection\|exclusion\|anti-restriction\|toxin` |
| unknown function | `hypothetical\|unknown\|uncharacteri[sz]ed\|DUF` |
| other | no match |

Rules are tried top to bottom and the first match wins. `connector` is checked first, so
"head-tail connector" is caught neither by the broad `head` pattern nor by `tail`; tests lock
this order in. RNA genes skip these rules (§4.2). Groups are non-capturing (`(?:...)`) because
pandas warns on capture groups in `str.contains`.

The table extends the original draft with `amidase`, `head`, `topoisomerase`,
`methyltransferase`, `replication`, `ssDNA`, `activator`, `anti.?terminat` and `transcription`
(`prohead` and `anti-sigma` were dropped as already covered). Measured on the demo genomes,
products left in `other` fell from 25 to 16 (lambda, 73 CDS), 10 to 5 (T7, 60) and 59 to 40
(T4, 278). Known debatable case: "inhibitor of host transcription" lands in transcription
regulation rather than host takeover.

### 8.2 Gene calling

```python
finder = pyrodigal.GeneFinder(meta=True)
genes = finder.find_genes(record.sequence.encode())
```

- `meta=True` because single-genome training needs ≥ 20 kb and gives poor results on small phages.
- Output rows: `start=g.begin`, `end=g.end`, `strand=g.strand`, `aa=g.translate(translation_table=11).rstrip("*")`, `product="hypothetical protein"`, `id=f"{record.id}_{n:04d}"`.
- Then set `annotation_origin="predicted"` and show a banner.

### 8.3 Metrics

| Metric | Formula |
|---|---|
| GC content | `(G + C) / (A + C + G + T)` (N excluded) |
| GC skew | `(G − C) / (G + C)`, 0 when `G + C = 0` |
| Windows | Window `w` (default 1,000), step `w // 2`, computed with cumulative sums (O(n)) |
| Coding density | Bases covered by ≥ 1 CDS / genome length (overlaps counted once) |
| % hypothetical | CDS with category `unknown function` / total CDS |
| Strand bias | CDS on + strand / total CDS |
| Mean CDS length | Mean `length_aa` |

### 8.4 Hallmark check

| Hallmark | Pattern on `product` | Expected in |
|---|---|---|
| Terminase large subunit | `terminase.*large\|large.*terminase\|TerL` | All tailed phages |
| Portal | `portal` | All tailed phages |
| Major capsid | `major (head\|capsid)\|MCP` | All tailed phages |
| Endolysin | `endolysin\|lysozyme\|amidase` | Lytic release |
| Holin | `holin` | Lytic release |
| Integrase | `integrase` | Temperate phages |
| Repressor | `repressor\|CI protein` | Temperate phages |

Integrase and repressor are shown as **informational**, not pass/fail. The tab states that missing hallmarks often mean incomplete annotation rather than a missing gene.

### 8.5 Rotation

```
new_pos = ((old_pos - offset - 1) % L) + 1      # offset = new_start - 1
```

`new_start` must be in `1..L`; `new_start = 1` returns the record unchanged. Each feature's real
parts (`models.feature_parts`) are moved with the formula; if the moved feature now crosses the
new origin it is stored unwrapped (`end + L`, `wraps_origin=True`, §4.4), and a feature that no
longer crosses it becomes a normal row. The sequence becomes `seq[offset:] + seq[:offset]`.
Reverse-strand start genes are rotated to their `start` coordinate. **Invariant:** for
`2 ≤ k ≤ L`, `rotate(rotate(r, k), L - k + 2)` returns the original coordinates (tested). The
`k = 1` case is excluded because `L - k + 2 = L + 1` is out of range; it is covered by the
identity rule instead.

## 9. Output specifications

### 9.1 CSV (`{record_id}_genes.csv`)

Columns in order: `id, type, start, end, strand, length_aa, product, category, category_source, phrog`. Strand written as `+`/`-`. UTF-8, comma separated, header row.

### 9.2 Protein FASTA (`{record_id}_proteins.faa`)

```
>{id} {product} [{start}..{end}({strand})] [category={category}]
MSKLV...  (60 characters per line)
```

CDS only. Records with empty `aa` are skipped.

### 9.3 GenBank (`{record_id}_phageviz.gbk`)

Built with Biopython `SeqRecord` and `SeqFeature`. Each CDS has `locus_tag`, `product`, `translation`, `note="category: …"`, and `inference="ab initio prediction:Pyrodigal"` for predicted genes. Topology and rotation offset are added to the `COMMENT`.

### 9.4 Figures

`fig.to_image(format="svg" | "png", width=1600, height=500, scale=2)` via `kaleido`. If kaleido isn't available, fall back to Plotly's built-in camera button and show a note.

## 10. UI state and caching

| Item | Mechanism | Key |
|---|---|---|
| Parsed records from upload | `st.cache_data` | file bytes hash |
| NCBI fetch / search | `st.cache_data` with TTL | normalised accession / term |
| PHROG lookup table | `st.cache_resource` | none (loaded once) |
| Demo records | `st.cache_data` | demo name |
| GC windows | `st.cache_data` | `(record_id, rotation_offset, window)` |
| Selected record index | `st.session_state["record_idx"]` | — |
| Selected gene | `st.session_state["gene_id"]` from `plotly_chart` `on_select` | — |
| NCBI search results | `st.session_state["ncbi_results"]` | — |
| Rotation start gene | `st.session_state["start_gene"]` | — |

When a new genome loads, `gene_id` and `start_gene` are cleared so stale selections don't carry over.

## 11. Edge cases

| Case | Handling |
|---|---|
| CDS crossing the origin | `CompoundLocation`: use `parts`, store unwrapped with `end > L` and `wraps_origin=True` (§4.4), draw two bars |
| Joined CDS (ribosomal frameshift, e.g. lambda tail assembly gene) | Use outer bounds for the map, keep the `translation` qualifier as given |
| Missing `locus_tag` | Generate `{record_id}_{n:04d}` |
| Duplicate IDs | Add `_2`, `_3` suffix |
| Lowercase or soft-masked sequence | Upper-case it |
| Sequence missing from GenBank (CON record) | `gbwithparts` for NCBI; error for uploads |
| GFF3 without `##FASTA` | Error with explanation |
| Windows line endings, BOM | Strip on decode |
| Accession without version | Accepted; NCBI returns the latest version, which is shown |
| Accession for a protein (`YP_`, `WP_`) | Matches regex but esummary on `nuccore` fails → `AccessionNotFound` with hint "this looks like a protein accession" |
| Record is a bacterial chromosome | Blocked by the length check |
| Genome with 0 genes called | Show map with sequence only and a warning |
| Window larger than genome | Clamp to `length // 10` |
| Extremely high GC or AT | No special handling; reported in summary |

## 12. Testing

### 12.1 Fixtures (`tests/data/`)

- Real genomes are read from `data/demo/` (lambda, T7, T4) rather than copied here.
- `make_fixtures.py` regenerates the synthetic files below with fixed seeds. Designed ORFs are
  built from a fixed codon table and their expected protein is stored in `/note="expected_aa=…"`,
  so tests don't depend on a translator to know the right answer.
  - `pharokka_sample.gbk`: Pharokka-style qualifiers, missing and duplicate locus tags, GTG
    starts (complete and 5'-partial), tRNA without product
  - `wraps_origin.gbk`: two circular records with origin-crossing ORFs (+ strand ATG start,
    − strand GTG start), plus a feature ending on the last base
- Malformed inputs (truncated GenBank, binary, empty, PNG, corrupt gzip, CON record, multi-record
  files) are built inside the tests from the files above.
- M4: `sample.gff3` with `##FASTA`, `sample.fasta`, `multi_record.fasta`, GFF3 without FASTA
- M2: `ncbi_esummary.xml`, `ncbi_efetch.gbk` (saved responses for mocking)

### 12.2 Test plan

| Module | Tests |
|---|---|
| parsers | CDS count matches source; coordinates match known genes (lambda gene `A` (terminase large subunit) at 711–2636, + strand); format detection; every malformed file raises `PhageVizError` |
| ncbi | Regex accepts and rejects the expected cases; `fetch_genbank` and `summarize` with mocked `Entrez` (via `monkeypatch`); every error path in §7.6 |
| ncbi (network, opt-in) | Live fetch of `NC_001416` returns a record of 48,502 bp |
| gene_calling | Recall ≥ 90% against lambda RefSeq |
| annotation | Keyword table cases; connector precedence; hallmark detection on lambda |
| metrics | GC of `"GGCC"` = 1.0, skew of `"GGGC"` = 0.5; window count; coding density ≤ 1 |
| transforms | Round-trip rotation invariant; split of origin-crossing features |
| exporters | CSV round-trips through pandas; FASTA parses with `SeqIO`; GenBank re-parses with the same CDS count |
| app | `streamlit.testing.v1.AppTest`: app loads demo without exceptions; NCBI tab shows setup message without secrets |

Network tests use `@pytest.mark.network`, excluded by default in `pytest.ini`.

### 12.3 CI

GitHub Actions on push and pull request: Python 3.10 and 3.12 matrix, `pip install -r requirements.txt`, `pytest --cov=src -m "not network"`.

## 13. Deployment

| Target | Notes |
|---|---|
| Streamlit Community Cloud | Main target. Add `ncbi.email` and `ncbi.api_key` in the app's Secrets settings. About 1 GB RAM. |
| Hugging Face Spaces | Docker SDK, secrets as Space secrets |
| Docker | `python:3.12-slim`, `EXPOSE 8501`, `CMD streamlit run app.py --server.address=0.0.0.0` |

`.streamlit/config.toml`:

```toml
[server]
maxUploadSize = 50

[theme]
base = "light"
```

## 14. Security and privacy

- Uploads are read into memory with `getvalue()` and never written to disk.
- No uploaded content is executed or passed to a shell.
- Logs contain only error types and record IDs, never sequences or filenames.
- The NCBI email and API key live in Streamlit secrets; `secrets.toml` is in `.gitignore`.
- Only accessions and search terms go to NCBI. Uploaded sequences are never sent to any outside service.
- HTML in record titles or product names is escaped before display (no `unsafe_allow_html` on user data).

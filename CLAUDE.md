# CLAUDE.md

Guidance for Claude when working in this repository.

## Project

**PhageViz** is a Streamlit app for viewing bacteriophage genomes. It loads a genome (demo file, upload, or NCBI accession), turns it into one data model, and shows an interactive gene map, GC/skew tracks, a gene table, summary stats and a hallmark-gene check. Results are exported as CSV, FASTA, GenBank and SVG/PNG.

There is **no backend API and no database**. Everything runs in memory in one Streamlit process. The only external service is NCBI E-utilities.

## Source of truth

Read these before starting any feature:

- `docs/PRD.md`: what to build and why (user stories, priorities, milestones)
- `docs/SPEC.md`: how to build it (schemas, module interfaces, formulas, error types, edge cases)

If the code and SPEC disagree, ask before changing either. If you make a design decision the SPEC doesn't cover, add it to SPEC.md in the same change.

## Commands

```bash
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements-dev.txt     # runtime pins + pytest/ruff; installs src/ in editable mode
NCBI_EMAIL=you@example.com python scripts/fetch_demo.py   # one-off: download demo genomes
streamlit run app.py                     # run the app (http://localhost:8501)
pytest -q                                # tests, network tests skipped
pytest -m network                        # live NCBI tests (needs secrets)
pytest --cov=src --cov-report=term-missing
ruff check . && ruff format .            # lint and format
```

## Architecture rules

```
app.py              → UI only: widgets, tabs, st.cache_*, session_state, error display
src/phageviz/       → all logic, pure Python. NEVER import streamlit here.
  models.py         → GenomeRecord dataclass + features schema
  parsers.py        → GenBank, GFF3, FASTA → list[GenomeRecord]
  ncbi.py           → Bio.Entrez search / summary / fetch
  gene_calling.py   → pyrodigal
  annotation.py     → PHROG categories, keyword rules, hallmarks
  metrics.py        → GC, skew, summary stats
  transforms.py     → genome rotation
  exporters.py      → CSV, protein FASTA, GenBank (return bytes)
  plots.py          → Plotly figures (return go.Figure)
  errors.py         → PhageVizError and subclasses
data/demo/          → bundled lambda, T7, T4 GenBank files
data/phrog_annot.tsv
tests/              → mirrors src/ ; fixtures in tests/data/
```

- Every input path ends in a `GenomeRecord`. Plotting and export code must not care where data came from.
- Caching (`@st.cache_data`, `@st.cache_resource`) lives only in `app.py`, wrapping pure functions from `src/`.
- Modules raise subclasses of `PhageVizError` with a user-facing message. `app.py` catches them and calls `st.error`. Never let a traceback reach the UI.

## Critical conventions

**Coordinates: internal is 1-based, inclusive, `start <= end`.**
- Biopython `FeatureLocation` is 0-based half-open → `start + 1`, `end`
- GFF3 and pyrodigal are already 1-based inclusive
- Slicing: `seq[start - 1 : end]`
- Origin-crossing features (circular only): store unwrapped, `end = last base after origin + L`, `wraps_origin=True` (SPEC §4.4). Never store `start > end` anywhere. Use `models.feature_parts()` to get real parts.

Off-by-one errors are the most likely bug in this project. Add a test with a known gene coordinate whenever you touch parsing or rotation (lambda gene `A`: 711–2636, + strand; lambda genome length 48,502 bp).

**Features DataFrame** must match the schema in SPEC §4.2 exactly (column names, dtypes, defaults). Use `models.validate_features(df)` after every parser.

**Categories** are the 10 PHROG categories, spelled exactly as in SPEC §4.3. Priority: Pharokka `function` → PHROG lookup → keyword rules. Always set `category_source`.

**Genetic code** defaults to table 11.

## NCBI

- Use `Bio.Entrez` only, never raw `requests` to eutils.
- Email and API key come from `st.secrets["ncbi"]`, passed into `ncbi.configure()` from `app.py`. Never hardcode them, never print or log them.
- Always call `summarize()` and check length before `fetch_genbank()`.
- `efetch` uses `db="nuccore"`, `rettype="gbwithparts"`, `retmode="text"`.
- Cache raw GenBank text, not parsed objects.
- Unit tests mock `Entrez` with `monkeypatch` and saved responses in `tests/data/`. Live calls only in tests marked `@pytest.mark.network`.

## Code style

- Python 3.10+, type hints on all public functions, `from __future__ import annotations`
- Google-style docstrings on public functions
- `ruff` for lint and format, line length 100
- Prefer vectorised numpy/pandas over Python loops for anything per-base or per-gene
- Plotly: one trace per category, never one trace per gene
- Keep functions small and pure; no global state in `src/`

## Testing

- Write or update tests in the same change as the code.
- Every parser must pass: CDS count matches the source file, a known gene's coordinates are exact, malformed input raises `PhageVizError`.
- Rotation must pass the round-trip invariant (SPEC §8.5).
- `streamlit.testing.v1.AppTest` smoke test: app loads the demo without exceptions.
- Target ≥ 80% coverage on `src/`.
- Run `pytest -q` before saying a task is done, and report the result.

## Workflow

Build in PRD milestone order. Finish and test one milestone before starting the next.

1. **M1:** models, GenBank parser, demo genomes, genome map, gene table, click-to-inspect
2. **M2:** NCBI search and fetch, caching, error handling, secrets
3. **M3:** GC/skew tracks, summary tab, hallmarks, category filter
4. **M4:** FASTA + pyrodigal, GFF3 parser, exports
5. **M5:** rotation, GenBank export, CI, Dockerfile, deployment, README screenshots

For each task:
1. Read the relevant PRD stories and SPEC sections.
2. Propose a short plan if the change touches more than two files.
3. Implement logic in `src/` with tests first, then wire it into `app.py`.
4. Run tests and lint.
5. Update docs if behaviour changed (SPEC for design, README for user-facing features).
6. Commit with a conventional message, e.g. `feat(ncbi): fetch genome by accession`.

## Don't

- Don't import `streamlit` in `src/`.
- Don't write uploaded files to disk or log sequences, filenames or secrets.
- Don't send uploaded sequences to any external service.
- Don't use `unsafe_allow_html=True` with record titles, product names or any user data.
- Don't add a database, auth, or a separate backend; they are out of scope (PRD §2).
- Don't add heavy dependencies (pyhmmer, PHROG HMM databases, external binaries) without asking; they are v2.
- Don't commit `.streamlit/secrets.toml`, `.venv/`, or large data files (> 5 MB).
- Don't silently change coordinates, schema columns or category names; they are contracts.

## Domain notes

- Phage genomes are 10 kb–~700 kb; the app rejects > 1 Mb and < 1 kb.
- Many RefSeq phage genomes are annotated "linear" even when the DNA circularises; don't infer topology.
- Pharokka qualifier names vary by version; parse defensively and fall back to keyword rules.
- pyrodigal uses `meta=True` because single-genome mode needs ≥ 20 kb of training sequence.
- Keyword categories and predicted genes are approximate; the UI must label them as such.

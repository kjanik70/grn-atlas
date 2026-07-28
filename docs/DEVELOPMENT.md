# Development & bootstrap

Practical guide to running, testing, and rebuilding GRN Atlas. For *what the tool does*
and the roadmap see [`ROADMAP.md`](../ROADMAP.md); to add a species see
[`ONBOARDING_SPECIES.md`](./ONBOARDING_SPECIES.md).

## Run it

```bash
# backend (FastAPI) — serves the API + reads backend/data/grn.sqlite3
cd backend && ../venv/bin/python -m uvicorn main:app --host 0.0.0.0 --port 8000
# frontend (Vite dev server)
npm install && npm run dev            # http://localhost:3001
```

The backend needs `backend/data/grn.sqlite3` (gitignored, ~420 MB). Build it from the
committed source caches (see below).

## Test

```bash
venv/bin/python -m pytest backend -q      # backend: unit + DB-invariant + API-contract
npx vitest run                            # frontend
npx vite build                            # production build sanity
npx oxlint src/...                        # lint
```

## Rebuild the database

`build_db.py` deletes and rebuilds `grn.sqlite3` from local caches in `backend/data/`
(no network needed) in ~10 s:

```bash
venv/bin/python backend/scripts/build_db.py
```

Everything the loaders add (sequence context, motif hits, pathways, traits, curated
symbols) is glob-loaded from committed `backend/data/*.json[.gz]` caches, so a fresh
clone reproduces the live DB. Targeted loaders (`load_seqctx.py`, `load_pathways.py`,
`load_traits.py`, `load_curated_symbols.py`) update the live DB in place without a full
rebuild.

## Compute dependencies (only for regenerating derived data)

These are **not** needed to run the app (the derived caches are committed), only to
re-fetch/re-derive them:

- **kallisto** (expression + dsRNA transcript stores). Install a linux binary under
  `tools/kallisto/` (gitignored):
  ```bash
  curl -sL https://github.com/pachterlab/kallisto/releases/download/v0.50.1/kallisto_linux-v0.50.1.tar.gz \
    | tar xz -C tools
  ```
- **BLAST+** (curated petunia symbols via homology; regulator mapping). `tblastn`/
  `makeblastdb` under `BLAST_BIN` (default `/tmp/blastwork/ncbi-blast-2.17.0+/bin`).
- Working files (FASTA, indexes, FASTQ) live under `backend/data/expr/` and `tools/`,
  both gitignored; only the resulting JSON/`.fasta.gz` caches are committed.

Regeneration scripts (all offline-cache-producing): `fetch_seqctx.py`, `motif_scan.py`,
`fetch_expression.py`, `fetch_pathways.py`, `fetch_traits.py`, `fetch_curated_symbols.py`,
`check_source_freshness.py` — driven by `backend/scripts/species_config.py`.

## Data-source currency

`GET /api/v1/provenance/freshness` (backed by `check_source_freshness.py`) reports each
source's loaded vs latest version. See the provenance manifest at `GET /api/v1/provenance`.

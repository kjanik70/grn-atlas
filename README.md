# GRN Atlas

Interactive multi-species **gene-regulatory-network atlas**: explore regulatory networks,
sequence/binding context, expression, pathways, traits, cross-species conservation,
predicted perturbations, and in-silico **dsRNA / RNAi design** — with every predicted or
inferred value clearly labelled distinct from measured data.

React + Cytoscape.js frontend · FastAPI + SQLite backend.

Species: **human, mouse, arabidopsis, tomato, petunia** (dahlia onboarding prepared).
Layers vary by species — see the live matrix at `GET /api/v1/species`.

---

## Quick start (a fresh clone)

Prereqs: **Python 3.10+**, **Node 18+**, `git`, network access (for the source fetch).

```bash
git clone <this-repo> grn-atlas && cd grn-atlas

# 1. Backend deps
python3 -m venv venv
venv/bin/pip install -r backend/requirements.txt

# 2. Fetch the source data (this repo does NOT redistribute third-party data), then build the DB
venv/bin/python backend/scripts/fetch_sources.py --tier light   # pulls sources into backend/data/ (network, minutes)
venv/bin/python backend/scripts/build_db.py                      # writes backend/data/grn.sqlite3 (gitignored)

# 3. Run the API (http://localhost:8000, docs at /docs)
cd backend && ../venv/bin/python -m uvicorn main:app --port 8000

# 4. In another shell: run the UI (http://localhost:3001, proxies /api to :8000)
npm install && npm run dev
```

Or with the Makefile: `make setup && make fetch && make db`, then `make backend` and
(elsewhere) `make frontend`.

> **Data is not committed.** Third-party data (each under its own upstream licence — see
> LICENSE) is fetched from source by `fetch_sources.py`; the ~420 MB SQLite DB is then built
> locally by `build_db.py`. The **core atlas** (genes, interactions, coordinates, orthologs,
> GO, pathways, traits, curated symbols, sequence context) comes from the `core`/`light`
> tiers. The **heavy layers** — expression (kallisto over public RNA-seq) and predicted
> binding (motif scans over multi-GB genomes) — are optional, need kallisto/BLAST+, and take
> much longer; `build_db` loads whatever caches are present, so a clone always yields a
> working atlas and those layers light up once regenerated.

## Tests

```bash
venv/bin/pip install -r backend/requirements-dev.txt
venv/bin/python -m pytest backend -q     # backend: unit + DB-invariant + API-contract
npm run test                             # frontend (vitest)
```

## Regenerating the heavy layers (optional)

The expression and predicted-binding layers are fetched separately because they need the
compute tools (kallisto, BLAST+) and hours of processing over large downloads:
`fetch_expression.py <species>` and `motif_scan.py <species> <genome>` (both driven by
`species_config.py`). See **[docs/DEVELOPMENT.md](docs/DEVELOPMENT.md)** for the bootstrap,
and **[docs/ONBOARDING_SPECIES.md](docs/ONBOARDING_SPECIES.md)** to add a
species.

## Docs
- **[ROADMAP.md](ROADMAP.md)** — living source of truth: capabilities, honest boundaries,
  plan, and a dated iteration log.
- **[docs/DEVELOPMENT.md](docs/DEVELOPMENT.md)** — run, test, rebuild, compute-dep bootstrap.
- **[docs/ONBOARDING_SPECIES.md](docs/ONBOARDING_SPECIES.md)** — add a new species.

## Data provenance & citations
Every integrated source (TRRUST, PlantRegMap, PLAZA, OMA, JASPAR, Plant Reactome, GWAS
Catalog, UniProt, mygene, …) is listed with version + DOI in the machine-readable manifest
at `GET /api/v1/provenance` (BibTeX at `/api/v1/citations.bib`), and a data-currency audit
is at `/api/v1/provenance/freshness`. **Each source keeps its own upstream licence** —
consult the manifest before redistributing derived data.

## Guiding principle
Never present predicted/inferred/curated data as measured. Inferred edges
(`Inferred:Arabidopsis` / `Inferred:Expression`), predicted binding sites (`JASPAR_scan`),
inferred gene labels, and homology-mapped symbols (`UniProt:homology`) are all flagged as
such in the API and UI.

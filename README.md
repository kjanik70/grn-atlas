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
>
> **This `light` build is not the full atlas.** Two core inputs (the measured Arabidopsis
> network + ATRM) are not auto-fetched, and the expression/binding layers are optional and
> tool-heavy. See **[Full data setup, caveats & quality checks](#full-data-setup-caveats--quality-checks)**
> below for the complete, verified setup.

## Tests

```bash
venv/bin/pip install -r backend/requirements-dev.txt
venv/bin/python -m pytest backend -q     # backend: unit + DB-invariant + API-contract
npm run test                             # frontend (vitest)
```

## Full data setup, caveats & quality checks

The `light` quick-start above gives a working atlas, but **not every layer is fully
automatic.** Here is exactly what each step provides, what needs manual work, and how to
confirm the build is complete.

### What each tier provides

| Tier / step | Command | Provides | Auto? | Needs |
|---|---|---|---|---|
| core | `fetch_sources.py --tier core` | genes, **human** network (TRRUST), coords, orthologs, GO | mostly | network |
| light | `fetch_sources.py --tier light` | + pathways, traits, sequence-context windows, curated UniProt symbols | mostly | network; BLAST+ for petunia symbols |
| manual core | *(see below)* | **measured Arabidopsis network** + its tomato/petunia projection | **no** | manual download |
| heavy | `fetch_expression.py`, `motif_scan.py` | expression + predicted binding | **no** | kallisto / BLAST+, hours, GBs |

`build_db.py` glob-loads whatever caches are present and **skips missing inputs gracefully**
(printing `(skip) …`), so a partial fetch always yields a working — if reduced — atlas.

### ⚠️ Caveats (know these before relying on a fresh clone)

1. **Two core inputs are NOT auto-fetched** — their upstreams are unreliable or need
   reshaping, so `fetch_sources.py` only prints guidance:
   - `backend/data/regulation_arabidopsis.tsv` — the **measured Arabidopsis TF→target
     network**. Without it you lose the Arabidopsis edges **and** the inferred tomato/petunia
     edges projected from them (a large share of the plant networks).
   - `backend/data/atrm_regulations.tsv` — ATRM literature-curated direction labels (refine
     Arabidopsis edge signs). Optional; the atlas works without it.
2. **Heavy layers are optional and slow** — expression (kallisto over dozens of public
   RNA-seq runs = hours) and predicted binding (motif scans over multi-GB genomes). A basic
   clone has neither; the dsRNA/expression/binding features light up once regenerated.
3. **`petunia` curated symbols need BLAST+** (`fetch_curated_symbols.py petunia` homology-maps
   real names). Skipped automatically if BLAST+ isn't on `BLAST_BIN`/PATH.
4. **A full from-scratch fetch has not been certified end-to-end** — it hits several live
   sources; expect occasional retries. Each fetcher is the same one that produced the shipped
   data, and graceful degradation is tested, but plan to spot-check (see quality checks).

### Ensuring the manual core files

- **`regulation_arabidopsis.tsv`** — a tab-separated file with **no header**, one edge per
  line, exactly four columns: `TF_locus  target_locus  activation|repression  confidence`
  (AGI ids, e.g. `AT1G01060`; confidence 0–1). Produce it from PlantRegMap's Arabidopsis
  regulation data (https://plantregmap.gao-lab.org/, TF→target / FunTFBS) reduced to those
  four columns.
- **`atrm_regulations.tsv`** *(optional)* — tab-separated **with a header row**; ≥5 columns
  where col 1 = TF locus, col 2 = target locus, col 5 = direction label `A` / `R` / `D`.
  Source: ATRM (http://atrm.cbi.pku.edu.cn/). Skip if unavailable.

Place both in `backend/data/`, then re-run `build_db.py`.

### Regenerating the heavy layers (optional)

Install the compute tools once (see **[docs/DEVELOPMENT.md](docs/DEVELOPMENT.md)** for the
kallisto/BLAST+ bootstrap), then per plant species:

```bash
# expression (kallisto index of the species CDS + a curated RNA-seq panel in species_config)
EXPR_SUBSAMPLE=3000000 venv/bin/python backend/scripts/fetch_expression.py petunia
# predicted binding (JASPAR-plant PWM scan over the genome)
venv/bin/python backend/scripts/motif_scan.py petunia /path/to/genome.fa
venv/bin/python backend/scripts/load_seqctx.py petunia         # load motifs+hits into the DB
```

Both are driven by `backend/scripts/species_config.py` (assembly, URLs, RNA-seq panel).
To add a species see **[docs/ONBOARDING_SPECIES.md](docs/ONBOARDING_SPECIES.md)**.

### Quality checks — confirm the build is complete & correct

```bash
# 1. Referential-integrity + sanity invariants over the built DB
venv/bin/python -m pytest backend/tests/test_db_invariants.py -q

# 2. Per-species layer coverage (network / orthologs / binding / expression / pathways / traits)
curl -s localhost:8000/api/v1/species | python3 -m json.tool

# 3. Source-currency audit (loaded vs latest upstream version)
curl -s localhost:8000/api/v1/provenance/freshness | python3 -m json.tool
```

A **complete** build (all tiers + manual core + heavy layers) should report roughly:
`~50,800` genes · human `~4,859` edges · arabidopsis `~91,844` edges · tomato measured
`12,719` + inferred `~197,618`. `build_db.py`'s own summary prints these counts — compare
them, and use `/api/v1/species` to see which layers are populated vs empty. If a layer is
unexpectedly empty, its source file wasn't fetched (check the `(skip)` lines from `build_db`).

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

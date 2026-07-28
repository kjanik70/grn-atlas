# GRN Atlas

Interactive multi-species gene-regulatory-network atlas: explore regulatory networks,
sequence/binding context, expression, pathways, traits, cross-species conservation,
predicted perturbations, and in-silico **dsRNA / RNAi design** — with every predicted
or inferred value clearly labelled distinct from measured data.

React + Cytoscape.js frontend, FastAPI + SQLite backend.

## Species & data layers
Human, mouse, arabidopsis, tomato, petunia (dahlia onboarding prepared). Layers vary by
species — see the live matrix at `GET /api/v1/species`. Sources (TRRUST, PlantRegMap,
PLAZA, OMA, JASPAR, Plant Reactome, GWAS Catalog, UniProt, …) are cited in the provenance
manifest (`GET /api/v1/provenance`, BibTeX at `/api/v1/citations.bib`).

## Docs
- **[ROADMAP.md](../ROADMAP.md)** — the living source of truth: capabilities, honest
  boundaries, the plan, and a dated iteration log.
- **[DEVELOPMENT.md](./DEVELOPMENT.md)** — run, test, rebuild the DB, compute-dep bootstrap.
- **[ONBOARDING_SPECIES.md](./ONBOARDING_SPECIES.md)** — add a new species (e.g. Dahlia).

## Quick start
```bash
cd backend && ../venv/bin/python -m uvicorn main:app --port 8000   # API (needs grn.sqlite3)
npm install && npm run dev                                          # UI on :3001
```
(If `backend/data/grn.sqlite3` is missing, build it: `venv/bin/python backend/scripts/build_db.py`.)

## Principle
Never present predicted/inferred/curated data as measured. Inferred edges
(`Inferred:Arabidopsis` / `Inferred:Expression`), predicted binding sites (`JASPAR_scan`),
inferred labels, and homology-mapped symbols are all flagged as such.

# Onboarding a new species (e.g. Dahlia)

This runbook makes adding a species turnkey. It is the concrete path we're preparing for
the **Dahlia** collaboration (Alex/Zach — incoming genome + RNA-seq + ~400-cultivar GWAS;
see the `grn-atlas-dahlia-collaboration` memory). Every added layer follows the project
rule: **predicted/inferred data is labelled distinct from measured** (see the honesty memory).

Central config: `backend/scripts/species_config.py` — add one entry per species
(assembly, reference URLs, `id_style`). A `dahlia` placeholder is already stubbed there.

Check readiness / progress any time: `GET /api/v1/species` shows the per-species layer
matrix (network / orthologs / binding / expression / pathways / traits); empty layers are
the remaining onboarding steps.

## Steps (each independent; do what the available data supports)

1. **Genes + genome/annotation.** Load gene models (id, symbol, description, coords) from
   the release's GFF into `genes` + `gene_locations`. Set the assembly tag.

2. **Sequence context + predicted TF binding** (needs genome FASTA + GFF):
   - promoter/gene-body windows: adapt `fetch_arabidopsis_seqctx.py` (identity crosswalk +
     windows) or the generic path; load with a `load_*_seqctx.py`.
   - JASPAR-plant PWM scan over target promoters: adapt `scan_arabidopsis_motifs.py`
     (symbol/`tf_motif_map` mapping) → `motif_hits` (`tier=JASPAR_scan`).
   - Register the assembly in `_ASSEMBLY_OF` (backend/main.py) → motif enrichment +
     sequence-context export light up.

3. **Expression + co-expression** (needs a CDS FASTA + a curated SRA panel):
   - build a kallisto index from the CDS; curate a tissue panel (aim for contrast:
     leaf/root/flower/fruit…, prefer same-study contrasts to limit batch effects).
   - copy `fetch_tomato_expression.py` (it reuses the shared ENA-streaming + kallisto
     helpers in `fetch_petunia_expression.py`) → `expression_<species>.json.gz`.
   - No backend change: `expression.get_matrix(species)` + the `/expression` and
     `/coexpression` endpoints auto-select by gene species.

4. **Orthologs** (enables cross-species conservation + Arabidopsis-network projection):
   load ortholog pairs (OMA/PLAZA BHIF, or genus phylogeny/bait-capture for Dahlia) into
   `orthologs`. Projecting the Arabidopsis network onto the new species yields
   `Inferred:Arabidopsis` edges (labelled inferred).

5. **Pathways** (optional): if the species is in Plant Reactome, extend `fetch_pathways.py`
   `SPECIES`; else project via orthologs.

6. **Traits / GWAS** (species-agnostic, ready now): produce a TSV with `gene_id` + `trait`
   (+ optional `pmid`, `source`) and run:
   ```
   python backend/scripts/ingest_trait_table.py dahlia zach_gwas.tsv --source "Dahlia GWAS (Zach 20XX)"
   ```
   Gene ids must match atlas ids for the species; non-matching rows are skipped and
   reported. `/trait_enrichment` then works for that species automatically.

## Dahlia specifics (fill in when data lands)
- Genome/annotation: from the G3 paper / NCBI BioProject (bioRxiv link pending from Alex).
- Expression: the paper's RNA-seq on SRA → step 3.
- GWAS: Zach's ~400 low-pass cultivar runs → gene-mapped trait table → step 6.
- Phylogeny: genus bait-capture (~1000 CompCos genes) → orthology context (step 4).
- Prioritise the **anthocyanin/floral-pigmentation** pathway (the shared scientific interest).

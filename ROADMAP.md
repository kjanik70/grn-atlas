# GRN Atlas — Capabilities & Roadmap (living document)

> This is the single, continuously-updated map of what GRN Atlas can do, where it
> falls short, and what we're building next. **Update it every iteration:** when a
> capability ships, move it up; when a gap closes, strike it; append to the
> Iteration Log. The working loop is: _build → test → document here → find gaps →
> plan → repeat._

Last updated: 2026-07-26 · Baseline: backend 72 tests, frontend 5 tests, build clean.

---

## 1. What we can do today (capabilities → problems solved)

### Regulatory structure
- **Who regulates gene X / what does TF Y target** — neighborhood + network view,
  filterable by direction, confidence, evidence source, and measured-vs-inferred;
  PubMed links per edge.
- **Path from A → B** — pathfinding (BFS over signed edges).
- **Master regulators / core circuit** — Organism tab (hub ranking, core-circuit graph).

### Gene-set interpretation
- **Functional theme of a set** — GO enrichment (hypergeometric + BH FDR).
- **Which TFs drive a set** — motif enrichment over the scanned-promoter background,
  now for tomato, petunia, and **arabidopsis** (e.g. AN2 top-enriched in petunia
  flavonoid promoters q=4e-5; BPC1 targets enrich the BPC1 motif in arabidopsis q=1.9e-83).
- **Actionable coordinates** — export of signed edges + confidence + genomic coords
  + promoter windows + predicted binding sites (tomato/petunia).

### Expression & co-expression (petunia)
- **Per-tissue expression profile** — `GET /api/v1/expression/{gene_id}`: TPM across a
  29-sample P. axillaris RNA-seq panel (vegetative + floral/pigmentation tissues),
  quantified with kallisto vs PLAZA pax CDS (Peaxi162 IDs). Predicted/shallow (subsampled).
- **Predicted co-expression** — `POST /api/v1/coexpression`: Pearson on log2(TPM+1),
  labeled `Inferred:Expression`, undirected (not causal, not measured regulation).
  `tf_only` restricts partners to candidate TF regulators. Shown in the gene detail panel.
- Verified: AN2 peaks in flower/corolla/petal-limb; co-expresses with petal-identity
  genes (PI/AP3); 3 of AN2's 18 network targets are independently co-expressed —
  expression corroborating a subset of projected edges.

### Perturbation prediction
- **Predict downstream effects of a TF knockout/over-expression** — `POST /api/v1/perturb`
  propagates signs along the network (activation +1 / repression −1 × intervention sign),
  confidence-weighted and depth-damped. Returns predicted up/down/unknown per reachable
  gene with the path. Unsigned edges → "unknown"; inferred routes flagged. Drives the
  Intervention Designer. Verified: AN2 knockout → its anthocyanin targets predicted down.

### Comparative / evolutionary
- **Edge conservation across species** — `/api/v1/conservation` joins orthologs +
  both networks (JAF13→CHS conserved petunia↔tomato; AN2→CHS diverged).
- **Synteny / orthologs** — genome view, ideograms, ribbons.

### Practical / trust
- **Gene identification despite messy annotation** — BLAST regulator ID + synonym search.
- **Cite & reproduce** — provenance manifest + BibTeX + versioned methods in exports.
- **Share** — permalinks; exports (PNG/SVG/GraphML/JSON/CSV/TSV); collections.
- **Teach** — student scaffolding (examples, glossary, inferred-edge explainers).

### Species & data
- 5 species: human, mouse, arabidopsis, tomato, petunia.
- Measured edges: TRRUST v2 (human), PlantRegMap FunTFBS (tomato, arabidopsis).
- Inferred edges: Arabidopsis network projected onto tomato/petunia (`Inferred:Arabidopsis`).
- Sequence context: promoter/gene-body windows + JASPAR-scan binding sites (tomato/petunia).
- Orthology/coords: OMA, PLAZA Dicots 4.5, DNA Zoo *P. axillaris* Hi-C.

---

## 2. Honest boundaries (what bounds rigor today)

- **Expression exists for petunia only, and is shallow** — 29 subsampled samples, one
  species; use for relative/co-expression signal, not absolute levels. Other species
  still have no expression axis.
- **No dynamics** — the cascade/intervention view is a toy, not a quantitative model.
- **Petunia edges are inferred** — hypotheses, not evidence; no measured petunia GRN.
- **No accessibility (ATAC), PPI/complexes, or phenotype/QTL linkage.**
- **Data currency** — PLAZA 4.5 is 2018; TRRUST v2 is older.

**Guiding principle:** never fabricate scientific data. Inferred/predicted/curated
data must always be labeled distinct from measured (inferred edges, JASPAR_scan
sites, BLAST_curated symbols). New inferred layers inherit this rule.

---

## 3. Expansion plan — items 1–6 (value-ranked)

Numbered by research value. **Execution order differs** (feasibility-first): the
data-free item ships first, then the expression linchpin, then the rest.

| # | Feature | Unlocks | New data needed | Effort |
|---|---------|---------|-----------------|--------|
| 1 | **Expression integration** | condition/tissue-specificity, co-expression, expression-weighted TFBS | RNA-seq atlases | High |
| 2 | **Network inference from expression** (GENIE3/ARACNe-style) | a data-driven petunia network (not just projected) | uses #1 data | High |
| 3 | **Perturbation modeling upgrade** | "predict effect of knocking out AN2" via signed-path propagation | none | Medium |
| 4 | **Base-resolution binding for human/arabidopsis** | motif enrichment + seq context for best-measured species | JASPAR/ReMap ChIP | Medium |
| 5 | **Broader enrichment + trait linkage** | KEGG/MapMan pathways; QTL/GWAS→gene | KEGG, GWAS/QTL | Medium |
| 6 | **Taxonomic scope / freshness** | wheat/cotton homoeologs; refreshed PLAZA/TRRUST | newer releases | Medium |

### Execution sequence (feasibility-first)
1. **#3 perturbation** — no external data; uses existing signed edges. *Start here.*
2. **#1 + #2 expression + inference** — the linchpin; do together (inference consumes
   the loaded expression). Removes the static + petunia-only-inferred ceiling.
3. **#4 base-resolution binding** — extends the strongest existing analysis to human/arabidopsis.
4. **#5 enrichment + trait linkage.**
5. **#6 scope/freshness.**

### Per-item design sketch
- **#3:** replace toy cascade with signed-path propagation over `interactions`
  (sign = activation/repression product along path, with sign-flip on repression;
  damped by confidence & depth). Endpoint `POST /api/v1/perturb` (KO/OE list →
  predicted up/down per reachable gene, with path evidence). Honest labeling:
  output is "predicted direction," qualitative, gated by inferred-edge inclusion.
- **#1:** new `expression(gene_id, species, sample, tpm)` + `samples(...)` tables from
  public RNA-seq (tomato PRJNA980935 / expression atlases; Arabidopsis; petunia if
  available). Endpoints: per-gene expression profile; expression-weighted edge/motif
  views. Cache to committed JSON like other fetchers; runtime stays offline.
- **#2:** `infer_network.py` — GENIE3-style tree importance (or correlation fallback)
  per species from #1's matrix → `Inferred:Expression` edges, labeled distinct from
  both measured and Arabidopsis-projected. Compare against projection in the UI.
- **#4:** extend `motifs`/`motif_hits` + crosswalk to human/arabidopsis using
  JASPAR/ReMap; reuse existing motif-enrichment + seq-context machinery.
- **#5:** `go_annotations`-style tables for KEGG/MapMan; QTL/GWAS→gene mapping table;
  reuse hypergeometric+BH; new "trait" lookup endpoint.
- **#6:** refresh fetchers to current PLAZA/TRRUST; add wheat/cotton with homoeolog
  handling.

---

## 4. Known gaps / backlog (revisit each iteration)
- ~~Perturbation model is a toy~~ ✅ shipped (#3): signed-path propagation, honest unknown/inferred labels.
- ~~Static network — no time/condition axis~~ ✅ shipped (#1) for petunia: 29-sample expression profiles.
- ~~Petunia has no data-derived network~~ ✅ shipped (#2): co-expression inference (`Inferred:Expression`).
  Follow-ups: extend expression to tomato/arabidopsis; deepen sampling; upgrade correlation → tree-based (GENIE3).
- ~~Sequence layer absent for arabidopsis~~ ✅ shipped (#4, plant side): TAIR10 JASPAR scan (95k sites).
  Human base-resolution binding (ReMap/JASPAR vertebrate) still pending — larger genome + peak ingest.
- Only GO enrichment; no pathway/trait ontologies (addressed by #5).
- Stale releases; narrow taxon set (addressed by #6).
- Older docs (START_HERE/PROJECT_SUMMARY etc.) predate recent features — consolidate.

## 5. Iteration Log
- **2026-07-26** — Established this roadmap + baseline (61 backend / 5 frontend green).
  Prior shipped: provenance/citations, cross-species conservation, motif enrichment.
- **2026-07-26** — Shipped **#3 perturbation**: `/api/v1/perturb` signed-path propagation
  replacing the toy cascade; rewired the Intervention Designer to it. +3 API tests
  (64 backend / 5 frontend green). Verified AN2 KO → anthocyanin targets down.
  Next: **#1 + #2** — expression integration + network inference (the linchpin). This
  needs external RNA-seq; first step is a data-availability check + a fetch script
  (cache to committed JSON like other fetchers, runtime stays offline).
- **2026-07-26** — Shipped **#1 + #2 (petunia)**. Reference unblocked via PLAZA pax CDS
  (SGN was down). Built `fetch_petunia_expression.py`: streams subsampled reads for a
  curated 29-sample panel from ENA, kallisto-quantifies to Peaxi162 TPM. New
  `expression.py` + `/api/v1/expression/{id}` (#1) and `/api/v1/coexpression` (#2,
  `Inferred:Expression`, undirected). `ExpressionPanel` in the gene detail view. +7 tests
  (71 backend / 5 frontend). Verified: AN2 pigmented-tissue-specific; co-expresses PI/AP3;
  corroborates 3/18 projected AN2 targets.
- **2026-07-26** — Shipped **#4 (plant side)**: Arabidopsis base-resolution TF binding.
  Reused the plant seqctx/scan machinery — PLAZA ath GFF → TAIR10 promoter windows;
  JASPAR-plant PWM scan (symbol-mapped, 346 TFs) → 95,132 predicted sites; targeted DB
  loader (no full rebuild); enabled arabidopsis in `_ASSEMBLY_OF` so motif enrichment +
  sequence-context now work for it. +1 DB-invariant test (72 backend). Verified: BPC1
  targets enrich BPC1 motif (q=1.9e-83), paralogs BPC5/6 co-enrich.
  Next: **human** base-resolution binding (ReMap/JASPAR vertebrate), then **#5**
  (KEGG/MapMan + trait linkage); also extend expression to tomato/arabidopsis.

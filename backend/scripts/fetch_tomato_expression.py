"""Build a tomato (S. lycopersicum) expression matrix keyed to Solyc gene IDs.

Same pipeline as petunia (reuses its ENA-streaming + kallisto helpers): PLAZA Dicots
4.5 sly CDS as the reference (headers carry the versioned Solyc id, matching the atlas),
a curated public RNA-seq panel spanning tissues (leaf, root, stem, flower, bud, fruit,
apex, cotyledon), subsampled and pseudo-aligned to per-gene TPM.

Shallow/predicted by design — for relative + co-expression signal, not absolute levels.
Output: backend/data/expression_tomato.json.gz
"""
import gzip
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fetch_petunia_expression as P  # noqa: E402  (reuse ENA + kallisto helpers)

DATA = HERE.parent / "data"
IDX = DATA / "expr" / "sly.idx"
OUT = DATA / "expression_tomato.json.gz"

# Curated panel: run -> (tissue, study). Spans distinct organs for co-expression contrast.
PANEL = {
    "DRR016684": ("leaf", "PRJDB3892"),
    "DRR016687": ("leaf", "PRJDB3892"),
    "DRR016686": ("fruit", "PRJDB3892"),
    "DRR177588": ("fruit_green", "PRJDB8570"),
    "DRR256501": ("fruit_8mm", "PRJDB10790"),
    "DRR092919": ("root", "PRJDB5790"),
    "DRR128403": ("root", "PRJDB8390"),
    "DRR092901": ("stem", "PRJDB5790"),
    "DRR092918": ("stem", "PRJDB5790"),
    "DRR092898": ("flower_closed", "PRJDB5790"),
    "DRR111122": ("flower_open", "PRJDB7574"),
    "DRR092914": ("bud_3_4mm", "PRJDB5790"),
    "DRR092915": ("bud_2mm", "PRJDB5790"),
    "DRR1004608": ("apex", "PRJDB11748"),
    "DRR1004609": ("apex", "PRJDB11748"),
    "DRR271999": ("cotyledon", "PRJDB11160"),
    "DRR272000": ("cotyledon", "PRJDB11160"),
    "DRR256503": ("ovary_anthesis", "PRJDB10790"),
    "DRR256504": ("ovary_anthesis", "PRJDB10790"),
    "DRR092900": ("leaf_mature", "PRJDB5790"),
}


def main():
    if not IDX.exists():
        sys.exit(f"missing kallisto index {IDX}")
    runs = list(PANEL)
    per_sample = {}
    for i, run in enumerate(runs, 1):
        print(f"[{i}/{len(runs)}] {run} ({PANEL[run][0]})", flush=True)
        ab = P.quant_run(run, IDX)
        if ab is None:
            continue
        gene_tpm = {}
        with open(ab) as fh:
            next(fh)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                g = P.gene_of(parts[0])
                gene_tpm[g] = gene_tpm.get(g, 0.0) + float(parts[4])
        per_sample[run] = gene_tpm

    ok = [r for r in runs if r in per_sample]
    genes = sorted({g for tpm in per_sample.values() for g in tpm})
    matrix = {g: [round(per_sample[r].get(g, 0.0), 3) for r in ok] for g in genes}
    out = {
        "meta": {"assembly": "ITAG2.3 (Solyc IDs)",
                 "reference": "PLAZA Dicots 4.5 sly CDS (selected transcript)",
                 "quantifier": "kallisto 0.50.1", "source": "ENA/SRA public RNA-seq, subsampled",
                 "n_reads_subsample": P.N_READS, "unit": "TPM",
                 "label": "predicted (shallow subsampled quantification)",
                 "note": "Proof-of-concept expression; relative/co-expression signal, not absolute."},
        "samples": [{"run": r, "tissue": PANEL[r][0], "study": PANEL[r][1]} for r in ok],
        "genes": matrix,
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(out, fh)
    print(f"\nwrote {OUT}  ({len(ok)} samples, {len(genes)} genes)", flush=True)


if __name__ == "__main__":
    main()

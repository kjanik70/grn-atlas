"""Build an Arabidopsis expression matrix keyed to AGI (AT…) gene IDs.

Same pipeline as petunia/tomato (reuses the shared ENA + kallisto helpers): PLAZA
Dicots 4.5 ath CDS reference (headers ">AT1G01010.1 | AT1G01010" -> gene_of gives the
AGI, matching the atlas), a curated public RNA-seq panel spanning vegetative shoot,
inflorescence (flower), root, and seedling.

Shallow/predicted (subsampled) — relative + co-expression signal, not absolute levels.
Output: backend/data/expression_arabidopsis.json.gz
"""
import gzip
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fetch_petunia_expression as P  # noqa: E402

DATA = HERE.parent / "data"
IDX = DATA / "expr" / "ath.idx"
OUT = DATA / "expression_arabidopsis.json.gz"

# Curated panel. PRJDB3784 gives a same-study vegetative-shoot vs inflorescence contrast
# (minimal batch effect); roots + seedlings add organ diversity.
PANEL = {
    "DRR031752": ("vegetative_shoot", "PRJDB3784"),
    "DRR031753": ("vegetative_shoot", "PRJDB3784"),
    "DRR031754": ("vegetative_shoot", "PRJDB3784"),
    "DRR031755": ("vegetative_shoot", "PRJDB3784"),
    "DRR031756": ("vegetative_shoot", "PRJDB3784"),
    "DRR031757": ("vegetative_shoot", "PRJDB3784"),
    "DRR031758": ("inflorescence", "PRJDB3784"),
    "DRR031759": ("inflorescence", "PRJDB3784"),
    "DRR031760": ("inflorescence", "PRJDB3784"),
    "DRR031761": ("inflorescence", "PRJDB3784"),
    "DRR031762": ("inflorescence", "PRJDB3784"),
    "DRR031763": ("inflorescence", "PRJDB3784"),
    "DRR016112": ("root", "PRJDB1593"),
    "DRR016113": ("root", "PRJDB1593"),
    "DRR070501": ("root", "PRJDB5141"),
    "DRR032000": ("seedling", "PRJDB3217"),
    "DRR032003": ("seedling", "PRJDB3217"),
    "DRR032004": ("seedling", "PRJDB3217"),
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
        "meta": {"assembly": "TAIR10", "reference": "PLAZA Dicots 4.5 ath CDS (selected transcript)",
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

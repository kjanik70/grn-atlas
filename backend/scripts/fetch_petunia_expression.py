"""Build a Petunia axillaris expression matrix keyed to Peaxi162 gene IDs.

Petunia has no ready-made Peaxi162-keyed expression atlas, so we build one:
  1. reference   = PLAZA Dicots 4.5 pax CDS (headers are Peaxi162 gene IDs)
  2. reads       = a curated panel of public P. axillaris RNA-seq runs (ENA/SRA),
                   spanning vegetative + floral/pigmentation tissues
  3. quantify    = kallisto pseudo-alignment -> TPM per transcript -> per gene

Reads are *subsampled* (first N reads streamed from ENA) to keep this tractable and
offline-friendly; the resulting matrix is a shallow, honest proof-of-concept — good
for co-expression signal on well-expressed genes (e.g. the anthocyanin pathway), not
a definitive atlas. Output is cached to a committed JSON so runtime stays offline.

Output: backend/data/expression_petunia.json.gz
  {meta, samples:[{run,tissue,study,n_reads}], genes:{gene_id:[tpm per sample]}}

Re-runnable: per-run kallisto output is cached under data/expr/quant/<run>/ and skipped
if already present.
"""
import gzip
import json
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
EXPR = DATA / "expr"
QUANT = EXPR / "quant"
IDX = EXPR / "pax.idx"
KALLISTO = HERE.parent.parent / "tools" / "kallisto" / "kallisto"
OUT = DATA / "expression_petunia.json.gz"

N_READS = int(os.environ.get("EXPR_SUBSAMPLE", "3000000"))
ASSEMBLY = "Peaxi162v1.6.2"

# Curated panel: run -> (tissue label, study). Chosen to span vegetative tissues
# (low pigment) and floral/pigmentation tissues (high pigment) for co-expression.
PANEL = {
    # PRJNA261953 tissue atlas (paired)
    "SRR1585615": ("apical_shoot", "PRJNA261953"),
    "SRR1585635": ("flower", "PRJNA261953"),
    "SRR1585830": ("seedling", "PRJNA261953"),
    "SRR1585954": ("callus", "PRJNA261953"),
    "SRR1585955": ("trichome", "PRJNA261953"),
    # PRJNA524676 corolla spatial (paired) -- pigmented petal tissue
    "SRR8644905": ("corolla_lobes", "PRJNA524676"),
    "SRR8644913": ("corolla_lobes", "PRJNA524676"),
    "SRR8644915": ("corolla_lobes", "PRJNA524676"),
    "SRR8644906": ("corolla_tube", "PRJNA524676"),
    "SRR8644910": ("corolla_tube", "PRJNA524676"),
    "SRR8644908": ("corolla_tube", "PRJNA524676"),
    "SRR8644907": ("corolla_lobes_tz", "PRJNA524676"),
    "SRR8644904": ("corolla_tube_tz", "PRJNA524676"),
    # PRJNA674380 petal limb (single) -- strongly pigmented
    "SRR12998769": ("petal_limb", "PRJNA674380"),
    "SRR12998762": ("petal_limb", "PRJNA674380"),
    "SRR12998768": ("petal_limb", "PRJNA674380"),
    "SRR12998755": ("petal_limb", "PRJNA674380"),
    "SRR12998767": ("petal_limb", "PRJNA674380"),
    "SRR12998766": ("petal_limb", "PRJNA674380"),
    # PRJNA1267051 bud development series (single)
    "SRR33679195": ("bud_9wk", "PRJNA1267051"),
    "SRR33679192": ("bud_12wk", "PRJNA1267051"),
    "SRR33679191": ("bud_18wk", "PRJNA1267051"),
    "SRR33679194": ("large_bud_12wk", "PRJNA1267051"),
    "SRR33679193": ("larger_bud_18wk", "PRJNA1267051"),
    "SRR33679190": ("young_bud_18wk", "PRJNA1267051"),
    # PRJNA533335 style length series (single) -- non-petal floral organ
    "SRR8930520": ("style_small", "PRJNA533335"),
    "SRR8930518": ("style_long", "PRJNA533335"),
    "SRR8930523": ("style_medium", "PRJNA533335"),
    "SRR8930526": ("style_medium", "PRJNA533335"),
}

UA = {"User-Agent": "grn-atlas/1.0"}


def ena_fastq_urls(run):
    """Return list of https fastq URLs for a run via the ENA filereport API."""
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + run +
           "&result=read_run&fields=fastq_ftp&format=tsv")
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        lines = r.read().decode().splitlines()
    if len(lines) < 2 or not lines[1].strip():
        return []
    ftp = lines[1].split("\t")[-1]
    return ["https://" + p for p in ftp.split(";") if p]


def stream_subsample(url, dest, n_reads):
    """Stream first n_reads reads (4*n lines) of a gzipped fastq from ENA to dest."""
    n_lines = n_reads * 4
    # curl | zcat | head : curl stops early on SIGPIPE once head has enough.
    cmd = f"curl -s --max-time 1200 '{url}' | zcat 2>/dev/null | head -n {n_lines} > '{dest}'"
    subprocess.run(["bash", "-c", cmd], check=False)
    return dest.exists() and dest.stat().st_size > 0


def quant_run(run):
    """Download (subsampled) + kallisto quant one run. Returns abundance.tsv path."""
    outdir = QUANT / run
    ab = outdir / "abundance.tsv"
    if ab.exists():
        return ab
    outdir.mkdir(parents=True, exist_ok=True)
    urls = ena_fastq_urls(run)
    if not urls:
        print(f"  [{run}] no fastq urls", flush=True)
        return None

    fastqs = []
    paired = len(urls) >= 2 and any(u.endswith("_1.fastq.gz") for u in urls)
    try:
        if paired:
            u1 = next(u for u in urls if u.endswith("_1.fastq.gz"))
            u2 = next(u for u in urls if u.endswith("_2.fastq.gz"))
            f1, f2 = EXPR / f"{run}_1.fq", EXPR / f"{run}_2.fq"
            if not (stream_subsample(u1, f1, N_READS) and stream_subsample(u2, f2, N_READS)):
                print(f"  [{run}] download failed", flush=True)
                return None
            fastqs = [str(f1), str(f2)]
            kcmd = [str(KALLISTO), "quant", "-i", str(IDX), "-o", str(outdir),
                    "-t", "2"] + fastqs
        else:
            u = next((u for u in urls if u.endswith(".fastq.gz")), urls[0])
            f = EXPR / f"{run}.fq"
            if not stream_subsample(u, f, N_READS):
                print(f"  [{run}] download failed", flush=True)
                return None
            fastqs = [str(f)]
            kcmd = [str(KALLISTO), "quant", "-i", str(IDX), "-o", str(outdir),
                    "-t", "2", "--single", "-l", "200", "-s", "20"] + fastqs
        r = subprocess.run(kcmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  [{run}] kallisto failed: {r.stderr[-300:]}", flush=True)
            return None
        return ab if ab.exists() else None
    finally:
        for f in fastqs:
            try:
                os.remove(f)
            except OSError:
                pass


def gene_of(target_id):
    """Peaxi162Scf00000g00013.1 -> Peaxi162Scf00000g00013"""
    return target_id.rsplit(".", 1)[0]


def main():
    if not IDX.exists():
        sys.exit(f"missing kallisto index {IDX} (run kallisto index first)")

    runs = list(PANEL)
    per_sample_tpm = {}   # run -> {gene: tpm}
    for i, run in enumerate(runs, 1):
        print(f"[{i}/{len(runs)}] {run} ({PANEL[run][0]})", flush=True)
        ab = quant_run(run)
        if ab is None:
            continue
        gene_tpm = {}
        with open(ab) as fh:
            next(fh)  # header: target_id length eff_length est_counts tpm
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                g = gene_of(parts[0])
                gene_tpm[g] = gene_tpm.get(g, 0.0) + float(parts[4])
        per_sample_tpm[run] = gene_tpm

    ok_runs = [r for r in runs if r in per_sample_tpm]
    all_genes = sorted({g for tpm in per_sample_tpm.values() for g in tpm})
    matrix = {g: [round(per_sample_tpm[r].get(g, 0.0), 3) for r in ok_runs] for g in all_genes}

    out = {
        "meta": {
            "assembly": ASSEMBLY,
            "reference": "PLAZA Dicots 4.5 pax CDS (selected transcript)",
            "quantifier": "kallisto 0.50.1",
            "source": "ENA/SRA public RNA-seq, subsampled",
            "n_reads_subsample": N_READS,
            "unit": "TPM",
            "label": "predicted (shallow subsampled quantification)",
            "note": "Proof-of-concept expression built from subsampled public reads; "
                    "use for relative/co-expression signal, not absolute quantification.",
        },
        "samples": [{"run": r, "tissue": PANEL[r][0], "study": PANEL[r][1]} for r in ok_runs],
        "genes": matrix,
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(out, fh)
    print(f"\nwrote {OUT}  ({len(ok_runs)} samples, {len(all_genes)} genes)", flush=True)


if __name__ == "__main__":
    main()

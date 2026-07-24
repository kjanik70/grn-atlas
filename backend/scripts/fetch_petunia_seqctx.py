"""
Petunia sequence-context windows (Path B, generalized to petunia).

Petunia's atlas gene IDs are already its assembly IDs (Peaxi162…), so no version
crosswalk is needed — just an identity crosswalk. Windows are computed on the
ORIGINAL P. axillaris v1.6.2 scaffolds (the coordinate system of the scaffold
FASTA, so promoter sequences are directly fetchable), NOT the atlas's lifted
1–7 chromosome display coords. Assembly is tagged accordingly.

Outputs (append alongside tomato's; loaded by build_db):
  - gene_id_crosswalk_petunia.json.gz : identity, assembly Peaxi162v1.6.2
  - gene_windows_petunia.json.gz      : promoter + gene_body on scaffolds (BED0)

Usage: python backend/scripts/fetch_petunia_seqctx.py
"""
import gzip
import json
import sqlite3
import urllib.request
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
DB_PATH = DATA_DIR / "grn.sqlite3"
CROSSWALK = DATA_DIR / "gene_id_crosswalk_petunia.json.gz"
WINDOWS = DATA_DIR / "gene_windows_petunia.json.gz"

PAX_GFF = ("https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04_5/GFF/pax/"
           "annotation.selected_transcript.all_features.pax.gff3.gz")
ASSEMBLY = "Peaxi162v1.6.2"
PROMOTER_UP, PROMOTER_DOWN = 2000, 500
UA = {"User-Agent": "grn-atlas-build/1.0"}


def atlas_petunia():
    conn = sqlite3.connect(DB_PATH)
    ids = {r[0] for r in conn.execute("SELECT id FROM genes WHERE species='petunia'")}
    conn.close()
    return ids


def promoter(start, end, strand):
    if strand >= 0:
        tss0 = start - 1
        return max(0, tss0 - PROMOTER_UP), tss0 + PROMOTER_DOWN
    return max(0, end - PROMOTER_DOWN), end + PROMOTER_UP


def main():
    keep = atlas_petunia()
    print(f"Atlas petunia genes: {len(keep)}")
    print("Downloading PLAZA pax GFF…")
    req = urllib.request.Request(PAX_GFF, headers=UA)
    with urllib.request.urlopen(req, timeout=180) as resp, gzip.GzipFile(fileobj=resp) as gz:
        lines = gz.read().decode("utf-8", "replace").splitlines()

    crosswalk, windows = [], []
    n = 0
    for line in lines:
        if not line or line.startswith("#"):
            continue
        c = line.split("\t")
        if len(c) != 9 or c[2] != "gene":
            continue
        attrs = dict(f.split("=", 1) for f in c[8].split(";") if "=" in f)
        gid = attrs.get("gene_id") or attrs.get("ID")
        if gid not in keep:
            continue
        n += 1
        scaf, start, end = c[0], int(c[3]), int(c[4])
        strand = 1 if c[6] == "+" else -1 if c[6] == "-" else 0
        crosswalk.append({"species": "petunia", "atlas_gene_id": gid, "ext_gene_id": gid,
                          "ext_assembly": ASSEMBLY, "relation": "1:1"})
        ps, pe = promoter(start, end, strand)
        windows.append({"ext_gene_id": gid, "assembly": ASSEMBLY, "window_type": "promoter",
                        "chromosome": scaf, "start": ps, "end": pe, "strand": strand})
        windows.append({"ext_gene_id": gid, "assembly": ASSEMBLY, "window_type": "gene_body",
                        "chromosome": scaf, "start": start - 1, "end": end, "strand": strand})

    with gzip.open(CROSSWALK, "wt", encoding="utf-8") as f:
        json.dump(crosswalk, f)
    with gzip.open(WINDOWS, "wt", encoding="utf-8") as f:
        json.dump(windows, f)
    print(f"matched {n} petunia genes")
    print(f"Wrote {CROSSWALK} ({len(crosswalk)} rows)")
    print(f"Wrote {WINDOWS} ({len(windows)} rows: promoter + gene_body)")


if __name__ == "__main__":
    main()

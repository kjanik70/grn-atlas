"""
Petunia motif scan (Path B, generalized). Same JASPAR PWM scan as tomato, on the
P. axillaris v1.6.2 scaffold assembly, over petunia promoter windows.

NB: every petunia regulatory edge is INFERRED (projected from Arabidopsis), so
these sites are predicted TF binding sites in petunia promoters supporting
*inferred* edges — a doubly-computational (but real-sequence) hypothesis. Labeled
tier='JASPAR_scan'. Especially aimed at flavonoid/floral genes (MYB/bHLH incl.
the petunia-derived MA0054 myb.Ph3 matrix).

Reuses the PWM machinery from scan_tomato_motifs.py.

Usage: python backend/scripts/scan_petunia_motifs.py /tmp/pax.fa
"""
import gzip
import json
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from scan_tomato_motifs import load_jaspar, encode, scan, SCALE, PVAL, COMP  # noqa: E402

DATA_DIR = Path(__file__).parent.parent / "data"
DB_PATH = DATA_DIR / "grn.sqlite3"
MOTIFS_JSON = DATA_DIR / "motifs_petunia.json"
HITS_JSON = DATA_DIR / "motif_hits_petunia.json.gz"
ASSEMBLY = "Peaxi162v1.6.2"


def main(fasta_path):
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    jaspar = load_jaspar()
    name_to_mid = {}
    for mid, (name, *_ ) in jaspar.items():
        for tok in re.split(r"[\/():]", name):
            tok = tok.strip().upper()
            if tok:
                name_to_mid.setdefault(tok, mid)

    tf_mats, tf_symbol = defaultdict(set), {}
    seq_map = DATA_DIR / "tf_motif_map_petunia.json"
    smap = json.loads(seq_map.read_text()) if seq_map.exists() else None
    for r in conn.execute("""SELECT DISTINCT g.id, g.symbol, g.synonyms FROM genes g
                             JOIN interactions i ON i.source_id=g.id
                             WHERE g.species='petunia'"""):
        if smap is not None:
            mats = set(smap.get(r["id"], []))
        else:
            cands = {(r["symbol"] or "").upper()}
            if r["synonyms"]:
                cands |= {s.upper() for s in r["synonyms"].split("; ")}
            mats = {name_to_mid[c] for c in cands if c in name_to_mid}
        if mats:
            tf_mats[r["id"]] = mats
            syns = r["synonyms"].split("; ") if r["synonyms"] else []
            tf_symbol[r["id"]] = syns[0] if syns else r["symbol"]
    print(f"petunia TFs mapped to a JASPAR motif: {len(tf_mats)}"
          + (" (sequence map)" if smap is not None else " (symbol match)"))

    # petunia promoter windows (identity crosswalk; keyed by atlas id)
    prom = {r["ext_gene_id"]: r for r in conn.execute(
        "SELECT ext_gene_id, chromosome, start, end, strand FROM gene_windows "
        "WHERE window_type='promoter' AND assembly=?", (ASSEMBLY,))}

    scan_jobs = defaultdict(set)   # target -> {tf}
    for r in conn.execute("SELECT source_id, target_id FROM interactions WHERE source_id LIKE 'Peaxi%'"):
        if r["source_id"] in tf_mats and r["target_id"] in prom:
            scan_jobs[r["target_id"]].add(r["source_id"])
    print(f"target promoters to scan: {len(scan_jobs)}")

    by_chrom = defaultdict(list)
    for ext in scan_jobs:
        by_chrom[prom[ext]["chromosome"]].append(ext)

    seqs = {}
    want = set(by_chrom)
    cur_name, cur = None, []

    def flush(name, chunks):
        if name is None or name not in want:
            return
        s = "".join(chunks).upper()
        for ext in by_chrom[name]:
            w = prom[ext]
            sub = s[w["start"]:w["end"]]
            if w["strand"] < 0:
                sub = sub.translate(COMP)[::-1]
            seqs[ext] = sub

    print("Reading petunia scaffold FASTA…")
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush(cur_name, cur)
                cur_name, cur = line[1:].split()[0], []
            else:
                cur.append(line.strip())
    flush(cur_name, cur)
    print(f"promoter sequences extracted: {len(seqs)}")

    motifs_out, hits = {}, []
    for ext, tfs in scan_jobs.items():
        seq = seqs.get(ext)
        if not seq:
            continue
        si, si_rc = encode(seq), encode(seq.translate(COMP)[::-1])
        w = prom[ext]
        for tf in tfs:
            for mid in tf_mats[tf]:
                name, lom, ints, thr, pval_of = jaspar[mid]
                L = ints.shape[0]
                syn_id = f"{mid}|{tf}"
                for strand, arr in ((1, si), (-1, si_rc)):
                    for off, sc in scan(arr, ints, thr, pval_of):
                        gs = (w["start"] + off) if strand == 1 else (w["end"] - (off + L))
                        eff = (w["strand"] or 1) * strand
                        hits.append({
                            "ext_gene_id": ext, "motif_id": syn_id, "assembly": ASSEMBLY,
                            "window_type": "promoter", "chromosome": w["chromosome"],
                            "start": gs, "end": gs + L, "strand": 1 if eff >= 0 else -1,
                            "score": round(sc / SCALE, 3), "p_value": float(pval_of.get(sc, PVAL)),
                            "tier": "JASPAR_scan", "site_confidence": 0.5,
                        })
                        motifs_out[syn_id] = {"motif_id": syn_id, "source": "JASPAR2024",
                                              "jaspar_id": mid, "tf_gene_id": tf,
                                              "tf_symbol": tf_symbol.get(tf)}

    MOTIFS_JSON.write_text(json.dumps(list(motifs_out.values()), indent=1))
    with gzip.open(HITS_JSON, "wt", encoding="utf-8") as f:
        json.dump(hits, f)
    print(f"Wrote {MOTIFS_JSON} ({len(motifs_out)} TF-motif rows)")
    print(f"Wrote {HITS_JSON} ({len(hits)} binding sites, p<{PVAL})")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/tmp/pax.fa")

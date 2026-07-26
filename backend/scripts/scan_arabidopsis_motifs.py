"""Arabidopsis motif scan (#4). Same JASPAR plant PWM scan as tomato/petunia, over
TAIR10 promoter windows, edge-driven (scan a target's promoter for the PWMs of its
regulators). JASPAR's plant collection is Arabidopsis-centric, so symbol matching
gives good TF->motif coverage without BLAST.

Sites are predicted (tier='JASPAR_scan'), labelled distinct from measured data.

Usage: python backend/scripts/scan_arabidopsis_motifs.py /tmp/ath.fa
"""
import gzip
import json
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np  # noqa: F401  (used indirectly via imported scan machinery)

sys.path.insert(0, str(Path(__file__).parent))
from scan_tomato_motifs import load_jaspar, encode, scan, SCALE, PVAL, COMP  # noqa: E402

DATA_DIR = Path(__file__).parent.parent / "data"
DB_PATH = DATA_DIR / "grn.sqlite3"
MOTIFS_JSON = DATA_DIR / "motifs_arabidopsis.json"
HITS_JSON = DATA_DIR / "motif_hits_arabidopsis.json.gz"
ASSEMBLY = "TAIR10"


def main(fasta_path):
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row

    jaspar = load_jaspar()
    name_to_mid = {}
    for mid, (name, *_) in jaspar.items():
        for tok in re.split(r"[\/():]", name):
            tok = tok.strip().upper()
            if tok:
                name_to_mid.setdefault(tok, mid)

    # arabidopsis TF -> matrices via symbol/synonym (or sequence map if present)
    tf_mats = defaultdict(set)
    tf_symbol = {}
    seq_map = DATA_DIR / "tf_motif_map_arabidopsis.json"
    smap = json.loads(seq_map.read_text()) if seq_map.exists() else None
    for r in conn.execute("""SELECT DISTINCT g.id, g.symbol, g.synonyms FROM genes g
                             JOIN interactions i ON i.source_id=g.id
                             WHERE g.species='arabidopsis'"""):
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
    print(f"arabidopsis TFs mapped to a JASPAR motif: {len(tf_mats)}"
          + (" (sequence map)" if smap is not None else " (symbol match)"))

    atlas2ext = {r["atlas_gene_id"]: r["ext_gene_id"] for r in conn.execute(
        "SELECT atlas_gene_id, ext_gene_id FROM gene_id_crosswalk WHERE species='arabidopsis'")}
    prom = {r["ext_gene_id"]: r for r in conn.execute(
        "SELECT ext_gene_id, chromosome, start, end, strand FROM gene_windows "
        "WHERE window_type='promoter' AND assembly=?", (ASSEMBLY,))}

    scan_jobs = defaultdict(set)   # ext_target -> set(tf_id)
    for r in conn.execute("""SELECT i.source_id, i.target_id FROM interactions i
                             JOIN genes t ON t.id=i.target_id
                             WHERE t.species='arabidopsis'"""):
        tf, tgt = r["source_id"], r["target_id"]
        ext = atlas2ext.get(tgt)
        if tf in tf_mats and ext in prom:
            scan_jobs[ext].add(tf)
    print(f"target promoters to scan: {len(scan_jobs)}")

    by_chrom = defaultdict(list)
    for ext in scan_jobs:
        by_chrom[str(prom[ext]["chromosome"])].append(ext)
    want_chroms = set(by_chrom)
    seqs = {}

    def flush(name, chunks):
        if name is None or name not in want_chroms:
            return
        s = "".join(chunks).upper()
        for ext in by_chrom[name]:
            w = prom[ext]
            sub = s[w["start"]:w["end"]]
            if w["strand"] < 0:
                sub = sub.translate(COMP)[::-1]
            seqs[ext] = sub

    print("Reading TAIR10 FASTA…")
    cur_name, cur = None, []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush(cur_name, cur)
                cur_name, cur = line[1:].split()[0], []
            else:
                cur.append(line.strip())
    flush(cur_name, cur)
    print(f"promoter sequences extracted: {len(seqs)}")

    motifs_out, hits_out = {}, []
    for ext, tfs in scan_jobs.items():
        seq = seqs.get(ext)
        if not seq:
            continue
        si = encode(seq)
        si_rc = encode(seq.translate(COMP)[::-1])
        w = prom[ext]
        for tf in tfs:
            for mid in tf_mats[tf]:
                name, lom, ints, thr, pval_of = jaspar[mid]
                L = ints.shape[0]
                syn_id = f"{mid}|{tf}"
                for strand, arr in ((1, si), (-1, si_rc)):
                    for off, sc in scan(arr, ints, thr, pval_of):
                        if strand == 1:
                            gs = w["start"] + off
                        else:
                            gs = w["end"] - (off + L)
                        eff_strand = (w["strand"] or 1) * strand
                        hits_out.append({
                            "ext_gene_id": ext, "motif_id": syn_id, "assembly": ASSEMBLY,
                            "window_type": "promoter", "chromosome": str(w["chromosome"]),
                            "start": gs, "end": gs + L, "strand": 1 if eff_strand >= 0 else -1,
                            "score": round(sc / SCALE, 3),
                            "p_value": float(pval_of.get(sc, PVAL)),
                            "tier": "JASPAR_scan", "site_confidence": 0.5,
                        })
                        motifs_out[syn_id] = {
                            "motif_id": syn_id, "source": "JASPAR2024", "jaspar_id": mid,
                            "tf_gene_id": tf, "tf_symbol": tf_symbol.get(tf),
                        }

    MOTIFS_JSON.write_text(json.dumps(list(motifs_out.values()), indent=1))
    with gzip.open(HITS_JSON, "wt", encoding="utf-8") as f:
        json.dump(hits_out, f)
    print(f"Wrote {MOTIFS_JSON} ({len(motifs_out)} TF-motif rows)")
    print(f"Wrote {HITS_JSON} ({len(hits_out)} binding sites, p<{PVAL})")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/tmp/ath.fa")

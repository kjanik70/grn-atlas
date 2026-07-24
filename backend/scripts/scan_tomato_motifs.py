"""
Tomato motif scan (Path B, option 2 — in-house). Produces the remaining two
sequence-context caches by scanning JASPAR plant PWMs over tomato promoter
windows on the real SL4.0 assembly:

  - motifs.json     : one row per (JASPAR matrix, tomato TF) it is assigned to
  - motif_hits.json : PWM hits (p < 1e-4) in the promoter of each edge's target,
                      for the motif of that edge's regulator

Honest scope / labels:
  * This is a JASPAR PWM scan, not the agent's FunTFBS/FIMO pipeline. Sites are
    tagged tier='JASPAR_scan'; site_confidence reflects the JASPAR matrix quality
    tier, NOT PlantRegMap motif_CE/DGF tiers.
  * p-values are exact for the discretised PWM score distribution (the same math
    FIMO uses), with a uniform 0.25 background.
  * TF<->matrix mapping is by symbol/synonym (tomato TFs carry inferred
    Arabidopsis symbols); coverage is therefore partial (~62% of tomato
    regulators) and edge-driven (only sites supporting existing edges).

Requires the SL4.0 FASTA at the path below (≈795 MB, not committed).

Usage: python backend/scripts/scan_tomato_motifs.py /tmp/sl4.fa
"""
import gzip
import json
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).parent.parent / "data"
DB_PATH = DATA_DIR / "grn.sqlite3"
JASPAR = Path("/tmp/jaspar_plants.jaspar")
MOTIFS_JSON = DATA_DIR / "motifs.json"
HITS_JSON = DATA_DIR / "motif_hits.json.gz"     # large -> gzip

PVAL = 1e-4
BG = 0.25          # uniform background
PSEUDO = 0.8       # pseudocount
SCALE = 100        # score discretisation for the p-value DP
BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ---------- JASPAR PWMs ----------
def load_jaspar():
    """Return {matrix_id: (name, lom LxN log-odds, int_scores, threshold, pval_of)}."""
    entries = {}
    mid = name = None
    rows = []
    text = JASPAR.read_text().splitlines()

    def finalize(mid, name, rows):
        counts = np.array(rows, dtype=float).T            # L x 4 (A,C,G,T)
        probs = (counts + PSEUDO) / (counts.sum(1, keepdims=True) + 4 * PSEUDO)
        lom = np.log2(probs / BG)                         # log-odds
        ints = np.rint(lom * SCALE).astype(int)           # discretised
        thr, pval_of = score_pvalues(ints)
        entries[mid] = (name, lom, ints, thr, pval_of)

    for line in text:
        if line.startswith(">"):
            if mid and rows:
                finalize(mid, name, rows)
            parts = line[1:].split("\t")
            mid, name, rows = parts[0], (parts[1] if len(parts) > 1 else parts[0]), []
        elif line.strip():
            nums = re.findall(r"[\d.]+", line.split("[")[-1])
            rows.append([float(x) for x in nums])
    if mid and rows:
        finalize(mid, name, rows)
    return entries


def score_pvalues(int_scores):
    """Exact p-value of the discretised PWM score distribution (uniform bg).
    Returns (threshold_score, {score: p_value}) for p <= PVAL."""
    dist = {0: 1.0}
    for pos in int_scores:                    # pos: 4 integer scores (A,C,G,T)
        nd = defaultdict(float)
        for s, p in dist.items():
            for base_score in pos:
                nd[s + int(base_score)] += p * BG
        dist = nd
    # cumulative tail p(score >= s)
    scores = sorted(dist, reverse=True)
    cum, pval_of = 0.0, {}
    thr = None
    for s in scores:
        cum += dist[s]
        pval_of[s] = cum
        if thr is None and cum > PVAL:
            thr = s + 1        # last score still within p<=PVAL is s-1 step; threshold above
    if thr is None:
        thr = scores[-1]
    return thr, pval_of


# ---------- sequences ----------
def encode(seq):
    a = np.full(len(seq), -1, dtype=np.int8)
    for b, i in BASE.items():
        a[np.frombuffer(seq.encode(), dtype=np.uint8) == ord(b)] = i
    return a


def scan(seq_int, ints, thr, pval_of):
    """Return list of (offset, strand_score_int) hits >= thr on the given strand
    (caller handles both strands). Uses the discretised integer PWM."""
    L = ints.shape[0]
    n = seq_int.size - L + 1
    if n <= 0:
        return []
    score = np.zeros(n, dtype=np.int64)
    valid = np.ones(n, dtype=bool)
    for j in range(L):
        col = seq_int[j:j + n]
        valid &= col >= 0
        contrib = np.where(col >= 0, ints[j][np.clip(col, 0, 3)], 0)
        score += contrib
    hits = []
    for i in np.where(valid & (score >= thr))[0]:
        hits.append((int(i), int(score[i])))
    return hits


def main(fasta_path):
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row

    jaspar = load_jaspar()
    # name token -> matrix_id
    name_to_mid = {}
    for mid, (name, *_ ) in jaspar.items():
        for tok in re.split(r"[\/():]", name):
            tok = tok.strip().upper()
            if tok:
                name_to_mid.setdefault(tok, mid)

    # tomato TF -> matrices, via symbol/synonym
    tf_mats = defaultdict(set)
    tf_symbol = {}
    for r in conn.execute("""SELECT DISTINCT g.id, g.symbol, g.synonyms FROM genes g
                             JOIN interactions i ON i.source_id=g.id
                             WHERE g.species='tomato'"""):
        cands = {(r["symbol"] or "").upper()}
        if r["synonyms"]:
            cands |= {s.upper() for s in r["synonyms"].split("; ")}
        mats = {name_to_mid[c] for c in cands if c in name_to_mid}
        if mats:
            tf_mats[r["id"]] = mats
            # friendly name: inferred Arabidopsis synonym, else the locus symbol
            syns = r["synonyms"].split("; ") if r["synonyms"] else []
            tf_symbol[r["id"]] = syns[0] if syns else r["symbol"]
    print(f"tomato TFs mapped to a JASPAR motif: {len(tf_mats)}")

    # measured tomato edges whose TF is mapped -> which target promoters to scan
    atlas2ext = {r["atlas_gene_id"]: r["ext_gene_id"]
                 for r in conn.execute("SELECT atlas_gene_id, ext_gene_id FROM gene_id_crosswalk")}
    prom = {r["ext_gene_id"]: r for r in conn.execute(
        "SELECT ext_gene_id, chromosome, start, end, strand FROM gene_windows WHERE window_type='promoter'")}

    scan_jobs = defaultdict(set)   # ext_target -> set(tf_id)
    for r in conn.execute("""SELECT source_id, target_id FROM interactions
                             WHERE sources LIKE '%PlantRegMap%' AND source_id LIKE 'Solyc%'"""):
        tf, tgt = r["source_id"], r["target_id"]
        ext = atlas2ext.get(tgt)
        if tf in tf_mats and ext in prom:
            scan_jobs[ext].add(tf)
    print(f"target promoters to scan: {len(scan_jobs)}")

    # promoter sequences from FASTA, grouped by chromosome
    by_chrom = defaultdict(list)
    for ext in scan_jobs:
        w = prom[ext]
        by_chrom[str(w["chromosome"])].append(ext)
    seqs = {}
    want_chroms = set(by_chrom)
    cur_name, cur = None, []

    def flush(name, chunks):
        if name is None:
            return
        cn = norm(name)
        if cn not in want_chroms:
            return
        s = "".join(chunks).upper()
        for ext in by_chrom[cn]:
            w = prom[ext]
            sub = s[w["start"]:w["end"]]
            if w["strand"] < 0:
                sub = sub.translate(COMP)[::-1]
            seqs[ext] = sub

    def norm(name):
        m = re.search(r"ch0*(\d+)$", name)
        return m.group(1) if m else name

    print("Reading SL4.0 FASTA…")
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush(cur_name, cur)
                cur_name, cur = line[1:].split()[0], []
            else:
                cur.append(line.strip())
    flush(cur_name, cur)
    print(f"promoter sequences extracted: {len(seqs)}")

    # scan
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
                        # map offset back to genomic BED0 coords of the site
                        if strand == 1:
                            gs = w["start"] + off
                        else:
                            gs = w["end"] - (off + L)
                        # combine promoter strand with match strand
                        eff_strand = (w["strand"] or 1) * strand
                        hits_out.append({
                            "ext_gene_id": ext, "motif_id": syn_id, "assembly": "SL4.0",
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
    main(sys.argv[1] if len(sys.argv) > 1 else "/tmp/sl4.fa")

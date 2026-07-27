"""In-silico dsRNA / RNAi analysis (predicted gene silencing).

A dsRNA is diced into ~21 nt siRNAs that can silence any transcript sharing a
sufficiently long exact match (either strand). This module predicts, for a dsRNA:
  - on-target coverage of an intended gene,
  - off-target genes (other transcripts sharing >=k-nt matches), ranked,
  - a specificity score,
and, in design mode, the most specific window within a target transcript.

STRONGLY PREDICTED, NOT MEASURED: exact k-mer matching is a specificity heuristic.
Real RNAi knockdown also depends on dicing efficiency, delivery/SIGS uptake, target
accessibility, and (in plants) transitivity/amplification — none modelled here.
Everything is labelled "predicted silencing".

Pure functions (no FastAPI / no DB) so they are unit-testable in isolation.
"""
import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def kmers(seq: str, k: int = 21):
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        w = seq[i:i + k]
        if "N" not in w:
            yield i, w


def query_kmers(dsrna: str, k: int = 21) -> set:
    """All k-mers of the dsRNA AND its reverse complement (both strands are diced)."""
    s = dsrna.upper()
    out = {w for _, w in kmers(s, k)}
    out |= {w for _, w in kmers(revcomp(s), k)}
    return out


def gene_of(header_token: str) -> str:
    return header_token.rsplit(".", 1)[0]


def load_transcripts(path: Path) -> Dict[str, str]:
    """gene_id -> concatenated transcript sequence(s) (isoforms joined by a k-blocker).
    Reads a (gzipped) FASTA whose headers start with a transcript id."""
    genes: Dict[str, List[str]] = defaultdict(list)
    opener = gzip.open if str(path).endswith(".gz") else open
    gid, buf = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if gid:
                    genes[gid].append("".join(buf))
                gid = gene_of(line[1:].split()[0])
                buf = []
            else:
                buf.append(line.strip())
        if gid:
            genes[gid].append("".join(buf))
    # join isoforms with a run of N so no k-mer spans two isoforms
    return {g: ("N" * 21).join(seqs).upper() for g, seqs in genes.items()}


def scan(dsrna: str, transcripts: Dict[str, str], k: int = 21,
         target_gene: Optional[str] = None) -> dict:
    """Predict silencing of `dsrna` across the transcriptome.

    Returns on-target sites (if target_gene given), ranked off-target genes, and a
    specificity score = on_target_sites / (on_target_sites + total_offtarget_sites).
    A 'site' is one k-mer position in a gene matching a dsRNA k-mer (either strand).
    """
    qk = query_kmers(dsrna, k)
    per_gene = {}
    for gid, seq in transcripts.items():
        n = sum(1 for _, w in kmers(seq, k) if w in qk)
        if n:
            per_gene[gid] = n

    on = per_gene.get(target_gene, 0) if target_gene else 0
    offs = [{"gene_id": g, "sites": n} for g, n in per_gene.items() if g != target_gene]
    offs.sort(key=lambda d: d["sites"], reverse=True)
    off_total = sum(o["sites"] for o in offs)
    denom = on + off_total
    specificity = round(on / denom, 4) if denom else 0.0
    return {
        "k": k,
        "dsrna_length": len(dsrna),
        "n_sirnas": len(qk),
        "on_target_gene": target_gene,
        "on_target_sites": on,
        "off_target_gene_count": len(offs),
        "off_target_total_sites": off_total,
        "off_targets": offs,
        "specificity": specificity,
        "silenced_genes": ([target_gene] if on else []) + [o["gene_id"] for o in offs],
    }


def design(target_gene: str, transcripts: Dict[str, str], k: int = 21,
           window: int = 250, step: int = 25) -> dict:
    """Pick the most specific dsRNA window within a target transcript: the window
    whose k-mers hit the fewest OTHER genes. Returns the window sequence + its
    off-target profile."""
    tseq = transcripts.get(target_gene)
    if not tseq or len(tseq) < k:
        return {"error": f"no transcript for {target_gene}"}
    # first, the target's own k-mers, and which other genes each one appears in
    tk = {w: i for i, w in kmers(tseq, k)}
    off_by_kmer: Dict[str, set] = {w: set() for w in tk}
    for gid, seq in transcripts.items():
        if gid == target_gene:
            continue
        for _, w in kmers(seq, k):
            if w in off_by_kmer:
                off_by_kmer[w].add(gid)
    # slide a window over the target; score = number of distinct off-target genes hit
    best = None
    L = len(tseq)
    win = min(window, L)
    for start in range(0, max(1, L - win + 1), step):
        sub = tseq[start:start + win]
        hit = set()
        for _, w in kmers(sub, k):
            hit |= off_by_kmer.get(w, set())
        if best is None or len(hit) < best["off_target_gene_count"]:
            best = {"start": start, "end": start + win, "sequence": sub,
                    "off_target_gene_count": len(hit),
                    "off_target_genes": sorted(hit)[:50]}
    best["target_gene"] = target_gene
    best["window"] = win
    best["note"] = "Predicted most-specific window (fewest off-target genes); verify with scan()."
    return best


_cache: Dict[str, Optional[Dict[str, str]]] = {}


def get_transcripts(species: str, data_dir: Path) -> Optional[Dict[str, str]]:
    if species not in _cache:
        path = data_dir / f"transcripts_{species}.fasta.gz"
        _cache[species] = load_transcripts(path) if path.exists() else None
    return _cache[species]

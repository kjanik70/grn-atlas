"""Targeted loader for the Arabidopsis sequence-context layer (#4) into grn.sqlite3
WITHOUT a full DB rebuild. Idempotent: clears existing arabidopsis/TAIR10 rows first.

Loads whatever caches are present:
  gene_id_crosswalk_arabidopsis.json.gz, gene_windows_arabidopsis.json.gz  (always)
  motifs_arabidopsis.json, motif_hits_arabidopsis.json.gz                  (if non-empty)

Run once before scanning (crosswalk + windows), then again after scanning (motifs + hits).
"""
import gzip
import json
import sqlite3
from pathlib import Path

DATA = Path(__file__).parent.parent / "data"
DB = DATA / "grn.sqlite3"


def _load_json(path):
    if not path.exists():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as f:
        return json.load(f)


def main():
    conn = sqlite3.connect(DB)
    cw = _load_json(DATA / "gene_id_crosswalk_arabidopsis.json.gz")
    win = _load_json(DATA / "gene_windows_arabidopsis.json.gz")
    conn.execute("DELETE FROM gene_id_crosswalk WHERE species='arabidopsis'")
    conn.executemany(
        "INSERT OR REPLACE INTO gene_id_crosswalk "
        "(species, atlas_gene_id, ext_gene_id, ext_assembly, relation) VALUES (?,?,?,?,?)",
        [(r["species"], r["atlas_gene_id"], r["ext_gene_id"], r["ext_assembly"], r["relation"]) for r in cw])
    conn.execute("DELETE FROM gene_windows WHERE assembly='TAIR10'")
    conn.executemany(
        "INSERT OR REPLACE INTO gene_windows "
        "(ext_gene_id, assembly, window_type, chromosome, start, end, strand) VALUES (?,?,?,?,?,?,?)",
        [(r["ext_gene_id"], r["assembly"], r["window_type"], r["chromosome"],
          r["start"], r["end"], r["strand"]) for r in win])
    print(f"loaded crosswalk={len(cw)} windows={len(win)}")

    motifs = _load_json(DATA / "motifs_arabidopsis.json")
    hits = _load_json(DATA / "motif_hits_arabidopsis.json.gz")
    if motifs:
        conn.execute("DELETE FROM motifs WHERE motif_id LIKE '%|AT%'")
        conn.executemany(
            "INSERT OR REPLACE INTO motifs (motif_id, source, jaspar_id, tf_gene_id, tf_symbol) "
            "VALUES (?,?,?,?,?)",
            [(m["motif_id"], m["source"], m.get("jaspar_id"), m.get("tf_gene_id"), m.get("tf_symbol"))
             for m in motifs])
    if hits:
        conn.execute("DELETE FROM motif_hits WHERE assembly='TAIR10'")
        conn.executemany(
            "INSERT INTO motif_hits (ext_gene_id, motif_id, assembly, window_type, chromosome, "
            "start, end, strand, score, p_value, tier, site_confidence) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            [(h["ext_gene_id"], h["motif_id"], h["assembly"], h["window_type"], h["chromosome"],
              h["start"], h["end"], h["strand"], h["score"], h["p_value"], h["tier"], h["site_confidence"])
             for h in hits])
    print(f"loaded motifs={len(motifs)} hits={len(hits)}")
    conn.commit()
    conn.close()


if __name__ == "__main__":
    main()

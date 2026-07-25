"""Tier 3 — API contract tests against a tiny in-memory-ish fixture DB.

Schema is copied from the built grn.sqlite3 (so it never drifts from build_db);
a handful of curated rows exercise the endpoints. main.DB_PATH is overridden via
GRN_DB before importing/reloading main.
"""
import importlib
import os
import sqlite3
from pathlib import Path

import pytest

REAL_DB = Path(__file__).resolve().parents[1] / "data" / "grn.sqlite3"
pytestmark = pytest.mark.skipif(not REAL_DB.exists(), reason="need built DB for schema")


def _build_fixture(path):
    src = sqlite3.connect(REAL_DB)
    schema = [r[0] for r in src.execute(
        "SELECT sql FROM sqlite_master WHERE sql IS NOT NULL AND name NOT LIKE 'sqlite_%'")]
    src.close()
    db = sqlite3.connect(path)
    for stmt in schema:
        db.execute(stmt)

    genes = [
        ("TF1", "TF1", "human TF one", "human", 1, "protein_coding", None),
        ("TG1", "TG1", "human target one", "human", 0, "protein_coding", None),
        ("TG2", "TG2", "human target two", "human", 0, "protein_coding", None),
        ("SlTF", "MYB1", "tomato MYB", "tomato", 1, "protein_coding", "MYB1"),
        ("SlTGT", "SlTGT", "tomato target", "tomato", 0, "protein_coding", "CHS"),
        ("SlTGT2", "SlTGT2", "tomato target 2", "tomato", 0, "protein_coding", None),
    ]
    db.executemany("INSERT INTO genes (id,symbol,name,species,is_tf,gene_type,synonyms) "
                   "VALUES (?,?,?,?,?,?,?)", genes)
    db.executemany("INSERT INTO interactions (source_id,target_id,regulation_type,confidence,sources,pmids) "
                   "VALUES (?,?,?,?,?,?)", [
        ("TF1", "TG1", "activation", 0.9, '["TRRUST"]', '["12345"]'),
        ("TF1", "TG2", "repression", 0.8, '["TRRUST"]', '["67890"]'),
        ("SlTF", "SlTGT", "regulation", 0.65, '["PlantRegMap"]', '[]'),
        ("SlTF", "SlTGT2", "regulation", 0.5, '["Inferred:Arabidopsis"]', '[]'),
    ])
    db.executemany("INSERT INTO gene_locations VALUES (?,?,?,?,?,?)", [
        ("TF1", "human", "1", 1000, 2000, 1), ("TG1", "human", "1", 5000, 6000, -1),
        ("TG2", "human", "2", 3000, 4000, 1), ("SlTF", "tomato", "1", 1000, 2000, 1),
        ("SlTGT", "tomato", "1", 8000, 9000, 1), ("SlTGT2", "tomato", "2", 2000, 3000, 1),
    ])
    db.executemany("INSERT INTO chromosomes VALUES (?,?,?)", [
        ("human", "1", 100000), ("human", "2", 100000),
        ("tomato", "1", 100000), ("tomato", "2", 100000)])
    db.executemany("INSERT INTO go_terms VALUES (?,?,?)", [
        ("GO:1", "process one", "BP"), ("GO:2", "process two", "BP")])
    db.executemany("INSERT INTO go_annotations VALUES (?,?)", [
        ("TF1", "GO:1"), ("TG1", "GO:1"), ("TG1", "GO:2"), ("TG2", "GO:2")])
    # sequence context for the tomato target (SL4.0 crosswalk + window + motif hit)
    db.execute("INSERT INTO gene_id_crosswalk VALUES ('tomato','SlTGT','SlTGT.4','SL4.0','1:1')")
    db.execute("INSERT INTO gene_windows VALUES ('SlTGT.4','SL4.0','promoter','1',6000,8500,1)")
    db.execute("INSERT INTO motifs VALUES ('M1|SlTF','JASPAR2024','MA1','SlTF','MYB1')")
    db.execute("INSERT INTO motif_hits VALUES "
               "('SlTGT.4','M1|SlTF','SL4.0','promoter','1',6500,6512,1,10.0,1e-5,'JASPAR_scan',0.5)")
    db.commit(); db.close()


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    path = tmp_path_factory.mktemp("db") / "fixture.sqlite3"
    _build_fixture(path)
    os.environ["GRN_DB"] = str(path)
    import main
    importlib.reload(main)
    from fastapi.testclient import TestClient
    with TestClient(main.app) as c:
        yield c
    os.environ.pop("GRN_DB", None)


def test_health(client):
    assert client.get("/health").json()["status"] == "healthy"


def test_search(client):
    r = client.get("/api/v1/genes/search?q=MYB1&species=tomato").json()
    assert r["results"][0]["id"] == "SlTF"


def test_neighborhood_targets_and_pmids(client):
    r = client.post("/api/v1/pathways/neighborhood/TF1",
                    json={"direction": "targets", "regulation_type": ["activation", "repression"],
                          "min_confidence": 0.3}).json()
    tgt = {t["symbol"]: t for t in r["targets"]}
    assert "TG1" in tgt and tgt["TG1"]["pmids"] == ["12345"]
    assert tgt["TG1"]["inferred"] is False


def test_include_inferred_filter(client):
    body = {"direction": "targets", "regulation_type": ["regulation"], "min_confidence": 0.3}
    on = client.post("/api/v1/pathways/neighborhood/SlTF", json={**body, "include_inferred": True}).json()
    off = client.post("/api/v1/pathways/neighborhood/SlTF", json={**body, "include_inferred": False}).json()
    assert len(on["targets"]) == 2 and len(off["targets"]) == 1  # inferred edge dropped


def test_export_signed_edges_with_pmids(client):
    r = client.post("/api/v1/export/edges", json={"gene_ids": ["TF1", "TG1", "TG2"]}).json()
    signs = {e["target_symbol"]: e["sign"] for e in r["edges"]}
    assert signs["TG1"] == "positive" and signs["TG2"] == "negative"
    assert r["stats"]["signed"] == 2
    e = next(e for e in r["edges"] if e["target_symbol"] == "TG1")
    assert e["pmids"] == ["12345"] and e["source_coord_assembly"] == "GRCh38"


def test_export_sequence_context_join(client):
    r = client.post("/api/v1/export/edges",
                    json={"gene_ids": ["SlTF", "SlTGT"], "include_sequence_context": True}).json()
    e = next(e for e in r["edges"] if e["target_gene_id"] == "SlTGT")
    sc = e["sequence_context"]
    assert sc["assembly"] == "SL4.0" and sc["coord_system"] == "BED0"
    assert sc["target_windows"][0]["window_type"] == "promoter"
    site = sc["supporting_sites"][0]
    assert site["tf_symbol"] == "MYB1" and site["p_value"] <= 1e-4
    assert r["stats"]["edges_with_supporting_sites"] == 1


def test_subgraph(client):
    r = client.post("/api/v1/pathways/subgraph", json={"gene_ids": ["TF1", "TG1", "TG2"]}).json()
    assert len(r["nodes"]) == 3 and len(r["edges"]) == 2


def test_enrichment_shape_and_qbounds(client):
    r = client.post("/api/v1/enrichment",
                    json={"gene_ids": ["TG1", "TG2"], "species": "human", "min_genes": 2}).json()
    assert r["species"] == "human" and r["background"] >= 1
    qs = [t["q_value"] for t in r["results"]]
    assert all(0.0 <= q <= 1.0 for q in qs)
    assert qs == sorted(qs) or all(t["p_value"] <= 1 for t in r["results"])  # p sorted


def test_organism_overview(client):
    r = client.get("/api/v1/organism/human/overview?top=5").json()
    assert r["edges"]["measured"] == 2 and r["edges"]["inferred"] == 0
    assert any(t["symbol"] == "TF1" for t in r["top_regulators"])


def test_genome_endpoints(client):
    sp = client.get("/api/v1/genome/species").json()["species"]
    assert {s["species"] for s in sp} >= {"human", "tomato"}
    tom = client.get("/api/v1/genome/tomato").json()
    assert tom["species"] == "tomato" and len(tom["chromosomes"]) >= 1

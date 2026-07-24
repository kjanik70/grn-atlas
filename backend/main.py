"""
GRN Atlas FastAPI Backend
Complete example implementation with all required endpoints
"""

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel
from typing import List, Optional, Dict, Any
from datetime import datetime
from pathlib import Path as FilePath
from collections import defaultdict
import json
import math
import logging
import sqlite3

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="GRN Atlas API",
    description="Gene Regulatory Network visualization backend",
    version="1.0.0"
)

# ============= CORS Configuration =============
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://localhost:3001", "http://localhost:8000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ============= Pydantic Models =============

class Gene(BaseModel):
    id: str
    symbol: str
    name: str
    species: str
    ensembl_id: Optional[str] = None
    is_tf: bool = False
    gene_type: Optional[str] = None
    # Inferred alternative names (e.g. Arabidopsis ortholog symbols for
    # tomato/petunia). Approximate — surfaced separately from the real symbol.
    synonyms: Optional[List[str]] = None
    
    class Config:
        json_schema_extra = {
            "example": {
                "id": "ENSG00000141510",
                "symbol": "TP53",
                "name": "Tumor protein 53",
                "species": "human",
                "ensembl_id": "ENSG00000141510",
                "is_tf": True,
                "gene_type": "protein_coding"
            }
        }

class GeneInteraction(BaseModel):
    id: str
    symbol: str
    name: str
    species: str
    is_tf: bool
    confidence: float
    regulation_type: str  # 'activation', 'repression', 'unknown'
    source_databases: List[str]
    pmids: List[str] = []
    inferred: bool = False  # True when projected from another species' network

class NetworkData(BaseModel):
    gene: Gene
    regulators: List[GeneInteraction]
    targets: List[GeneInteraction]
    stats: Dict[str, int]

class PathGene(BaseModel):
    id: str
    symbol: str
    name: str

class Path(BaseModel):
    genes: List[PathGene]
    regulation_types: List[str]
    confidences: List[float]
    sources: List[List[str]]
    overall_confidence: float

class CascadeEffect(BaseModel):
    id: str
    symbol: str
    level: int
    direction: str  # 'up', 'down'
    magnitude: float
    confidence: float

class CascadeResult(BaseModel):
    cascade: List[CascadeEffect]
    average_confidence: float
    affected_genes: int

class Intervention(BaseModel):
    tf_id: str
    direction: str  # 'up', 'down'
    magnitude: float

class CascadeRequest(BaseModel):
    target_gene_id: str
    interventions: List[Intervention]
    depth: int = 3
    return_nodes: bool = True


class PathFindingRequest(BaseModel):
    source_gene_id: str
    target_symbol: str
    max_depth: int = 3
    limit: int = 20
    min_confidence: float = 0.3
    regulation_type: List[str] = ["activation", "repression", "regulation"]
    include_inferred: bool = True

class NeighborhoodRequest(BaseModel):
    max_depth: int = 1
    direction: str = "both"
    regulation_type: List[str] = ["activation", "repression", "regulation"]
    min_confidence: float = 0.3
    include_inferred: bool = True

# ============= Database Service =============
# Backed by a local SQLite database built from the full TRRUST v2 human
# TF-target corpus (https://www.grnpedia.org/trrust/), with gene names
# enriched from mygene.info at build time. See backend/scripts/build_db.py
# and backend/scripts/fetch_gene_names.py. No network access at runtime.

DB_PATH = FilePath(__file__).parent / "data" / "grn.sqlite3"


class GeneDatabase:
    """SQLite-backed gene/interaction lookups"""

    def __init__(self, db_path: FilePath):
        if not db_path.exists():
            logger.info("Database not found, building from local TRRUST data...")
            from scripts.build_db import build
            build()
        self.conn = sqlite3.connect(db_path, check_same_thread=False)
        self.conn.row_factory = sqlite3.Row

    def _row_to_gene(self, row) -> Gene:
        keys = row.keys()
        raw_syn = row["synonyms"] if "synonyms" in keys else None
        return Gene(
            id=row["id"],
            symbol=row["symbol"],
            name=row["name"],
            species=row["species"],
            is_tf=bool(row["is_tf"]),
            gene_type=row["gene_type"],
            synonyms=[s for s in raw_syn.split("; ") if s] if raw_syn else None,
        )

    def search_genes(self, query: str, limit: int = 10, species: Optional[str] = None) -> List[Gene]:
        """Search for genes by symbol or name"""
        sql = "SELECT * FROM genes WHERE (symbol LIKE ? OR name LIKE ? OR synonyms LIKE ?)"
        params: List[Any] = [f"%{query}%", f"%{query}%", f"%{query}%"]
        if species:
            sql += " AND species = ?"
            params.append(species)
        sql += " ORDER BY (symbol = ? COLLATE NOCASE) DESC, LENGTH(symbol) ASC LIMIT ?"
        params.extend([query, limit])
        rows = self.conn.execute(sql, params).fetchall()
        return [self._row_to_gene(r) for r in rows]

    def get_gene(self, gene_id: str) -> Optional[Gene]:
        """Get gene by ID"""
        row = self.conn.execute("SELECT * FROM genes WHERE id = ?", (gene_id,)).fetchone()
        return self._row_to_gene(row) if row else None

    def find_gene_by_symbol_species(self, symbol: str, species: str) -> Optional[Gene]:
        """Find a gene by symbol in a specific species"""
        row = self.conn.execute(
            "SELECT * FROM genes WHERE symbol = ? COLLATE NOCASE AND species = ?",
            (symbol, species)
        ).fetchone()
        return self._row_to_gene(row) if row else None

    @staticmethod
    def _row_to_interaction(row) -> GeneInteraction:
        sources = json.loads(row["sources"])
        pmids = json.loads(row["pmids"]) if "pmids" in row.keys() and row["pmids"] else []
        return GeneInteraction(
            id=row["id"], symbol=row["symbol"], name=row["name"], species=row["species"],
            is_tf=bool(row["is_tf"]), confidence=row["confidence"],
            regulation_type=row["regulation_type"], source_databases=sources,
            pmids=pmids, inferred=any(s.startswith("Inferred") for s in sources),
        )

    def get_regulators(self, gene_id: str, min_confidence: float = 0.0,
                       include_inferred: bool = True) -> List[GeneInteraction]:
        """Get regulators of a gene"""
        sql = """
            SELECT g.*, i.regulation_type, i.confidence, i.sources, i.pmids
            FROM interactions i JOIN genes g ON g.id = i.source_id
            WHERE i.target_id = ? AND i.confidence >= ?
        """
        if not include_inferred:
            sql += " AND i.sources NOT LIKE '%Inferred%'"
        rows = self.conn.execute(sql, (gene_id, min_confidence)).fetchall()
        return [self._row_to_interaction(r) for r in rows]

    def get_targets(self, gene_id: str, min_confidence: float = 0.0,
                    include_inferred: bool = True) -> List[GeneInteraction]:
        """Get targets of a gene"""
        sql = """
            SELECT g.*, i.regulation_type, i.confidence, i.sources, i.pmids
            FROM interactions i JOIN genes g ON g.id = i.target_id
            WHERE i.source_id = ? AND i.confidence >= ?
        """
        if not include_inferred:
            sql += " AND i.sources NOT LIKE '%Inferred%'"
        rows = self.conn.execute(sql, (gene_id, min_confidence)).fetchall()
        return [self._row_to_interaction(r) for r in rows]

# Initialize database
db = GeneDatabase(DB_PATH)

# ============= Root Endpoint =============

@app.get("/")
async def root():
    """Health check and API info"""
    return {
        "name": "GRN Atlas API",
        "version": "1.0.0",
        "status": "running",
        "timestamp": datetime.utcnow().isoformat()
    }

# ============= Gene Search Endpoints =============

@app.get("/api/v1/genes/search", response_model=Dict[str, List[Gene]])
async def search_genes(
    q: str = Query(..., min_length=1, description="Search query"),
    limit: int = Query(10, le=50, description="Maximum results"),
    species: Optional[str] = Query(None, description="Filter by species")
):
    """
    Search for genes by symbol or name
    
    - **q**: Search query (gene symbol or name)
    - **limit**: Maximum number of results (1-50)
    - **species**: Optional species filter
    """
    try:
        results = db.search_genes(q, limit, species)
        return {"results": results}
    except Exception as e:
        logger.error(f"Search error: {e}")
        raise HTTPException(status_code=500, detail="Search failed")

@app.get("/api/v1/genes/{gene_id}", response_model=Gene)
async def get_gene(gene_id: str):
    """
    Get gene details by Ensembl ID
    """
    gene = db.get_gene(gene_id)
    if not gene:
        raise HTTPException(status_code=404, detail="Gene not found")
    return gene

@app.get("/api/v1/genes/symbol/{symbol}", response_model=Gene)
async def get_gene_by_symbol(symbol: str):
    """
    Get gene details by symbol
    """
    results = db.search_genes(symbol, limit=1)
    if not results:
        raise HTTPException(status_code=404, detail="Gene not found")
    return results[0]

# ============= Pathway Endpoints =============

@app.post("/api/v1/pathways/neighborhood/{gene_id}", response_model=NetworkData)
async def get_neighborhood(gene_id: str, request: NeighborhoodRequest = NeighborhoodRequest()):
    """
    Get regulatory neighborhood around a gene

    - **gene_id**: Target gene Ensembl ID
    - **max_depth**: Maximum network hops (1-5)
    - **direction**: 'both', 'regulators', or 'targets'
    - **regulation_type**: Filter by regulation type
    - **min_confidence**: Minimum confidence score (0.3-0.9)
    """
    gene = db.get_gene(gene_id)
    if not gene:
        raise HTTPException(status_code=404, detail="Gene not found")

    # Get regulators and targets
    regulators = db.get_regulators(gene_id, request.min_confidence, request.include_inferred) if request.direction in ["both", "regulators"] else []
    targets = db.get_targets(gene_id, request.min_confidence, request.include_inferred) if request.direction in ["both", "targets"] else []

    # Filter by regulation type
    regulators = [r for r in regulators if r.regulation_type in request.regulation_type]
    targets = [t for t in targets if t.regulation_type in request.regulation_type]
    
    return NetworkData(
        gene=gene,
        regulators=regulators,
        targets=targets,
        stats={
            "regulators": len(regulators),
            "targets": len(targets),
            "paths": 0  # Can be calculated from path finding algorithm
        }
    )
@app.post("/api/v1/pathways/pathfinding", response_model=Dict[str, List[Path]])
async def find_paths(request: PathFindingRequest):
    """
    Find regulatory paths between two genes using BFS algorithm
    """
    source = db.get_gene(request.source_gene_id)
    target_results = db.search_genes(request.target_symbol, limit=1)
    
    if not source or not target_results:
        raise HTTPException(status_code=404, detail="Gene not found")
    
    target_gene = target_results[0]
    
    paths = []
    queue = [(request.source_gene_id, [source], [], [], [], {request.source_gene_id})]
    max_queue = 50000

    while queue and len(paths) < request.limit and len(queue) < max_queue:
        current_id, current_path, regulations, confidences, edge_sources, path_visited = queue.pop(0)

        if current_id == target_gene.id:
            path_genes = [PathGene(id=g.id, symbol=g.symbol, name=g.name) for g in current_path]
            paths.append(Path(
                genes=path_genes,
                regulation_types=regulations,
                confidences=confidences,
                sources=edge_sources,
                overall_confidence=sum(confidences) / len(confidences) if confidences else 0.0
            ))
            continue

        if len(current_path) <= request.max_depth:
            targets = db.get_targets(current_id, request.min_confidence, request.include_inferred)
            for target in targets:
                if target.id not in path_visited and target.regulation_type in request.regulation_type:
                    target_gene_obj = db.get_gene(target.id)
                    if target_gene_obj:
                        queue.append((
                            target.id,
                            current_path + [target_gene_obj],
                            regulations + [target.regulation_type],
                            confidences + [target.confidence],
                            edge_sources + [target.source_databases],
                            path_visited | {target.id}
                        ))

    paths.sort(key=lambda p: (-p.overall_confidence, len(p.genes)))
    return {"paths": paths}

@app.post("/api/v1/pathway/predict-cascade", response_model=CascadeResult)
async def predict_cascade(request: CascadeRequest):
    """
    Predict cascade effects of regulatory interventions

    Simulates regulatory cascade using simple propagation model
    """
    if not request.interventions:
        raise HTTPException(status_code=400, detail="At least one intervention required")

    gene = db.get_gene(request.target_gene_id)
    if not gene:
        raise HTTPException(status_code=404, detail="Gene not found")

    cascade_effects = []

    # Simple cascade simulation - in production, use ODE/boolean network model
    targets = db.get_targets(request.target_gene_id, min_confidence=0.5)
    for i, target in enumerate(targets[:5]):  # Limit to first 5 targets
        # Calculate cascade magnitude based on interventions
        magnitude = 1.0
        for intervention in request.interventions:
            if intervention.direction == "up":
                magnitude *= intervention.magnitude
            else:
                magnitude *= (2.0 - intervention.magnitude)
        
        cascade_effects.append(CascadeEffect(
            id=target.id,
            symbol=target.symbol,
            level=1,
            direction="up" if magnitude > 1.0 else "down",
            magnitude=abs(magnitude),
            confidence=target.confidence * 0.95  # Slightly reduced confidence for predictions
        ))
    
    return CascadeResult(
        cascade=cascade_effects,
        average_confidence=sum(e.confidence for e in cascade_effects) / len(cascade_effects) if cascade_effects else 0.0,
        affected_genes=len(cascade_effects)
    )

# ============= Orthology Endpoints =============

@app.get("/api/v1/genes/orthology/{gene_id}")
async def get_orthology(
    gene_id: str,
    species: Optional[str] = Query(None, description="Comma-separated species list")
):
    """
    Get orthologous genes and their regulatory networks across species
    """
    gene = db.get_gene(gene_id)
    if not gene:
        raise HTTPException(status_code=404, detail="Gene not found")
    
    # Parse species list
    target_species = species.split(",") if species else ["human", "arabidopsis", "rice"]
    
    result = {}
    for sp in target_species:
        if sp == gene.species:
            match = gene
        else:
            match = db.find_gene_by_symbol_species(gene.symbol, sp)

        if not match:
            result[sp] = {"found": False, "ortholog_symbol": gene.symbol, "regulators": [], "targets": []}
            continue

        regulators = db.get_regulators(match.id, min_confidence=0.5)
        targets = db.get_targets(match.id, min_confidence=0.5)

        result[sp] = {
            "found": True,
            "ortholog_symbol": match.symbol,
            "regulators": [
                {
                    "id": r.id,
                    "symbol": r.symbol,
                    "is_tf": r.is_tf,
                    "confidence": r.confidence,
                    "regulation_type": r.regulation_type
                }
                for r in regulators
            ],
            "targets": [
                {
                    "id": t.id,
                    "symbol": t.symbol,
                    "is_tf": t.is_tf,
                    "confidence": t.confidence,
                    "regulation_type": t.regulation_type
                }
                for t in targets
            ]
        }

    return result

# ============= Gene-set analysis: subgraph + GO enrichment =============

class SubgraphRequest(BaseModel):
    gene_ids: List[str]
    min_confidence: float = 0.0
    include_inferred: bool = True


@app.post("/api/v1/pathways/subgraph")
async def get_subgraph(request: SubgraphRequest):
    """Induced sub-network: genes in the set + interactions among them."""
    ids = list(dict.fromkeys(request.gene_ids))
    if not ids:
        return {"nodes": [], "edges": []}
    placeholders = ",".join("?" * len(ids))
    node_rows = db.conn.execute(
        f"SELECT id, symbol, name, species, is_tf FROM genes WHERE id IN ({placeholders})", ids
    ).fetchall()
    nodes = [{"id": r["id"], "symbol": r["symbol"], "name": r["name"],
              "species": r["species"], "is_tf": bool(r["is_tf"])} for r in node_rows]
    known = {n["id"] for n in nodes}
    edge_sql = (
        f"SELECT source_id, target_id, regulation_type, confidence, sources FROM interactions "
        f"WHERE source_id IN ({placeholders}) AND target_id IN ({placeholders}) AND confidence >= ?"
    )
    params = ids + ids + [request.min_confidence]
    edges = []
    for r in db.conn.execute(edge_sql, params).fetchall():
        if r["source_id"] not in known or r["target_id"] not in known:
            continue
        sources = json.loads(r["sources"])
        inferred = any(s.startswith("Inferred") for s in sources)
        if inferred and not request.include_inferred:
            continue
        edges.append({"source": r["source_id"], "target": r["target_id"],
                      "regulation_type": r["regulation_type"], "confidence": r["confidence"],
                      "source_databases": sources, "inferred": inferred})
    return {"nodes": nodes, "edges": edges}


class ExportRequest(BaseModel):
    """Export regulatory edges with sequence-fetch context."""
    gene_ids: List[str]
    min_confidence: float = 0.0
    include_inferred: bool = True
    signed_only: bool = False        # keep only activation/repression edges
    promoter_upstream: int = 2000    # bp 5' of TSS
    promoter_downstream: int = 500   # bp 3' of TSS
    format: str = "json"             # json | tsv


# activation/repression -> the positive/negative sign; everything else is unsigned.
_SIGN = {"activation": "positive", "repression": "negative"}


def _promoter_window(loc, upstream: int, downstream: int, chrom_len: Optional[int]):
    """Derive TSS and a strand-aware promoter window from a gene locus."""
    if not loc:
        return None, None, None
    strand = loc["strand"] or 0
    tss = loc["end"] if strand < 0 else loc["start"]
    if strand < 0:
        ws, we = tss - downstream, tss + upstream
    else:
        ws, we = tss - upstream, tss + downstream
    ws = max(0, ws)
    if chrom_len:
        we = min(chrom_len, we)
    return tss, ws, we


@app.post("/api/v1/export/edges")
async def export_edges(request: ExportRequest):
    """Regulatory edges annotated with sign, confidence, provenance, genomic
    coordinates for both partners, and derived promoter windows — i.e. everything
    needed to fetch promoter sequence downstream. Sequences themselves are not
    served (no assembly FASTA loaded); this emits window *coordinates*."""
    ids = list(dict.fromkeys(request.gene_ids))
    if not ids:
        return {"edges": [], "stats": {"edges": 0}, "params": request.dict()}
    ph = ",".join("?" * len(ids))

    genes = {
        r["id"]: r for r in db.conn.execute(
            f"SELECT id, symbol, name, species, is_tf FROM genes WHERE id IN ({ph})", ids
        ).fetchall()
    }
    locs = {
        r["gene_id"]: r for r in db.conn.execute(
            f"SELECT gene_id, species, chromosome, start, end, strand "
            f"FROM gene_locations WHERE gene_id IN ({ph})", ids
        ).fetchall()
    }
    chrom_len = {
        (r["species"], r["chromosome"]): r["length"]
        for r in db.conn.execute("SELECT species, chromosome, length FROM chromosomes").fetchall()
    }

    def side(gene_id, prefix):
        g, loc = genes.get(gene_id), locs.get(gene_id)
        clen = chrom_len.get((loc["species"], loc["chromosome"])) if loc else None
        tss, ws, we = _promoter_window(loc, request.promoter_upstream, request.promoter_downstream, clen)
        return {
            f"{prefix}_gene_id": gene_id,
            f"{prefix}_symbol": g["symbol"] if g else None,
            f"{prefix}_species": g["species"] if g else None,
            f"{prefix}_is_tf": bool(g["is_tf"]) if g else None,
            f"{prefix}_chromosome": loc["chromosome"] if loc else None,
            f"{prefix}_start": loc["start"] if loc else None,
            f"{prefix}_end": loc["end"] if loc else None,
            f"{prefix}_strand": loc["strand"] if loc else None,
            f"{prefix}_tss": tss,
            f"{prefix}_promoter_start": ws,
            f"{prefix}_promoter_end": we,
        }

    rows = db.conn.execute(
        f"SELECT source_id, target_id, regulation_type, confidence, sources, pmids "
        f"FROM interactions WHERE source_id IN ({ph}) AND target_id IN ({ph}) AND confidence >= ?",
        ids + ids + [request.min_confidence],
    ).fetchall()

    edges, complete = [], 0
    for r in rows:
        sources = json.loads(r["sources"])
        inferred = any(s.startswith("Inferred") for s in sources)
        if inferred and not request.include_inferred:
            continue
        sign = _SIGN.get(r["regulation_type"], "unsigned")
        if request.signed_only and sign == "unsigned":
            continue
        edge = {**side(r["source_id"], "source"), **side(r["target_id"], "target"),
                "regulation_type": r["regulation_type"], "sign": sign,
                "confidence": r["confidence"], "sources": sources,
                "pmids": json.loads(r["pmids"]) if r["pmids"] else [],
                "inferred": inferred}
        if edge["source_promoter_start"] is not None and edge["target_promoter_start"] is not None:
            complete += 1
        edges.append(edge)

    stats = {
        "edges": len(edges),
        "edges_with_complete_coordinates": complete,
        "signed": sum(1 for e in edges if e["sign"] != "unsigned"),
        "unsigned": sum(1 for e in edges if e["sign"] == "unsigned"),
        "inferred": sum(1 for e in edges if e["inferred"]),
    }

    if request.format == "tsv":
        cols = [k for k in (edges[0].keys() if edges else [])]
        lines = ["\t".join(cols)]
        for e in edges:
            lines.append("\t".join(
                ";".join(map(str, e[c])) if isinstance(e[c], list) else
                ("" if e[c] is None else str(e[c])) for c in cols))
        return PlainTextResponse("\n".join(lines), media_type="text/tab-separated-values")

    return {"edges": edges, "stats": stats, "params": request.dict()}


class EnrichmentRequest(BaseModel):
    gene_ids: List[str]
    species: Optional[str] = None
    max_terms: int = 40
    min_genes: int = 2


# Lazily-built per-species GO index: {species: (N, term_k, gene_terms)}.
_go_index: Dict[str, Any] = {}


def _go_index_for(species: str):
    if species not in _go_index:
        gene_terms: Dict[str, set] = {}
        term_k: Dict[str, int] = defaultdict(int)
        rows = db.conn.execute(
            "SELECT a.gene_id, a.go_id FROM go_annotations a JOIN genes g ON g.id = a.gene_id "
            "WHERE g.species = ?", (species,)
        ).fetchall()
        for gid, go_id in rows:
            gene_terms.setdefault(gid, set()).add(go_id)
        for terms in gene_terms.values():
            for t in terms:
                term_k[t] += 1
        _go_index[species] = (len(gene_terms), dict(term_k), gene_terms)
    return _go_index[species]


def _log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _hypergeom_sf(k: int, n: int, K: int, N: int) -> float:
    """P(X >= k) for drawing n from N with K successes (overrepresentation)."""
    logCNn = _log_choose(N, n)
    total = 0.0
    hi = min(n, K)
    for i in range(k, hi + 1):
        lp = _log_choose(K, i) + _log_choose(N - K, n - i) - logCNn
        if lp > -700:
            total += math.exp(lp)
    return min(total, 1.0)


@app.post("/api/v1/enrichment")
async def enrichment(request: EnrichmentRequest):
    """GO-term overrepresentation for a gene set (hypergeometric + BH FDR)."""
    ids = list(dict.fromkeys(request.gene_ids))
    species = request.species
    if not species and ids:
        row = db.conn.execute("SELECT species FROM genes WHERE id = ?", (ids[0],)).fetchone()
        species = row["species"] if row else None
    if not species:
        raise HTTPException(status_code=400, detail="Could not determine species")

    N, term_k, gene_terms = _go_index_for(species)
    if N == 0:
        return {"species": species, "background": 0, "study": 0, "results": []}

    study = [g for g in ids if g in gene_terms]
    n = len(study)
    study_k: Dict[str, int] = defaultdict(int)
    for g in study:
        for t in gene_terms[g]:
            study_k[t] += 1

    tested = []
    for go_id, k in study_k.items():
        if k < request.min_genes:
            continue
        K = term_k.get(go_id, 0)
        p = _hypergeom_sf(k, n, K, N)
        tested.append((go_id, k, K, p))

    # Benjamini–Hochberg FDR.
    tested.sort(key=lambda x: x[3])
    m = len(tested)
    results = []
    prev_q = 1.0
    for rank in range(m - 1, -1, -1):
        go_id, k, K, p = tested[rank]
        q = min(prev_q, p * m / (rank + 1))
        prev_q = q
        term = db.conn.execute("SELECT name, namespace FROM go_terms WHERE go_id = ?", (go_id,)).fetchone()
        results.append({
            "go_id": go_id, "name": term["name"] if term else go_id,
            "namespace": term["namespace"] if term else "",
            "study_count": k, "background_count": K, "p_value": p, "q_value": q,
        })
    results.sort(key=lambda r: r["p_value"])
    return {"species": species, "background": N, "study": n,
            "results": results[:request.max_terms]}


# ============= Genome / Synteny Endpoints =============

def _chromosome_sort_key(name: str):
    """Order chromosomes numerically (1,2,...) then alphabetically (X, Y, ...)."""
    return (0, int(name)) if name.isdigit() else (1, name)


@app.get("/api/v1/genome/species")
async def genome_species():
    """List species that have genome coordinate data, with their chromosomes."""
    cur = db.conn.execute
    species_rows = cur(
        "SELECT DISTINCT species FROM chromosomes ORDER BY species"
    ).fetchall()
    result = []
    for (species,) in species_rows:
        chroms = cur(
            """
            SELECT c.chromosome, c.length,
                   (SELECT COUNT(*) FROM gene_locations g
                    WHERE g.species = c.species AND g.chromosome = c.chromosome) AS gene_count
            FROM chromosomes c WHERE c.species = ?
            """,
            (species,)
        ).fetchall()
        chroms = sorted(
            [{"name": r["chromosome"], "length": r["length"], "gene_count": r["gene_count"]}
             for r in chroms],
            key=lambda c: _chromosome_sort_key(c["name"])
        )
        result.append({
            "species": species,
            "chromosomes": chroms,
            "gene_count": sum(c["gene_count"] for c in chroms),
        })
    return {"species": result}


@app.get("/api/v1/genome/orthologs")
async def genome_orthologs(
    species_a: str = Query(..., description="First species"),
    species_b: str = Query(..., description="Second species"),
):
    """Ortholog pairs between two species, joined to both genes' loci."""
    cur = db.conn.execute
    # Orthologs may be stored in either direction; normalize to (a -> b).
    rows = cur(
        """
        SELECT gene_a, gene_b, species_a, species_b, rel_type, score FROM orthologs
        WHERE (species_a = ? AND species_b = ?) OR (species_a = ? AND species_b = ?)
        """,
        (species_a, species_b, species_b, species_a)
    ).fetchall()

    def locus(gene_id):
        r = cur(
            "SELECT l.chromosome, l.start, l.end, g.symbol, g.is_tf "
            "FROM gene_locations l JOIN genes g ON g.id = l.gene_id WHERE l.gene_id = ?",
            (gene_id,)
        ).fetchone()
        if not r:
            return None
        return {"gene_id": gene_id, "symbol": r["symbol"], "chromosome": r["chromosome"],
                "start": r["start"], "end": r["end"], "is_tf": bool(r["is_tf"])}

    pairs = []
    for r in rows:
        # Orient so 'a' matches the requested species_a.
        if r["species_a"] == species_a:
            ga, gb = r["gene_a"], r["gene_b"]
        else:
            ga, gb = r["gene_b"], r["gene_a"]
        la, lb = locus(ga), locus(gb)
        if not la or not lb:
            continue
        pairs.append({
            "symbol": la["symbol"], "rel_type": r["rel_type"], "score": r["score"],
            "a": la, "b": lb,
        })
    return {"species_a": species_a, "species_b": species_b, "pairs": pairs}


@app.get("/api/v1/genome/{species}")
async def genome_detail(species: str):
    """Get all chromosomes for a species with their positioned genes."""
    cur = db.conn.execute
    chrom_rows = cur(
        "SELECT chromosome, length FROM chromosomes WHERE species = ?", (species,)
    ).fetchall()
    if not chrom_rows:
        raise HTTPException(status_code=404, detail="No genome data for species")

    gene_rows = cur(
        """
        SELECT l.gene_id, l.chromosome, l.start, l.end, l.strand,
               g.symbol, g.is_tf
        FROM gene_locations l JOIN genes g ON g.id = l.gene_id
        WHERE l.species = ?
        """,
        (species,)
    ).fetchall()

    by_chrom: Dict[str, list] = {r["chromosome"]: [] for r in chrom_rows}
    for r in gene_rows:
        by_chrom.setdefault(r["chromosome"], []).append({
            "id": r["gene_id"], "symbol": r["symbol"],
            "start": r["start"], "end": r["end"], "strand": r["strand"],
            "is_tf": bool(r["is_tf"]),
        })

    chromosomes = sorted(
        [{"name": r["chromosome"], "length": r["length"],
          "genes": sorted(by_chrom.get(r["chromosome"], []), key=lambda g: g["start"])}
         for r in chrom_rows],
        key=lambda c: _chromosome_sort_key(c["name"])
    )
    return {"species": species, "chromosomes": chromosomes}


# ============= Statistics Endpoints =============

@app.get("/api/v1/stats")
async def get_stats():
    """Get overall database statistics from live data"""
    cur = db.conn.execute
    total_genes = cur("SELECT COUNT(*) FROM genes").fetchone()[0]
    total_interactions = cur("SELECT COUNT(*) FROM interactions").fetchone()[0]
    species_list = [r[0] for r in cur("SELECT DISTINCT species FROM genes ORDER BY species").fetchall()]
    return {
        "species": len(species_list),
        "species_list": species_list,
        "genes": total_genes,
        "interactions": total_interactions,
        "databases": ["TRRUST", "PlantRegMap"],
        "version": "1.0.0"
    }

@app.get("/api/v1/stats/species/{species}")
async def get_species_stats(species: str):
    """Get species-specific statistics from live data"""
    cur = db.conn.execute
    genes = cur("SELECT COUNT(*) FROM genes WHERE species = ?", (species,)).fetchone()[0]
    if genes == 0:
        raise HTTPException(status_code=404, detail="Species not found")
    tfs = cur("SELECT COUNT(*) FROM genes WHERE species = ? AND is_tf = 1", (species,)).fetchone()[0]
    gene_ids = [r[0] for r in cur("SELECT id FROM genes WHERE species = ?", (species,)).fetchall()]
    placeholders = ",".join("?" * len(gene_ids))
    interactions = cur(
        f"SELECT COUNT(*) FROM interactions WHERE source_id IN ({placeholders})",
        gene_ids
    ).fetchone()[0]
    return {
        "species": species,
        "genes": genes,
        "transcription_factors": tfs,
        "interactions": interactions
    }

# ============= Error Handlers =============

@app.exception_handler(HTTPException)
async def http_exception_handler(request, exc):
    return {
        "error": exc.detail,
        "status_code": exc.status_code
    }

# ============= Startup Events =============

@app.on_event("startup")
async def startup_event():
    logger.info("GRN Atlas API starting up...")
    logger.info("Database initialized with mock data")
    logger.info("CORS enabled for frontend development")

@app.on_event("shutdown")
async def shutdown_event():
    logger.info("GRN Atlas API shutting down...")

# ============= Health Check Endpoint =============

@app.get("/health")
async def health_check():
    """Health check endpoint for monitoring"""
    return {
        "status": "healthy",
        "timestamp": datetime.utcnow().isoformat(),
        "version": "1.0.0"
    }

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        app,
        host="0.0.0.0",
        port=8000,
        log_level="info"
    )

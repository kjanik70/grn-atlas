"""
GRN Atlas FastAPI Backend
Complete example implementation with all required endpoints
"""

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse, JSONResponse

import provenance
import expression
import rnai
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
from datetime import datetime
from pathlib import Path as FilePath
from collections import defaultdict
import json
import math
import os
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
    # Friendliest display label + whether it is an inferred (ortholog) name rather
    # than a native symbol. Populated by friendly_label(); never fabricated.
    label: Optional[str] = None
    label_inferred: bool = False
    # provenance when `symbol` is a curated real name (e.g. 'UniProt', 'UniProt:homology')
    symbol_source: Optional[str] = None

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

from genelabels import friendly_label  # noqa: E402  (pure, unit-tested separately)


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
    label_inferred: bool = False   # True when `symbol` is an inferred ortholog name
    symbol_source: Optional[str] = None

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

class PerturbInterv(BaseModel):
    gene_id: str
    action: str = "ko"  # 'ko' (knock-out / down) or 'oe' (over-express / up)

class PerturbRequest(BaseModel):
    interventions: List[PerturbInterv]
    depth: int = 4
    min_confidence: float = 0.0
    include_inferred: bool = True
    min_effect: float = 0.05  # prune predicted effects weaker than this
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

DB_PATH = FilePath(os.environ.get("GRN_DB") or (FilePath(__file__).parent / "data" / "grn.sqlite3"))
DATA_DIR = DB_PATH.parent


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
        syns = [s for s in raw_syn.split("; ") if s] if raw_syn else None
        label, inferred = friendly_label(row["symbol"], row["id"], syns)
        sym_src = row["symbol_source"] if "symbol_source" in keys else None
        return Gene(
            id=row["id"],
            symbol=row["symbol"],
            name=row["name"],
            species=row["species"],
            is_tf=bool(row["is_tf"]),
            gene_type=row["gene_type"],
            synonyms=syns,
            label=label,
            label_inferred=inferred,
            symbol_source=sym_src,
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
        keys = row.keys()
        raw_syn = row["synonyms"] if "synonyms" in keys else None
        syns = [s for s in raw_syn.split("; ") if s] if raw_syn else None
        label, label_inf = friendly_label(row["symbol"], row["id"], syns)
        return GeneInteraction(
            id=row["id"], symbol=label, name=row["name"], species=row["species"],
            is_tf=bool(row["is_tf"]), confidence=row["confidence"],
            regulation_type=row["regulation_type"], source_databases=sources,
            pmids=pmids, inferred=any(s.startswith("Inferred") for s in sources),
            label_inferred=label_inf,
            symbol_source=row["symbol_source"] if "symbol_source" in keys else None,
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
            path_genes = [PathGene(id=g.id, symbol=(getattr(g, "label", None) or g.symbol), name=g.name)
                          for g in current_path]
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

_EDGE_SIGN = {"activation": 1, "repression": -1}
_PERTURB_DECAY = 0.7  # per-level attenuation of predicted effect magnitude


@app.post("/api/v1/perturb")
async def perturb(request: PerturbRequest):
    """Predict the qualitative direction of change of downstream genes after a
    set of TF perturbations, by propagating signs along the regulatory network.

    This is a *predicted*, qualitative model (signed-path propagation over measured/
    inferred edges), NOT a quantitative simulation. Effect sign = product of edge
    signs (activation +1, repression -1) times the intervention sign (ko -1, oe +1);
    magnitude is a confidence-weighted, depth-damped heuristic. An edge with an
    unsigned regulation type makes the downstream direction 'unknown'.
    """
    if not request.interventions:
        raise HTTPException(status_code=400, detail="At least one intervention required")

    seeds = {}
    for iv in request.interventions:
        g = db.get_gene(iv.gene_id)
        if not g:
            raise HTTPException(status_code=404, detail=f"Gene not found: {iv.gene_id}")
        seeds[iv.gene_id] = (-1 if iv.action == "ko" else 1, g.symbol)

    # best[gene] = dict(magnitude, sign, unknown, level, path, uses_inferred)
    best: Dict[str, Dict[str, Any]] = {}
    depth = max(1, min(request.depth, 6))
    # frontier: (gene_id, sign, magnitude, level, path_symbols, uses_inferred, unknown)
    frontier = [(gid, s, 1.0, 0, [sym], False, False) for gid, (s, sym) in seeds.items()]

    while frontier:
        gid, sign, mag, level, path, inf, unknown = frontier.pop()
        if level >= depth:
            continue
        for t in db.get_targets(gid, min_confidence=request.min_confidence,
                                include_inferred=request.include_inferred):
            esign = _EDGE_SIGN.get(t.regulation_type, 0)
            n_unknown = unknown or esign == 0
            n_sign = sign * (esign if esign else 1)
            n_mag = mag * max(t.confidence, 0.01) * _PERTURB_DECAY
            n_inf = inf or any(str(s).startswith("Inferred") for s in t.source_databases)
            if n_mag < request.min_effect:
                continue
            prev = best.get(t.id)
            if t.id not in seeds and (prev is None or n_mag > prev["magnitude"]):
                best[t.id] = {"symbol": t.symbol, "magnitude": round(n_mag, 4),
                              "sign": n_sign, "unknown": n_unknown, "level": level + 1,
                              "path": path + [t.symbol], "uses_inferred": n_inf}
            # keep expanding as long as this route is the strongest seen for t
            if prev is None or n_mag > prev["magnitude"]:
                frontier.append((t.id, n_sign, n_mag, level + 1,
                                 path + [t.symbol], n_inf, n_unknown))

    def direction(e):
        if e["unknown"]:
            return "unknown"
        return "up" if e["sign"] > 0 else "down"

    effects = [{"gene_id": gid, "symbol": e["symbol"], "predicted_direction": direction(e),
                "magnitude": e["magnitude"], "level": e["level"], "path": e["path"],
                "uses_inferred": e["uses_inferred"]}
               for gid, e in best.items()]
    effects.sort(key=lambda e: e["magnitude"], reverse=True)

    return {
        "interventions": [{"gene_id": gid, "symbol": sym,
                           "action": "ko" if s < 0 else "oe"}
                          for gid, (s, sym) in seeds.items()],
        "effects": effects,
        "stats": {
            "affected": len(effects),
            "up": sum(1 for e in effects if e["predicted_direction"] == "up"),
            "down": sum(1 for e in effects if e["predicted_direction"] == "down"),
            "unknown": sum(1 for e in effects if e["predicted_direction"] == "unknown"),
            "uses_inferred": any(e["uses_inferred"] for e in effects),
        },
        "note": "Predicted qualitative directions from signed-path propagation, "
                "not a quantitative simulation.",
    }


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

# ============= Provenance & citations =============

@app.get("/api/v1/provenance")
async def get_provenance():
    """Machine-readable data-source versions, citations, and analysis methods —
    for reproducibility and to cite the atlas in a publication."""
    return provenance.manifest()


@app.get("/api/v1/provenance/freshness")
async def get_freshness():
    """Data-currency audit: each source's loaded version vs the latest available
    release, so users can judge how current the atlas data is (#6)."""
    return provenance.freshness()


@app.get("/api/v1/citations.bib")
async def get_citations():
    """BibTeX for every data source the atlas integrates."""
    return PlainTextResponse(provenance.bibtex(), media_type="application/x-bibtex")


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
        f"SELECT id, symbol, name, species, is_tf, synonyms FROM genes WHERE id IN ({placeholders})", ids
    ).fetchall()
    nodes = [{"id": r["id"], "symbol": r["symbol"], "name": r["name"],
              "species": r["species"], "is_tf": bool(r["is_tf"]),
              "synonyms": r["synonyms"].split("; ") if r["synonyms"] else []} for r in node_rows]
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
    promoter_upstream: int = 2000    # bp 5' of TSS (derived-window fallback)
    promoter_downstream: int = 500   # bp 3' of TSS
    include_sequence_context: bool = False   # attach ingested windows + motif sites
    window_types: List[str] = ["promoter"]   # promoter | gene_body | atac
    max_site_pvalue: float = 1e-4
    format: str = "json"             # json | tsv


# Assembly the atlas gene_locations coordinates live on, per species (for tagging).
_ATLAS_ASSEMBLY = {
    "human": "GRCh38", "mouse": "GRCm39", "arabidopsis": "TAIR10",
    "tomato": "SL2.50", "petunia": "Peaxi162_HiC",
}


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
            f"{prefix}_coord_assembly": _ATLAS_ASSEMBLY.get(g["species"]) if g else None,
            f"{prefix}_coord_system": "GFF1",
        }

    # Sequence context (Path B): batch-load ingested windows + motif sites for
    # the target genes, via the annotation-version crosswalk. Tables always
    # exist but may be empty (no ingestion bundle yet) -> no context attached.
    atlas2ext, xrelation, windows_by_ext, hits_by_ext = {}, {}, {}, {}
    if request.include_sequence_context:
        for x in db.conn.execute(
            f"SELECT atlas_gene_id, ext_gene_id, ext_assembly, relation "
            f"FROM gene_id_crosswalk WHERE atlas_gene_id IN ({ph})", ids
        ).fetchall():
            atlas2ext[x["atlas_gene_id"]] = (x["ext_gene_id"], x["ext_assembly"])
            xrelation[x["atlas_gene_id"]] = x["relation"]
        ext_ids = [e for (e, _) in atlas2ext.values()]
        wt = request.window_types or []
        if ext_ids and wt:
            eph, wtph = ",".join("?" * len(ext_ids)), ",".join("?" * len(wt))
            for w in db.conn.execute(
                f"SELECT ext_gene_id, window_type, chromosome, start, end, strand "
                f"FROM gene_windows WHERE ext_gene_id IN ({eph}) AND window_type IN ({wtph})",
                ext_ids + wt
            ).fetchall():
                windows_by_ext.setdefault(w["ext_gene_id"], []).append(
                    {"window_type": w["window_type"], "chromosome": w["chromosome"],
                     "start": w["start"], "end": w["end"], "strand": w["strand"]})
            for h in db.conn.execute(
                f"SELECT h.ext_gene_id, h.motif_id, h.window_type, h.chromosome, h.start, h.end, "
                f"h.strand, h.score, h.p_value, h.tier, h.site_confidence, m.tf_gene_id, m.tf_symbol "
                f"FROM motif_hits h JOIN motifs m ON m.motif_id = h.motif_id "
                f"WHERE h.ext_gene_id IN ({eph}) AND h.window_type IN ({wtph}) AND h.p_value <= ?",
                ext_ids + wt + [request.max_site_pvalue]
            ).fetchall():
                hits_by_ext.setdefault(h["ext_gene_id"], []).append(dict(h))

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
        # Attach ingested windows + the motif sites that support THIS edge
        # (sites in the target's promoter bound by a motif of the source TF).
        if request.include_sequence_context:
            ext = atlas2ext.get(r["target_id"])
            if ext:
                ext_id, assembly = ext
                # Prefer a real curated symbol (e.g. AN2, != locus id); otherwise
                # keep the scan-time synonym label (e.g. HYH, PIF4).
                tf_label = (edge["source_symbol"]
                            if edge["source_symbol"] and edge["source_symbol"] != r["source_id"]
                            else h["tf_symbol"])
                sites = [
                    {"motif_id": h["motif_id"], "tf_symbol": tf_label,
                     "window_type": h["window_type"], "chromosome": h["chromosome"],
                     "start": h["start"], "end": h["end"], "strand": h["strand"],
                     "score": h["score"], "p_value": h["p_value"],
                     "tier": h["tier"], "site_confidence": h["site_confidence"]}
                    for h in hits_by_ext.get(ext_id, []) if h["tf_gene_id"] == r["source_id"]
                ]
                edge["sequence_context"] = {
                    "assembly": assembly, "coord_system": "BED0",
                    "target_ext_gene_id": ext_id,
                    "crosswalk_relation": xrelation.get(r["target_id"]),
                    "target_windows": windows_by_ext.get(ext_id, []),
                    "supporting_sites": sites,
                }
        edges.append(edge)

    stats = {
        "edges": len(edges),
        "edges_with_complete_coordinates": complete,
        "signed": sum(1 for e in edges if e["sign"] != "unsigned"),
        "unsigned": sum(1 for e in edges if e["sign"] == "unsigned"),
        "inferred": sum(1 for e in edges if e["inferred"]),
        "edges_with_sequence_context": sum(1 for e in edges if e.get("sequence_context")),
        "edges_with_supporting_sites": sum(
            1 for e in edges if e.get("sequence_context", {}).get("supporting_sites")),
    }

    prov = provenance.manifest()

    if request.format == "tsv":
        # Flat table only; the nested sequence_context is JSON-only.
        cols = [k for k in (edges[0].keys() if edges else []) if k != "sequence_context"]
        header = [f"# GRN Atlas export v{prov['atlas_version']} — generated {prov['generated']}",
                  f"# promoter window: {prov['methods']['promoter_window']}",
                  f"# inferred edges: {prov['methods']['inferred_edges']}",
                  f"# motif sites: {prov['methods']['motif_scan']}",
                  "# full provenance + citations: GET /api/v1/provenance , /api/v1/citations.bib"]
        lines = header + ["\t".join(cols)]
        for e in edges:
            lines.append("\t".join(
                ";".join(map(str, e[c])) if isinstance(e[c], list) else
                ("" if e[c] is None else str(e[c])) for c in cols))
        return PlainTextResponse("\n".join(lines), media_type="text/tab-separated-values")

    return {"edges": edges, "stats": stats, "params": request.dict(), "provenance": prov}


class ConservationRequest(BaseModel):
    gene_ids: List[str]
    species_b: str                       # compare against this species
    min_confidence: float = 0.0
    include_inferred: bool = True


@app.post("/api/v1/conservation")
async def conservation(request: ConservationRequest):
    """Cross-species conservation of regulatory edges: for each edge among the
    given genes (species A), is the corresponding edge (via orthologs) present in
    species B? Joins the ortholog map with both species' networks."""
    ids = list(dict.fromkeys(request.gene_ids))
    if not ids:
        return {"species_b": request.species_b, "edges": [], "stats": {}}
    ph = ",".join("?" * len(ids))
    inferred_a = "" if request.include_inferred else " AND i.sources NOT LIKE '%Inferred%'"

    a_edges = db.conn.execute(
        f"SELECT source_id, target_id, regulation_type, confidence FROM interactions i "
        f"WHERE source_id IN ({ph}) AND target_id IN ({ph}) AND confidence >= ?{inferred_a}",
        ids + ids + [request.min_confidence]).fetchall()
    if not a_edges:
        return {"species_b": request.species_b, "edges": [], "stats": {"edges": 0, "conserved": 0}}

    # orthologs of the involved genes in species B (either orientation)
    orth = defaultdict(set)
    rows = db.conn.execute(
        f"SELECT gene_a, gene_b, species_a, species_b FROM orthologs "
        f"WHERE (gene_a IN ({ph}) AND species_b = ?) OR (gene_b IN ({ph}) AND species_a = ?)",
        ids + [request.species_b] + ids + [request.species_b]).fetchall()
    for r in rows:
        if r["species_b"] == request.species_b:
            orth[r["gene_a"]].add(r["gene_b"])
        else:
            orth[r["gene_b"]].add(r["gene_a"])

    # B-side edges among all candidate orthologs
    b_ids = sorted({g for s in orth.values() for g in s})
    b_edge = {}
    if b_ids:
        bph = ",".join("?" * len(b_ids))
        inferred_b = "" if request.include_inferred else " AND i.sources NOT LIKE '%Inferred%'"
        for r in db.conn.execute(
            f"SELECT source_id, target_id, regulation_type, confidence FROM interactions i "
            f"WHERE source_id IN ({bph}) AND target_id IN ({bph}){inferred_b}", b_ids + b_ids).fetchall():
            b_edge[(r["source_id"], r["target_id"])] = (r["regulation_type"], r["confidence"])

    def sym(gene_id):
        r = db.conn.execute("SELECT symbol FROM genes WHERE id = ?", (gene_id,)).fetchone()
        return r["symbol"] if r else gene_id

    out, n_cons = [], 0
    for e in a_edges:
        matches = []
        for a in orth.get(e["source_id"], ()):
            for b in orth.get(e["target_id"], ()):
                if (a, b) in b_edge:
                    reg, conf = b_edge[(a, b)]
                    matches.append({"source_ortholog": a, "source_ortholog_symbol": sym(a),
                                    "target_ortholog": b, "target_ortholog_symbol": sym(b),
                                    "regulation_type": reg, "confidence": conf})
        conserved = bool(matches)
        n_cons += conserved
        out.append({
            "source_gene_id": e["source_id"], "source_symbol": sym(e["source_id"]),
            "target_gene_id": e["target_id"], "target_symbol": sym(e["target_id"]),
            "regulation_type": e["regulation_type"], "confidence": e["confidence"],
            "conserved": conserved,
            "source_has_ortholog": bool(orth.get(e["source_id"])),
            "target_has_ortholog": bool(orth.get(e["target_id"])),
            "b_edges": matches,
        })
    return {"species_b": request.species_b, "edges": out,
            "stats": {"edges": len(out), "conserved": n_cons,
                      "both_orthologs": sum(1 for e in out if e["source_has_ortholog"] and e["target_has_ortholog"])}}


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


_pathway_index: Dict[str, Any] = {}


def _pathway_index_for(species: str):
    """(N annotated genes, pathway_k background counts, gene_pathways) for a species."""
    if species not in _pathway_index:
        gene_pathways: Dict[str, set] = {}
        pw_k: Dict[str, int] = defaultdict(int)
        try:
            rows = db.conn.execute(
                "SELECT a.gene_id, a.pathway_id FROM pathway_annotations a "
                "JOIN genes g ON g.id = a.gene_id WHERE g.species = ?", (species,)).fetchall()
        except sqlite3.OperationalError:
            rows = []  # pathway tables not loaded
        for gid, pid in rows:
            gene_pathways.setdefault(gid, set()).add(pid)
        for pids in gene_pathways.values():
            for p in pids:
                pw_k[p] += 1
        _pathway_index[species] = (len(gene_pathways), dict(pw_k), gene_pathways)
    return _pathway_index[species]


@app.post("/api/v1/pathway_enrichment")
async def pathway_enrichment(request: EnrichmentRequest):
    """Reactome pathway over-representation for a gene set (hypergeometric + BH FDR).

    Curated pathway membership (Plant Reactome). Currently plant species
    (arabidopsis, tomato); returns an empty result + note otherwise.
    """
    ids = list(dict.fromkeys(request.gene_ids))
    species = request.species
    if not species and ids:
        row = db.conn.execute("SELECT species FROM genes WHERE id = ?", (ids[0],)).fetchone()
        species = row["species"] if row else None

    N, pw_k, gene_pathways = _pathway_index_for(species)
    if N == 0:
        return {"species": species, "background": 0, "study": 0, "results": [],
                "note": "pathway annotations available for arabidopsis and tomato"}

    study = [g for g in ids if g in gene_pathways]
    n = len(study)
    study_k: Dict[str, int] = defaultdict(int)
    for g in study:
        for p in gene_pathways[g]:
            study_k[p] += 1

    tested = []
    for pid, k in study_k.items():
        if k < request.min_genes:
            continue
        tested.append((pid, k, pw_k.get(pid, 0), _hypergeom_sf(k, n, pw_k.get(pid, 0), N)))
    tested.sort(key=lambda x: x[3])
    m = len(tested)
    results, prev_q = [], 1.0
    for rank in range(m - 1, -1, -1):
        pid, k, K, p = tested[rank]
        q = min(prev_q, p * m / (rank + 1))
        prev_q = q
        term = db.conn.execute("SELECT name, source FROM pathways WHERE pathway_id = ?", (pid,)).fetchone()
        results.append({"pathway_id": pid, "name": term["name"] if term else pid,
                        "source": term["source"] if term else "",
                        "study_count": k, "background_count": K, "p_value": p, "q_value": q})
    results.sort(key=lambda r: r["p_value"])
    return {"species": species, "background": N, "study": n,
            "results": results[:request.max_terms]}


_trait_index: Dict[str, Any] = {}


def _trait_index_for(species: str):
    """(N annotated genes, trait_k background, gene_traits) — GWAS trait associations."""
    if species not in _trait_index:
        gene_traits: Dict[str, set] = {}
        trait_k: Dict[str, int] = defaultdict(int)
        try:
            rows = db.conn.execute(
                "SELECT a.gene_id, a.trait FROM trait_associations a "
                "JOIN genes g ON g.id = a.gene_id WHERE g.species = ?", (species,)).fetchall()
        except sqlite3.OperationalError:
            rows = []
        for gid, trait in rows:
            gene_traits.setdefault(gid, set()).add(trait)
        for traits in gene_traits.values():
            for t in traits:
                trait_k[t] += 1
        _trait_index[species] = (len(gene_traits), dict(trait_k), gene_traits)
    return _trait_index[species]


@app.get("/api/v1/traits/{gene_id}")
async def gene_traits(gene_id: str):
    """GWAS Catalog trait associations for a gene (statistical, not mechanistic)."""
    try:
        rows = db.conn.execute(
            "SELECT trait, pubmed_id, source FROM trait_associations WHERE gene_id = ? "
            "ORDER BY trait", (gene_id,)).fetchall()
    except sqlite3.OperationalError:
        rows = []
    return {"gene_id": gene_id,
            "traits": [{"trait": r["trait"], "pubmed_id": r["pubmed_id"], "source": r["source"]}
                       for r in rows],
            "note": "GWAS associations (SNP→mapped gene→trait), not mechanistic regulation."}


@app.post("/api/v1/trait_enrichment")
async def trait_enrichment(request: EnrichmentRequest):
    """GWAS trait over-representation for a gene set (hypergeometric + BH FDR).

    Answers "which phenotypes are this set / this regulon's targets enriched for".
    Human only (GWAS Catalog); returns [] + note otherwise.
    """
    ids = list(dict.fromkeys(request.gene_ids))
    species = request.species
    if not species and ids:
        row = db.conn.execute("SELECT species FROM genes WHERE id = ?", (ids[0],)).fetchone()
        species = row["species"] if row else None

    N, trait_k, gene_traits = _trait_index_for(species)
    if N == 0:
        have = [r[0] for r in db.conn.execute(
            "SELECT DISTINCT g.species FROM trait_associations a JOIN genes g ON g.id=a.gene_id")]
        return {"species": species, "background": 0, "study": 0, "results": [],
                "note": f"trait associations available for: {', '.join(sorted(have)) or 'none'}"}

    study = [g for g in ids if g in gene_traits]
    n = len(study)
    study_k: Dict[str, int] = defaultdict(int)
    for g in study:
        for t in gene_traits[g]:
            study_k[t] += 1

    tested = []
    for trait, k in study_k.items():
        if k < request.min_genes:
            continue
        tested.append((trait, k, trait_k.get(trait, 0), _hypergeom_sf(k, n, trait_k.get(trait, 0), N)))
    tested.sort(key=lambda x: x[3])
    m = len(tested)
    results, prev_q = [], 1.0
    for rank in range(m - 1, -1, -1):
        trait, k, K, p = tested[rank]
        q = min(prev_q, p * m / (rank + 1))
        prev_q = q
        results.append({"trait": trait, "study_count": k, "background_count": K,
                        "p_value": p, "q_value": q})
    results.sort(key=lambda r: r["p_value"])
    return {"species": species, "background": N, "study": n,
            "results": results[:request.max_terms]}


# ---- Motif enrichment: which TFs' predicted binding sites are over-represented
# in a gene set's promoters, vs the scanned-promoter background ----

_ASSEMBLY_OF = {"tomato": "SL4.0", "petunia": "Peaxi162v1.6.2", "arabidopsis": "TAIR10"}
_motif_index: Dict[str, Any] = {}


def _motif_index_for(species: str):
    """(N, tf_bg_count, gene_tfs, atlas2ext, tf_symbol) for a species with a scan."""
    if species not in _motif_index:
        assembly = _ASSEMBLY_OF.get(species)
        gene_tfs: Dict[str, set] = defaultdict(set)      # ext gene -> {tf_gene_id}
        tf_symbol: Dict[str, str] = {}
        if assembly:
            for r in db.conn.execute(
                "SELECT DISTINCT h.ext_gene_id, m.tf_gene_id, g.symbol "
                "FROM motif_hits h JOIN motifs m ON m.motif_id = h.motif_id "
                "JOIN genes g ON g.id = m.tf_gene_id WHERE h.assembly = ?", (assembly,)
            ).fetchall():
                gene_tfs[r["ext_gene_id"]].add(r["tf_gene_id"])
                tf_symbol[r["tf_gene_id"]] = r["symbol"]
        tf_bg = defaultdict(int)
        for tfs in gene_tfs.values():
            for tf in tfs:
                tf_bg[tf] += 1
        atlas2ext = {r["atlas_gene_id"]: r["ext_gene_id"]
                     for r in db.conn.execute("SELECT atlas_gene_id, ext_gene_id FROM gene_id_crosswalk")}
        _motif_index[species] = (len(gene_tfs), dict(tf_bg), gene_tfs, atlas2ext, tf_symbol)
    return _motif_index[species]


@app.post("/api/v1/motif_enrichment")
async def motif_enrichment(request: EnrichmentRequest):
    """Hypergeometric over-representation of each TF's predicted binding sites in
    the promoters of a gene set, vs the scanned-promoter background (BH FDR)."""
    ids = list(dict.fromkeys(request.gene_ids))
    species = request.species
    if not species and ids:
        row = db.conn.execute("SELECT species FROM genes WHERE id = ?", (ids[0],)).fetchone()
        species = row["species"] if row else None
    if species not in _ASSEMBLY_OF:
        return {"species": species, "background": 0, "study": 0, "results": [],
                "note": "motif scan available for tomato and petunia only"}

    N, tf_bg, gene_tfs, atlas2ext, tf_symbol = _motif_index_for(species)
    if N == 0:
        return {"species": species, "background": 0, "study": 0, "results": []}

    study = [atlas2ext.get(g) for g in ids]
    study = [e for e in study if e in gene_tfs]
    n = len(study)
    study_k = defaultdict(int)
    for e in study:
        for tf in gene_tfs[e]:
            study_k[tf] += 1

    tested = []
    for tf, k in study_k.items():
        if k < request.min_genes:
            continue
        tested.append((tf, k, tf_bg.get(tf, 0), _hypergeom_sf(k, n, tf_bg.get(tf, 0), N)))
    tested.sort(key=lambda x: x[3])
    m = len(tested)
    results, prev_q = [], 1.0
    for rank in range(m - 1, -1, -1):
        tf, k, K, p = tested[rank]
        q = min(prev_q, p * m / (rank + 1))
        prev_q = q
        results.append({"tf_gene_id": tf, "tf_symbol": tf_symbol.get(tf, tf),
                        "study_count": k, "background_count": K, "p_value": p, "q_value": q})
    results.sort(key=lambda r: r["p_value"])
    return {"species": species, "background": N, "study": n,
            "results": results[:request.max_terms]}


# ============= Expression + co-expression (petunia) =============

class CoexpRequest(BaseModel):
    gene_id: str
    top: int = 25
    min_abs_r: float = 0.7
    min_expr: float = 5.0
    tf_only: bool = False  # restrict partners to transcription factors (regulator candidates)


def _gene_meta(gene_ids):
    """friendly label + symbol + is_tf for a list of gene ids, in one query."""
    if not gene_ids:
        return {}
    ph = ",".join("?" * len(gene_ids))
    out = {}
    for r in db.conn.execute(
            f"SELECT id, symbol, is_tf, synonyms FROM genes WHERE id IN ({ph})", list(gene_ids)):
        syns = [s for s in (r["synonyms"] or "").split("; ") if s]
        label, inferred = friendly_label(r["symbol"], r["id"], syns)
        out[r["id"]] = {"symbol": label, "is_tf": bool(r["is_tf"]), "label_inferred": inferred}
    return out


def _species_of(gene_id: str) -> Optional[str]:
    row = db.conn.execute("SELECT species FROM genes WHERE id = ?", (gene_id,)).fetchone()
    return row["species"] if row else None


@app.get("/api/v1/expression/{gene_id}")
async def gene_expression(gene_id: str):
    """Per-sample TPM profile for a gene across its species' RNA-seq panel (#1).

    Predicted from subsampled public reads (kallisto vs PLAZA CDS); relative,
    not absolute. Available where an expression panel is built (petunia, tomato).
    """
    species = _species_of(gene_id)
    mx = expression.get_matrix(species) if species else None
    if mx is None:
        return {"gene_id": gene_id, "available": False,
                "note": f"no expression panel for {species or 'this gene'}"}
    prof = mx.profile(gene_id)
    if prof is None:
        return {"gene_id": gene_id, "available": False,
                "note": "no expression for this gene in the panel"}
    meta = _gene_meta([gene_id]).get(gene_id, {})
    return {"available": True, "species": species, "symbol": meta.get("symbol"),
            "is_tf": meta.get("is_tf"), "matrix_meta": mx.meta, **prof}


@app.post("/api/v1/coexpression")
async def coexpression(request: CoexpRequest):
    """Predicted co-expression partners of a gene across its species' panel (#2).

    Pearson correlation on log2(TPM+1). This is an inferred, UNDIRECTED association
    (labelled Inferred:Expression) — not measured regulation and not a causal
    direction. With tf_only=True, restricts partners to TFs (candidate regulators).
    """
    species = _species_of(request.gene_id)
    mx = expression.get_matrix(species) if species else None
    if mx is None:
        return {"gene_id": request.gene_id, "available": False,
                "results": [], "note": f"no expression panel for {species or 'this gene'}"}

    candidates = None
    if request.tf_only:
        candidates = [r["id"] for r in db.conn.execute(
            "SELECT id FROM genes WHERE is_tf = 1 AND species = ?", (species,)).fetchall()]
    hits = mx.coexpressed(request.gene_id, top=request.top, min_abs_r=request.min_abs_r,
                          min_expr=request.min_expr, candidates=candidates)
    meta = _gene_meta([h["gene_id"] for h in hits])
    for h in hits:
        m = meta.get(h["gene_id"], {})
        h["symbol"] = m.get("symbol", h["gene_id"])
        h["is_tf"] = m.get("is_tf", False)
        h["source"] = "Inferred:Expression"
    return {"gene_id": request.gene_id, "available": True,
            "n_samples": mx.n, "results": hits,
            "note": "Predicted co-expression (undirected), not measured regulation."}


# ============= dsRNA / RNAi design + off-target analysis =============

def clean_dsrna(seq: str, k: int) -> str:
    """Validate a pasted dsRNA via rnai.validate_dsrna; surface errors as HTTP 400."""
    try:
        return rnai.validate_dsrna(seq, k)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


class DsRnaRequest(BaseModel):
    sequence: Optional[str] = None       # the dsRNA to test (analyze mode)
    target_gene_id: Optional[str] = None  # intended target (design mode if no sequence)
    species: Optional[str] = None
    k: int = Field(21, ge=15, le=28)     # siRNA length
    max_off_targets: int = Field(50, ge=1, le=500)
    design_window: int = Field(250, ge=40, le=1000)
    predict_effect: bool = True          # summarise downstream effect of the silenced genes


@app.post("/api/v1/dsrna")
async def dsrna_analysis(request: DsRnaRequest):
    """Predict which genes a dsRNA would silence (on-target + off-target) and, in
    design mode, propose the most specific dsRNA window for a target gene.

    PREDICTED, NOT MEASURED: exact k-mer (siRNA) matching is a specificity heuristic;
    real RNAi knockdown also depends on dicing, delivery/SIGS uptake, target
    accessibility, and plant transitivity/amplification (not modelled).
    """
    species = request.species
    if not species and request.target_gene_id:
        species = _species_of(request.target_gene_id)
    if not species:
        raise HTTPException(status_code=400, detail="species (or a target_gene_id) required")

    transcripts = rnai.get_transcripts(species, DATA_DIR)
    if transcripts is None:
        return {"species": species, "available": False,
                "note": f"no transcript store for {species} (add transcripts_{species}.fasta.gz)"}

    design = None
    dsrna = clean_dsrna(request.sequence, request.k) if request.sequence else None
    if not dsrna:
        if not request.target_gene_id:
            raise HTTPException(status_code=400, detail="provide a sequence or a target_gene_id")
        design = rnai.design(request.target_gene_id, transcripts, k=request.k,
                             window=request.design_window)
        if "error" in design:
            raise HTTPException(status_code=404, detail=design["error"])
        dsrna = design["sequence"]

    result = rnai.scan(dsrna, transcripts, k=request.k, target_gene=request.target_gene_id)

    # annotate genes with symbol + tissue expression context
    shown = result["off_targets"][:request.max_off_targets]
    ids = [g for g in [request.target_gene_id] if g] + [o["gene_id"] for o in shown]
    meta = _gene_meta(ids)
    emx = expression.get_matrix(species) if species in set(expression.species_with_expression()) else None
    for o in shown:
        m = meta.get(o["gene_id"], {})
        o["symbol"] = m.get("symbol", o["gene_id"])
        o["is_tf"] = m.get("is_tf", False)
        o["label_inferred"] = m.get("label_inferred", False)
        prof = emx.profile(o["gene_id"]) if emx else None
        o["mean_tpm"] = prof["mean_tpm"] if prof else None
    result["off_targets"] = shown

    if request.target_gene_id:
        tm = meta.get(request.target_gene_id, {})
        tprof = emx.profile(request.target_gene_id) if emx else None
        result["on_target"] = {"gene_id": request.target_gene_id, "symbol": tm.get("symbol"),
                               "label_inferred": tm.get("label_inferred", False),
                               "sites": result["on_target_sites"],
                               "mean_tpm": tprof["mean_tpm"] if tprof else None}

    # optional: predicted downstream effect of the silenced set (feeds the perturb model)
    effect = None
    # off-targets are scanned against the full transcriptome; only genes present in the
    # network can be propagated, so filter before perturbation.
    silenced = result["silenced_genes"][:10]
    if silenced:
        ph = ",".join("?" * len(silenced))
        in_net = {r[0] for r in db.conn.execute(
            f"SELECT id FROM genes WHERE id IN ({ph})", silenced)}
        silenced = [g for g in silenced if g in in_net]
    if request.predict_effect and silenced:
        pr = await perturb(PerturbRequest(
            interventions=[PerturbInterv(gene_id=g, action="ko") for g in silenced], depth=3))
        effect = {"affected": pr["stats"]["affected"], "up": pr["stats"]["up"],
                  "down": pr["stats"]["down"], "unknown": pr["stats"]["unknown"],
                  "top": pr["effects"][:8]}

    return {"species": species, "available": True, "mode": "design" if design else "analyze",
            "design": design, "predicted_effect": effect, **result,
            "note": "Predicted silencing from exact siRNA k-mer matches (both strands); "
                    "not a guarantee of knockdown. Feed 'silenced_genes' to /perturb for "
                    "the full downstream cascade."}


_MAX_SCREEN_GENES = 300


class DsRnaScreenRequest(BaseModel):
    gene_ids: Optional[List[str]] = None
    pathway_id: Optional[str] = None     # screen every gene in a pathway
    species: Optional[str] = None
    k: int = Field(21, ge=15, le=28)
    design_window: int = Field(250, ge=40, le=1000)
    predict_effect: bool = True


@app.post("/api/v1/dsrna/screen")
async def dsrna_screen(request: DsRnaScreenRequest):
    """Batch dsRNA-designability screen across a gene set or a whole pathway: for each
    gene, the off-target burden of its most-specific window (one transcriptome pass),
    ranked cleanest-first — so you can choose the best RNAi target(s) to alter a pathway.
    Optionally reports the predicted downstream effect of silencing the whole set.
    """
    species = request.species
    genes = list(request.gene_ids or [])
    if request.pathway_id:
        rows = db.conn.execute(
            "SELECT a.gene_id FROM pathway_annotations a JOIN genes g ON g.id=a.gene_id "
            "WHERE a.pathway_id=?" + ("" if not species else " AND g.species=?"),
            (request.pathway_id, species) if species else (request.pathway_id,)).fetchall()
        genes += [r["gene_id"] for r in rows]
    genes = list(dict.fromkeys(genes))
    if not genes:
        raise HTTPException(status_code=400, detail="provide gene_ids or a pathway_id")
    if len(genes) > _MAX_SCREEN_GENES:
        raise HTTPException(status_code=400,
                            detail=f"too many genes to screen ({len(genes)}; max {_MAX_SCREEN_GENES})")
    if not species:
        species = _species_of(genes[0])

    transcripts = rnai.get_transcripts(species, DATA_DIR)
    if transcripts is None:
        return {"species": species, "available": False, "results": [],
                "note": f"no transcript store for {species}"}

    ranked = rnai.screen(genes, transcripts, k=request.k, window=request.design_window)
    meta = _gene_meta([r["gene_id"] for r in ranked])
    emx = expression.get_matrix(species) if species in set(expression.species_with_expression()) else None
    for r in ranked:
        m = meta.get(r["gene_id"], {})
        r["symbol"] = m.get("symbol", r["gene_id"])
        r["label_inferred"] = m.get("label_inferred", False)
        prof = emx.profile(r["gene_id"]) if emx else None
        r["mean_tpm"] = prof["mean_tpm"] if prof else None

    effect = None
    if request.predict_effect and ranked:
        pr = await perturb(PerturbRequest(
            interventions=[PerturbInterv(gene_id=r["gene_id"], action="ko") for r in ranked[:15]],
            depth=3))
        effect = {"affected": pr["stats"]["affected"], "up": pr["stats"]["up"],
                  "down": pr["stats"]["down"], "unknown": pr["stats"]["unknown"],
                  "top": pr["effects"][:8]}

    return {"species": species, "available": True, "n_genes": len(ranked),
            "designable": sum(1 for r in ranked if r["designable"]),
            "results": ranked, "predicted_effect": effect,
            "note": "Predicted dsRNA designability (fewest off-target genes in the best "
                    "window). Verify a chosen gene with /dsrna design mode."}


# ============= Organism overview =============

@app.get("/api/v1/organism/{species}/overview")
async def organism_overview(
    species: str,
    top: int = Query(25, le=100),
    min_confidence: float = Query(0.0),
    include_inferred: bool = Query(True),
):
    """Whole-organism summary: gene/coverage counts, edge counts split by
    evidence (measured vs inferred), and the top regulators by out-degree —
    the entry point for a network too large to render whole."""
    cur = db.conn.execute
    gene_count = cur("SELECT COUNT(*) FROM genes WHERE species = ?", (species,)).fetchone()[0]
    if gene_count == 0:
        raise HTTPException(status_code=404, detail="Unknown species")
    located = cur("SELECT COUNT(*) FROM gene_locations WHERE species = ?", (species,)).fetchone()[0]
    tf_count = cur("SELECT COUNT(*) FROM genes WHERE species = ? AND is_tf = 1", (species,)).fetchone()[0]

    edge_rows = cur(
        """
        SELECT CASE WHEN i.sources LIKE '%Inferred%' THEN 'inferred' ELSE 'measured' END AS kind,
               COUNT(*) n
        FROM interactions i JOIN genes g ON g.id = i.source_id
        WHERE g.species = ? GROUP BY kind
        """, (species,)
    ).fetchall()
    edges = {r["kind"]: r["n"] for r in edge_rows}
    edges = {"measured": edges.get("measured", 0), "inferred": edges.get("inferred", 0)}
    edges["total"] = edges["measured"] + edges["inferred"]

    inferred_clause = "" if include_inferred else " AND i.sources NOT LIKE '%Inferred%'"
    top_rows = cur(
        f"""
        SELECT i.source_id AS id, g.symbol, g.name, g.is_tf, g.synonyms, COUNT(*) AS out_degree
        FROM interactions i JOIN genes g ON g.id = i.source_id
        WHERE g.species = ? AND i.confidence >= ?{inferred_clause}
        GROUP BY i.source_id ORDER BY out_degree DESC LIMIT ?
        """, (species, min_confidence, top)
    ).fetchall()
    top_regulators = [
        {"id": r["id"], "symbol": r["symbol"], "name": r["name"],
         "is_tf": bool(r["is_tf"]), "out_degree": r["out_degree"],
         "synonyms": r["synonyms"].split("; ") if r["synonyms"] else []}
        for r in top_rows
    ]

    return {
        "species": species,
        "genes": gene_count,
        "transcription_factors": tf_count,
        "genes_with_coordinates": located,
        "regulators": cur(
            "SELECT COUNT(DISTINCT source_id) FROM interactions i JOIN genes g ON g.id = i.source_id "
            "WHERE g.species = ?", (species,)).fetchone()[0],
        "edges": edges,
        "top_regulators": top_regulators,
    }


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

@app.get("/api/v1/species")
async def species_capabilities():
    """Per-species capability matrix: which data layers are populated (network,
    orthologs, binding sites, expression, pathways, traits). This is the
    onboarding-readiness view — a new species (e.g. dahlia) appears here and fills
    in as its data lands."""
    cur = db.conn.execute
    species = [r[0] for r in cur("SELECT DISTINCT species FROM genes ORDER BY species")]
    expr_species = set(expression.species_with_expression())
    rows = []
    for sp in species:
        assembly = _ASSEMBLY_OF.get(sp)
        genes = cur("SELECT COUNT(*) FROM genes WHERE species=?", (sp,)).fetchone()[0]
        measured = cur("SELECT COUNT(*) FROM interactions i JOIN genes t ON t.id=i.target_id "
                       "WHERE t.species=? AND i.sources NOT LIKE '%Inferred%'", (sp,)).fetchone()[0]
        inferred = cur("SELECT COUNT(*) FROM interactions i JOIN genes t ON t.id=i.target_id "
                       "WHERE t.species=? AND i.sources LIKE '%Inferred%'", (sp,)).fetchone()[0]
        orthologs = cur("SELECT COUNT(*) FROM orthologs o JOIN genes g ON g.id=o.gene_a "
                        "WHERE g.species=?", (sp,)).fetchone()[0]
        binding = cur("SELECT COUNT(*) FROM motif_hits WHERE assembly=?", (assembly,)).fetchone()[0] if assembly else 0
        pathways = cur("SELECT COUNT(*) FROM pathway_annotations a JOIN genes g ON g.id=a.gene_id "
                       "WHERE g.species=?", (sp,)).fetchone()[0]
        traits = cur("SELECT COUNT(*) FROM trait_associations a JOIN genes g ON g.id=a.gene_id "
                     "WHERE g.species=?", (sp,)).fetchone()[0]
        emx = expression.get_matrix(sp) if sp in expr_species else None
        rows.append({
            "species": sp, "assembly": assembly, "genes": genes,
            "layers": {
                "network": {"measured_edges": measured, "inferred_edges": inferred},
                "orthologs": orthologs,
                "binding_sites": binding,
                "expression_samples": emx.n if emx else 0,
                "pathway_annotations": pathways,
                "trait_associations": traits,
            },
        })
    return {"species": rows,
            "note": "Layer counts reflect currently-loaded data; empty layers are "
                    "onboarding opportunities."}


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
    return JSONResponse(status_code=exc.status_code,
                        content={"error": exc.detail, "status_code": exc.status_code})

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

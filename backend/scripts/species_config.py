"""Central per-species onboarding config (single source of truth for the generic
ingestion scripts and the onboarding runbook).

Adding a new species = add an entry here + drop in its reference files, then run the
generic fetchers (fetch_seqctx.py, scan_motifs.py, fetch_expression.py) and loaders.
The `dahlia` entry is a prepared placeholder for the incoming Alex/Zach data — fill in
the URLs/assembly once the NCBI BioProject / genome release is available.
"""

# Each entry: assembly tag, PLAZA species code (if in PLAZA), reference URLs, and the
# id_style used to map reference transcript ids -> atlas gene ids:
#   "strip_isoform" : Peaxi162…​.1  -> gene id (rsplit '.')            (petunia/tomato)
#   "agi"           : AT1G01010.1   -> AT1G01010                       (arabidopsis)
# expression_panel is a scripts-local reference to the curated SRA panel (kept in the
# per-species fetch script for now; listed here for discoverability).
SPECIES = {
    "petunia": {
        "status": "loaded",
        "assembly": "Peaxi162v1.6.2",
        "plaza_code": "pax",
        "gff_url": "https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04_5/GFF/pax/"
                   "annotation.selected_transcript.all_features.pax.gff3.gz",
        "cds_url": "https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04_5/Fasta/"
                   "cds.selected_transcript.pax.fasta.gz",
        "genome_url": None,
        "id_style": "strip_isoform",
        "expression_panel": "fetch_petunia_expression.py",
    },
    "tomato": {
        "status": "loaded",
        "assembly": "SL4.0",
        "plaza_code": "sly",
        "gff_url": "https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04_5/GFF/sly/"
                   "annotation.selected_transcript.all_features.sly.gff3.gz",
        "cds_url": "https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04_5/Fasta/"
                   "cds.selected_transcript.sly.fasta.gz",
        "genome_url": None,
        "id_style": "strip_isoform",
        "expression_panel": "fetch_tomato_expression.py",
    },
    "arabidopsis": {
        "status": "loaded",
        "assembly": "TAIR10",
        "plaza_code": "ath",
        "gff_url": "https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04_5/GFF/ath/"
                   "annotation.selected_transcript.all_features.ath.gff3.gz",
        "cds_url": "https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04_5/Fasta/"
                   "cds.selected_transcript.ath.fasta.gz",
        "genome_url": "https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04_5/Genomes/ath.con.gz",
        "id_style": "agi",
        "expression_panel": "fetch_arabidopsis_expression.py",
    },
    # ---- Prepared placeholder for the incoming Dahlia collaboration data ----
    # Fill assembly + URLs once the G3 / NCBI BioProject genome+annotation is released
    # (Alex/Zach, 2026). Dahlia (Dahlia pinnata/variabilis, taxid 42159) is a strong fit
    # for the anthocyanin/floral-pigmentation focus.
    "dahlia": {
        "status": "pending",
        "assembly": None,
        "plaza_code": None,            # not in PLAZA; use the BioProject genome+GFF+CDS
        "gff_url": None,
        "cds_url": None,
        "genome_url": None,
        "id_style": "strip_isoform",   # confirm against the released annotation
        "expression_panel": None,      # build from the paper's RNA-seq on SRA
        "notes": "Incoming: genome+GFF+RNA-seq (G3/BioProject) + ~400 GWAS cultivar runs "
                 "(Zach) for trait mapping. See memory grn-atlas-dahlia-collaboration.",
    },
}


def get(species):
    return SPECIES.get(species)

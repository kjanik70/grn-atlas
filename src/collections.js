// Curated starter gene sets, surfaced as one-click entry points into the
// gene-set Analyze panel. Petunia is the classic model for anthocyanin
// (flower pigmentation) regulation, so these lean floral/pigmentation.
//
// Caveat: gene IDs are drawn from the atlas's functional descriptions. The
// flavonoid *structural* genes are reliably identified; the lineage-specific
// master regulators (AN2/AN1/AN11) are NOT cleanly annotated here and are
// omitted rather than guessed — see the planned BLAST-based identification.

export const COLLECTIONS = [
  {
    id: 'petunia-anthocyanin',
    name: 'Petunia — anthocyanin pathway',
    species: 'petunia',
    description: 'Flavonoid/anthocyanin biosynthesis structural genes (CHS→UF3GT).',
    geneIds: [
      'Peaxi162Scf00047g01225', // Chalcone synthase A
      'Peaxi162Scf00164g00313', // Chalcone synthase B
      'Peaxi162Scf00080g01317', // Chalcone isomerase
      'Peaxi162Scf00328g01214', // F3H (naringenin 3-dioxygenase)
      'Peaxi162Scf00150g00218', // F3'5'H
      'Peaxi162Scf00032g00067', // DFR-like 1
      'Peaxi162Scf00329g00024', // DFR-like 1 (paralog)
      'Peaxi162Scf00427g00022', // Anthocyanidin 3-O-glucosyltransferase
    ],
  },
  {
    id: 'petunia-floral-mads',
    name: 'Petunia — floral MADS-box',
    species: 'petunia',
    description: 'MADS-box / AGAMOUS-like floral-development regulators.',
    geneIds: [
      'Peaxi162Scf00016g00324', // AGAMOUS-like AGL62
      'Peaxi162Scf00274g00457', // AGAMOUS-like 30
      'Peaxi162Scf00008g01613', // MADS-box TF 50
      'Peaxi162Scf00013g00725', // MADS-box TF 18
    ],
  },
];

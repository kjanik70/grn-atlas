// Pick the friendliest short label for a gene. Tomato/petunia genes have no
// native symbol (their `symbol` is just the locus id), so fall back to the
// inferred Arabidopsis-ortholog synonym, flagged as inferred so it is never
// mistaken for a curated symbol.
export function geneLabel(gene) {
  if (!gene) return { label: '', inferred: false, id: '' };
  const id = gene.id || gene.symbol || '';
  if (gene.symbol && gene.symbol !== id) {
    return { label: gene.symbol, inferred: false, id };
  }
  if (gene.synonyms && gene.synonyms.length > 0) {
    return { label: gene.synonyms[0], inferred: true, id };
  }
  return { label: gene.symbol || id, inferred: false, id };
}

import React, { useState, useEffect, useCallback } from 'react';
import { analysisAPI } from '../services/apiService';
import SubgraphGraph from './SubgraphGraph';
import '../styles/GeneSetPanel.css';

const NS_LABEL = { BP: 'process', CC: 'component', MF: 'function', '': '' };

function fmtP(p) {
  if (p === 0) return '0';
  if (p < 0.001) return p.toExponential(1);
  return p.toFixed(3);
}

// Compute out-degree hubs from subgraph edges.
function topHubs(nodes, edges, k = 5) {
  const out = {};
  edges.forEach((e) => { out[e.source] = (out[e.source] || 0) + 1; });
  const byId = Object.fromEntries(nodes.map((n) => [n.id, n]));
  return Object.entries(out)
    .sort((a, b) => b[1] - a[1])
    .slice(0, k)
    .map(([id, deg]) => ({ symbol: byId[id]?.symbol || id, deg }));
}

export default function GeneSetPanel({ open, onClose, initialGeneIds, species, includeInferred }) {
  const [text, setText] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [subgraph, setSubgraph] = useState(null);
  const [enrichment, setEnrichment] = useState(null);
  const [lastIds, setLastIds] = useState([]);

  // Download the sequence-context export (signed edges + coords + promoter windows).
  const exportEdges = async () => {
    const res = await fetch('/api/v1/export/edges', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        gene_ids: lastIds, include_inferred: includeInferred, format: 'tsv',
        promoter_upstream: 2000, promoter_downstream: 500,
      }),
    });
    const text = await res.text();
    const url = URL.createObjectURL(new Blob([text], { type: 'text/tab-separated-values' }));
    const a = document.createElement('a');
    a.href = url; a.download = `grn_edges_export.tsv`; a.click();
    URL.revokeObjectURL(url);
  };

  const analyze = useCallback(async (geneIds, sp) => {
    if (!geneIds || geneIds.length < 2) {
      setError('Provide at least 2 genes.');
      return;
    }
    setLoading(true);
    setError(null);
    try {
      const [sg, enr] = await Promise.all([
        analysisAPI.subgraph(geneIds, { includeInferred }),
        analysisAPI.enrich(geneIds, sp),
      ]);
      setSubgraph(sg);
      setEnrichment(enr);
      setLastIds(geneIds);
    } catch (e) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  }, [includeInferred]);

  // Auto-run when opened from the network ("analyze this network").
  useEffect(() => {
    if (open && initialGeneIds && initialGeneIds.length) {
      setText(initialGeneIds.join(', '));
      analyze(initialGeneIds, species);
    }
  }, [open, initialGeneIds, species, analyze]);

  if (!open) return null;

  // Resolve pasted tokens (symbols or ids, any species) to gene ids, then analyze.
  const runFromText = async () => {
    const tokens = text.split(/[\s,]+/).map((t) => t.trim()).filter(Boolean);
    if (tokens.length < 2) { setError('Provide at least 2 genes.'); return; }
    setLoading(true);
    setError(null);
    try {
      const resolved = await Promise.all(tokens.map(async (tok) => {
        const r = await fetch(`/api/v1/genes/search?q=${encodeURIComponent(tok)}&limit=1`);
        const d = await r.json();
        return d.results?.[0] || null;
      }));
      const hits = resolved.filter(Boolean);
      if (hits.length < 2) { setError('Could not resolve those genes.'); setLoading(false); return; }
      const counts = {};
      hits.forEach((h) => { counts[h.species] = (counts[h.species] || 0) + 1; });
      const sp = Object.entries(counts).sort((a, b) => b[1] - a[1])[0][0];
      await analyze(hits.map((h) => h.id), sp);
    } catch (e) {
      setError(e.message);
      setLoading(false);
    }
  };

  const hubs = subgraph ? topHubs(subgraph.nodes, subgraph.edges) : [];

  return (
    <div className="gs-overlay" onClick={onClose}>
      <div className="gs-modal" onClick={(e) => e.stopPropagation()}>
        <div className="gs-header">
          <h2>Gene-set analysis</h2>
          <button className="gs-close" onClick={onClose} aria-label="Close">×</button>
        </div>

        <p className="gs-hint">
          Paste gene IDs (human symbols, or locus IDs like AT1G01060 / Solyc… / Peaxi…),
          then Analyze. Species: <strong>{species || 'auto'}</strong>.
        </p>
        <textarea
          className="gs-input"
          rows={3}
          value={text}
          placeholder="TP53, MYC, CDKN1A, BAX, MDM2 …"
          onChange={(e) => setText(e.target.value)}
        />
        <button className="gs-run" onClick={runFromText} disabled={loading}>
          {loading ? 'Analyzing…' : 'Analyze'}
        </button>

        {error && <div className="gs-error">{error}</div>}

        {subgraph && (
          <div className="gs-section">
            <h3>Induced network</h3>
            <p className="gs-metrics">
              {subgraph.nodes.length} genes · {subgraph.edges.length} interactions among them
            </p>
            <button className="gs-export" onClick={exportEdges}
              title="Signed edges + confidence + genomic coordinates + promoter windows (TSS −2000/+500)">
              ⤓ Export edges + promoter windows (TSV)
            </button>
            {hubs.length > 0 && (
              <div className="gs-hubs">
                <span className="gs-label">Top regulators:</span>{' '}
                {hubs.map((h) => (
                  <span key={h.symbol} className="gs-hub">{h.symbol} <em>({h.deg})</em></span>
                ))}
              </div>
            )}
            {subgraph.edges.length > 0 && (
              <SubgraphGraph nodes={subgraph.nodes} edges={subgraph.edges} />
            )}
          </div>
        )}

        {enrichment && (
          <div className="gs-section">
            <h3>GO enrichment</h3>
            <p className="gs-metrics">
              {enrichment.study} of the genes annotated · background {enrichment.background} ·{' '}
              {enrichment.results.length} enriched terms (FDR-adjusted)
            </p>
            {enrichment.results.length === 0 ? (
              <p className="gs-metrics">No significantly enriched terms.</p>
            ) : (
              <table className="gs-table">
                <thead>
                  <tr><th>GO term</th><th></th><th>genes</th><th>q-value</th></tr>
                </thead>
                <tbody>
                  {enrichment.results.map((r) => (
                    <tr key={r.go_id}>
                      <td>
                        <a href={`https://amigo.geneontology.org/amigo/term/${r.go_id}`}
                           target="_blank" rel="noopener noreferrer">{r.name}</a>
                      </td>
                      <td className="gs-ns">{NS_LABEL[r.namespace] ?? r.namespace}</td>
                      <td className="gs-num">{r.study_count}/{r.background_count}</td>
                      <td className="gs-num">{fmtP(r.q_value)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            )}
          </div>
        )}
      </div>
    </div>
  );
}

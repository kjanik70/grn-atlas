import React, { useState, useEffect, useCallback } from 'react';
import { analysisAPI } from '../services/apiService';
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

  const run = useCallback(async (geneIds) => {
    if (!geneIds || geneIds.length < 2) {
      setError('Provide at least 2 genes.');
      return;
    }
    setLoading(true);
    setError(null);
    try {
      const [sg, enr] = await Promise.all([
        analysisAPI.subgraph(geneIds, { includeInferred }),
        analysisAPI.enrich(geneIds, species),
      ]);
      setSubgraph(sg);
      setEnrichment(enr);
    } catch (e) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  }, [species, includeInferred]);

  // Auto-run when opened from the network ("analyze this network").
  useEffect(() => {
    if (open && initialGeneIds && initialGeneIds.length) {
      setText(initialGeneIds.join(', '));
      run(initialGeneIds);
    }
  }, [open, initialGeneIds, run]);

  if (!open) return null;

  const runFromText = () => {
    const ids = text.split(/[\s,]+/).map((t) => t.trim()).filter(Boolean);
    run(ids);
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
            {hubs.length > 0 && (
              <div className="gs-hubs">
                <span className="gs-label">Top regulators:</span>{' '}
                {hubs.map((h) => (
                  <span key={h.symbol} className="gs-hub">{h.symbol} <em>({h.deg})</em></span>
                ))}
              </div>
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

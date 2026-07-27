import React, { useState } from 'react';
import { analysisAPI } from '../services/apiService';
import '../styles/GeneSetPanel.css';

// dsRNA / RNAi design + off-target analysis. Everything shown is PREDICTED silencing
// (exact siRNA k-mer matching), not a guarantee of knockdown.
export default function DsRnaPanel({ open, onClose }) {
  const [seq, setSeq] = useState('');
  const [target, setTarget] = useState('');
  const [species, setSpecies] = useState('petunia');
  const [setText, setSetText] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [res, setRes] = useState(null);
  const [screen, setScreen] = useState(null);

  if (!open) return null;

  const runScreen = async () => {
    const ids = setText.split(/[\s,]+/).map((t) => t.trim()).filter(Boolean);
    if (ids.length < 1) { setError('Provide gene ids to screen.'); return; }
    setLoading(true); setError(null); setScreen(null);
    try {
      const r = await analysisAPI.dsrnaScreen(ids, species);
      if (r.available === false) setError(r.note); else setScreen(r);
    } catch (e) { setError(e.message); } finally { setLoading(false); }
  };

  const run = async () => {
    if (!seq.trim() && !target.trim()) { setError('Provide a dsRNA sequence or a target gene id.'); return; }
    setLoading(true); setError(null); setRes(null);
    try {
      const r = await analysisAPI.dsrna({
        sequence: seq.trim() || null, targetGeneId: target.trim() || null, species,
      });
      if (r.available === false) setError(r.note || 'No transcript store for this species.');
      else setRes(r);
    } catch (e) { setError(e.message); } finally { setLoading(false); }
  };

  const pe = res?.predicted_effect;
  return (
    <div className="gs-overlay" onClick={onClose}>
      <div className="gs-modal" onClick={(e) => e.stopPropagation()}>
        <div className="gs-header">
          <h2>dsRNA / RNAi design</h2>
          <button className="gs-close" onClick={onClose} aria-label="Close">×</button>
        </div>
        <p className="gs-hint">
          Predict which genes a dsRNA would silence (on- and off-target), or leave the
          sequence blank and give a <strong>target gene</strong> to design the most specific
          window. Predicted from siRNA k-mer matches — not a guarantee of knockdown.
        </p>

        <div className="gs-cons-controls">
          <label className="gs-label">Species</label>
          <select value={species} onChange={(e) => setSpecies(e.target.value)}>
            <option value="petunia">petunia</option>
            <option value="tomato">tomato</option>
            <option value="arabidopsis">arabidopsis</option>
            <option value="dahlia">dahlia</option>
          </select>
          <input className="gs-input" style={{ width: 'auto', flex: 1 }} placeholder="target gene id (e.g. Peaxi162Scf00118g00310 = AN2)"
            value={target} onChange={(e) => setTarget(e.target.value)} />
        </div>
        <textarea className="gs-input" rows={3} value={seq}
          placeholder="dsRNA sequence (optional — blank + target gene = design mode)"
          onChange={(e) => setSeq(e.target.value)} />
        <button className="gs-run" onClick={run} disabled={loading}>
          {loading ? 'Analyzing…' : (seq.trim() ? 'Analyze dsRNA' : 'Design dsRNA')}
        </button>
        {error && <div className="gs-error">{error}</div>}

        <div className="gs-section">
          <h3>Screen a gene set / pathway</h3>
          <p className="gs-hint">Rank genes by how cleanly a specific dsRNA can be designed (fewest off-targets).</p>
          <textarea className="gs-input" rows={2} value={setText}
            placeholder="gene ids to screen (e.g. anthocyanin pathway genes), space/comma separated"
            onChange={(e) => setSetText(e.target.value)} />
          <button className="gs-run" onClick={runScreen} disabled={loading}>
            {loading ? 'Screening…' : 'Screen set'}
          </button>
          {screen && (
            <>
              <p className="gs-metrics">
                {screen.designable}/{screen.n_genes} genes have a fully-specific window
                {screen.predicted_effect && ` · silencing all → ↓${screen.predicted_effect.down} down`}
              </p>
              <table className="gs-table">
                <thead><tr><th>Gene</th><th>best-window off</th><th>transcript off</th><th>mean TPM</th></tr></thead>
                <tbody>
                  {screen.results.map((r) => (
                    <tr key={r.gene_id}>
                      <td>{r.symbol}{r.designable && <span className="gs-cons-yes"> ✓</span>}</td>
                      <td className="gs-num">{r.best_window_off_targets}</td>
                      <td className="gs-num">{r.transcript_off_targets}</td>
                      <td className="gs-num">{r.mean_tpm != null ? r.mean_tpm : '—'}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </>
          )}
        </div>

        {res && (
          <>
            <div className="gs-section">
              <h3>{res.mode === 'design' ? 'Designed dsRNA' : 'dsRNA'}</h3>
              <p className="gs-metrics">
                {res.dsrna_length} bp · {res.n_sirnas} siRNAs ·{' '}
                specificity {(res.specificity * 100).toFixed(0)}% ·{' '}
                {res.off_target_gene_count} off-target gene{res.off_target_gene_count === 1 ? '' : 's'}
              </p>
              {res.on_target && (
                <p className="gs-metrics">
                  On-target <strong>{res.on_target.symbol || res.on_target.gene_id}</strong>:{' '}
                  {res.on_target.sites} sites{res.on_target.mean_tpm != null ? ` · mean ${res.on_target.mean_tpm} TPM` : ''}
                </p>
              )}
              {res.design && (
                <textarea className="gs-input" rows={3} readOnly value={res.design.sequence}
                  title={`window ${res.design.start}-${res.design.end}`} />
              )}
            </div>

            {res.off_targets.length > 0 && (
              <div className="gs-section">
                <h3>Predicted off-targets</h3>
                <table className="gs-table">
                  <thead><tr><th>Gene</th><th>siRNA sites</th><th>mean TPM</th></tr></thead>
                  <tbody>
                    {res.off_targets.map((o) => (
                      <tr key={o.gene_id}>
                        <td>{o.symbol}{o.is_tf && <span className="gs-ns"> TF</span>}</td>
                        <td className="gs-num">{o.sites}</td>
                        <td className="gs-num">{o.mean_tpm != null ? o.mean_tpm : '—'}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            )}

            {pe && (
              <div className="gs-section">
                <h3>Predicted downstream effect</h3>
                <p className="gs-metrics">
                  Knocking down the silenced gene(s): {pe.affected} affected ·{' '}
                  ↑{pe.up} ↓{pe.down} ?{pe.unknown}
                </p>
                {pe.top.length > 0 && (
                  <table className="gs-table">
                    <thead><tr><th>Gene</th><th>dir</th><th>mag</th></tr></thead>
                    <tbody>
                      {pe.top.map((e, i) => (
                        <tr key={i}>
                          <td>{e.symbol}</td>
                          <td>{e.predicted_direction === 'up' ? '↑' : e.predicted_direction === 'down' ? '↓' : '?'}</td>
                          <td className="gs-num">{e.magnitude.toFixed(2)}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                )}
              </div>
            )}
          </>
        )}
      </div>
    </div>
  );
}

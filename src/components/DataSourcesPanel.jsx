import React, { useEffect, useState } from 'react';
import '../styles/DataSourcesPanel.css';

const SOURCES = [
  {
    name: 'TRRUST v2',
    provides: 'Human TF–target regulatory interactions (literature-curated, with PubMed IDs).',
    cite: 'Han et al. (2018), Nucleic Acids Research 46:D380.',
    url: 'https://www.grnpedia.org/trrust/',
  },
  {
    name: 'PlantRegMap / PlantTFDB',
    provides: 'Arabidopsis and tomato TF–target regulation (FunTFBS functional binding sites).',
    cite: 'Tian et al. (2020), Nucleic Acids Research 48:D1104.',
    url: 'https://plantregmap.gao-lab.org/',
  },
  {
    name: 'OMA (Orthologous MAtrix)',
    provides: 'Genomic coordinates and cross-species orthologs for human, mouse, and Arabidopsis.',
    cite: 'Altenhoff et al. (2021), Nucleic Acids Research 49:D373.',
    url: 'https://omabrowser.org/',
  },
  {
    name: 'PLAZA Dicots 4.5',
    provides: 'Plant gene coordinates, synteny anchor points, BHIF orthology, gene descriptions (Arabidopsis, tomato, petunia).',
    cite: 'Van Bel et al. (2018), Nucleic Acids Research 46:D1190.',
    url: 'https://bioinformatics.psb.ugent.be/plaza/',
  },
  {
    name: 'DNA Zoo — Petunia axillaris Hi-C',
    provides: 'Chromosome-scale (7-chromosome) scaffolding of the P. axillaris v1.6.2 assembly.',
    cite: 'Dudchenko et al. (2017), Science 356:92; DNA Zoo Consortium.',
    url: 'https://www.dnazoo.org/',
  },
  {
    name: 'mygene.info',
    provides: 'Human/mouse gene names and identifiers.',
    cite: 'Xin et al. (2016), Genome Biology 17:91.',
    url: 'https://mygene.info/',
  },
];

const GLOSSARY = [
  ['Transcription factor (TF)', 'A protein that binds DNA to switch other genes on or off.'],
  ['Regulator / target', 'A regulator (TF) controls a target gene; edges point regulator → target.'],
  ['Activation / repression', 'Whether a regulator increases (activation) or decreases (repression) its target.'],
  ['Ortholog', 'The “same” gene in a different species, descended from a common ancestor.'],
  ['Synteny', 'Conserved gene order along chromosomes between species, shown as ribbons.'],
  ['Confidence', 'How well-supported an interaction is (evidence count / method); filterable.'],
  ['Inferred edge', 'A regulatory edge predicted for tomato/petunia by projecting the Arabidopsis network through orthology — a prediction, not a measurement.'],
];

export default function DataSourcesPanel({ open, onClose }) {
  const [stats, setStats] = useState(null);

  useEffect(() => {
    if (!open) return;
    fetch('/api/v1/stats').then((r) => r.json()).then(setStats).catch(() => {});
  }, [open]);

  if (!open) return null;

  return (
    <div className="ds-overlay" onClick={onClose}>
      <div className="ds-modal" onClick={(e) => e.stopPropagation()}>
        <div className="ds-header">
          <h2>Data sources &amp; citations</h2>
          <button className="ds-close" onClick={onClose} aria-label="Close">×</button>
        </div>

        {stats && (
          <p className="ds-stats">
            {stats.species} species · {stats.genes?.toLocaleString()} genes ·{' '}
            {stats.interactions?.toLocaleString()} interactions
          </p>
        )}

        <p className="ds-note">
          GRN Atlas combines measured regulation with <strong>inferred</strong> edges —
          the Arabidopsis network projected onto tomato and petunia through orthology.
          Inferred edges are shown dashed and labeled, and can be hidden with the
          “Include inferred edges” filter. They are predictions, not measurements.
        </p>

        <ul className="ds-list">
          {SOURCES.map((s) => (
            <li key={s.name} className="ds-item">
              <a className="ds-name" href={s.url} target="_blank" rel="noopener noreferrer">{s.name}</a>
              <div className="ds-provides">{s.provides}</div>
              <div className="ds-cite">{s.cite}</div>
            </li>
          ))}
        </ul>

        <h3 className="ds-subhead">Glossary</h3>
        <dl className="ds-glossary">
          {GLOSSARY.map(([term, def]) => (
            <div key={term} className="ds-term">
              <dt>{term}</dt>
              <dd>{def}</dd>
            </div>
          ))}
        </dl>
      </div>
    </div>
  );
}

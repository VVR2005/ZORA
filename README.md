<div align="center">

# ZORA

**Sequence & Mutation Analysis Workstation**

[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://python.org)
[![PySide6](https://img.shields.io/badge/GUI-PySide6-41cd52)](https://pypi.org/project/PySide6/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Status](https://img.shields.io/badge/status-beta-yellow)]()

*A desktop application for genomic sequence analysis, mutation detection, CRISPR guide design, molecular docking, and interactive visualization — all in one place.*

</div>

---

## Overview

ZORA brings together the core tools of computational genomics into a single, coherent desktop interface. Instead of switching between web servers, command-line tools, and spreadsheet editors, you can load sequences, detect and classify mutations, align them, design CRISPR guides, run molecular docking, and visualize results — all without leaving the application.

**Key capabilities:**

- **Sequence Analysis** — GC content, codon usage, ORF finding, melting temperature, molecular weight, motif search
- **Mutation Detection** — Base-level SNP/indel classification (silent, missense, nonsense, stop-loss), batch variant simulation
- **Pairwise Alignment** — Needleman-Wunsch (global) and Smith-Waterman (local) with threaded execution
- **CRISPR Design** — SpCas9 (NGG) guide RNA discovery with GC-content ranking
- **Molecular Docking** — AutoDock Vina integration with receptor/ligand preparation via Open Babel and RDKit Meeko
- **Visualization** — 7 plot types including GC content, conservation, mutation distributions, codon usage, and SNP density heatmaps
- **Project System** — Save/load/recover full analysis state with crash recovery
- **External Tool Integration** — One-click send to PyMOL, ChimeraX, and SnapGene
- **PDB Browser** — Import, fetch from RCSB, and interact with protein structures
- **Plugin System** — Extend functionality with Python plugins
- **Theming** — 6 dark/light themes with 8 contrast accent variants

---

## Quick Start

```bash
# Install dependencies
pip install pyside6 matplotlib numpy

# Run
python zora_main.py
```

*Optional dependencies for docking:* `rdkit`, `meeko`, AutoDock Vina, Open Babel.

---

## Screenshots

> *Screenshots coming soon.*

---

## Documentation

| Resource | Location |
|----------|----------|
| Test data guide | [`test_data/README_TEST_DATA.md`](test_data/README_TEST_DATA.md) |
| Example plugin | [`plugins/example_plugin.py`](plugins/example_plugin.py) |

---

## Features

### Sequence Analysis

Load FASTA, FASTQ, GenBank, CSV, or plain text. Compute GC content and skew, molecular weight, melting temperature, nucleotide frequencies, 3-frame translation, open reading frames, and per-amino-acid codon usage.

### Mutation Detection

Compare a reference and mutated sequence to classify every difference as silent, missense, nonsense, stop-loss, insertion, or deletion. Simulate individual mutations or generate batches of random variants with tunable SNP/insertion/deletion rates. Export results as CSV, JSON, or text reports.

### Pairwise Alignment

Global (Needleman-Wunsch) and local (Smith-Waterman) alignment with configurable scoring. Both run in background threads to keep the interface responsive.

### Visualization

Seven plot types rendered with matplotlib: sliding-window GC content, per-position mutation bar charts, grouped nucleotide frequencies, codon usage (horizontal bar + donut pie), mutation type pie, 1D SNP density heatmap, and multi-sequence conservation with standard deviation bands. Export as SVG or PNG.

### CRISPR Guide Design

Scan any sequence for SpCas9 NGG PAM sites, extract 20 bp guide sequences, and rank by GC content (optimal 40–60%). Displays the top 50 guides with position, strand, and score.

### Molecular Docking

Prepare receptor and ligand PDBQT files via Open Babel or MGLTools, generate ligands from SMILES via RDKit Meeko, configure the binding box manually or auto-detect from mutation positions, run AutoDock Vina with configurable exhaustiveness, and visualize results as an affinity bar chart. Open receptor + best pose directly in PyMOL.

### External Tool Integration

Send sequences to PyMOL, ChimeraX, or SnapGene with one click. Non-standard amino acids (U, O, B, Z, J) are automatically mapped to standard residues before transfer.

---

## Test Data

The [`test_data/`](test_data/) directory contains curated datasets for testing every major feature:

| File | Contents |
|------|----------|
| `reference_brca1.fasta` | BRCA1 exon 11: reference + mutated (550 bp) |
| `synthetic_snp_test.fasta` | WT + single SNP + multiple SNPs (300 bp) |
| `crispr_target.fasta` | EGFR exon 19 + TP53 exon 5 |
| `alignment_test.fasta` | Human vs Neanderthal mtDNA (400 bp) |
| `sample_1crn.pdb` | Crambin protein structure (46 aa) |

---

## Project Structure

```
zora_main.py            # Main application (run this)
zora/                   # Python package
├── __init__.py         # Public API
├── models.py           # Data structures, sequence utilities, alignment
├── project.py          # Project management, file I/O, NCBI fetch
├── docking.py          # AutoDock Vina integration
├── utils.py            # Qt widgets, theming, external tools, plugins
├── workers.py          # Background thread workers
└── core.py             # Backwards-compatible re-exports
plugins/                # User-extensible plugins
test_data/              # Test datasets
TEST/                   # Sample project directory
```

---

## Dependencies

| Package | Required | Purpose |
|---------|----------|---------|
| PySide6 | Yes | Qt GUI framework |
| matplotlib | Yes | Plotting and visualization |
| numpy | Yes | Numerical operations |
| rdkit | Docking | SMILES parsing and ligand preparation |
| meeko | Docking | PDBQT writing from RDKit molecules |
| Open Babel | Docking | PDBQT conversion |
| AutoDock Vina | Docking | Molecular docking |

*Set custom paths for external tools via environment variables:* `ZORA_VINA_PATH`, `ZORA_OBABEL_PATH`, `ZORA_MGLTOOLS_ROOT`.

---

## License

[MIT](LICENSE) © 2026 Veera Rahul

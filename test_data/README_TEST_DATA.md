# ZORA Test Datasets

## Files

| File | Contents | What to test |
|------|----------|-------------|
| `reference_brca1.fasta` | BRCA1 exon 11: reference + mutated (550bp) | Mutation detection, codon impact, classification |
| `synthetic_snp_test.fasta` | WT + single SNP + multiple SNPs (300bp) | SNP detection, distribution graphs |
| `crispr_target.fasta` | EGFR exon19 + TP53 exon5 | CRISPR sgRNA design |
| `alignment_test.fasta`      | Human vs Neanderthal mtDNA (400bp) | Alignment (Needleman-Wunsch / Smith-Waterman), conservation plot |
| `sample_1crn.pdb`           | Crambin protein PDB (46aa)         | PDB import, RCSB fetch |

## Workflows to Try

### 1. Basic Analysis
1. Open `reference_brca1.fasta` → click **Analyze** → see GC%, codon usage, ORFs
2. Click **RevComp** → reverse complement added to project

### 2. Mutation Detection
1. Open `reference_brca1.fasta` → both sequence entries load
2. Select Reference (seq 0) and Mutated (seq 1) → click **Detect Mutations**
3. Click mutation table rows → see AA change in annotation panel
4. Click **View Mutations in Sequence** → see bracketed mutations

### 3. Graphs
1. Load any two sequences → run mutation detection
2. Try all graph buttons: GC Plot, Mutation Distribution, Pie, Heatmap, Conservation

### 4. CRISPR
1. Load `crispr_target.fasta` → go to CRISPR tab → click **Design sgRNAs**

### 5. Alignment
1. Load `alignment_test.fasta` → select Human vs Neanderthal → click alignment button

### 6. Motif Search
1. **Tools → Find Motif / Pattern** or press Ctrl+F
2. Try: `ATG` (start codons), `TATA` (TATA box), or regex `GC[AT]GC`

### 7. Mutation Simulation
1. Load any sequence → go to Mutations tab
2. Set position and base → click **Simulate Mutation**
3. New mutated sequence appears in project tree

### 8. Batch Variants
1. Load any sequence → set count (e.g. 50) → **Generate Variants**

### 9. Integrations
1. Go to **Integrations** tab
2. Fetch from NCBI: enter `NM_007294` (BRCA1) → click **Fetch from NCBI**
3. Import PDB: click **Import PDB File** → select `sample_1crn.pdb`
4. Launch PyMOL/ChimeraX: click buttons (if installed)
5. Discover Plugins: click **Discover Plugins** to load from `plugins/` directory
6. Export FASTA for use in SnapGene, ChimeraX, etc.

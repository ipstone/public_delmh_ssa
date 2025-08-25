# DELMH-SSA: Deletion Microhomology Single-Strand Annealing Analysis

Analysis pipeline for detecting microhomology and single-strand annealing patterns in genomic deletions.

## System Requirements

**Software:**
- Python 3.7+
- pandas ≥1.3.0  
- biopython ≥1.78
- snakemake ≥6.0
- conda/mamba

**Tested on:**
- Ubuntu 18.04/20.04 LTS
- CentOS 7/8
- macOS 10.15+

**Hardware:**
- 8GB+ RAM recommended

## Installation

```bash
# Clone repository
git clone https://github.com/ipstone/public_delmh_ssa.git
cd public_delmh_ssa

# Create conda environment  
conda create -n delmh python=3.8 pandas biopython snakemake
conda activate delmh

# Download reference genome
wget ftp://ftp.broadinstitute.org/bundle/2.3/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz
```

## Demo

Run example analysis on breast560 dataset (275 mutations, 6 samples):

```bash
# Update reference path in smk/run_breast560_delmh_ssa.smk if needed
# Run analysis
snakemake -s smk/run_breast560_delmh_ssa.smk -c1

# View results
head output/breast560_delmh_ssa_found_finetune.tsv
```

## Usage

**Input format:** TSV with columns: `CHROM`, `POS`, `REF`, `ALT`, `SAMPLE_ID`, `NORMAL_ID`

**Run on your data:**
```bash
# 1. Place deletion data in input/ directory
# 2. Copy and modify smk/run_breast560_delmh_ssa.smk for your dataset
# 3. Update REF_FASTA path
# 4. Run pipeline
snakemake -s smk/your_analysis.smk -c1
```

**Key parameters:**
- `HOMOLOGY_CUTOFF`: Similarity threshold (default: 0.8)
- `REF_FASTA`: Reference genome path

**Output:** TSV with microhomology metrics including sequence patterns, lengths, and coordinates.

## Files

```
├── lib/         # Python analysis scripts
├── smk/         # Snakemake workflows  
├── input/       # Example deletion data
└── output/      # Analysis results
```
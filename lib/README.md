# Library Code for Deletion Microhomology and SSA Calculation

Core Python analysis library for calculating deletion microhomology and single-strand annealing (SSA) patterns. These scripts process genomic deletion data to identify microhomology sequences.

## Available Scripts

- **`calc_delmh_ssa_breast560_finetune.py`**: Main SSA/homeology calculation algorithm with configurable similarity thresholds (default: 80%)

## Algorithm Overview

The SSA calculation algorithm:
1. Analyzes deletion sequences for homology patterns
2. Reports longest homology matches along deleted sequences  
3. Uses configurable similarity cutoff (default: 80% homology)
4. Fine-tuned algorithm requires terminal nucleotides to be ATGC (not gaps)

## Input Format

TSV file with required columns:
- `CHROM`: Chromosome
- `POS`: Position 
- `REF`: Reference sequence
- `ALT`: Alternate sequence
- `SAMPLE_ID`: Sample identifier
- `NORMAL_ID`: Normal sample identifier

## Usage

**Via Snakemake (recommended):**
```bash
snakemake -s smk/run_breast560_delmh_ssa.smk -c1
```

**Direct usage:**
```bash
python lib/calc_delmh_ssa_breast560_finetune.py \
  --ref-fasta human_g1k_v37.fasta \
  --homology-cutoff 0.8 \
  -i input/breast560_deletions_for_delmh.tsv \
  -o output/breast560_delmh_ssa_found_finetune.tsv
```

## Output

TSV file with microhomology metrics including:
- Sequence patterns and lengths
- Genomic coordinates  
- Homology ratios and positions
- Up/downstream microhomology annotations

## Parameters

- `--homology-cutoff`: Similarity threshold (default: 0.8)
- `--ref-fasta`: Reference genome FASTA file path
- `-i`: Input deletion file
- `-o`: Output results file
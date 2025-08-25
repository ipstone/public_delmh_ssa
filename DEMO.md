# Demo Instructions

## Quick Start Demo

Run analysis on included breast560 example data.

### Prerequisites
- Completed installation (see INSTALL.md)
- Activated conda environment: `conda activate delmh`

### Demo Steps

```bash
# 1. Navigate to repository
cd public_delmh_ssa

# 2. Check input data
head input/breast560_deletions_for_delmh.tsv

# 3. Run analysis pipeline
snakemake -s smk/run_breast560_delmh_ssa.smk -c1
```

### Expected Output

**Output file:** `output/breast560_delmh_ssa_found_finetune.tsv`

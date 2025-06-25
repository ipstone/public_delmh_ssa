# Library Code for Deletion Microhomology and SSA Calculation

This folder contains the core 'library' code to calculate deletions with microhomology and homeology at the lower level. The generated data is used for higher level analysis later.

## Key Files

- **`calc_delmh.py`**: Traditional microhomology calculation with perfect matches
- **`calc_delmh_ssa_serena_finetune.py`**: Enhanced SSA/homeology calculation with configurable similarity thresholds
- **`dataprep_*.R`**: Data preparation scripts for different datasets

## Format Convention Differences

There is a convention difference between Serena's table indels characterization calculations:

### Serena's Format
- Uses the last unchanged base as the position
- Lists the unchanged base in both REF and ALT fields

### Standard Format (BRCA data EU etc.)
- Uses the deleted base as the position
- Uses the deleted base as the REF (and the bases after)
- Uses '-' as the ALT to signal the deletion

**Note**: These convention differences are handled in the code modifications. Serena's R code has been modified to incorporate this one base difference.

## Usage

The library files are typically called through the Snakemake workflows in the `smk/` directory, but can also be used directly:

```bash
# Traditional microhomology calculation
python lib/calc_delmh.py --ref-fasta <reference.fasta> -i <input.tsv> -o <output.tsv>

# SSA/homeology calculation with fine-tuning
python lib/calc_delmh_ssa_serena_finetune.py --ref-fasta <reference.fasta> --homology-cutoff 0.8 -i <input.tsv> -o <output.tsv>
```
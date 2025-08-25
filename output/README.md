# Output Directory

This directory contains the results from deletion microhomology and SSA analysis.

## Expected Output Files

When you run the analysis pipeline, the following files will be generated:

### Primary Results
- **`*_delmh_found.tsv`**: Traditional microhomology calculation results
- **`*_delmh_ssa_found_finetune.tsv`**: SSA/homeology calculation results with fine-tuning

### Dataset-Specific Files
- **`breast560_delmh_*.tsv`**: Results from breast560 cohort (560 breast cancer samples)
- **`pcawg_delmh_*.tsv`**: Results from PCAWG dataset analysis
- **`TCGA-WGS_delmh_*.tsv`**: Results from TCGA whole genome sequencing data

## Output Format

The TSV files contain the following key columns:
- **`up_mh`**, **`down_mh`**: Longest homology sequences found upstream and downstream
- **`delmh_ssa_len`**: Length of the longest homeology sequence
- **`del_len`**: Length of the deletion
- Sample and genomic coordinate information

## Note

Output files are not included in the public repository but will be generated when you run the analysis pipeline on your datasets.
# R Analysis Scripts

This folder contains higher level analysis code in R to understand the different definitions of deletion microhomology and homeology calculation, and their effects in different tumor genomic/mutation backgrounds such as BRCA1 and BRCA2 biallelic mutations vs controls.

## Key Files

### Library Functions
- **`lib_calculate_delmh_metrics.R`**: Statistical metrics calculation functions
- **`lib_load_serena_delmh_ssa_sv_data.R`**: Load Serena dataset with SSA calculations
- **`lib_load_serena_clinical_sig_data.R`**: Load Serena clinical and signature data
- **`lib_load_pcawg_tcga_ssa_data.R`**: Load PCAWG and TCGA datasets

### Figure Generation Scripts
- **`figures_serena_*.R`**: Generate figures for Serena dataset analysis
- **`figures_BOPP*.R`**: Generate figures for BOPP dataset analysis

## Analysis Focus

The R scripts perform statistical comparisons and visualizations focusing on:

1. **BRCA1/BRCA2 vs Control**: Deletion microhomology patterns in different genetic backgrounds
2. **HRDetect Analysis**: Homologous recombination deficiency scoring
3. **Deletion Length Distributions**: Analysis across different size categories
4. **Homeology Frequency**: Statistical comparisons of SSA patterns

## Usage

Run analysis scripts from the R/ directory:

```bash
# Calculate metrics
Rscript lib_calculate_delmh_metrics.R

# Generate figures
Rscript figures_BOPP.R
```
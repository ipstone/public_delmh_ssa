# R Analysis Scripts

This folder contains higher level analysis code in R to understand the different definitions of deletion microhomology and homeology calculation, and their effects in different tumor genomic/mutation backgrounds such as BRCA1 and BRCA2 biallelic mutations vs controls.

## Key Files

### Library Functions
- **`lib_calculate_delmh_metrics.R`**: Statistical metrics calculation functions
- **`lib_load_breast560_delmh_ssa_sv_data.R`**: Load breast560 cohort dataset with SSA calculations
- **`lib_load_breast560_clinical_sig_data.R`**: Load breast560 cohort clinical and signature data
- **`lib_load_pcawg_tcga_ssa_data.R`**: Load PCAWG and TCGA datasets

### Figure Generation Scripts
- **`figures_breast560_*.R`**: Generate figures for breast560 cohort dataset analysis
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
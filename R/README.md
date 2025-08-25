# R Analysis Scripts

Example analysis scripts for DELMH-SSA pipeline results. Demonstrates visualization and statistical analysis of deletion microhomology patterns.

## Available Scripts

### Example Analysis Scripts
- **`figures_breast560-dataset.R`**: Histogram analysis of deletion and microhomology length distributions
- **`figures_BOPP.R`**: Boxplot comparison analysis for BRCA2 vs wild-type samples

## Analysis Examples

**figures_breast560-dataset.R** demonstrates:
- Deletion length distribution histograms
- Microhomology length analysis (≥5bp cutoff)
- Basic summary statistics
- Automated plot generation and data export

**figures_BOPP.R** demonstrates:
- Comparative analysis between sample groups
- Statistical significance testing (Wilcoxon test)
- Clinical annotation integration
- Publication-ready boxplots

## Quick Start

```bash
# Run from repository root directory
cd public_delmh_ssa

# Generate histogram analysis
Rscript R/figures_breast560-dataset.R

# Generate comparison plots  
Rscript R/figures_BOPP.R
```

## Output

Scripts generate:
- PDF plots in `plots/` directory
- TSV data files alongside plots
- Console summary statistics

## Requirements

R packages: `data.table`, `dplyr`, `ggpubr`, `ggprism`

```r
install.packages(c("data.table", "dplyr", "ggpubr", "ggprism"))
```
# Deletion Microhomology and Single-Strand Annealing (SSA) Analysis

This repository contains the computational pipeline for analyzing deletion microhomology and single-strand annealing (SSA) patterns in genomic data. We use Serena's 560 breast cancer dataset, PCAWG, and TCGA whole genome sequencing data to explore ways to calculate potential deletions from SSA.

## Algorithm

- Go through the length of deletions, check if homology, keeping note of the homology length 
- In the end, report the longest homology (length) along the deleted sequence, and the exact matching sequence
- The constant `HOMOLOGY_CUTOFF = 0.8` determines the length of homeology to be reported

## Key Features

- **Traditional Microhomology**: Calculate perfect match microhomology at deletion breakpoints
- **Homeology Detection**: Detect imperfect homology sequences with configurable similarity thresholds (default: 80%)
- **Fine-tuned Algorithm**: Enhanced algorithm requiring terminal nucleotides to be ATGC (not gaps)
- **Multiple Datasets**: Support for PCAWG, TCGA-WGS, and custom cancer genomics datasets
- **Flexible Analysis**: Comprehensive R analysis pipeline for statistical comparisons

## Repository Structure

```
├── lib/                    # Core calculation libraries
│   ├── calc_delmh.py      # Traditional microhomology calculation
│   ├── calc_delmh_ssa_serena_finetune.py  # SSA/homeology calculation
│   └── dataprep_*.R       # Data preparation scripts
├── smk/                   # Snakemake workflow files
│   ├── run_pcawg_delmh_ssa.smk
│   └── run_TCGA-WGS_delmh_ssa.smk
├── R/                     # R analysis and visualization scripts
│   ├── lib_*.R           # Analysis utility functions
│   └── figures_*.R       # Figure generation scripts
├── input/                 # Input data files
└── output/               # Analysis results (empty in public repo)
```

## Installation

```bash
# Create conda environment
conda create -n delmh biopython pandas
conda activate delmh
```

## Configuration

Before running the analysis, you need to update the following paths in the configuration files:

1. **Reference Genome**: Update `REF_FASTA` path in Snakemake files (`smk/*.smk`) to point to your reference genome FASTA file
2. **Data Paths**: Update data paths in the R scripts and data preparation files to point to your datasets:
   - PCAWG data: Update `pcawg_data_path` in `lib/dataprep_prepare_pcawg_deletion.R`
   - TCGA data: Update `tcga_data_path` in `lib/dataprep_prepare_TCGA-WGS_deletion.R`
   - Clinical data: Update paths in R analysis scripts as needed

## Usage

### Running the Analysis Pipeline

Use Snakemake workflows to run the analysis:

```bash
# Run PCAWG analysis
snakemake -s smk/run_pcawg_delmh_ssa.smk -c2

# Run TCGA analysis
snakemake -s smk/run_TCGA-WGS_delmh_ssa.smk -c2

# Clean results for new runs
snakemake -s smk/run_pcawg_delmh_ssa.smk -c1 -R clean
```


### R Analysis

```bash
# Run statistical analysis (from R/ directory)
Rscript lib_calculate_delmh_metrics.R
Rscript figures_serena_01_s05_plot_grid-search_finetune_del-window_mh-cutoff.R
```

## Input Format

Input files should be in TSV format. The pipeline handles two different deletion format conventions:

1. **Serena's format**: Uses the last unchanged base as the position, listing the unchanged base in both REF and ALT
2. **Standard format**: Uses the deleted base as the position, with deleted bases as REF and '-' as ALT

## Output

- **Primary Output**: `output/delmh_ssa_found.tsv`
- **Fields**: `up_mh`, `down_mh` fields report the longest length of deleted sequence with homology at least the homology ratio (default 0.8)

## Key Parameters

- `HOMOLOGY_CUTOFF = 0.8`: Default homology ratio threshold for SSA calculation
- `REF_FASTA`: Reference genome path (typically human_g1k_v37.fasta)

## Datasets Supported

- **Serena Dataset**: 560 breast cancer samples with indel/deletion data
- **PCAWG**: Pan-Cancer Analysis of Whole Genomes
- **TCGA-WGS**: The Cancer Genome Atlas Whole Genome Sequencing data

## Contact

For questions or issues, please open an issue on GitHub.
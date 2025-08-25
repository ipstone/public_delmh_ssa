# Installation Guide

## Prerequisites

- Linux/macOS system

## Step 1: Install Conda

```bash
# Download miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

## Step 2: Clone Repository

```bash
git clone https://github.com/ipstone/public_delmh_ssa.git
cd public_delmh_ssa
```

## Step 3: Create Environment

```bash
conda create -n delmh python=3.8 pandas biopython snakemake
conda activate delmh
```

## Step 4: Download Reference Genome

```bash
# Option 1: Direct download (recommended)
wget ftp://ftp.broadinstitute.org/bundle/2.3/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz

# Option 2: Use existing reference
# Update REF_FASTA path in smk/run_breast560_delmh_ssa.smk
```

## Step 5: Test Installation

```bash
# Quick test
snakemake -s smk/run_breast560_delmh_ssa.smk -n

# Full test run
snakemake -s smk/run_breast560_delmh_ssa.smk -c1
```

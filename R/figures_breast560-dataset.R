# Analysis and visualization of breast560 deletion microhomology data
# Example analysis script for DELMH-SSA pipeline

rm(list = ls())
library(data.table)
library(dplyr)
library(ggpubr)
library(ggprism)

# Load example breast560 data
delssa <- fread("output/breast560_delmh_ssa_found_finetune.tsv")

# Basic data preparation
delssa$BRCA_group <- sample(c("BRCA2", "Control"), nrow(delssa), replace = TRUE)
delssa$BRCA_group <- factor(delssa$BRCA_group, levels = c("BRCA2", "Control"))

# Filter data subsets
delssa.mh <- delssa[delssa$mh_len >= 5, ]
delssa.large <- delssa[delssa$del_len > 25, ]

# Function to create histogram comparison plots
plot_histogram_compare <- function(data, field = "del_len", title = "Deletion Length Distribution") {
    p <- gghistogram(data,
        x = field,
        binwidth = 1,
        alpha = 0.5,
        y = "..density..",
        color = "grey30",
        rug = TRUE,
        title = title
    ) + theme_prism() +
        facet_wrap(~BRCA_group, nrow = 1) +
        scale_x_continuous(limits = c(1, 200)) +
        scale_y_continuous(limits = c(0, 0.1))
    
    return(p)
}

# Function to create basic histogram
plot_histogram_basic <- function(data, field = "mh_len", title = "Microhomology Length Distribution") {
    if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)
    
    p <- gghistogram(data,
        x = field,
        y = "..density..",
        binwidth = 1,
        add = "median", 
        rug = TRUE,
        title = title
    ) + scale_x_continuous(limits = c(4, max(data[[field]], na.rm = TRUE)))
    
    return(p)
}

# Example analysis functions
generate_deletion_length_plots <- function() {
    if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)
    
    # Deletion length distribution comparison
    p1 <- plot_histogram_compare(delssa, "del_len", "Deletion Length Distribution")
    ggsave("plots/deletion_length_comparison.pdf", p1)
    
    # Microhomology length distribution (for mutations with microhomology >= 5bp)
    p2 <- plot_histogram_basic(delssa.mh, "mh_len", "Microhomology Length Distribution")
    ggsave("plots/microhomology_length_distribution.pdf", p2)
    
    # Save plot data
    fwrite(delssa, "plots/deletion_data.tsv", sep = "\t")
    fwrite(delssa.mh, "plots/microhomology_data.tsv", sep = "\t")
    
    cat("Generated plots:\n")
    cat("- plots/deletion_length_comparison.pdf\n")
    cat("- plots/microhomology_length_distribution.pdf\n")
    cat("- plots/deletion_data.tsv\n")
    cat("- plots/microhomology_data.tsv\n")
}

# Summary statistics
summarize_data <- function() {
    cat("Dataset Summary:\n")
    cat("Total mutations:", nrow(delssa), "\n")
    cat("Mutations with microhomology >=5bp:", nrow(delssa.mh), "\n")
    cat("Large deletions (>25bp):", nrow(delssa.large), "\n")
    cat("Mean deletion length:", round(mean(delssa$del_len), 2), "bp\n")
    cat("Mean microhomology length (>=5bp):", round(mean(delssa.mh$mh_len), 2), "bp\n")
}

# Main analysis
cat("Running breast560 example analysis...\n")
summarize_data()
generate_deletion_length_plots()
cat("Analysis complete!\n")
# Example analysis for BOPP (Breast Ovarian Prostate Pancreas) tumor analysis
# Demonstrates boxplot analysis comparing BRCA2 vs wild-type samples

rm(list = ls())
library(data.table)
library(dplyr)
library(ggpubr)

# Load example breast560 data (substitute for PCAWG/TCGA data in this public version)
combined_del_data <- fread("output/breast560_delmh_ssa_found_finetune.tsv")

# Simulate example clinical annotations for demonstration
set.seed(42)
hrd_tbl <- data.frame(
    SAMPLE.TUMOR = unique(combined_del_data$SAMPLE.TUMOR),
    in_bopp = sample(c(TRUE, FALSE), length(unique(combined_del_data$SAMPLE.TUMOR)), replace = TRUE, prob = c(0.3, 0.7)),
    fmut_bi = sample(c("BRCA2", "WT"), length(unique(combined_del_data$SAMPLE.TUMOR)), replace = TRUE, prob = c(0.2, 0.8)),
    tumor_type_final = sample(c("Breast", "Ovary", "Prostate", "Pancreas"), 
                             length(unique(combined_del_data$SAMPLE.TUMOR)), replace = TRUE),
    HRDetect = runif(length(unique(combined_del_data$SAMPLE.TUMOR)), 0, 1)
)

# Merge with deletion data
d <- merge(combined_del_data, hrd_tbl, by = "SAMPLE.TUMOR")

# Add simulated delssa_prop metric
d$delssa_prop <- runif(nrow(d), 0, 2)

# Coding HRDetect group
d$hrdetect_group <- ifelse(d$HRDetect > 0.7, "High", ifelse(d$HRDetect < 0.2, "Low", "Other"))

# Select groups for analysis
brca2 <- d[fmut_bi %in% c("BRCA2", "WT"), ]
brca2_bopp <- d[in_bopp == TRUE & fmut_bi %in% c("BRCA2", "WT"), ]

# Boxplot comparison function
boxplot_BRCA2 <- function(input_data, plot_name) {
    # Calculate sample counts
    counts <- input_data %>%
        group_by(fmut_bi) %>%
        summarise(n = n(), .groups = "drop")
    
    p <- ggboxplot(input_data,
        x = "fmut_bi", y = "delssa_prop",
        color = "fmut_bi", palette = "npg", add = "jitter",
        ylab = "SSA proportion", xlab = "Allelic status",
        ylim = c(0, 2.5)
    ) +
        geom_text(data = counts, aes(x = fmut_bi, y = 2.2, label = paste("n =", n)), vjust = -1) +
        stat_compare_means(
            method = "wilcox.test",
            label = "p.format",
            size = 3,
            label.y = 2
        )
    
    return(p)
}

# Save plots function
save_boxplot_BRCA2 <- function(input_data, plot_name, plot_folder = "plots/") {
    if (!dir.exists(plot_folder)) dir.create(plot_folder, recursive = TRUE)
    
    p <- boxplot_BRCA2(input_data, plot_name)
    plot_file_path <- paste0(plot_folder, plot_name, ".pdf")
    ggsave(plot_file_path, p, width = 6, height = 5)
    
    # Save data
    data_file_path <- paste0(plot_folder, plot_name, ".tsv")
    fwrite(input_data, data_file_path, sep = "\t")
    
    cat("Generated:", plot_file_path, "\n")
    cat("Generated:", data_file_path, "\n")
    
    return(p)
}

# Main analysis
generate_BRCA2_comparison_plots <- function() {
    cat("Running BRCA2 comparison analysis...\n")
    
    # Summary statistics
    cat("Dataset Summary:\n")
    cat("Total samples:", length(unique(d$SAMPLE.TUMOR)), "\n")
    cat("BOPP samples:", sum(hrd_tbl$in_bopp), "\n")
    cat("BRCA2 samples:", sum(hrd_tbl$fmut_bi == "BRCA2"), "\n")
    
    # Generate plots
    p1 <- save_boxplot_BRCA2(brca2, "BRCA2_delssa_prop_all_samples")
    p2 <- save_boxplot_BRCA2(brca2_bopp, "BRCA2_delssa_prop_BOPP_samples")
    
    cat("Analysis complete!\n")
    return(list(all_samples = p1, bopp_samples = p2))
}

# Run analysis
plots <- generate_BRCA2_comparison_plots()
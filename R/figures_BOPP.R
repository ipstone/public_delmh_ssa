# Make plots for the BOPP tumors in PCAWG and TCGA data

rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)
library(dplyr)
library(ggpubr)

source("R/lib_load_pcawg_tcga_ssa_data.R")
# get the hrd_table with icgc_donor_id which will match with our deletion calculated data
hrd_tbl <- get_hrd_table_for_combined_samples()
hrd_tbl <- hrd_tbl[, .(icgc_donor_id, in_bopp, fmut_bi, tumor_type_final, pair, HRDetect)]
table(hrd_tbl$in_bopp)
table(hrd_tbl$tumor_type_final)
table(hrd_tbl$tumor_type_final, hrd_tbl$in_bopp)


pcawg <- fread("output/pcawg_combined_methods_result_delmh.tsv")
tcga <- fread("output/TCGA-WGS_combined_methods_result_delmh.tsv")
table(pcawg$SAMPLE.TUMOR %in% tcga$SAMPLE.TUMOR)
combined_del_data <- rbind(pcawg, tcga)
d <- merge(combined_del_data, hrd_tbl, by.x = "SAMPLE.TUMOR", by.y = "icgc_donor_id")
# Coding HRDetect group
d$hrdetect_group <- ifelse(d$HRDetect > 0.7, "High", ifelse(d$HRDetect < 0.2, "Low", "Other"))

# Select the group of patients to plot
brca2 <- d[fmut_bi %in% c("BRCA2", "WT"), ]
brca2_bopp <- d[in_bopp == TRUE & fmut_bi %in% c("BRCA2", "WT"), ]
table(brca2$fmut_bi)
table(brca2_bopp$fmut_bi)
table(brca2$tumor_type_final)
table(brca2_bopp$tumor_type_final)
table(d$tumor_type_final)
table(d$tumor_type_final, d$fmut_bi)

boxplot_BRCA2 <- function(input_data, plot_name) {
    brca2_pcawg <- input_data

    # Calculate the number of cases for each group
    counts <- brca2_pcawg %>%
        group_by(fmut_bi) %>%
        summarise(n = n(), .groups = "drop")

    p <- ggboxplot(brca2_pcawg,
        x = "fmut_bi", y = "delssa_prop",
        color = "fmut_bi", palette = "npg", add = "jitter",
        ylab = "delssa proportion", xlab = "Allelic status",
        ylim = c(0, 5)
    ) +
        geom_text(data = counts, aes(x = fmut_bi, y = 4.5, label = paste("n =", n)), vjust = -1) +
        stat_compare_means(
            method = "wilcox.test",
            label = "p.format",
            size = 3,
            label.y = 4
        )
    return(p)
}

save_boxplot_BRCA2 <- function(input_data, plot_name, plot_folder = "plots/BOPP_indels/") {
    mkdirp(plot_folder)
    p <- boxplot_BRCA2(input_data, plot_name)
    plot_file_path <- paste0(plot_folder, plot_name, ".pdf")
    ggsave(plot_file_path, p, width = 5, height = 5)
    return(p)
}

save_BRCA2_delssa_prop_plots <- function() {
    # Make plots in all pcawg, and then in BOPP tumors
    p.brca2 <- save_boxplot_BRCA2(brca2, "BRCA2_delssa_prop_combined")
    p.brca2_bopp <- save_boxplot_BRCA2(brca2_bopp, "BRCA2_delssa_prop_BOPP_in_combined")
    return(list(p.brca2, p.brca2_bopp))
}


boxplot_hrdetect <- function(pdata, plot_filename, plot_folder) {
    # Use the pdata, to make boxplot with jittering, and stats comparision by hrdetect_group category

    pt <- ggplot(pdata, aes(
        x = hrdetect_group,
        y = delssa_prop, group = hrdetect_group
    )) +
        # geom_violin(position = "identity") +
        # Instead of violin plot, make box plot with jitter
        ggplot2::geom_boxplot(position = "identity") +
        ggplot2::geom_jitter(position = position_jitter(0.2), size = 0.6) +
        ggpubr::stat_compare_means() +
        # ggforce::geom_sina(aes(x = BRCA_group, y = delmh_prop, group = BRCA_group),
        #     position = "identity",
        #     size = 0.6
        # ) +
        # Instead of using geom_sina, we can use boxplot + jitter
        ggpubr::theme_pubr() +
        labs(x = "", y = "Homeology rate", title = "Comparison short deletions between HRDetect high and low group") +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),
            text = element_text(size = 15)
        )

    # plot_folder <- "plots/grid_search_window_mh-cutoff/"
    # plot_folder <- "plots/figure_boxplot/"

    if (!dir.exists(plot_folder)) {
        dir.create(plot_folder, recursive = TRUE)
    }
    # Save plot and data
    plot_filename_path <- paste0(plot_folder, "/", plot_filename)
    ggsave(plot_filename_path, pt, width = 10, height = 10)

    data_filename <- gsub(".pdf", ".tsv", plot_filename_path)
    fwrite(pdata, data_filename)
}

d_hrdetect_highlow_noBRCA12 <- d[hrdetect_group != "Other" & !(fmut_bi %in% c("BRCA1", "BRCA2")), ]
d_hrdetect_highlow_noBRCA12.BOPP <- d_hrdetect_highlow_noBRCA12[in_bopp == TRUE, ]

save_figure_BOPP_s04_boxplot <- function() {
    boxplot_hrdetect(d_hrdetect_highlow_noBRCA12.BOPP,
        "BOPP_boxplot_delssa_homeology_hrdetect_highlow_no_BRCA12.pdf",
        plot_folder = "plots/Figures/BOPP_s04_boxplot_delssa_homeology_hrdetect_highlow_no_BRCA12/"
    )

    ## Include BRCA1 and BRCA2 tumors
    # boxplot_hrdetect(mh_ssa_hrdetect.highlow, "boxplot_delssa_homeology_hrdetect_highlow.pdf")
}
# save_figure_BOPP_s04_boxplot()

# Make the plot in Serena breast cancer dataset for the figures mentioned by Manisha and Erika
#
# -- for the moment, just plot all samples (biallelic BRCA2 vs control), may fine tune by tumor subtype later
#
# -- As the data are using Serena's data here, we will run this script within ssa_serena subfolder

rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)
library(dplyr)
library(ggpubr)
library(ggprism)


# Erika specified plot:

#----------------------------------------------------------------------------
# Figure -1: homoloyg/homeology plots
# <5bp, 5-25, >25bp deletion buckets with one panel of the three boxplots (each with BRCA2 vs control) with homology >=3bp and a second panel of the three boxplots (each with BRCA2 vs control) with homeology >=5bp with 20% mismatch.
## -- this can be done with:
#
# .       R/grid-search_plots/plot_grid-search_finetune_del-window_mh-cutoff.R
#----------------------------------------------------------------------------


# Figure -2: deletion length distribution and microhomology length distribution ----
#
# A. deletion length histogram for BRCA2 vs control (deletions with homeology) plotted as a proportion of total deletions (the same way YJ has done with the large deletions)
# B. homeology length histogram for BRCA2 plotted as YJ has done for the large deletions
# ---

source("R/lib_load_serena_delmh_ssa_sv_data.R")
# Load allelic info
allelic_info <- get_davies_biallelic_samples()

# delmh <- load_serena_delmh()
delssa <- load_serena_delssa_finetune()
delssa.brca2 <- delssa[delssa$SAMPLE.TUMOR %in%
    c(allelic_info$brca2_biallelic, allelic_info$control), ]
delssa.brca2.only <- delssa[delssa$SAMPLE.TUMOR %in% allelic_info$brca2_biallelic, ]

delssa.brca2$BRCA_group <- ifelse(delssa.brca2$SAMPLE.TUMOR %in% allelic_info$brca2_biallelic, "BRCA2", "Control")
delssa.brca2$BRCA_group <- factor(delssa.brca2$BRCA_group, levels = c("BRCA2", "Control"))

# Prepare data object for histogram plot
delssa.brca2.mh <- delssa.brca2[delssa.brca2$mh_len >= 5, ]
delssa.brca2.large <- delssa.brca2[delssa.brca2$del_len > 25, ]

plot_histogram_field_compare <- function(
    input_data,
    plot_name,
    plot_field = "del_len",
    plot_folder = "plots/Figures/Serena02_histogram_of_deletion_lengths_for_brca2_and_control/",
    plot_title = "Deletion length distribution") {
    #
    mkdirp(plot_folder)
    # Make histogram plots for the plot fiele, using BRCA_group to facet
    delssa.pdata <- input_data

    #----------------
    # Originally using gghistogram, but not showing what's needed such as starting at 1
    p <- gghistogram(delssa.pdata,
        x = plot_field,
        # binwidth = 0.2,
        binwidth = 1,
        alpha = 0.5,
        y = "..density..",
        # color = "BRCA_group",
        # fill = "BRCA_group",
        color = "grey30",
        # color = "black"
        # fill = "lightgrey",
        # add = "median",
        rug = TRUE,
        # set the mininum of x-axis to 1 # somehow this doesnt work
        # xlim = c(1, max(delssa.pdata[[plot_field]])),
        # xlim = c(1, 200),
        title = plot_title
        # ) + facet_wrap(~BRCA_group, nrow = 1, scales = "free_y") +
    ) + theme_prism() +
        facet_wrap(~BRCA_group, nrow = 1) +
        scale_x_continuous(limits = c(1, 200)) +
        scale_y_continuous(limits = c(0, 0.1))
    # scale_x_continuous(limits = c(1, max(delssa.pdata[[plot_field]])))
    # # xlim(1, max(delssa.pdata[[plot_field]])) +
    # # ylim(0, 0.04)

    #----------------
    # Use Yingjie's plot functions
    # p <- delssa.pdata %>% ggplot(aes(x = del_len, y = stat(density * width))) + # see https://stackoverflow.com/questions/22181132/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion-by-group
    #     geom_histogram(color = "black", binwidth = 0.2, alpha = 0.5) +
    #     theme_prism() +
    #     xlab("Deletion length (bp)") +
    #     ylab("Proportion of deletions") +
    #     ggtitle(plot_title) +
    #     facet_wrap(~BRCA_group, scales = "fixed")
    # xlim(1, max(delssa.pdata[[plot_field]]))

    #----------------
    # Debug/verify using seperate dataset for the plots
    pdata.brca2 <- delssa.pdata[BRCA_group == "BRCA2"]
    pdata.control <- delssa.pdata[BRCA_group == "Control"]
    # Make histogram for del_len, for its histogram density
    p.brca2 <- gghistogram(pdata.brca2,
        x = plot_field,
        binwidth = 1,
        y = "..density..",
        color = "grey30",
        alpha = 0.5,
        rug = TRUE,
        title = "BRCA2"
    ) + theme_prism() +
        scale_x_continuous(limits = c(1, 200)) +
        scale_y_continuous(limits = c(0, 0.1))

    p.control <- gghistogram(pdata.control,
        x = plot_field,
        binwidth = 1,
        y = "..density..",
        color = "grey30",
        alpha = 0.5,
        rug = TRUE,
        title = "Control"
    ) + theme_prism() +
        scale_x_continuous(limits = c(1, 200)) +
        scale_y_continuous(limits = c(0, 0.1))

    plot_file_path.brca2 <- paste0(plot_folder, plot_name, "_", plot_field, "_BRCA2.pdf")
    plot_file_path.control <- paste0(plot_folder, plot_name, "_", plot_field, "_Control.pdf")
    ggsave(plot_file_path.brca2, p.brca2)
    ggsave(plot_file_path.control, p.control)


    #----------------
    # Save the final plot and data
    plot_file_path <- paste0(
        plot_folder,
        plot_name, "_",
        plot_field, ".pdf"
    )
    plot_data_file_path <- paste0(
        plot_folder, plot_name, "_",
        plot_field, ".tsv"
    )
    ggsave(plot_file_path, p)
    fwrite(delssa.pdata, plot_data_file_path, sep = "\t")
}

plot_histogram_field <- function(
    input_data,
    plot_name,
    plot_field = "del_len",
    plot_folder = "plots/Figures/Serena02_histogram_of_deletion_lengths_for_brca2_and_control/",
    plot_title = "Deletion length distribution") {
    # Make histogram plots for the plot fiele, using BRCA_group to facet
    mkdirp(plot_folder)
    delssa.pdata <- input_data

    p <- gghistogram(delssa.pdata,
        x = plot_field,
        y = "..density..",
        binwidth = 1,
        # y = stat(density * width),
        add = "median", rug = TRUE,
        title = plot_title
    ) + scale_x_continuous(limits = c(4, max(delssa.pdata[[plot_field]], na.rm = TRUE)))
    # xlim(1, max(delssa.pdata[[plot_field]]))

    # + scale_x_continuous(limits = c(1, max(delssa.pdata[[plot_field]], na.rm = TRUE))) # Setting limits properly
    # + xlim(1, max(delssa.pdata[[plot_field]]))

    # p <- delssa.pdata %>% ggplot(aes(x = plot_field, y = stat(density * width))) + # see https://stackoverflow.com/questions/22181132/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion-by-group
    #     geom_histogram(color = "black", binwidth = 0.2, alpha = 0.5) +
    #     theme_prism() +
    #     xlab(plot_field) +
    #     ylab("Proportion of deletions") +
    #     ggtitle(plot_title)

    plot_file_path <- paste0(plot_folder, plot_name, "_", plot_field, ".pdf")
    plot_data_file_path <- paste0(plot_folder, plot_name, "_", plot_field, ".tsv")

    ggsave(plot_file_path, p)
    fwrite(delssa.pdata, plot_data_file_path, sep = "\t")
}

plot_histogram_field2 <- function(
    input_data,
    plot_name,
    plot_field = "del_len",
    plot_folder = "plots/Figures/Serena02_histogram_of_deletion_lengths_for_brca2_and_control/",
    plot_title = "Deletion length distribution") {
    # Make histogram plots for the plot field, using BRCA_group to facet
    mkdirp(plot_folder)
    delssa.pdata <- input_data
    p <- ggplot(delssa.pdata, aes_string(x = plot_field)) +
        geom_histogram(aes(y = ..density..), binwidth = 1) +
        scale_x_continuous(limits = c(1, max(delssa.pdata[[plot_field]], na.rm = TRUE))) +
        labs(title = plot_title)

    plot_file_path <- paste0(plot_folder, plot_name, "_", plot_field, ".pdf")
    plot_data_file_path <- paste0(plot_folder, plot_name, "_", plot_field, ".tsv")

    ggsave(plot_file_path, p)
    fwrite(delssa.pdata, plot_data_file_path, sep = "\t")
}

save_figure02_histogram <- function() {
    plot_histogram_field_compare(delssa.brca2,
        plot_name = "BRCA2_and_Control",
        plot_field = "del_len",
        plot_folder = "plots/Figures/Serena02_histogram_of_deletion_lengths_for_brca2_and_control/",
        plot_title = "Deletion length distribution for BRCA2 vs Control"
    )
}

save_figure03_mh_histogram <- function() {
    plot_histogram_field(delssa.brca2.only,
        plot_name = "BRCA2_",
        plot_field = "mh_len",
        plot_folder = "plots/Figures/Serena03_histogram_of_homeology_lengths_in_brca2/",
        plot_title = "Homeology length distribution for BRCA2"
    )
}

# Figure -4

# Figure -3: boxplot with HRDtect score categories ----
# box plot of homeologous deletions for HR detect high vs low (to combine with YJ's data - slide 15
#
# -- YJ's plot make HRDetect low and high group comparision without BRCA12 or /with BRCA12 samples
# .   the cutoff for HRDetect score is HRDetect-high: >0.7 HRDetect-low: < 0.2
#

# Loading functions to calculate the delssa prop rate
# source("R/lib_calculate_delmh_metrics.R")

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
            text = element_text(size = 25)
        )

    # plot_folder <- "plots/grid_search_window_mh-cutoff/"
    # plot_folder <- "plots/figure_boxplot/"

    if (!dir.exists(plot_folder)) {
        dir.create(plot_folder, recursive = TRUE)
    }
    # Save plot and data
    plot_filename_path <- paste0(plot_folder, "/", plot_filename)
    ggsave(plot_filename, pt, width = 10, height = 10)
    data_filename <- gsub(".pdf", ".tsv", plot_filename_path)
    fwrite(pdata, data_filename, sep = "\t")
}

mh_ssa <- calc_delmh_prop_cutoff(delssa, 25, 5)
# -- note here this calculation has a small difference with Erika's email about >=26bp

mh_ssa$delssa_prop <- mh_ssa$delmh_prop
hrdetect <- fread("input/Davies_HRDetect_2017/David_2017_HRDetect_score.csv")
hrdetect$Sample <- gsub("a", "", hrdetect$Sample)
mh_ssa_hrdetect <- merge(mh_ssa, hrdetect, by.x = "SAMPLE.TUMOR", by.y = "Sample")
mh_ssa_hrdetect$hrdetect_group <- ifelse(mh_ssa_hrdetect$predictorProb > 0.7, "High", ifelse(mh_ssa_hrdetect$predictorProb < 0.2, "Low", "Other"))

## Exclude BRCA1, BRCA2 tumors
brca12_tumor_samples <- c(allelic_info$brca1_biallelic, allelic_info$brca2_biallelic)

# for mh_ssa_hrdetect , further limit to hrdetect high and low case
table(mh_ssa_hrdetect$hrdetect_group)

mh_ssa_hrdetect.highlow <- mh_ssa_hrdetect[hrdetect_group != "Other", ]
mh_ssa_hrdetect.highlow.no_BRCA12 <- mh_ssa_hrdetect.highlow[!mh_ssa_hrdetect.highlow$SAMPLE.TUMOR %in% brca12_tumor_samples, ]

save_figure05_boxplot <- function() {
    boxplot_hrdetect(mh_ssa_hrdetect.highlow.no_BRCA12, "boxplot_delssa_homeology_hrdetect_highlow_no_BRCA12.pdf",
        plot_folder = "plots/Figures/Serena05_boxplot_delssa_homeology_hrdetect_highlow_no_BRCA12/"
    )

    ## Include BRCA1 and BRCA2 tumors
    # boxplot_hrdetect(mh_ssa_hrdetect.highlow, "boxplot_delssa_homeology_hrdetect_highlow.pdf")
}


# Load Serena repeat masker annotated data
repeat_masker <- fread("output/del_annotated_with_repeatmasker.csv")
table(repeat_masker$repClass)
table(repeat_masker$repClass == "")
table(repeat_masker$repName)
table(repeat_masker$repFamily)

convert_repeat_cat_table <- data.table(table(repeat_masker$repClass)) %>%
    arrange(-N) %>%
    # select the top 8 classes
    head(8)




repeat_masker$repClass_coded <- ifelse(repeat_masker$repClass %in% convert_repeat_cat_table$V1, repeat_masker$repClass, "Other")
# For repClass_coded that is "", set to "Non-Repeat-Region"
repeat_masker$repClass_coded[repeat_masker$repClass == ""] <- "Non-Repeat-Region"
table(repeat_masker$repClass_coded)

original_code <- function() {
    all_repeat_table <- data.table(table(repeat_masker$repClass_coded))
    all_repeat_table$pct_deletions <- all_repeat_table$N / sum(all_repeat_table$N)
    all_repeat_table <- all_repeat_table %>% arrange(desc(pct_deletions))
    all_repeat_table$V1 <- factor(all_repeat_table$V1, levels = all_repeat_table$V1, ordered = TRUE)


    p <- all_repeat_table %>% ggplot(aes(y = pct_deletions, x = V1)) +
        geom_bar(stat = "identity", position = position_dodge(), fill = "royalblue") +
        theme_prism(base_size = 10) +
        xlab("") +
        ylab("Proportion") +
        ggtitle("Repeat class around deletions with homeology") +
        scale_fill_viridis_d(name = "Repeat Class") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_x_discrete(drop = F)

    ggsave(p, filename = "plots/Figures/Serena04_repeat_region_plots/all_deletions_repClass_barplot.pdf")
    data_filename <- "plots/Figures/Serena04_repeat_region_plots/all_deletions_repClass_barplot.tsv"
    fwrite(all_repeat_table, data_filename, sep = "\t")
}

generate_repeat_region_barplot <- function(del_rep_data, plot_folder, plot_label) {
    del_rep_data$repClass_coded <- ifelse(del_rep_data$repClass %in% convert_repeat_cat_table$V1, del_rep_data$repClass, "Other")
    del_rep_data$repClass_coded[del_rep_data$repClass == ""] <- "Non-Repeat-Region"

    all_repeat_table <- data.table(table(del_rep_data$repClass_coded))
    all_repeat_table$pct_deletions <- all_repeat_table$N / sum(all_repeat_table$N)
    all_repeat_table <- all_repeat_table %>% arrange(desc(pct_deletions))
    all_repeat_table$V1 <- factor(all_repeat_table$V1, levels = all_repeat_table$V1, ordered = TRUE)

    p <- all_repeat_table %>% ggplot(aes(y = pct_deletions, x = V1)) +
        geom_bar(stat = "identity", position = position_dodge(), fill = "royalblue") +
        theme_prism(base_size = 10) +
        xlab("") +
        ylab("Proportion") +
        ggtitle("Repeat class around deletions with homeology") +
        # scale_fill_viridis_d(name = "Repeat Class") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_x_discrete(drop = F)

    ggsave(paste0(plot_folder, plot_label, "_repClass_barplot.pdf"), plot = p)
    fwrite(all_repeat_table, paste0(plot_folder, plot_label, "_repClass_barplot.tsv"), sep = "\t")
}

generate_repeat_region_piechart <- function(
    del_rep_data,
    plot_folder = "plots/Figures/Serena04_repeat_region_plots/",
    plot_label) {
    # Coding category for visualization
    del_rep_data$repClass_coded <- ifelse(
        del_rep_data$repClass %in% convert_repeat_cat_table$V1, del_rep_data$repClass,
        "Other"
    )

    del_rep_data$repClass_coded[del_rep_data$repClass == ""] <- "Non-Repeat-Region"
    del_rep_data$repeat_region <- ifelse(del_rep_data$repClass_coded == "Non-Repeat-Region", "Other region", "Repeat region")

    # Create piechart for del_rep_data$repeat_region
    # Count the number of occurrences for each category
    del_rep_counts <- table(del_rep_data$repeat_region)

    # Convert the table to a data frame for plotting
    del_rep_df <- as.data.frame(del_rep_counts)
    names(del_rep_df) <- c("repeat_region", "count")

    # Caclculate percentage for each V1 group
    del_rep_df$pct <- del_rep_df$count / sum(del_rep_df$count)
    # del_rep_df$plabel <- paste0(del_rep_df$repeat_region, " (", round(del_rep_df$pct * 100, 1), "% n=", del_rep_df$count, ")")
    del_rep_df$plabel <- paste0(del_rep_df$repeat_region, " (", round(del_rep_df$pct * 100, 1), "%)")

    # Create the pie chart
    p <- ggpie(del_rep_df, "count",
        label = "plabel",
        fill = "repeat_region",
        lab.font = 8,
        # color = "repeat_region",
        # palette = "npg"
    ) + # make the label font also of 14
        theme(
            # text = element_text(size = 14), # aAdjust this size as needed
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 20),
        )


    ggsave(paste0(plot_folder, plot_label, "_repRegion_piechat.pdf"),
        plot = p,
    )
    fwrite(del_rep_df, paste0(plot_folder, plot_label, "_repRegion_piechat.tsv"), sep = "\t")
}


ran_main <- function() {
    save_figure02_histogram()
    # # save_figure03_mh_histogram()
    # save_figure05_boxplot()

    generate_repeat_region_barplot(del_rep_data = repeat_masker, plot_folder = "plots/Figures/Serena04_repeat_region_plots/", plot_label = "all_deletions")
    #
    # deletion with homeology definitions
    repeat_mask.homeology <- repeat_masker[repeat_masker$mh_len.delssa >= 5, ]
    repeat_mask.nonhomeology <- repeat_masker[repeat_masker$mh_len.delssa < 5, ]
    #
    generate_repeat_region_barplot(del_rep_data = repeat_mask.homeology, plot_folder = "plots/Figures/Serena04_repeat_region_plots/", plot_label = "homeology_deletions")
    #
    generate_repeat_region_barplot(del_rep_data = repeat_mask.nonhomeology, plot_folder = "plots/Figures/Serena04_repeat_region_plots/", plot_label = "nonhomeology_deletions")
}


running_main <- function() {
    # For deletions > 25bps ----
    # table(repeat_masker$del_len.delmh == repeat_masker$indel.length)
    repeat_mask_25bp <- repeat_masker[repeat_masker$indel.length >= 25, ]
    repeat_mask_25bp.homeology <- repeat_mask_25bp[repeat_mask_25bp$mh_len.delssa >= 5, ]
    table(repeat_mask_25bp$repClass)

    del_rep_data <- repeat_mask_25bp
    table(del_rep_data$repClass)

    generate_repeat_region_piechart(repeat_mask_25bp, "plots/Figures/Serena04_repeat_region_plots/", "all-deletions_greater25bp")
    generate_repeat_region_piechart(repeat_mask_25bp.homeology, "plots/Figures/Serena04_repeat_region_plots/", "homeology_deletions_greater25bp")

    # For deletions 5-25bps ----
    repeat_mask_5_25bp <- repeat_masker[repeat_masker$indel.length >= 5 & repeat_masker$indel.length < 25, ]
    repeat_mask_5_25bp.homeology <- repeat_mask_5_25bp[repeat_mask_5_25bp$mh_len.delssa >= 5, ]
    table(repeat_mask_5_25bp$repClass)
    generate_repeat_region_piechart(repeat_mask_5_25bp, "plots/Figures/Serena04_repeat_region_plots/", "all-deletions_5-25bp")
    generate_repeat_region_piechart(repeat_mask_5_25bp.homeology, "plots/Figures/Serena04_repeat_region_plots/", "homeology_deletions_5-25bp")
}
# ran_main()
running_main()

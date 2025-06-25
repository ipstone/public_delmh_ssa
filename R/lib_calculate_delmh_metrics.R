#  the conventions for the function names del and the next word stands for what type of deletion calculations it is calculated from
# -- For example, delssa: mean the deletion usig ssa /homelogy method
# .               delserena: using the serena MH calculaton method




####################################################
# Encode the metrics for del
# Calculate the mean deletion length for each sample

calculate_del_metrics <- function(del) {
    # Calculate the deletion count for each sample
    del_count <- del %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(del_count = n())
    del_count$log2_del_count <- log2(del_count$del_count)

    # Calculate the mean deletion length
    del_mean_len <- del %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(mean_del_len = mean(del_len))

    # Calculate the mean delmh-ssa length for each sample
    delmh_ssa_mean_len <- del %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(mean_delmh_ssa_len = mean(delmh_ssa_len))

    # Calculate the percentage of deletions with delmh_ssa_len >=5
    delmh_ssa_len_5 <- del %>%
        filter(delmh_ssa_len >= 5) %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(delmh_ssa_len_5 = n())

    # Calculate the percentage of deletions with delmh_ssa_len >=3
    delmh_ssa_len_3 <- del %>%
        filter(delmh_ssa_len >= 3) %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(delmh_ssa_len_3 = n())

    # Calculate the percentage of deletions with delmh_ssa_len >=2
    delmh_ssa_len_2 <- del %>%
        filter(delmh_ssa_len >= 2) %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(delmh_ssa_len_2 = n())

    # Combined the above counts as del_metrics, fill the missing values with 0
    del_metrics <- del_count %>%
        inner_join(del_mean_len, by = "SAMPLE.TUMOR") %>%
        inner_join(delmh_ssa_mean_len, by = "SAMPLE.TUMOR") %>%
        left_join(delmh_ssa_len_5, by = "SAMPLE.TUMOR") %>%
        left_join(delmh_ssa_len_3, by = "SAMPLE.TUMOR") %>%
        left_join(delmh_ssa_len_2, by = "SAMPLE.TUMOR") %>%
        mutate(delmh_ssa_len_5 = ifelse(is.na(delmh_ssa_len_5), 0, delmh_ssa_len_5)) %>%
        mutate(delmh_ssa_len_5_pct = delmh_ssa_len_5 / del_count) %>%
        mutate(delmh_ssa_len_3 = ifelse(is.na(delmh_ssa_len_3), 0, delmh_ssa_len_3)) %>%
        mutate(delmh_ssa_len_3_pct = delmh_ssa_len_3 / del_count) %>%
        mutate(delmh_ssa_len_2 = ifelse(is.na(delmh_ssa_len_2), 0, delmh_ssa_len_2)) %>%
        mutate(delmh_ssa_len_2_pct = delmh_ssa_len_2 / del_count)


    # Code del_metrics ER status and biallelic status
    del_metrics$tumor_subtype <- ifelse(del_metrics$SAMPLE.TUMOR %in% serena$erPosSamples, "ER+", "Triple-")
    del_metrics$allelic_status <- ifelse(del_metrics$SAMPLE.TUMOR %in% allelic_info$biallelic, "biallelic", ifelse(del_metrics$SAMPLE.TUMOR %in% allelic_info$control, "control", "other"))

    # Add brca1_biallelic and brca2_biallelic status
    del_metrics$brca1_biallelic <- ifelse(del_metrics$SAMPLE.TUMOR %in% allelic_info$brca1_biallelic, "brca1_biallelic", "other")
    del_metrics$brca2_biallelic <- ifelse(del_metrics$SAMPLE.TUMOR %in% allelic_info$brca2_biallelic, "brca2_biallelic", "other")

    # Add davies del.mh.prop, del.mh.count data
    davies_sample$davies_total_del_count <- davies_sample$del.mh.count + davies_sample$del.rep.count + davies_sample$del.none.count

    del_metrics <- del_metrics %>%
        left_join(davies_sample[, .(
            sample,
            del.mh.prop,
            del.mh.count,
            davies_total_del_count
        )], by = c("SAMPLE.TUMOR" = "sample"))

    return(data.table(del_metrics))
}

# Create function to calculate_sv_metrics (based on calculate_del_metrics)
calculate_sv_metrics <- function(sv) {
    sv_metrics <- calculate_del_metrics(sv)

    # Calculate the mean delmh_sv_len for each sample
    delmh_sv_mean_len <- sv %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(mean_delmh_sv_len = mean(delmh_sv_len))

    # Calculate the percentage of delmh_sv_len with delmh_sv_len >=5
    delmh_sv_len_5 <- sv %>%
        filter(delmh_sv_len >= 5) %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(delmh_sv_len_5 = n())

    # Combine the above counts as sv_metrics, fill the missing values with 0
    sv_metrics <- sv_metrics %>%
        inner_join(delmh_sv_mean_len, by = "SAMPLE.TUMOR") %>%
        left_join(delmh_sv_len_5, by = "SAMPLE.TUMOR") %>%
        mutate(delmh_sv_len_5 = ifelse(is.na(delmh_sv_len_5), 0, delmh_sv_len_5)) %>%
        mutate(delmh_sv_len_5_pct = delmh_sv_len_5 / del_count)

    return(data.table(sv_metrics))
}

calc_delmh_prop_cutoff <- function(data, delen_cutoff, mhlen_cutoff) {
    # This code is good for both --delmh-- (most tradiational definition) AND
    # --delssa-- (allowing 80% match) definition
    #
    # Filter data
    filtered_data <- data %>%
        dplyr::filter(del_len >= delen_cutoff, mh_len >= mhlen_cutoff) %>%
        dplyr::group_by(SAMPLE.TUMOR)

    # Group by SAMPLE.TUMOR
    grouped_data <- data %>%
        dplyr::group_by(SAMPLE.TUMOR)

    # Calculate total number of deletions
    total_deletions <- grouped_data %>%
        dplyr::summarise(total = n())

    # Calculate number of deletions above cutoffs
    deletions_above_cutoffs <- filtered_data %>%
        dplyr::summarise(above_cutoff = n())

    # Join total_deletions and deletions_above_cutoffs, keep all rows of total deletions
    joined_data <- dplyr::left_join(total_deletions, deletions_above_cutoffs, by = "SAMPLE.TUMOR")
    # Fill in all the missed deletions_above_cutoffs as 0
    joined_data$above_cutoff[is.na(joined_data$above_cutoff)] <- 0

    # Calculate percentage
    joined_data$delmh_prop <- joined_data$above_cutoff / joined_data$total * 100
    joined_data$delen_cutoff <- delen_cutoff
    joined_data$mhlen_cutoff <- mhlen_cutoff
    joined_data <- data.table(joined_data)

    return(joined_data)
}

calc_delssa_strict_prop_cutoff <- function(
    data = merged,
    delen_cutoff,
    mhlen_cutoff) {
    #
    #  this is using the strict criterion to calculate the delssa rate

    # Filter data
    filtered_data <- data %>%
        dplyr::filter(
            del_len >= delen_cutoff,
            mh_len >= mhlen_cutoff
        ) %>%
        dplyr::group_by(SAMPLE.TUMOR)

    # Group by SAMPLE.TUMOR
    grouped_data <- data %>%
        # filter(del_len.delssa >= delen_cutoff) %>% # This is quite different with the bucket definition: the bucket was calcuating the rate within eacvh bucket, here is the rate of all the deletions
        dplyr::group_by(SAMPLE.TUMOR)

    # Calculate total number of deletions
    total_deletions <- grouped_data %>%
        dplyr::summarise(total = n())


    ## Adding the strict requirement to have '-' mismatch in the mh string
    deletions_above_cutoffs <- filtered_data %>%
        # Select rows where up_mh.delssa, or down_mh.delssa contains '-'
        # dplyr::filter(grepl("-", up_mh.delssa) | grepl("-", down_mh.delssa)) %>%
        dplyr::filter(grepl("-", finetune_mh)) %>%
        # -- this is an approximation ... most should come from down_mh.delssa
        dplyr::summarise(above_cutoff = n())


    # Join total_deletions and deletions_above_cutoffs, keep all rows of total deletions
    joined_data <- dplyr::left_join(total_deletions, deletions_above_cutoffs, by = "SAMPLE.TUMOR")
    # Fill in all the missed deletions_above_cutoffs as 0
    joined_data$above_cutoff[is.na(joined_data$above_cutoff)] <- 0

    # Calculate percentage
    joined_data$delmh_prop <- joined_data$above_cutoff / joined_data$total * 100
    joined_data$delen_cutoff <- delen_cutoff
    joined_data$mhlen_cutoff <- mhlen_cutoff
    joined_data <- data.table(joined_data)
    joined_data$sample <- joined_data$SAMPLE.TUMOR

    return(joined_data)
}

# Adding this function to make plots for dellen>=cutoff, and mhlen>=cutoff
calc_delserena_prop_cutoff <- function(
    data = merged,
    delen_cutoff,
    mhlen_cutoff) {
    #
    # This function calcs the Serena MH calculation method, using cutoff window approach
    # - the outoput of these function all use the same convention of delmh_prop, delen_cutoff, mhlen_cutoff in the output data.table
    #

    # Filter data
    filtered_data <- data %>%
        dplyr::filter(
            del_len >= delen_cutoff,
            mh_len >= mhlen_cutoff
        ) %>%
        dplyr::group_by(SAMPLE.TUMOR)

    # Group by SAMPLE.TUMOR
    grouped_data <- data %>%
        dplyr::group_by(SAMPLE.TUMOR)

    # Calculate total number of deletions
    total_deletions <- grouped_data %>%
        dplyr::summarise(total = n())

    deletions_above_cutoffs <- filtered_data %>%
        filter(classification == "Microhomology-mediated") %>% # only keep classification that's microhomology-mediated
        dplyr::summarise(above_cutoff = n())


    # Join total_deletions and deletions_above_cutoffs, keep all rows of total deletions
    joined_data <- dplyr::left_join(total_deletions, deletions_above_cutoffs, by = "SAMPLE.TUMOR")
    # Fill in all the missed deletions_above_cutoffs as 0
    joined_data$above_cutoff[is.na(joined_data$above_cutoff)] <- 0

    # Calculate percentage
    joined_data$delmh_prop <- joined_data$above_cutoff / joined_data$total * 100
    joined_data$delen_cutoff <- delen_cutoff
    joined_data$mhlen_cutoff <- mhlen_cutoff
    joined_data <- data.table(joined_data)
    joined_data$sample <- joined_data$SAMPLE.TUMOR

    return(joined_data)
}

calc_delserena_mh_prop <- function(
    data = merged,
    delen_cutoff,
    mhlen_cutoff) {
    #
    # This function calcs the Serena MH calculation method (default without any cutoff)
    # - the outoput of these function all use the same convention of delmh_prop, delen_cutoff, mhlen_cutoff in the output data.table
    #

    # Filter data
    filtered_data <- data %>%
        dplyr::group_by(SAMPLE.TUMOR)

    # Group by SAMPLE.TUMOR
    grouped_data <- data %>%
        dplyr::group_by(SAMPLE.TUMOR)

    # Calculate total number of deletions
    total_deletions <- grouped_data %>%
        dplyr::summarise(total = n())

    deletions_above_cutoffs <- filtered_data %>%
        filter(classification == "Microhomology-mediated") %>% # only keep classification that's microhomology-mediated
        dplyr::summarise(above_cutoff = n())


    # Join total_deletions and deletions_above_cutoffs, keep all rows of total deletions
    joined_data <- dplyr::left_join(total_deletions, deletions_above_cutoffs, by = "SAMPLE.TUMOR")
    # Fill in all the missed deletions_above_cutoffs as 0
    joined_data$above_cutoff[is.na(joined_data$above_cutoff)] <- 0

    # Calculate percentage
    joined_data$delmh_prop <- joined_data$above_cutoff / joined_data$total * 100
    joined_data <- data.table(joined_data)
    joined_data$sample <- joined_data$SAMPLE.TUMOR
    return(joined_data)
}

calc_del_basics <- function(del) {
    # Calculate the deletion count for each sample
    del_count <- del %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(del_count = n())
    # del_count$log2_del_count <- log2(del_count$del_count)

    # Calculate the mean deletion length
    del_mean_len <- del %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(mean_del_len = mean(del_len))

    # Combine the basics
    del_metrics <- del_count %>%
        inner_join(del_mean_len, by = "SAMPLE.TUMOR")

    return(data.table(del_metrics))
}

calculate_delmh_cutoff_metrics <- function(
    delmh_file,
    delssa_file,
    delserena_file,
    del_length_cutoff,
    mh_length_cutoff) {
    ##
    # Calculate and combine the metrics of short deletions from delmh, delssa, del_serena calculation using del_length, mh_length cutoffs, return a combined data.table containing all 3 types of MH_rate for each sample
    ##
    delmh_file <- "output/serena_delmh_found.tsv"
    delssa_file <- "output/serena_delmh_ssa_found_finetune.tsv"
    delserena_file <- "output/serena_deletion_classified.tsv"
    # TODO: -- this is the work in progress with the feedbacks from Manisha and Erika

    delmh <- load_delmh(delmh_file)
    delssa <- load_delssa_finetune(delssa_file) # Load func has some coding
    delserena <- load_delserena(delserena_file)


    # Calculate the deletion count for each sample
    del_count <- del %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(del_count = n())
    # del_count$log2_del_count <- log2(del_count$del_count)

    # Calculate the mean deletion length
    del_mean_len <- del %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(mean_del_len = mean(del_len))

    # calculate the mh length
    delmh$mh_len <- nchar(delmh$down_mh)
    ###########################################
    # UNFINISHED
    ###########################################


    # Calculate the percentage of deletions with delmh_ssa_len >=5
    delmh_ssa_len_5 <- del %>%
        filter(delmh_ssa_len >= 5) %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(delmh_ssa_len_5 = n())

    # Calculate the percentage of deletions with delmh_ssa_len >=3
    delmh_ssa_len_3 <- del %>%
        filter(delmh_ssa_len >= 3) %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(delmh_ssa_len_3 = n())

    # Calculate the percentage of deletions with delmh_ssa_len >=2
    delmh_ssa_len_2 <- del %>%
        filter(delmh_ssa_len >= 2) %>%
        group_by(SAMPLE.TUMOR) %>%
        summarize(delmh_ssa_len_2 = n())

    # Combined the above counts as del_metrics, fill the missing values with 0
    del_metrics <- del_count %>%
        inner_join(del_mean_len, by = "SAMPLE.TUMOR") %>%
        inner_join(delmh_ssa_mean_len, by = "SAMPLE.TUMOR") %>%
        left_join(delmh_ssa_len_5, by = "SAMPLE.TUMOR") %>%
        left_join(delmh_ssa_len_3, by = "SAMPLE.TUMOR") %>%
        left_join(delmh_ssa_len_2, by = "SAMPLE.TUMOR") %>%
        mutate(delmh_ssa_len_5 = ifelse(is.na(delmh_ssa_len_5), 0, delmh_ssa_len_5)) %>%
        mutate(delmh_ssa_len_5_pct = delmh_ssa_len_5 / del_count) %>%
        mutate(delmh_ssa_len_3 = ifelse(is.na(delmh_ssa_len_3), 0, delmh_ssa_len_3)) %>%
        mutate(delmh_ssa_len_3_pct = delmh_ssa_len_3 / del_count) %>%
        mutate(delmh_ssa_len_2 = ifelse(is.na(delmh_ssa_len_2), 0, delmh_ssa_len_2)) %>%
        mutate(delmh_ssa_len_2_pct = delmh_ssa_len_2 / del_count)


    # Code del_metrics ER status and biallelic status
    del_metrics$tumor_subtype <- ifelse(del_metrics$SAMPLE.TUMOR %in% serena$erPosSamples, "ER+", "Triple-")
    del_metrics$allelic_status <- ifelse(del_metrics$SAMPLE.TUMOR %in% allelic_info$biallelic, "biallelic", ifelse(del_metrics$SAMPLE.TUMOR %in% allelic_info$control, "control", "other"))

    # Add brca1_biallelic and brca2_biallelic status
    del_metrics$brca1_biallelic <- ifelse(del_metrics$SAMPLE.TUMOR %in% allelic_info$brca1_biallelic, "brca1_biallelic", "other")
    del_metrics$brca2_biallelic <- ifelse(del_metrics$SAMPLE.TUMOR %in% allelic_info$brca2_biallelic, "brca2_biallelic", "other")

    # Add davies del.mh.prop, del.mh.count data
    davies_sample$davies_total_del_count <- davies_sample$del.mh.count + davies_sample$del.rep.count + davies_sample$del.none.count

    del_metrics <- del_metrics %>%
        left_join(davies_sample[, .(
            sample,
            del.mh.prop,
            del.mh.count,
            davies_total_del_count
        )], by = c("SAMPLE.TUMOR" = "sample"))

    return(data.table(del_metrics))
}

# Loading existing biallelic information
source("R/lib_load_serena_clinical_sig_data.R")
source("R/lib_calculate_delmh_metrics.R")


load_davies_sample <- function() {
    sigs <- get_serena_signatures()
    serena <- get_serena_samples()
    allelic_info <- get_davies_biallelic_samples()

    # Use the clinical info from Davies' sample for ER status
    # Update this path to your Davies supplementary table location
    davies_sample <- fread("input/Davies_suppltable_1_sampleInfo.csv")
    names(davies_sample) <- ip_fixname(davies_sample)
    davies_sample$sample <- gsub("(a|b).*$", "", davies_sample$Sample)
    table(davies_sample$ERstatus)
    table(davies_sample$isUsedForEvaluation)
    table(davies_sample$ERstatus, davies_sample$isUsedForEvaluation)
    return(davies_sample)
}


load_del <- function() {
    # This is the initial function used for loading delssa daa
    # used in ssa_serena/lib/compare_delmh_ssa_biallelic_vs_control.R
    # Load the delmh-ssa results
    del <- fread("output/serena_delmh_ssa_found.tsv")
    # Lable del ER status and biallelic status
    del$tumor_subtype <- ifelse(del$SAMPLE.TUMOR %in% serena$erPosSamples, "ER+", "Triple-")
    del$allelic_status <- ifelse(del$SAMPLE.TUMOR %in% allelic_info$biallelic, "biallelic", ifelse(del$SAMPLE.TUMOR %in% allelic_info$control, "control", "other"))
    del$del_len <- nchar(del$REF)
    del$delmh_ssa_len <- nchar(del$down_mh)
    return(del)
}
# function to load delmh data
load_delmh <- function(input_file) {
    # Load delmh original calculation, using current delmh calc in the original 30gy paper
    # delmh <- fread("output/serena_delmh_found.tsv")
    delmh <- fread(input_file)

    delmh$del_len <- nchar(delmh$REF)
    delmh$method <- "delmh"
    delmh$sample_key <- paste(delmh$SAMPLE.TUMOR, delmh$CHROM, delmh$POS, delmh$REF, sep = "_")
    table(duplicated(delmh$sample_key))
    delmh$mh_found <- ifelse(delmh$down_found == "TRUE" | delmh$up_found == "TRUE", "TRUE", "FALSE")
    delmh$mh_len <- pmax(nchar(delmh$down_mh), nchar(delmh$up_mh))
    summary(delmh$mh_len)
    table(is.na(delmh$mh_len))
    table(delmh$mh_found)


    # Most of the mh sequence in down_mh is longer than up_mh
    table(nchar(delmh$up_mh) > nchar(delmh$down_mh))
    table(nchar(delmh$up_mh) < nchar(delmh$down_mh))
    return(delmh)
}

load_serena_delmh <- function() load_delmh("output/serena_delmh_found.tsv")

load_delssa_finetune <- function(input_file = "output/serena_delmh_ssa_found_finetune.tsv") {
    # Load/prepare delssa calculation where mistmatch is allowed up to 20% mismatch, so at least 80% match

    # Load delmh ssa calculation
    delssa <- fread(input_file)
    delssa$del_len <- nchar(delssa$REF)
    delssa$method <- "delssa_finetune"
    delssa$sample_key <- paste(delssa$SAMPLE.TUMOR, delssa$CHROM, delssa$POS, delssa$REF, sep = "_")
    table(duplicated(delssa$sample_key))

    delssa$mh_found <- ifelse(delssa$down_found == "TRUE" | delssa$up_found == "TRUE", "TRUE", "FALSE")
    # delssa$mh_found <- delssa$down_found
    delssa$mh_len <- pmax(nchar(delssa$down_mh), nchar(delssa$up_mh))
    # delssa$mh_len <- nchar(delssa$down_mh)

    #### FUTURE
    # -- Though likely the convention of deletion caller may prefer downstream calls but there are possibility that reverse is true -- unless this be fully ruled out,we are checking on both ends
    delssa$finetune_mh <- with(delssa, ifelse(nchar(down_mh) > nchar(up_mh), down_mh, up_mh))

    summary(delssa$mh_len)
    table(is.na(delssa$mh_len))
    table(delssa$mh_found)
    return(delssa)
}

load_serena_delssa_finetune <- function() load_delssa_finetune("output/serena_delmh_ssa_found_finetune.tsv")

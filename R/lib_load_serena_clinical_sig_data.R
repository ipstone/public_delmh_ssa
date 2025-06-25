# functions to load Serena/Davies' clinical data and siganture data

library(data.table)
library(magrittr)
library(ipfun)
library(dplyr)


get_serena_signatures <- function() {
    ## Signature and other info from Serena's paper
    serena_sigs <- fread("input/serena_nature_2016_breast_tumor_560/serena_nature_2016_560_breast_tumor_table21_signature_contribution.csv")
    serena_lst <- fread("input/serena_nature_2016_breast_tumor_560/serena_nature_2016_560_breast_tumor_table17_HRD_intermediary.csv")
    serena_sigs[, total := rowSums(.SD), .SDcols = 2:13]
    serena_sigs$sig3 <- serena_sigs[[4]] / serena_sigs$total
    # -- For these signature cutoff, for now just using sig3 >=0.3, LST>=20
    serena_sigs$sample <- gsub("(a|b).*$", "", serena_sigs[[1]])
    serena_lst$sample <- gsub("(a|b).*$", "", serena_lst$Sample)
    signatures <- merge(serena_sigs, serena_lst[, .(sample, LST)], by = "sample")
    return(signatures)
}

get_davies_biallelic_samples <- function() {
    sample <- fread("input/serena_nature_2016_breast_tumor_560/serena_nature_2016_560_breast_tumor_table1_clindata.csv")
    sample$sampleid <- gsub("(a|b).*$", "", sample$sample_name)
    # This is the cohort we are using for final testing 320 ER+/HER2-
    allSamples <- unique(sample$sampleid)

    # Return the biallelic/control samples list based on Davies HRDetect paper
    ## Loading mutation information to determine biallelic vs. control samples
    brca_mut <- fread("input/Davies_HRDetect_2017/Davies_suppltable_4a_mutation_details.csv")
    brca_mut$sample <- gsub("(a|b).*$", "", brca_mut$Sample)
    # -- use this to define the 77 cases with biallelic
    ddr_mut <- fread("input/Davies_HRDetect_2017/Davies_suppltable_4b_otherHRgene_mutation_details.csv")
    ddr_mut$sample <- gsub("(a|b).*$", "", ddr_mut$Sample)

    # -- use this to remote the other DDR gene (PALB2) mutations
    davies_sample <- fread("input/Davies_HRDetect_2017/Davies_suppltable_1_sampleInfo.csv")

    # Organize breast and OV samples list
    biallelic <- brca_mut[isBrcaMonoallelic == F]
    monoallelic <- ddr_mut[Gene == "PALB2"]
    signatures <- get_serena_signatures()
    dom_signatures <- signatures[signatures$LST >= 15 & signatures$sig3 >= 0.2]

    # Finalize sample categories
    control_samples <- allSamples[!(allSamples %in% biallelic$sample) &
        !(allSamples %in% monoallelic$sample) &
        !(allSamples %in% dom_signatures$sample)]

    biallelic_samples <- unique(biallelic$sample)

    # Further add BRCA1, BRCA2 biallelic samples
    brca1_biallelic <- brca_mut[isBrcaMonoallelic == F & Gene == "BRCA1"]
    brca2_biallelic <- brca_mut[isBrcaMonoallelic == F & Gene == "BRCA2"]
    brca1_biallelic_samples <- unique(brca1_biallelic$sample)
    brca2_biallelic_samples <- unique(brca2_biallelic$sample)

    fr <- list()
    fr$control <- control_samples
    fr$biallelic <- biallelic_samples
    fr$brca1_biallelic <- brca1_biallelic_samples
    fr$brca2_biallelic <- brca2_biallelic_samples
    return(fr)
}

get_serena_samples <- function() {
    # Loading Serena breast tumor type info
    sample <- fread("input/serena_nature_2016_breast_tumor_560/serena_nature_2016_560_breast_tumor_table1_clindata.csv")
    sample$sampleid <- gsub("(a|b).*$", "", sample$sample_name)

    # This is the cohort we are using for final testing 320 ER+/HER2-
    er_pos_sample <- sample[final.ER == "positive" & final.HER2 == "negative"]
    erPosSamples <- unique(er_pos_sample$sampleid)
    trp_neg_sample <- sample[final.ER == "negative" &
        final.HER2 == "negative" &
        final.PR == "negative"]
    trpNegSamples <- unique(trp_neg_sample$sampleid)

    rl <- list()
    rl$erPosSamples <- erPosSamples
    rl$trpNegSamples <- trpNegSamples
    return(rl)
}

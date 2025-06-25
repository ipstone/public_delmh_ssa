# Here are the function to load/proces PCAWG and TCGA data
library(data.table)

# Loading existing biallelic information
source("R/lib_load_serena_clinical_sig_data.R")
source("R/lib_calculate_delmh_metrics.R")


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

load_BOPP_delssa_finetune <- function() {
    # Combine the PCAWG and TCGA data
    pcawg <- load_delssa_finetune("output/pcawg_delmh_ssa_found_finetune.tsv")
    tcga <- load_delssa_finetune("output/TCGA-WGS_delmh_ssa_found_finetune.tsv")
    table(names(pcawg) == names(tcga))
    combined <- rbind(pcawg, tcga)
    return(combined)
    ## load_delssa_finetune("output/serena_delmh_ssa_found_finetune.tsv")
    ## -- previous code when working on Serena's
}

load_BOPP_delmh <- function() {
    # Combine the PCAWG and TCGA data
    pcawg <- load_delssa_finetune("output/pcawg_delmh_found.tsv")
    tcga <- load_delssa_finetune("output/TCGA-WGS_delmh_found.tsv")
    table(names(pcawg) == names(tcga))
    combined <- rbind(pcawg, tcga)
    return(combined)
    ## -- previous code when working on Serena's
}


# Get the tcga_hrd_table
get_hrd_table_for_TCGA_samples <- function() {
    nature_2020_keytable <- fread("input/nature2020_panCancerWGS_suppl_table1.txt")
    nature_2020_keytable <- nature_2020_keytable[, .(icgc_donor_id, submitted_specimen_id)]
    tcga_table <- nature_2020_keytable[grepl("TCGA", submitted_specimen_id), ]
    # str(tcga_table)
    # Classes ‘data.table’ and 'data.frame':  801 obs. of  2 variables:
    #  $ icgc_donor_id        : chr  "DO22145" "DO38966" "DO10172" "DO17323" ...
    #  $ submitted_specimen_id: chr  "TCGA-EZ-7264-01A" "TCGA-DJ-A2Q1-01A" "TCGA-A6-6141-01A" "TCGA-A3-3387-01A" ...

    # For tcga_table submitted_specimen_id, extract the 3rd part of ID, separated by '-'
    tcga_table$pair <- sapply(strsplit(tcga_table$submitted_specimen_id, "-"), function(x) x[3])


    # -- HRD table logic for processing PCAWG/TCGA data
    hrd_tbl <- fread("input/hrd-supp-table_ZC.tsv")
    hrd_tbl.simple <- hrd_tbl[, .(pair, dataset, in_bopp, fmut_bi, tumor_type_final)]

    tcga_table <- merge(tcga_table, hrd_tbl.simple, by = "pair")
    return(tcga_table)
}

get_hrd_table_for_combined_samples <- function() {
    # HRD status table from Marcin's paper
    hrd_tbl <- fread("input/hrd-supp-table_ZC.tsv")

    tcga_hrd_table <- get_hrd_table_for_TCGA_samples()
    tcga_hrd_table_select <- tcga_hrd_table[, .(icgc_donor_id, pair)]
    hrd_table <- merge(hrd_tbl, tcga_hrd_table_select, by = "pair", all.x = TRUE)

    fdn(hrd_table, "dono")
    table(hrd_table$dataset, is.na(hrd_table$icgc_donor_id))
    # in hrd_table, for dataset==PCAWG , assign icgc_donor_id field with value from  pair field
    hrd_table$icgc_donor_id[hrd_table$dataset == "PCAWG"] <- hrd_table$pair[hrd_table$dataset == "PCAWG"]
    table(hrd_table$dataset, is.na(hrd_table$icgc_donor_id))
    hrd_table_final <- hrd_table[!is.na(icgc_donor_id)]
    table(hrd_table_final$dataset, is.na(hrd_table_final$icgc_donor_id))
    return(hrd_table_final)
}

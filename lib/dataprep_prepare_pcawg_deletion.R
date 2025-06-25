# Use the pcawg consensus deletion calls
rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)
library(dplyr)

# Load PCAWG data - update this path to your PCAWG data location
# Download from: https://dcc.icgc.org/releases/PCAWG
pcawg_data_path <- "path/to/pcawg/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz"
d <- fread(pcawg_data_path)

# The convention of the ICGC data is similar to the Jorge lab's data, the REF is the deleted sequence, and ALT is designated as '-'

convert_pcawg_dels_for_delmh <- function(maf_data) {
    # subset for indels data
    d <- maf_data[Variant_Type == "DEL"] # only deletions
    d$Project <- d$Project_Code
    d$Sample <- d$Donor_ID
    d$ID <- d$Hugo_Symbol
    d$Genome <- "GRCh37"
    d$assembly_version <- "GRCh37"
    d$mut_type <- "DEL"
    d$SAMPLE.TUMOR <- d$Donor_ID
    d$SAMPLE.NORMAL <- d$Matched_Norm_Sample_Barcode
    d$CHROM <- d$Chromosome
    d$POS <- d$Start_position
    d$REF <- d$Reference_Allele
    d$ALT <- d$Tumor_Seq_Allele2
    d$Type <- "SOMATIC"
    d$variantCaller <- d$i_Callers
    d$total_read_count <- d$t_alt_count + d$t_ref_count
    d$mutant_allele_read_count <- d$t_alt_count

    # Processing duplicated samples from
    d$label <- paste(d$Sample,
        d$mut_type,
        d$Chromosome,
        d$Start_position,
        d$End_position,
        d$Reference_Allele,
        d$Tumor_Seq_Allele2,
        sep = "-"
    )

    # One way is to only include verified sequencing result
    d2 <- d[!duplicated(d$label)]

    # Create the specific columns to be returned for delmh calcualtion
    d.prepared <- d2[, .(
        SAMPLE.TUMOR,
        SAMPLE.NORMAL,
        CHROM,
        POS,
        REF,
        ALT,
        mut_type,
        variantCaller,
        assembly_version,
        Project,
        total_read_count,
        mutant_allele_read_count
    )]
    return(d.prepared)
}

d.prepared <- convert_pcawg_dels_for_delmh(d)
fwrite(d.prepared, "input/pcawg_deletions_for_delmh.tsv", sep = "\t")

# Use snakemake pipeline to run delmh/ssa calculation for longest partial matching homology
# Example workflow for breast560 cohort dataset - uses subset of data for demonstration

REF_FASTA = "/data/riazlab/reference/GATK_bundle/2.3/human_g1k_v37.fasta"
# Specify the level of homology allowed
HOMOLOGY_CUTOFF = 0.8


rule all:
    input:
        "output/breast560_delmh_ssa_found_finetune.tsv"

rule delmh_ssa:
    input:
        "input/breast560_deletions_for_delmh.tsv"
    output:
        "output/breast560_delmh_ssa_found_finetune.tsv"
    shell:
        "python lib/calc_delmh_ssa_breast560_finetune.py --ref-fasta {REF_FASTA} --homology-cutoff {HOMOLOGY_CUTOFF} -i {input} -o {output}"

rule clean:
    shell:
        "rm {output}"
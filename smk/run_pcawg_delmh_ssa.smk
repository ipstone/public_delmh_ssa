# Use snakemake pipeline to run delmh/ssa cacluation for longest partial matching homology

# Set path to your reference genome FASTA file
REF_FASTA = "path/to/your/reference_genome.fasta"
# Specify the level of homology allowed
HOMOLOGY_CUTOFF = 0.8
# 0.9 is a lot more stingent, for better specificity (resulting mostly perfectly matched homology sequence)


rule all:
    input:
        "output/pcawg_delmh_found.tsv",
        "output/pcawg_delmh_ssa_found_finetune.tsv"

rule delmh_mh:
    input:
        "input/pcawg_deletions_for_delmh.tsv"
    output:
        "output/pcawg_delmh_found.tsv"
    shell:
        "python lib/calc_delmh.py --ref-fasta {REF_FASTA} -i {input} -o {output}"

rule delmh_ssa:
    input:
        "input/pcawg_deletions_for_delmh.tsv"
    output:
        "output/pcawg_delmh_ssa_found_finetune.tsv"
    shell:
        "python lib/calc_delmh_ssa_breast560_finetune.py --ref-fasta {REF_FASTA} --homology-cutoff {HOMOLOGY_CUTOFF} -i {input} -o {output}"

rule clean:
    shell:
        "rm {output}"

#############################################
# Bison nasal microbiome pipeline
# (SNS/DNS directory structure)
#############################################

# Define runs for each depth
SNS_RUNS = ["run1", "run2"]
DNS_RUNS = ["run1", "run2"]

rule all:
    input:
        "results/alpha_diversity.png",
        "results/ordination_pcoa.png",
        "results/deseq2_results.csv"

#############################################
# Merge wf-16s outputs and build metadata
#############################################
rule merge_and_metadata:
    input:
        abund = expand("data/SNS/{run}/abundance_table_species.tsv", run=SNS_RUNS) +
                expand("data/DNS/{run}/abundance_table_species.tsv", run=DNS_RUNS),

        stats = expand("data/SNS/{run}/alignment-stats.tsv", run=SNS_RUNS) +
                expand("data/DNS/{run}/alignment-stats.tsv", run=DNS_RUNS)

    output:
        abund = "results/merged_abundance.tsv",
        stats = "results/merged_stats.tsv",
        meta  = "results/metadata.csv"

    conda:
        "workflow/envs/r.yaml"

    script:
        "workflow/scripts/merge_and_build_metadata.R"


#############################################
# Downstream analysis (paired + batch-aware)
#############################################
rule analysis:
    input:
        abund = "results/merged_abundance.tsv",
        stats = "results/merged_stats.tsv",
        meta  = "results/metadata.csv"

    output:
        alpha = "results/alpha_diversity.png",
        beta  = "results/ordination_pcoa.png",
        deseq = "results/deseq2_results.csv"

    conda:
        "workflow/envs/r.yaml"

    script:
        "workflow/scripts/analysis.R"

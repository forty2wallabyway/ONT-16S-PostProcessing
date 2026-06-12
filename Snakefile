#############################################
# Multi-run wf-16s integration pipeline
#############################################

RUNS = ["run1", "run2", "run3"]

rule all:
    input:
        "results/merged_abundance.tsv",
        "results/merged_stats.tsv",
        "results/alpha_diversity.png",
        "results/ordination_pcoa.png",
        "results/deseq2_results.csv"

# ---------------------------
# Merge multiple wf-16s outputs
# ---------------------------
rule merge_tables:
    input:
        abund = expand("data/{run}/abundance_table_species.tsv", run=RUNS),
        stats = expand("data/{run}/alignment-stats.tsv", run=RUNS)
    output:
        abund = "results/merged_abundance.tsv",
        stats = "results/merged_stats.tsv"
    conda:
        "workflow/envs/r.yaml"
    script:
        "workflow/scripts/merge_tables.R"

# ---------------------------
# Downstream analysis
# ---------------------------
rule analysis:
    input:
        abund = "results/merged_abundance.tsv",
        stats = "results/merged_stats.tsv",
        meta  = "data/metadata.csv"
    output:
        alpha = "results/alpha_diversity.png",
        beta  = "results/ordination_pcoa.png",
        deseq = "results/deseq2_results.csv"
    conda:
        "workflow/envs/r.yaml"
    script:
        "workflow/scripts/analysis.R"

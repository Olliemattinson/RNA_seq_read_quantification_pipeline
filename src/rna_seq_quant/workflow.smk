import os

from rna_seq_quant.functions.sum_transcript_to_gene import sum_transcript_to_gene
from rna_seq_quant.functions.merge_quants import merge_quants_for_exp

configfile: "data_config.yaml"


ENV_DIR = "envs"

_GENOME = config["genome"]

rule all:
    input:
        expand(
            "data/quants/{experiment}_total_gene_quant_tpm.txt",
            experiment=config["experiment"],
        ),
        expand(
            "data/quants/{experiment}_total_gene_quant_counts.txt",
            experiment=config["experiment"],
        ),
        expand(
            "data/diff_exp/{experiment}_DESeq2_LRT_results.csv",
            experiment=config["experiment"],
        ),
        expand(
            "fastqc/{experiment}_{reads}/{experiment}_{reads}_2_fastqc.html",
            experiment=config["experiment"],
            reads=config["reads"],
        ),
        expand("data/tpm_boxplots/{experiment}", experiment=config["experiment"]),
        expand("multiqc/{experiment}.html", experiment=config["experiment"]),
    wildcard_constraints:
        experiment="[^_]+_[^_]+",
        reads="[^_]+_[^_]+",


rule fastqc_primary_qc:
    input:
        reads1="data/reads/{experiment}_{reads}_1.fastq.gz",
        reads2="data/reads/{experiment}_{reads}_2.fastq.gz",
    output:
        "fastqc/{experiment}_{reads}/{experiment}_{reads}_2_fastqc.html",
    params:
        output_stem="fastqc/{experiment}_{reads}",
    threads: 10
    conda:
        os.path.join(ENV_DIR, "fastqc_env.yaml")
    wildcard_constraints:
        experiment="[^_]+_[^_]+",
        reads="[^_]+_[^_]+",
    shell:
        "mkdir -p {params.output_stem};"
        "fastqc {input.reads1} {input.reads2} --threads {threads} -o {params.output_stem}"


rule multiqc_combined_qc:
    input:
        expand(
            "fastqc/{experiment}_{reads}/{experiment}_{reads}_2_fastqc.html",
            experiment=config["experiment"],
            reads=config["reads"],
        )
    output:
        "multiqc/{experiment}.html",
    params:
        input_pattern="fastqc/*",
        output_dir="multiqc",
    conda:
        os.path.join(ENV_DIR, "multiqc_env.yaml")
    shell:
        "multiqc {params.input_pattern} -o {params.output_dir}"


rule salmon_index_transcriptome:
    input:
        "data/transcriptomes/{transcriptome}_transcript.fa",
    output:
        "data/transcriptomes/{transcriptome}_transcriptome_index",
    conda:
        os.path.join(ENV_DIR, "salmon_env.yaml")
    shell:
        "salmon index -t {input} -i {output}"


rule salmon_quantify_reads:
    input:
        reads1="data/reads/{experiment}_{reads}_1.fastq.gz",
        reads2="data/reads/{experiment}_{reads}_2.fastq.gz",
        transcriptome_index_dir=f"data/transcriptomes/{_GENOME}_transcriptome_index",
    output:
        "data/quants/{experiment}_{reads}_quant/quant.sf",
    params:
        output_dir="data/quants/{experiment}_{reads}_quant",
    threads: 10
    conda:
        os.path.join(ENV_DIR, "salmon_env.yaml")
    wildcard_constraints:
        experiment="[^_]+_[^_]+",
        reads="[^_]+_[^_]+",
    shell:
        "salmon quant -i {input.transcriptome_index_dir} -l A -1 {input.reads1} -2 "
        "{input.reads2} -p {threads} --validateMappings -o {params.output_dir}"


rule merge:
    input:
        expand(
            "data/quants/{experiment}_{reads}_quant/quant.sf",
            experiment=config["experiment"],
            reads=config["reads"],
        ),
    wildcard_constraints:
        experiment="[^_]+_[^_]+",
        reads="[^_]+_[^_]+",
    params:
        output_prefix="data/quants/{experiment}_total_transcript_quant_",
    output:
        tpm="data/quants/{experiment}_total_transcript_quant_tpm.txt",
        counts="data/quants/{experiment}_total_transcript_quant_counts.txt",
    run:
        merge_quants(
            input,
            params.output_prefix,
        )



rule sum_transcript_to_gene:
    input:
        tpm="data/quants/{experiment}_total_transcript_quant_tpm.txt",
        counts="data/quants/{experiment}_total_transcript_quant_counts.txt",
        annotation_info=f"data/annotation_info/{_GENOME}_annotation_info.tsv",
    output:
        tpm="data/quants/{experiment}_total_gene_quant_tpm.txt",
        counts="data/quants/{experiment}_total_gene_quant_counts.txt",
    run:
        sum_transcript_to_gene(
            input.counts,
            input.annotation_info,
            output.counts,
        )
        sum_transcript_to_gene(
            input.tpm,
            input.annotation_info,
            output.tpm,
        )


rule deseq2_diff_exp_analysis:
    input:
        "data/quants/{experiment}_total_gene_quant_counts.txt",
    conda:
        os.path.join(ENV_DIR, "deseq2_env.yaml")
    params:
        output_dir="data/diff_exp",
    output:
        "data/diff_exp/{experiment}_DESeq2_LRT_results.csv",
    shell:
        "Rscript diff_exp_analysis.R {input} {params.output_dir} {wildcards.experiment}"


rule tpm_boxplots:
    input:
        tpm="data/quants/{experiment}_total_gene_quant_tpm.txt",
        gene_list="data/Genes_of_interest.txt",
    conda:
        "envs/RNA_seq_read_quant_env.yaml"
    output:
        directory("data/tpm_boxplots/{experiment}"),
    shell:
        "mkdir -p data/tpm_boxplots/{wildcards.experiment};"
        "Rscript tpm_boxplot_generator.R {input.tpm} {input.gene_list} {output}"

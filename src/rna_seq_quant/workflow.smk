import os

from rna_seq_quant.functions.sum_transcript_to_gene import sum_transcript_to_gene
from rna_seq_quant.functions.merge_quants import merge_quants_for_exp

configfile: "data_config.yaml"


ENV_DIR = "envs"

_GENOME = config["genome"]

reads_prefix = "data/reads/{experiment}_{reads}"
reads1 = f"{reads_prefix}_1.fastq.gz"
reads2 = f"{reads_prefix}_2.fastq.gz"
multiqc_file = "multiqc/{experiment}.html"
deseq_results = "data/diff_exp/{experiment}_DESeq2_LRT_results.csv"
tpm_boxplots_dir = "data/tpm_boxplots/{experiment}"

rule all:
    input:
        expand(multiqc_file, experiment=config["experiment"]),
        expand(deseq_results, experiment=config["experiment"]),
        expand(tpm_boxplots_dir, experiment=config["experiment"]),
    wildcard_constraints:
        experiment="[^_]+_[^_]+",
        reads="[^_]+_[^_]+",


rule fastqc_primary_qc:
    input:
        reads1=reads1,
        reads2=reads2,
    output:
        "fastqc/{experiment}_{reads}/{experiment}_{reads}_fastqc.html",
    params:
        output_stem="fastqc/{experiment}_{reads}",
    threads: 10
    conda:
        os.path.join(ENV_DIR, "fastqc_env.yaml")
    shell:
        "mkdir -p {params.output_stem};"
        "fastqc {input.reads1} {input.reads2} --threads {threads} -o {params.output_stem}"


rule multiqc_combined_qc:
    input:
        expand(
            rules.fastqc_primary_qc.output,
            experiment=config["experiment"],
            reads=config["reads"],
        )
    output:
        multiqc_file,
    params:
        input_pattern="fastqc/*",
        output_dir="multiqc",
    conda:
        os.path.join(ENV_DIR, "multiqc_env.yaml")
    shell:
        "multiqc {params.input_pattern} -o {params.output_dir}"


rule salmon_index_transcriptome:
    input:
        f"data/transcriptomes/{_GENOME}_transcript.fa",
    output:
        directory(f"data/transcriptomes/{_GENOME}_transcriptome_index"),
    conda:
        os.path.join(ENV_DIR, "salmon_env.yaml")
    shell:
        "salmon index -t {input} -i {output}"


rule salmon_quantify_reads:
    input:
        reads1=reads1,
        reads2=reads2,
        transcriptome_index_dir=rules.salmon_index_transcriptome.output,
    output:
        "data/quants/{experiment}_{reads}_quant/quant.sf",
    params:
        output_dir="data/quants/{experiment}_{reads}_quant",
    threads: 10
    conda:
        os.path.join(ENV_DIR, "salmon_env.yaml")
    shell:
        "salmon quant -i {input.transcriptome_index_dir} -l A -1 {input.reads1} -2 "
        "{input.reads2} -p {threads} --validateMappings -o {params.output_dir}"


rule merge:
    input:
        expand(
            rules.salmon_quantify_reads.output,
            experiment=config["experiment"],
            reads=config["reads"],
        ),
    output:
        tpm="data/quants/{experiment}_total_transcript_quant_tpm.txt",
        counts="data/quants/{experiment}_total_transcript_quant_counts.txt",
    params:
        output_prefix="data/quants/{experiment}_total_transcript_quant_",
    run:
        merge_quants(
            input,
            params.output_prefix,
        )



rule sum_transcript_to_gene:
    input:
        tpm=rules.merge.output.tpm,
        counts=rules.merge.output.counts,
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
        rules.sum_transcript_to_gene.output.counts,
    output:
        deseq_results,
    conda:
        os.path.join(ENV_DIR, "deseq2_env.yaml")
    params:
        output_dir="data/diff_exp",
    shell:
        "Rscript diff_exp_analysis.R {input} {params.output_dir} {wildcards.experiment}"


rule tpm_boxplots:
    input:
        tpm=rules.sum_transcript_to_gene.output.tpm,
        gene_list="data/Genes_of_interest.txt",
    output:
        directory(tpm_boxplots_dir),
    conda:
        "envs/RNA_seq_read_quant_env.yaml"
    shell:
        "mkdir -p data/tpm_boxplots/{wildcards.experiment};"
        "Rscript tpm_boxplot_generator.R {input.tpm} {input.gene_list} {output}"

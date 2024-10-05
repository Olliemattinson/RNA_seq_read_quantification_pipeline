import os


configfile: "RNA_seq_read_quant_config.yaml"


ENV_DIR = "envs"


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
        directory(
            expand("data/tpm_boxplots/{experiment}", experiment=config["experiment"])
        ),
    wildcard_constraints:
        experiment="[^_]+_[^_]+",
        reads="[^_]+_[^_]+",


rule fastqc_primary_qc:
    input:
        #reads1='data/reads/{experiment}_{reads}_1.fastq',
        #reads2='data/reads/{experiment}_{reads}_2.fastq'
        reads1="data/reads/{experiment}_{reads}_1.fq.gz",
        reads2="data/reads/{experiment}_{reads}_2.fq.gz",
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
        "fastqc/{experiment}_{reads}/{experiment}_{reads}_2_fastqc.html",
    output:
        "multiqc/{experiment}.html",  # COMPLETE RULE + ADD TO RULE ALL
    conda:
        os.path.join(ENV_DIR, "multiqc_env.yaml")


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
        #reads1='data/reads/{experiment}_{reads}_1.fastq',
        #reads2='data/reads/{experiment}_{reads}_2.fastq',
        reads1="data/reads/{experiment}_{reads}_1.fq.gz",
        reads2="data/reads/{experiment}_{reads}_2.fq.gz",
        transcriptome_index_dir=expand(
            "data/transcriptomes/{transcriptome}_transcriptome_index",
            transcriptome=config["transcriptome"],
        ),
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
    conda:
        "envs/RNA_seq_read_quant_env.yaml"
    wildcard_constraints:
        experiment="[^_]+_[^_]+",
        reads="[^_]+_[^_]+",
    params:
        experiment_plus_reads=expand(
            "{experiment}_{reads}",
            experiment=config["experiment"],
            reads=config["reads"],
        ),
        path_prefix="data/quants/",
        path_suffix="_quant/quant.sf",
        output_prefix="data/quants/{experiment}_total_transcript_quant_",
    output:
        tpm="data/quants/{experiment}_total_transcript_quant_tpm.txt",
        counts="data/quants/{experiment}_total_transcript_quant_counts.txt",
    shell:
        "python RNA_seq_read_quant_merge_quants.py {params.path_prefix} {params.path_suffix} "
        "{params.output_prefix} {params.experiment_plus_reads}"


rule sum_transcript_to_gene:
    input:
        tpm="data/quants/{experiment}_total_transcript_quant_tpm.txt",
        counts="data/quants/{experiment}_total_transcript_quant_counts.txt",
    conda:
        "envs/RNA_seq_read_quant_env.yaml"
    params:
        transcriptome=config["transcriptome"],
    output:
        tpm="data/quants/{experiment}_total_gene_quant_tpm.txt",
        counts="data/quants/{experiment}_total_gene_quant_counts.txt",
    shell:
        "python RNA_seq_read_quant_sum_transcript_to_gene.py {params.transcriptome} "
        "{output.tpm} {input.tpm};"
        "python RNA_seq_read_quant_sum_transcript_to_gene.py {params.transcriptome} "
        "{output.counts} {input.counts};"


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
        "Rscript RNA_seq_read_quant_DESeq2.R {input} {params.output_dir} {wildcards.experiment}"


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
        "Rscript RNA_seq_read_quant_tpm_boxplot_generator.R {input.tpm} {input.gene_list} {output}"

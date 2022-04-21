configfile: 'RNA_seq_read_quant_config.yaml'

rule all:
    input:
        expand('data/quants/{experiment}_total_gene_quant.txt',experiment=config['experiment']),
        expand('fastqc/{experiment}_{reads}/{experiment}_{reads}_2_fastqc.html',experiment=config['experiment'],reads=config['reads'])
    wildcard_constraints:
        experiment='[^_]+_[^_]+',
        reads='[^_]+_[^_]+'

rule read_quality_report:
    input:
        reads1='data/reads/{experiment}_{reads}_1.fastq',
        reads2='data/reads/{experiment}_{reads}_2.fastq'
    threads: 10
    conda:
        'envs/RNA_seq_read_quant_env.yaml'
    wildcard_constraints:
        experiment='[^_]+_[^_]+',
        reads='[^_]+_[^_]+'
    params:
        output_stem='fastqc/{experiment}_{reads}'
    output:
        'fastqc/{experiment}_{reads}/{experiment}_{reads}_2_fastqc.html'
    shell:
        'mkdir -p {params.output_stem};'
        'fastqc {input.reads1} {input.reads2} --threads {threads} -o {params.output_stem}'

rule index_transcriptome:
    input:
        'data/transcriptomes/{transcriptome}_transcript.fa'
    conda:
        'envs/RNA_seq_read_quant_env.yaml'
    output:
        'data/transcriptomes/{transcriptome}_transcriptome_index'
    shell:
        'salmon index -t {input} -i {output}'

rule quantify_reads:
    input:
        reads1='data/reads/{experiment}_{reads}_1.fastq',
        reads2='data/reads/{experiment}_{reads}_2.fastq',
        transcriptome_index_dir=expand('data/transcriptomes/{transcriptome}_transcriptome_index',transcriptome=config['transcriptome'])
    threads: 10
    conda:
        'envs/RNA_seq_read_quant_env.yaml'
    wildcard_constraints:
        experiment='[^_]+_[^_]+',
        reads='[^_]+_[^_]+'
    params:
        output_dir='data/quants/{experiment}_{reads}_quant'
    output:
        'data/quants/{experiment}_{reads}_quant/quant.sf'
    shell:
        'salmon quant -i {input.transcriptome_index_dir} -l A -1 {input.reads1} -2 '
        '{input.reads2} -p {threads} --validateMappings -o {params.output_dir}'

rule merge:
    input:
        expand('data/quants/{experiment}_{reads}_quant/quant.sf',experiment=config['experiment'],reads=config['reads'])
    conda:
        'envs/RNA_seq_read_quant_env.yaml'
    wildcard_constraints:
        experiment='[^_]+_[^_]+',
        reads='[^_]+_[^_]+'
    params:
        experiment_plus_reads=expand('{experiment}_{reads}',experiment=config['experiment'],reads=config['reads']),
        path_prefix='data/quants/',
        path_suffix='_quant/quant.sf'
    output:
        'data/quants/{experiment}_total_transcript_quant.txt'
    shell:
        'python RNA_seq_read_quant_merge_quants.py {params.path_prefix} {params.path_suffix} '
        '{output} {params.experiment_plus_reads}'

rule sum_transcript_to_gene:
    input:
        'data/quants/{experiment}_total_transcript_quant.txt'
    conda:
        'envs/RNA_seq_read_quant_env.yaml'
    params:
        transcriptome=config['transcriptome']
    output:
        'data/quants/{experiment}_total_gene_quant.txt'
    shell:
        'python RNA_seq_read_quant_sum_transcript_to_gene.py {params.transcriptome} '
        '{output} {input}'

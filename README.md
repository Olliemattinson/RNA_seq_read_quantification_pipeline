# RNA_seq_read_quantification_pipeline

Snakemake pipeline for:
- Checking RNA-seq read quality using Fastqc
- Quantifying RNA-seq reads against a transcriptome using Salmon, and outputting both transcript- and gene-level quantification files in both counts and tpm.
- Carrying out differential expression analysis using DESeq2, comparing all conditions in the experiment.

#### To run pipeline:
- Place transcriptome file ```<transcriptome>_transcript.fa``` in ```data/transcriptomes/```
- Place paired read files ```<experiment>_<reads>_1.fastq``` and ```<experiment>_<reads>_2.fastq``` in ```data/reads/```
- Activate a conda environment with snakemake installed
- Run the following command:
```
snakemake -s RNA_seq_read_quant_snakefile.smk --cores 10 --use-conda
```
#### Output:
- Fastqc
  - Read quality reports in: ```fastqc/<experiment>_<reads>/```
- Salmon
  - Transcript-level quantification file (tpm) at: ```data/quants/<experiment>_total_transcript_quant_tpm.txt```
  - Transcript-level quantification file (counts) at: ```data/quants/<experiment>_total_transcript_quant_counts.txt```
  - Gene-level quantification file (tpm) at: ```data/quants/<experiment>_total_gene_quant_tpm.txt```
  - Gene-level quantification file (counts) at: ```data/quants/<experiment>_total_gene_quant_counts.txt```
- DESeq2
  - Likelihood ratio test across all conditions at ```data/diff_exp/<experiment>_DESeq2_LRT_results.csv```
  - Pairwise differential expression comparisons at ```data/diff_exp/<experiment>_DESeq2_<condition1>_vs_<condition2>.csv


#### Pipeline:

![plot](pipeline.svg)

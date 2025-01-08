# RNA_seq_read_quantification_pipeline

Snakemake pipeline for:
- Checking RNA-seq read quality using Fastqc
- Quantifying RNA-seq reads against a transcriptome using Salmon, and outputting both transcript- and gene-level quantification files in both counts and tpm.
- Carrying out differential expression analysis using DESeq2, comparing all conditions in the experiment.

## Installation
Create the `RNA_seq_read_quant_env` conda environment
```
conda env create -f environment.yaml
```
This will install all the 3rd party dependencies, as well as installing the local `rna-seq-quant` package in editable mode
Now activate the conda environment
```
conda activate RNA_seq_read_quant_env
```

## Usage

#### To run pipeline:
- Place transcriptome file ```<transcriptome>_transcript.fa``` in ```data/transcriptomes/```
- Place paired read files ```<experiment>_<reads>_1.fastq``` and ```<experiment>_<reads>_2.fastq``` in ```data/reads/```
- Place tab-delimited text file ```Genes_of_interest.txt``` with 'Gene_accession' and 'Gene_name' column, containing any genes that you want a boxplot of tpm expression levels for, in ```data```
- Activate a conda environment with snakemake installed
- Run the following command:
```
snakemake -s workflow.smk --cores 10 --use-conda
```
#### Output:
- Fastqc
  - Read quality reports in: ```fastqc/<experiment>_<reads>/```
- Salmon
  - Transcript-level quantification file (tpm) at: ```data/quants/<experiment>_total_transcript_quant_tpm.txt```
  - Transcript-level quantification file (counts) at: ```data/quants/<experiment>_total_transcript_quant_counts.txt```
  - Gene-level quantification file (tpm) at: ```data/quants/<experiment>_total_gene_quant_tpm.txt```
  - Gene-level quantification file (counts) at: ```data/quants/<experiment>_total_gene_quant_counts.txt```
  - Tpm boxplots for genes of interest in ```data/tpm_boxplots/<experiment>```
    - Note: currently the section of the pipeline that produces tpm boxplots will only run if the ```data/tpm_boxplots/<experiment>``` directory does not already exist
- DESeq2
  - Likelihood ratio test across all conditions at ```data/diff_exp/<experiment>_DESeq2_LRT_results.csv```
  - Pairwise differential expression comparisons at ```data/diff_exp/<experiment>_DESeq2_<condition1>_vs_<condition2>_IF.csv``` and ```data/diff_exp/<experiment>_DESeq2_<condition1>_vs_<condition2>_MC.csv```
    - IF = Using DESeq2 'Independent Filtering' method to automatically filter out low read-count genes in order to optimise the number of genes with adjusted p-values lower than significance threshold alpha after correcting for multiple comparisons.
    - MC = Using own 'Manual Cut-off' (Independent Filtering turned off) method to choose a cut-off count level for reads to be excluded from analysis.
    - Pairwise differential expression comparison files are the standard DESeq2 output files, but with added columns:
      - ```diff_expression``` = condition (if either) in pairwise comparison that shows significantly higher expression level
      - ```DE_BFC_#_comparisons``` = condition (if either) in parwise comparison that shows significantly higher expression level after correcting for multiple comparisons using Bonferroni correction. # ranges between 2 and the maximum number of pairwise comparisons between conditions in the dataset (i.e. 5 conditions = max 10 comparisons). Thus, you can choose the relevant number of comparisons, and thus the appropriate significance level, depending on how many between-condition comparisons you are interested in.


## Pipeline:

![plot](pipeline.svg)

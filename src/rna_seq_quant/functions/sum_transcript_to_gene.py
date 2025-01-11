import pandas as pd


def sum_transcript_to_gene(input_csv: str, output_csv: str, annotation_info: str):
    """
    Use annotation_info file to map transcript names to gene names, and sum transcript tpms and counts
    to gene.

    Args:
        input_csv (str): input csv containing transcript tpms and counts
        output_csv (str): output csv containing gene tpms and counts
        annotation_info (str): annotation info tsv that maps transcripts to genes
    """
    # Import dataframes
    input_df = pd.read_csv(input_csv, sep="\t")
    gene_transcript_mapping = pd.read_csv(
        annotation_info,
        sep="\t",
        usecols=["locusName", "transcriptName"],
    )

    # Merge dataframes to add gene names to expression dataframe
    merged_df = pd.merge(
        input_df,
        gene_transcript_mapping,
        left_on="Transcript_name",
        right_on="transcriptName",
        how="inner",
    )

    # Sum transcript tpm and counts to gene
    sum_columns = [
        col
        for col in merged_df.columns
        if col.endswith("transcript_counts") or col.endswith("transcript_TPM")
    ]
    grouped_df = merged_df.groupby(("locusName"))[sum_columns].sum().reset_index()

    # Rename transcript columns to gene columns
    grouped_df.columns = [
        col.replace("transcript_counts", "gene_counts").replace(
            "transcript_TPM", "gene_TPM"
        )
        for col in grouped_df.columns
    ]

    grouped_df.to_csv(output_csv, sep="\t", index=False)

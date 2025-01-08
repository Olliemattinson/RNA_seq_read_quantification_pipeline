import os
import pandas as pd
import sys
import copy


def _get_experiment_and_reads_from_file_path(file_path: str) -> tuple[str, str]:
    """
    Extract experiment and reads names from quantified reads directory name.

    Args:
        file_path (str): Path to quantified reads directory. Must be in the form path/to/{experiment}_{reads}_quant/quant.sf.
    Returns:
        experiment (str): Experiment name extracted from quantified reads directory.
        reads (str): Reads name extracted from quantified reads directory.
    """
    dir_name = os.path.basename(os.path.dirname(file_path))
    dir_name_parts = dir_name.split("_")
    if len(dir_name_parts) != 5 or dir_name_parts[4] != "quant":
        raise ValueError(
            f"Invalid quantified reads directory name: {dir_name}. Should be in the form '{{experiment}}_{{reads}}_quant'"
        )
    publication, species, condition, accession, _ = dir_name_parts
    experiment = f"{publication}_{species}"
    reads = f"{condition}_{accession}"

    return experiment, reads


def _get_condition_and_accession_from_reads_name(reads_name: str) -> tuple[str, str]:
    """
    Extract condition and accession from reads name.

    Args:
        reads_name (str): Must be in the form {condition}_{accession}.
    Returns:
        condition (str): Condition name extracted from reads name.
        accession (str): Accession extracted from reads name.
    """
    reads_name_parts = reads_name.split("_")
    if len(reads_name_parts) != 2:
        raise ValueError(
            f"Invalid reads name: {reads_name}. Should be in the form '{{condition}}_{{accession}}'"
        )
    condition, accession = reads_name_parts

    return condition, accession


def merge_quants_for_exp(input_quants: list[str], output_prefix: str):
    """
    Takes a list of quant.sf files, one for each set of reads. Outputs two combined TSV files.
    One contains transcript abundances in TPM, including average TPM columns for each condition.
    One contains transcript abundances in raw read counts.

    Args:
        input_quants (list[str]): Input quant.sf files
        output_prefix (str): Prefix for output files
    """
    conditions = set()
    df_tpm_merged = pd.DataFrame(columns=["Transcript_name"])
    df_counts_merged = pd.DataFrame(columns=["Transcript_name"])
    for input_csv_path in input_quants:
        # Extract reads and condition names
        _, reads = _get_experiment_and_reads_from_file_path(input_csv_path)
        condition, _ = _get_condition_and_accession_from_reads_name(reads)
        conditions.add(condition)

        df = pd.read_csv(input_csv_path)

        # Create TPM dataframe. Drop and rename columns. Merge.
        df_tpm = df.drop(["Length", "EffectiveLength", "NumReads"], axis=1)
        df_tpm = df_tpm.rename(
            columns={"Name": "Transcript_name", "TPM": f"{reads}_transcript_TPM"}
        )
        df_tpm_merged = pd.merge(
            df_tpm_merged, df_tpm, on="Transcript_name", how="outer"
        )
        df_tpm_merged.fillna(0.0, inplace=True)

        # Create counts dataframe. Drop and rename columns. Merge.
        df_counts = df.drop(["Length", "EffectiveLength", "TPM"], axis=1)
        count_col_name = f"{reads}_transcript_counts"
        df_counts = df_counts.rename(
            columns={
                "Name": "Transcript_name",
                "NumReads": count_col_name,
            }
        )
        df_counts_merged = pd.merge(
            df_counts_merged, df_counts, on="Transcript_name", how="outer"
        )
        df_counts_merged.fillna(0, inplace=True)
        df_counts_merged[count_col_name] = df_counts_merged[count_col_name].astype(int)

    # For TPM dataframe, create mean column per condition
    for condition in conditions:
        condition_cols = [
            col for col in df_tpm_merged.columns if col.startswith(condition)
        ]
        df_tpm_merged[f"{condition}_mean_transcript_TPM"] = df_tpm_merged[
            condition_cols
        ].mean(numeric_only=True, axis=1)

    # Write dataframes to CSV
    df_tpm_merged.to_csv(f"{output_prefix}_tpm.txt", sep="\t", mode="w", index=False)
    df_counts_merged.to_csv(
        f"{output_prefix}_counts.txt", sep="\t", mode="w", index=False
    )

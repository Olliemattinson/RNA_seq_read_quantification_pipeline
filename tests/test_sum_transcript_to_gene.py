import os
import tempfile
import unittest

import pandas as pd

from rna_seq_quant.functions.sum_transcript_to_gene import sum_transcript_to_gene

class Test_Sum_Transcript_To_Gene(unittest.TestCase):
    input_counts_df = pd.DataFrame(
        {
            "Transcript_name": ["transcript_1.1", "transcript_1.2", "transcript_2.1"],
            "conditionA_accession1_transcript_counts": [50, 100, 200],
            "conditionA_accession2_transcript_counts": [0, 200, 200],
        }
    )
    input_tpm_df = pd.DataFrame(
        {
            "Transcript_name": ["transcript_1.1", "transcript_1.2", "transcript_2.1"],
            "conditionA_accession1_transcript_TPM": [10.0, 20.0, 40.0],
            "conditionA_accession2_transcript_TPM": [0.0, 40.0, 40.0],
        }
    )
    annotation_info_df = pd.DataFrame(
        {
            "locusName": ["gene_1", "gene_1", "gene_2", "gene_3"],
            "transcriptName": ["transcript_1.1", "transcript_1.2", "transcript_2.1", "transcript_3.1"]
        }
    )
    expected_output_counts_df = pd.DataFrame(
        {
            "locusName": ["gene_1", "gene_2"],
            "conditionA_accession1_gene_counts": [150, 200],
            "conditionA_accession2_gene_counts": [200, 200],
        }
    )
    expected_output_tpm_df = pd.DataFrame(
        {
            "locusName": ["gene_1", "gene_2"],
            "conditionA_accession1_gene_TPM": [30.0, 40.0],
            "conditionA_accession2_gene_TPM": [40.0, 40.0],
        }
    )

    def test_sum_transcript_to_gene(self):
        for unit in ("counts", "tpm"):
            with tempfile.TemporaryDirectory() as tempdir:
                input_counts_csv = os.path.join(tempdir, f"input_{unit}.csv")
                self.input_counts_df.to_csv(input_counts_csv, sep="\t", index=False)
                annotation_info_csv = os.path.join(tempdir, "annotation_info.csv")
                self.annotation_info_df.to_csv(annotation_info_csv, sep="\t", index=False)
                output_counts_csv = os.path.join(tempdir, f"output_{unit}.csv")
                sum_transcript_to_gene(input_counts_csv, output_counts_csv, annotation_info_csv)
                pd.testing.assert_frame_equal(
                    self.expected_output_counts_df,
                    pd.read_csv(output_counts_csv, sep="\t"),
                )
            
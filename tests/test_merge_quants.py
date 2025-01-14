import unittest
import os
import tempfile
import pandas as pd

from rna_seq_quant.functions.merge_quants import merge_quants_for_exp


class Test_Merge_Quants(unittest.TestCase):
    input_quant_1 = {
        "Name": ["transcript_1", "transcript_2", "transcript_3"],
        "Length": [100, 100, 100],
        "EffectiveLength": [100, 100, 100],
        "NumReads": [50, 100, 200],
        "TPM": [10.0, 20.0, 40.0],
    }
    input_quant_2 = {
        "Name": ["transcript_2", "transcript_3"],
        "Length": [100, 100],
        "EffectiveLength": [100, 100],
        "NumReads": [200, 200],
        "TPM": [40.0, 40.0],
    }
    input_quant_3 = {
        "Name": ["transcript_1", "transcript_2", "transcript_3"],
        "Length": [100, 100, 100],
        "EffectiveLength": [100, 100, 100],
        "NumReads": [200, 100, 50],
        "TPM": [40.0, 20.0, 10.0],
    }
    input_quant_4 = {
        "Name": ["transcript_1", "transcript_2"],
        "Length": [100, 100],
        "EffectiveLength": [100, 100],
        "NumReads": [200, 100],
        "TPM": [40.0, 40.0],
    }
    expected_merged_quants_tpm = {
        "Transcript_name": ["transcript_1", "transcript_2", "transcript_3"],
        "conditionA_accession1_transcript_TPM": [10.0, 20.0, 40.0],
        "conditionA_accession2_transcript_TPM": [0.0, 40.0, 40.0],
        "conditionB_accession1_transcript_TPM": [40.0, 20.0, 10.0],
        "conditionB_accession2_transcript_TPM": [40.0, 40.0, 0.0],
        "conditionA_mean_transcript_TPM": [5.0, 30.0, 40.0],
        "conditionB_mean_transcript_TPM": [40.0, 30.0, 5.0],
    }
    expected_merged_quants_counts = {
        "Transcript_name": ["transcript_1", "transcript_2", "transcript_3"],
        "conditionA_accession1_transcript_counts": [50, 100, 200],
        "conditionA_accession2_transcript_counts": [0, 200, 200],
        "conditionB_accession1_transcript_counts": [200, 100, 50],
        "conditionB_accession2_transcript_counts": [200, 100, 0],
    }

    def test_merge_quants_for_exp(self):
        with tempfile.TemporaryDirectory() as tempdir:
            # Write input quants to temp csvs
            csv_list = []
            for dict, dir_name, file_name in (
                (
                    self.input_quant_1,
                    "publication_species_conditionA_accession1_quant",
                    "input_quant_1.csv",
                ),
                (
                    self.input_quant_2,
                    "publication_species_conditionA_accession2_quant",
                    "input_quant_2.csv",
                ),
                (
                    self.input_quant_3,
                    "publication_species_conditionB_accession1_quant",
                    "input_quant_3.csv",
                ),
                (
                    self.input_quant_4,
                    "publication_species_conditionB_accession2_quant",
                    "input_quant_4.csv",
                ),
            ):
                df = pd.DataFrame(dict)
                dir_path = os.path.join(tempdir, dir_name)
                os.makedirs(dir_path, exist_ok=True)
                file_path = os.path.join(dir_path, file_name)
                df.to_csv(file_path, index=False)
                csv_list.append(file_path)

            # Get output file paths
            output_prefix = os.path.join(tempdir, "output")
            output_tpm_file = f"{output_prefix}_tpm.txt"
            output_counts_file = f"{output_prefix}_counts.txt"

            # Run function
            merge_quants_for_exp(csv_list, output_prefix)

            # Read in output csvs
            expected_df_tpm_merged = pd.read_csv(output_tpm_file, sep="\t")
            expected_df_counts_merged = pd.read_csv(output_counts_file, sep="\t")

        # Assert output csvs are as expected
        pd.testing.assert_frame_equal(
            expected_df_tpm_merged,
            pd.DataFrame(self.expected_merged_quants_tpm),
            check_like=True,
        )
        pd.testing.assert_frame_equal(
            expected_df_counts_merged,
            pd.DataFrame(self.expected_merged_quants_counts),
            check_like=True,
        )

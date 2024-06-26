import unittest
import os
import tempfile
import pandas as pd

from functions.merge_quants import merge_quants_for_exp

class Test_Merge_Quants(unittest.TestCase):
    input_quant_1 = {
        "Name": ['transcript_1', 'transcript_2', 'transcript_3'],
        "Length": [100, 100, 100],
        "Effective_Length": [100, 100, 100],
        "NumReads": [50, 100, 200],
        "TPM": [10, 20, 40],
    }
    input_quant_2 = {
        "Name": ['transcript_2', 'transcript_3'],
        "Length": [100, 100],
        "Effective_Length": [100, 100],
        "NumReads": [200, 200],
        "TPM": [40, 40],
    }
    input_quant_3 = {
        "Name": ['transcript_1', 'transcript_2', 'transcript_3'],
        "Length": [100, 100, 100],
        "Effective_Length": [100, 100, 100],
        "NumReads": [200, 100, 50],
        "TPM": [40, 20, 10],
    }
    input_quant_4 = {
        "Name": ['transcript_1', 'transcript_2'],
        "Length": [100, 100],
        "Effective_Length": [100, 100],
        "NumReads": [200, 100],
        "TPM": [40, 40],
    }
    expected_merged_quants_tpm = {
        "Transcript_name": ['transcript_1', 'transcript_2', 'transcript_3'],
        "A_accession1_transcript_TPM": [10, 20, 40],
        "A_accession2_transcript_TPM": [0, 40, 40],
        "B_accession1_transcript_TPM": [40, 20, 10],
        "B_accession2_transcript_TPM": [40, 40, 0],
        "A_mean_transcript_TPM": [5, 30, 40],
        "B_mean_transcript_TPM": [40, 30, 5],
    }
    expected_merged_quants_counts = {
        "Transcript_name": ['transcript_1', 'transcript_2', 'transcript_3'],
        "A_accession1_transcript_counts": [50, 100, 200],
        "A_accession2_transcript_counts": [0, 200, 200],
        "B_accession1_transcript_counts": [200, 100, 50],
        "B_accession2_transcript_counts": [200, 100],
    }


    def test_merge_quants_for_exp(self):

        with tempfile.TemporaryDirectory() as tempdir:
            # Write input quants to temp csvs
            csv_list = []
            for dict, file_name in (
                (self.input_quant_1, "input_quant_1.csv"),
                (self.input_quant_2, "input_quant_2.csv"),
                (self.input_quant_3, "input_quant_3.csv"),
                (self.input_quant_4, "input_quant_4.csv"),
            ):
                dict.to_csv(file_name)
                csv_list.append(file_name)

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
        self.assertEqual(expected_df_tpm_merged, self.expected_merged_quants_tpm)
        self.assertEqual(expected_df_counts_merged, self.expected_merged_quants_counts)

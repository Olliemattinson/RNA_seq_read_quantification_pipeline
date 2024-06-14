import os
import numpy as np
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
    if len(dir_name_parts) != 3 or dir_name_parts[2] != "quant":
        raise ValueError(f"Invalid quantified reads directory name: {dir_name}. Should be in the form '{{experiment}}_{{reads}}_quant'")
    experiment, reads, _ = dir_name_parts

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
        raise ValueError(f"Invalid reads name: {reads_name}. Should be in the form '{{condition}}_{{accession}}'")
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
        df_tpm = df.drop(['Length','EffectiveLength','NumReads'], axis=1)
        df_tpm = df_tpm.rename(columns={'Name':'Transcript_name', 'TPM': f"{reads}_transcript_TPM"})
        df_tpm_merged = pd.merge(df_tpm_merged, df_tpm, on="Transcript_name", how="outer")

        # Create counts dataframe. Drop and rename columns. Merge.
        df_counts = df.drop(['Length', 'EffectiveLength', 'TPM'])
        df_counts = df_counts.rename(columns={'Name': 'Transcript_name', 'NumReads': f'{reads}_transcript_counts'})
        df_counts_merged = pd.merge(df_counts_merged, df_counts, on="Transcript_name", how="outer")

    # For TPM dataframe, create mean column per condition
    for condition in conditions:
        condition_cols = [col for col in df_tpm_merged.columns if col.startswith(condition)]
        df_tpm_merged[f"{condition}_mean_transcript_TPM"]=df_tpm_merged[condition_cols].mean(numeric_only=True, axis=1)
    
    # Write dataframes to CSV
    df_tpm_merged.to_csv(f"{output_prefix}_tpm.txt", sep='\t', mode='w', index=False)
    df_counts_merged.to_csv(f"{output_prefix}_couts.txt", sep='\t', mode="w", index=False)


# TODO: Test and then remove old code below



### Parse command line args
# Retrieve reads_name,file,condition and store in compound data structure in the format:
    # dict={condition1={reads_name1:file,reads_name2:file,...},condition2={reads_name1:file,reads_name2:file,...}}
# Copy dict to form two separate dictionaries, dictTPM and dictCounts, to be worked on separately
n=len(sys.argv)
dict={}
path_prefix=sys.argv[1]
path_suffix=sys.argv[2]
output_prefix=sys.argv[3]
for i in range(4,n):
    experiment_plus_reads=sys.argv[i]
    condition=experiment_plus_reads.split('_')[2]
    dict[condition]={}
for i in range(4,n):
    experiment_plus_reads=sys.argv[i]
    file=path_prefix+experiment_plus_reads+path_suffix
    reads_name='_'.join(experiment_plus_reads.split('_')[2:])
    condition=experiment_plus_reads.split('_')[2]
    dict[condition][reads_name]=file
dictTPM=copy.deepcopy(dict)
dictCounts=copy.deepcopy(dict)

### Generate TPM dataframe, including mean column for each condition
# Create and process dataframe from each file and store in compound structure in the format:
    # dict={condition1={reads_name1:df,reads_name2:df,...},condition2={reads_name1:df,reads_name2:df,...}}
# Create a merged dataframe for each condition with added mean_transcript_TPM column
# Create a final merged dataframe of all conditions, containing mean for each condition (dfFinal)
dfFinalTPM = pd.DataFrame(columns=['Transcript_name'])
for condition in dictTPM:
    for reads_name in dictTPM[condition]:
        file=dictTPM[condition][reads_name]
        dfTPM=pd.read_csv(file,sep='\t')
        dfTPM=dfTPM.drop(['Length','EffectiveLength','NumReads'],axis=1)
        dfTPM=dfTPM.rename(columns={'Name':'Transcript_name'})
        dfTPM=dfTPM.rename(columns={'TPM':reads_name+'_transcript_TPM'})
        dictTPM[condition][reads_name]=dfTPM
    dictTPM[condition]['Merged']=pd.DataFrame(columns=['Transcript_name'])
    for reads_name in dictTPM[condition]:
        dfTPM=dictTPM[condition][reads_name]
        if reads_name != 'Merged':
            dictTPM[condition]['Merged']=pd.merge(dictTPM[condition]['Merged'],dfTPM,on='Transcript_name',how='outer')
    for reads_name in dictTPM[condition]:
        dfTPM=dictTPM[condition][reads_name]
        if reads_name == 'Merged':
            dfTPM[condition+'_mean_transcript_TPM']=dfTPM.mean(numeric_only=True,axis=1)
            dfFinalTPM = pd.merge(dfFinalTPM,dfTPM,on='Transcript_name',how='outer')

### Generate counts dataframe, NOT including mean column for each condition
# Same as above for TPM, except:
    # Use dictCounts rather than dictTPM (same thing, but dictTPM has been modified)
    # Remove TPM column rather than NumReads column
    # Remove line for generating mean for each condition
dfFinalCounts = pd.DataFrame(columns=['Transcript_name'])
for condition in dictCounts:
    for reads_name in dictCounts[condition]:
        file=dictCounts[condition][reads_name]
        dfCounts=pd.read_csv(file,sep='\t')
        dfCounts=dfCounts.drop(['Length','EffectiveLength','TPM'],axis=1)
        dfCounts=dfCounts.rename(columns={'Name':'Transcript_name'})
        dfCounts=dfCounts.rename(columns={'NumReads':reads_name+'_transcript_counts'})
        dictCounts[condition][reads_name]=dfCounts
    dictCounts[condition]['Merged']=pd.DataFrame(columns=['Transcript_name'])
    for reads_name in dictCounts[condition]:
        dfCounts=dictCounts[condition][reads_name]
        if reads_name != 'Merged':
            dictCounts[condition]['Merged']=pd.merge(dictCounts[condition]['Merged'],dfCounts,on='Transcript_name',how='outer')
    for reads_name in dictCounts[condition]:
        dfCounts=dictCounts[condition][reads_name]
        if reads_name == 'Merged':
            #dfCounts[condition+'_mean_transcript_Counts']=dfCounts.mean(numeric_only=True,axis=1)
            dfFinalCounts = pd.merge(dfFinalCounts,dfCounts,on='Transcript_name',how='outer')


# Write to CSV
dfFinalTPM.to_csv(output_prefix+'tpm.txt',sep='\t',mode='w',index=False)
dfFinalCounts.to_csv(output_prefix+'counts.txt',sep='\t',mode='w',index=False)



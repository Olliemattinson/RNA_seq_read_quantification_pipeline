import numpy as np
import pandas as pd
import sys
import copy

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

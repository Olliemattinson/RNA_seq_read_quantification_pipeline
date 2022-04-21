import numpy as np
import pandas as pd
import sys

# Parse command line args
# Retrieve reads_name,file,condition and store in compound data structure in the format:
    # dict={condition1={reads_name1:file,reads_name2:file,...},condition2={reads_name1:file,reads_name2:file,...}}
n=len(sys.argv)
dict={}
path_prefix=sys.argv[1]
path_suffix=sys.argv[2]
output=sys.argv[3]
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

# Create and process dataframe from each file and store in compound structure in the format:
    # dict={condition1={reads_name1:df,reads_name2:df,...},condition2={reads_name1:df,reads_name2:df,...}}
# Create a merged dataframe for each condition with added mean_transcript_TPM column
# Create a final merged dataframe of all conditions, containing mean for each condition (dfFinal)
dfFinal = pd.DataFrame(columns=['Transcript_name'])
for condition in dict:
    for reads_name in dict[condition]:
        file=dict[condition][reads_name]
        df=pd.read_csv(file,sep='\t')
        df=df.drop(['Length','EffectiveLength','NumReads'],axis=1)
        df=df.rename(columns={'Name':'Transcript_name'})
        df=df.rename(columns={'TPM':reads_name+'_transcript_TPM'})
        dict[condition][reads_name]=df
    dict[condition]['Merged']=pd.DataFrame(columns=['Transcript_name'])
    for reads_name in dict[condition]:
        df=dict[condition][reads_name]
        if reads_name != 'Merged':
            dict[condition]['Merged']=pd.merge(dict[condition]['Merged'],df,on='Transcript_name',how='outer')
    for reads_name in dict[condition]:
        df=dict[condition][reads_name]
        if reads_name == 'Merged':
            df[condition+'_mean_transcript_TPM']=df.mean(numeric_only=True,axis=1)
            dfFinal = pd.merge(dfFinal,df,on='Transcript_name',how='outer')

# Write to CSV
dfFinal.to_csv(output,sep='\t',mode='w',index=False)

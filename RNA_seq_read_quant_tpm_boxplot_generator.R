library(ggplot2)
library(tidyverse)

args=commandArgs(trailingOnly=TRUE)
tpmFile = args[1]
goiFile = args[2]
outputDir = args[3]

# Import tpm dataframe
tpmDF <- read.csv(file=tpmFile,sep='\t',header=TRUE)

# Import genes_of_interest dataframe
goiDF <- read.csv(file=goiFile,sep='\t',header=TRUE)

# Remove mean columns
tpmDF <- tpmDF %>% select(-contains('mean'))

# Convert DF to long format
tpmDFLong <- tpmDF %>%
  gather(Condition, tpm, 2:ncol(tpmDF), na.rm = TRUE)

# Format run column to just condition
tpmDFLong$Condition <- sapply(strsplit(tpmDFLong$Condition,'_'),'[',1)

# Plot boxplot for each gene of interest
for (gene in goiDF$Gene_name)
{
  geneDF <- tpmDFLong[tpmDFLong$Gene_name==gene,]
  name <- goiDF$Name[goiDF$Gene_name==gene]
  graph_title <- paste(gene,name,sep='_')
  file_name <- paste(graph_title,'.png',sep='')
  file_path <- paste(outputDir,file_name,sep='/')
  p <- ggplot(geneDF,aes(x=Condition,y=tpm))+
    geom_boxplot()+
    labs(title=graph_title)+
    theme_bw()
  ggsave(file_path,p,width=1.5,height=1,units="in",scale=5)
}



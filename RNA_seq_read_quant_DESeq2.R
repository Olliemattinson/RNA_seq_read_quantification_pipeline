library("DESeq2")
library("tidyverse")

# Command line args
args=commandArgs(trailingOnly=TRUE)
countsFile = args[1]
outputDir = args[2]
experiment = args[3]

# Import counts file
countsDF <- read.table(file=countsFile,sep='\t',header=TRUE)
# Convert 'Gene' column to row names
countsDFRowNames <- countsDF[,1]
countsDF <- countsDF[,-1]
rownames(countsDF) <- countsDFRowNames
print(head(countsDF))

# Generate metadata file from countsDF
runs <- colnames(countsDF)
metadataDF <- data.frame(run=runs)
metadataDF$condition <- sapply(strsplit(metadataDF$run,'_'),'[',1)
conditions <- unique(metadataDF$condition)
print(metadataDF)

# Provide counts data and metadata file to DDSeq
countsMatrix <- round(as.matrix(countsDF))
dds <- DESeqDataSetFromMatrix(countData=countsMatrix,
                              colData=metadataDF,
                              design=~condition)

# Define summary statistics function
    # Prints summary statistics after application of results filters
summary_statistics <- function(res, filter) {
  number_DE_genes <- sum(res$padj < alpha, na.rm = TRUE)
  print(sprintf("Number of DE expressed genes for %s is:", filter))
  print(number_DE_genes)
  print(sprintf("Summary for results using %s:", filter))
  print(summary(res))
}

### Differential expression analysis
# Define minimum baseMean (mean of all normalised counts across all samples) for inclusion
cutoff<-2
# Define p-value threshold (alpha)
alpha<-0.05

# Likelihood ratio test to determine which genes show differential expression between at least 2 conditions
ddsLRT <- DESeq(dds,test='LRT',reduced=~1)
resLRT <- results(ddsLRT)
ddsLRTFiltered <- ddsLRT[resLRT$baseMean>cutoff, ]
resLRTFiltered <- results(ddsLRTFiltered, independentFiltering = FALSE)
write.csv(as.data.frame(resLRTFiltered),file=sprintf('%s/%s_DESeq2_LRT_results.csv',outputDir,experiment))

# Comparisons between all condition combinations, with multiple comparison correction
dds <- DESeq(dds)
for (condition1 in conditions)
{
  for (condition2 in conditions)
  {
    if (condition1!=condition2)
    {
      comparison=paste(condition1,'vs',condition2,sep='_')
      res <- results(dds,pAdjustMethod='BH',contrast=c('condition',condition1,condition2),alpha=alpha)
      ddsFiltered <- dds[res$baseMean>cutoff,]
      resFiltered <- results(ddsFiltered,independentFiltering=FALSE)
      resFiltered$diff_expression <- with(resFiltered, ifelse(resFiltered$log2FoldChange>0 & resFiltered$padj<alpha, condition1,ifelse(resFiltered$log2FoldChange<0 & resFiltered$padj<alpha, condition2,'')))
      write.csv(as.data.frame(resFiltered),file=sprintf('%s/%s_DESeq2_%s_results.csv',outputDir,experiment,comparison))
    }
  }
}





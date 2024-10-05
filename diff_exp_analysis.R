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
conditionN <- length(conditions)
maxComparisonN <- (((conditionN^2)+conditionN)/2)-conditionN
print(metadataDF)

# Provide counts data and metadata file to DDSeq
countsMatrix <- round(as.matrix(countsDF))
dds <- DESeqDataSetFromMatrix(countData=countsMatrix,
                              colData=metadataDF,
                              design=~condition)

# Define summary statistics function
    # Prints summary statistics after application of results filters
summary_statistics <- function(condition1, condition2, res, filter) {
  number_DE_genes <- sum(res$padj < alpha, na.rm = TRUE)
  print(sprintf("Number of DE expressed genes for %s vs. %s using %s is:", condition1, condition2, filter))
  print(number_DE_genes)
  print(sprintf("Summary for %s vs. %s. using %s:", condition1, condition2, filter))
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
write.csv(as.data.frame(resLRT),file=sprintf('%s/%s_DESeq2_LRT_results.csv',outputDir,experiment))

# Comparisons between all condition combinations, with multiple comparison correction
ddsPairwise <- DESeq(dds)

for (i in seq(1,length(conditions)-1))
{
  for (j in seq(i+1,length(conditions)))
  {
    condition1=conditions[i]
    condition2=conditions[j]
    comparison=paste(condition1,'vs',condition2,sep='_')
    
    # Results using DESeq2 Independent Filtering (IF)
      # i.e. automatically remove genes with baseMean below automatically determined threshold 
      # to maximise number of higher expression genes that reach significance threshold alpha 
      # after multiple testing
    resPairwiseIF <- results(ddsPairwise,pAdjustMethod='BH',contrast=c('condition',condition1,condition2),alpha=alpha)
    summary_statistics(condition1, condition2, resPairwiseIF,'Independent Filtering')
    # Add new columns to indicate direction of significant, differential expression for single up to maximum comparisons (appropriate corrected using Bonferroni)
    resPairwiseIF$diff_expression <- with(resPairwiseIF, ifelse(resPairwiseIF$log2FoldChange>0 & resPairwiseIF$padj<alpha, condition1,ifelse(resPairwiseIF$log2FoldChange<0 & resPairwiseIF$padj<alpha, condition2,'')))
    if (maxComparisonN>1)
    {
      for (comparisonN in seq(2,maxComparisonN))
      {
        alphaBF<-0.05/comparisonN
        colName<-paste('DE_BFC',comparisonN,'comparisons',sep='_')
        resPairwiseIF[,colName] <- with(resPairwiseIF, ifelse(resPairwiseIF$log2FoldChange>0 & resPairwiseIF$padj<alphaBF, condition1,ifelse(resPairwiseIF$log2FoldChange<0 & resPairwiseIF$padj<alphaBF, condition2,'')))
      }
    }
    # Clean up DF
    resPairwiseIF <- resPairwiseIF[!is.na(resPairwiseIF$log2FoldChange),]
    resPairwiseIF <- resPairwiseIF[!is.na(resPairwiseIF$padj),]
    # Order DF by log2FoldChange
    resPairwiseIF$log2FoldChange <- as.numeric(resPairwiseIF$log2FoldChange)
    resPairwiseIFOrdered <- resPairwiseIF[order(abs(resPairwiseIF$log2FoldChange),decreasing=TRUE),]
    # Write output file
    write.csv(as.data.frame(resPairwiseIFOrdered),file=sprintf('%s/%s_DESeq2_%s_results_IF.csv',outputDir,experiment,comparison))
    
    # Results using Manual, predefined baseMean Cut-off (MC)
    resPairwiseMC <- results(ddsPairwise,pAdjustMethod='BH',contrast=c('condition',condition1,condition2),independentFiltering=FALSE,alpha=alpha)
    ddsPairwiseMCFiltered <- ddsPairwise[resPairwiseMC$baseMean > cutoff,]
    resPairwiseMCFiltered <- results(ddsPairwiseMCFiltered,pAdjustMethod='BH',contrast=c('condition',condition1,condition2),independentFiltering=FALSE)
    summary_statistics(condition1, condition2, resPairwiseMCFiltered,'Manual Cut-off')
    # Add new columns to indicate direction of significant, differential expression for single up to maximum comparisons (appropriate corrected using Bonferroni)
    resPairwiseMCFiltered$diff_expression <- with(resPairwiseMCFiltered, ifelse(resPairwiseMCFiltered$log2FoldChange>0 & resPairwiseMCFiltered$padj<alpha, condition1,ifelse(resPairwiseMCFiltered$log2FoldChange<0 & resPairwiseMCFiltered$padj<alpha, condition2,'')))
    if (maxComparisonN>1)
    {
      for (comparisonN in seq(2,maxComparisonN))
      {
        alphaBF<-0.05/comparisonN
        colName<-paste('DE_BFC',comparisonN,'comparisons',sep='_')
        resPairwiseMCFiltered[,colName] <- with(resPairwiseMCFiltered, ifelse(resPairwiseMCFiltered$log2FoldChange>0 & resPairwiseMCFiltered$padj<alphaBF, condition1,ifelse(resPairwiseMCFiltered$log2FoldChange<0 & resPairwiseMCFiltered$padj<alphaBF, condition2,'')))
      }
    }
    # Clean up DF
    resPairedMCFiltered <- resPairwiseMCFiltered[!is.na(resPairwiseMCFiltered$padj),]
    # Order DF by log2FoldChange
    resPairwiseMCFiltered$log2FoldChange <- as.numeric(resPairwiseMCFiltered$log2FoldChange)
    resPairwiseMCFilteredOrdered <- resPairwiseMCFiltered[order(abs(resPairwiseMCFiltered$log2FoldChange),decreasing=TRUE),]
    # Write output file
    write.csv(as.data.frame(resPairwiseMCFilteredOrdered),file=sprintf('%s/%s_DESeq2_%s_results_MC.csv',outputDir,experiment,comparison))
  }
}




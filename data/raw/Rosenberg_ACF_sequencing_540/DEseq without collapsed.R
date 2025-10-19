# DEseq2 with replicates as individual samples


# DEseq on ACF raw data
library(useful)
library(DESeq2)
setwd("C:/Users/Administrator/Desktop/ACF paper/Rosenberg_ACF_sequencing_540")
genes <- read.csv("gene_metadata.csv", stringsAsFactors = F)
# read in counts for 2 replicates
raw_UMI1 <- read.csv("rawbcm1.csv", stringsAsFactors = F, check.names = F)
cleaned_UMI1 <- raw_UMI1[,3:42]
row.names(cleaned_UMI1) <- raw_UMI1$Target

raw_UMI2 <- read.csv("rawbcm2.csv", stringsAsFactors = F, check.names = F)
cleaned_UMI2 <- raw_UMI2[,3:42]
row.names(cleaned_UMI2) <- raw_UMI2$Target

coldata <- read.csv("sample_metadata.csv", stringsAsFactors = F, check.names = F)
row.names(coldata) <- coldata$Sample
coldata$Condition <- as.factor(coldata$Condition)
coldata$Tissue <- as.factor(coldata$Tissue)

# combine two sets
colnames(cleaned_UMI1) <- paste(colnames(cleaned_UMI1), "1", sep = "_")
colnames(cleaned_UMI2) <- paste(colnames(cleaned_UMI2), "2", sep = "_")
UMI <- cbind(cleaned_UMI1, cleaned_UMI2)
coldata <- rbind(coldata,coldata)
row.names(coldata)[1:40] <- colnames(cleaned_UMI1)
row.names(coldata)[41:80] <- colnames(cleaned_UMI2)

# filter UMI table, if a gene is expressed in < 50% samples, will drop. 329 genes kept
UMI <- UMI[rowSums(UMI>1) > 40,]

# work with epithelial samples first
coldata.e <- coldata[coldata$Tissue == 'Epithelial',]
dds.e <- DESeqDataSetFromMatrix(
  countData = UMI[,row.names(coldata.e)],
  colData = coldata.e,
  design = ~ Condition
)

dds.e<- DESeq(dds.e)
res.e <- results(dds.e)
res.e <- as.data.frame(res.e)
res.e$Gene <- mapvalues(row.names(res.e), from = genes$Target, to = genes$Gene, warn_missing = F )
res.e$Category <- mapvalues(row.names(res.e), from = genes$Target, to = genes$Category, warn_missing = F )
write.csv(res.e,"DESeq2_output.epi.non_collapsed.ACF_vs_Norm.csv")

# work with stromal samples
coldata.s <- coldata[coldata$Tissue == 'Stromal',]
dds.s <- DESeqDataSetFromMatrix(
  countData = UMI[,row.names(coldata.s)],
  colData = coldata.s,
  design = ~ Condition
)

dds.s<- DESeq(dds.s)
res.s <- results(dds.s)
res.s <- as.data.frame(res.s)
res.s$Gene <- mapvalues(row.names(res.s), from = genes$Target, to = genes$Gene, warn_missing = F )
res.s$Category <- mapvalues(row.names(res.s), from = genes$Target, to = genes$Category, warn_missing = F )

write.csv(res.s,"DESeq2_output.stromal.non_collapsed.ACF_vs_Norm.csv")

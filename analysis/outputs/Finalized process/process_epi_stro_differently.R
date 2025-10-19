# finalized DESeq2 process of ACF data, this script did 2 runs, on only epi or only stro

library(useful)
library(DESeq2)
library(plyr)
library(pheatmap)

# read in expression matrix with 2 replicates
setwd("C:/Users/Administrator/Desktop/ACF paper/Finalized process")
raw_UMI1 <- read.csv("rawbcm1.csv", stringsAsFactors = F, check.names = F)
cleaned_UMI1 <- raw_UMI1[,3:42]
row.names(cleaned_UMI1) <- raw_UMI1$Target

raw_UMI2 <- read.csv("rawbcm2.csv", stringsAsFactors = F, check.names = F)
cleaned_UMI2 <- raw_UMI2[,3:42]
row.names(cleaned_UMI2) <- raw_UMI2$Target

coldata <- read.csv("sample_metadata.csv", stringsAsFactors = F, check.names = F)
genes <- read.csv("gene_metadata.csv", stringsAsFactors = F, check.names = F)
row.names(coldata) <- coldata$Sample
coldata$Condition <- as.factor(coldata$Condition)
coldata$Tissue <- as.factor(coldata$Tissue)
coldata$Patient <- as.factor(coldata$Patient)
coldata$Age <- as.factor(coldata$Age)
coldata$group <- factor(paste(coldata$Tissue, coldata$Condition, sep = "_"))

# in coldata$group, we have EpithelialACF, EpithelialNormal, StromalACF, StromalNormal

# combine two sets
colnames(cleaned_UMI1) <- paste(colnames(cleaned_UMI1), "1", sep = "_")
colnames(cleaned_UMI2) <- paste(colnames(cleaned_UMI2), "2", sep = "_")
UMI <- cbind(cleaned_UMI1, cleaned_UMI2)
coldata <- rbind(coldata,coldata)
row.names(coldata)[1:40] <- colnames(cleaned_UMI1)
row.names(coldata)[41:80] <- colnames(cleaned_UMI2)


# DGE on epi
coldata.epi <- coldata[which(coldata$Tissue == 'Epithelial'),]
UMI.epi <- UMI[,row.names(coldata.epi)]
dds.epi <- DESeqDataSetFromMatrix(
  countData = UMI.epi,
  colData = coldata.epi,
  design = ~ Patient + Condition
)
dds.epi<- DESeq(dds.epi)
output.epi <- results(dds.epi)
output.epi$Gene <- mapvalues(row.names(output.epi), from = genes$Target, to = genes$Gene, warn_missing = F )

norm.epi.counts <- counts(dds.epi, normalized=TRUE)
lognorm.epi.counts <- log2(norm.epi.counts + 1)
write.csv(norm.epi.counts, "counts.epi.normalized.csv")
write.csv(lognorm.epi.counts, "counts.epi.lognormalized.csv")

merged.epi <- merge(colData(dds.epi), t(lognorm.epi.counts), by="row.names")
write.csv(merged.epi, "merged_counts.epi.lognormalized.csv")
write.csv(output.epi,"DESeq2_output.epi_only.ACF_vs_Norm.csv")




# DGE on stro
coldata.stro <- coldata[which(coldata$Tissue == 'Stromal'),]
UMI.stro <- UMI[,row.names(coldata.stro)]
dds.stro <- DESeqDataSetFromMatrix(
  countData = UMI.stro,
  colData = coldata.stro,
  design = ~ Patient + Condition
)
dds.stro<- DESeq(dds.stro)
output.stro <- results(dds.stro)
output.stro$Gene <- mapvalues(row.names(output.stro), from = genes$Target, to = genes$Gene, warn_missing = F )

norm.stro.counts <- counts(dds.stro, normalized=TRUE)
lognorm.stro.counts <- log2(norm.stro.counts + 1)
write.csv(norm.stro.counts, "counts.stro.normalized.csv")
write.csv(lognorm.stro.counts, "counts.stro.lognormalized.csv")
merged.stro <- merge(colData(dds.stro), t(lognorm.stro.counts), by="row.names")
write.csv(merged.stro, "merged_counts.stro.lognormalized.csv")

write.csv(output.stro,"DESeq2_output.stro_only.ACF_vs_Norm.csv")

# calculate normalized expression matrix for each
vsd.epi <- varianceStabilizingTransformation(dds.epi)
norm.epi <- assay(vsd.epi)
write.csv(norm.epi, "VST_normalized_count.epi_only.csv")
vsd.stro <- varianceStabilizingTransformation(dds.stro)
norm.stro <- assay(vsd.stro)
write.csv(norm.stro, "VST_normalized_count.stro_only.csv")

# gene-gene correlation
norm.epi <- norm.epi[order(row.names(norm.epi), decreasing = F),]
norm.stro <- norm.stro[order(row.names(norm.stro), decreasing = F),]
corr.epi <- cor(t(norm.epi))
write.csv(corr.epi, 'gene_corr.epi.csv')
corr.stro <- cor(t(norm.stro))
write.csv(corr.stro, 'gene_corr.stro.csv')
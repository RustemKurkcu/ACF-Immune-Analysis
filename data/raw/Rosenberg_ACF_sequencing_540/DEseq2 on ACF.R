
install.packages("htmltools")
library(htmltools)
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")


# DEseq on ACF raw data
library(useful)
library(DESeq2)
setwd("C:/Users/Administrator/Desktop/ACF paper/Rosenberg_ACF_sequencing_540")
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

# reaplicate 1 process
# study only Epithelial data, compare epithelial vs normal
coldata.epi <- coldata[which(coldata$Tissue == 'Epithelial'),]
count_tab.epi <- DESeqDataSetFromMatrix(countData = cleaned_UMI1[,which(colnames(cleaned_UMI1) %in% coldata.epi$Sample)], colData = coldata.epi, design = ~Condition)
count_tab.epi <- DESeq(count_tab.epi)
count_tab.epi.output <- results(count_tab.epi)
write.csv(count_tab.epi.output,'DESeq2_output.epithelial.rep1.ACF_vs_Norm.csv')

# study only Stromal data, compare epithelial vs normal
coldata.s <- coldata[which(coldata$Tissue == "Stromal"),]
count_tab.s <- DESeqDataSetFromMatrix(countData = cleaned_UMI1[,which(colnames(cleaned_UMI1) %in% coldata.s$Sample)], colData = coldata.s, design = ~Condition)
count_tab.s <- DESeq(count_tab.s)
count_tab.s.output <- results(count_tab.s)
write.csv(count_tab.s.output,'DESeq2_output.stromal.rep1.ACF_vs_Norm.csv')

# reaplicate 2 process
# study only Epithelial data, compare epithelial vs normal
coldata.epi <- coldata[which(coldata$Tissue == 'Epithelial'),]
count_tab.epi <- DESeqDataSetFromMatrix(countData = cleaned_UMI2[,which(colnames(cleaned_UMI2) %in% coldata.epi$Sample)], colData = coldata.epi, design = ~Condition)
count_tab.epi <- DESeq(count_tab.epi)
count_tab.epi.output <- results(count_tab.epi)
write.csv(count_tab.epi.output,'DESeq2_output.epithelial.rep2.ACF_vs_Norm.csv')

# study only Stromal data, compare epithelial vs normal
coldata.s <- coldata[which(coldata$Tissue == "Stromal"),]
count_tab.s <- DESeqDataSetFromMatrix(countData = cleaned_UMI2[,which(colnames(cleaned_UMI2) %in% coldata.s$Sample)], colData = coldata.s, design = ~Condition)
count_tab.s <- DESeq(count_tab.s)
count_tab.s.output <- results(count_tab.s)
write.csv(count_tab.s.output,'DESeq2_output.stromal.rep2.ACF_vs_Norm.csv')
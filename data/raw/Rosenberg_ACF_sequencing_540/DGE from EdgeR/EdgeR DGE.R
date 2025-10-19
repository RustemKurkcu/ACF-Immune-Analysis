# use EdgeR for analysis
# setwd("C:/Users/Administrator/Desktop/ACF paper/Rosenberg_ACF_sequencing_540/DGE from EdgeR")

library(edgeR)
library(plyr)

# read in raw expression matrix from replicate 1
raw_UMI1 <- read.csv("rawbcm1.csv", stringsAsFactors = F, check.names = F)
cleaned_UMI1 <- raw_UMI1[,3:42]
row.names(cleaned_UMI1) <- raw_UMI1$Target

# read in raw expression matrix from replicate 2
raw_UMI2 <- read.csv("rawbcm2.csv", stringsAsFactors = F, check.names = F)
cleaned_UMI2 <- raw_UMI2[,3:42]
row.names(cleaned_UMI2) <- raw_UMI2$Target

# read in gene metadata from ThermoFisher gene panel
# read in sample information, indicating patientID, tissue type and condition
genes <- read.csv("gene_metadata.csv", stringsAsFactors = F)
samples <- read.csv('sample_metadata.csv', stringsAsFactors = F)

# pool two replicates
colnames(cleaned_UMI1) <- paste(colnames(cleaned_UMI1), "1", sep = "_")
colnames(cleaned_UMI2) <- paste(colnames(cleaned_UMI2), "2", sep = "_")
UMI <- cbind(cleaned_UMI1, cleaned_UMI2)
samples <- rbind(samples, samples)
samples$Sample[1:40] <- paste(samples$Sample[1:40], "1", sep = "_")
samples$Sample[41:80] <- paste(samples$Sample[41:80], "2", sep = "_")

# filter UMI table, if a gene is expressed in < 50% samples, will drop, 329 genes kept
UMI <- UMI[rowSums(UMI>1) > 40,]
genes_kept <- raw_UMI1[,1:2]
genes_kept <- genes_kept[genes_kept$Target %in% row.names(UMI),]

# process epithelial samples
samples.e <- samples[samples$Tissue == 'Epithelial',]
UMI.e <- UMI[,samples.e$Sample]
y.e <- DGEList(counts = as.matrix(UMI.e), genes=genes_kept)
y.e$samples$lib.size <- colSums(y$counts)

# make design matrix
patient.e <- as.factor(samples.e$Patient)
condition.e <- as.factor(samples.e$Condition)
design.e <- model.matrix(~patient.e+condition.e)
rownames(design.e) <- colnames(y.e)

y.e <- calcNormFactors(y.e)
y.e <- estimateDisp(y.e, design.e, robust = T)
y.e$common.dispersion
plotBCV(y.e)
fit.e <- glmFit(y.e, design.e)
lrt.e <- glmLRT(fit.e)

# get FDR
top <- topTags(lrt.s, n=Inf)
top <- top$table
output.e <- lrt.e$table
output.e$Gene <- mapvalues(row.names(output.e), from = genes$Target, to = genes$Gene, warn_missing = F )
output.e$Category <- mapvalues(row.names(output.e), from = genes$Target, to = genes$Category, warn_missing = F )
output.e$FDR <- mapvalues(row.names(output.e), from = top$Target, to = top$FDR)
write.csv(output.e, "EdgeR_DGE.epithelial.fdr.csv")

# process stromal samples
samples.s <- samples[samples$Tissue == 'Stromal',]
UMI.s <- UMI[,samples.s$Sample]
y.s <- DGEList(counts = as.matrix(UMI.s), genes=genes_kept)
y.s$samples$lib.size <- colSums(y$counts)

# make design matrix
patient.s <- as.factor(samples.s$Patient)
condition.s <- as.factor(samples.s$Condition)
design.s <- model.matrix(~patient.s+condition.s)
rownames(design.s) <- colnames(y.s)

y.s <- calcNormFactors(y.s)
y.s <- estimateDisp(y.s, design.s, robust = T)
y.s$common.dispersion
plotBCV(y.s)
fit.s <- glmFit(y.s, design.s)
lrt.s <- glmLRT(fit.s)

# get FDR
top <- topTags(lrt.s, n=Inf)
top <- top$table
output.s <- lrt.s$table
output.s$Gene <- mapvalues(row.names(output.s), from = genes$Target, to = genes$Gene, warn_missing = F )
output.s$Category <- mapvalues(row.names(output.s), from = genes$Target, to = genes$Category, warn_missing = F )
output.s$FDR <- mapvalues(row.names(output.s), from = top$Target, to = top$FDR)
write.csv(output.s, "EdgeR_DGE.stromal.fdr.csv")

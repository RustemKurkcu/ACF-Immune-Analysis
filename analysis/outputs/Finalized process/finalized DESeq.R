# finalized DESeq2 process of ACF data
library(tidyr)
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
dds.all <- DESeqDataSetFromMatrix(
  countData = UMI,
  colData = coldata,
  design = ~ Patient + group
)

# get normalized counts
vsd <- varianceStabilizingTransformation(dds.all)
norm.all <- assay(vsd)
write.csv(norm.all, "VST_normalized_count.all.csv")

# PCA with all genes
pcadata.allgenes <- plotPCA(vsd, intgroup=c("Patient", "Condition", "Tissue"))
pca_output.allgenes <- pcadata.allgenes$data
write.csv(pca_output.allgenes, "PCA_output.all_sample.all_genes.csv")

# PCA with top 50 genes
pcadata.topgenes <- plotPCA(vsd, intgroup=c("Patient", "Condition", "Tissue"), ntop = 50)
pca_output.topgenes <- pcadata.topgenes$data
write.csv(pca_output.topgenes, "PCA_output.all_sample.50_genes.csv")

# run DGE
dds.all<- DESeq(dds.all)

# grab epithelial ACF vs epithelial Normal
output.epi <- results(dds.all, contrast=c("group","EpithelialACF","EpithelialNormal"))
output.str <- results(dds.all, contrast=c("group","StromalACF","StromalNormal"))
output.epi <- as.data.frame(output.epi)
output.str <- as.data.frame(output.str) 

output.epi$Gene <- mapvalues(row.names(output.epi), from = genes$Target, to = genes$Gene, warn_missing = F )
output.epi$Category <- mapvalues(row.names(output.epi), from = genes$Target, to = genes$Category, warn_missing = F )
output.str$Gene <- mapvalues(row.names(output.str), from = genes$Target, to = genes$Gene, warn_missing = F )
output.str$Category <- mapvalues(row.names(output.str), from = genes$Target, to = genes$Category, warn_missing = F )
write.csv(output.epi,"DESeq2_output.epi.ACF_vs_Norm.csv")
write.csv(output.str,"DESeq2_output.str.ACF_vs_Norm.csv")

# try another normalized counts
norm.all.counts <- counts(dds.all, normalized=TRUE)
lognorm.all.counts <- log2(norm.all.counts + 1)
write.csv(norm.all.counts, "counts.normalized.csv")
write.csv(lognorm.all.counts, "counts.lognormalized.csv")

# find detection rate
det.tab <- as.data.frame(matrix(nrow = nrow(UMI), ncol = 4))
colnames(det.tab) <- paste("det in", unique(coldata$group), sep = " ")

coldata.Epi.A <- coldata[which(coldata$group == "EpithelialACF"),]
UMI.Epi.A <- UMI[,row.names(coldata.Epi.A)]
det.tab$`det in EpithelialACF` <- rowSums(UMI.Epi.A != 0)

coldata.Str.A <- coldata[which(coldata$group == "StromalACF"),]
UMI.Str.A <- UMI[,row.names(coldata.Str.A)]
det.tab$`det in StromalACF` <- rowSums(UMI.Str.A != 0)

coldata.Epi.N <- coldata[which(coldata$group == "EpithelialNormal"),]
UMI.Epi.N <- UMI[,row.names(coldata.Epi.N)]
det.tab$`det in EpithelialNormal` <- rowSums(UMI.Epi.N != 0)

coldata.Str.N <- coldata[which(coldata$group == "StromalNormal"),]
UMI.Str.N <- UMI[,row.names(coldata.Str.N)]
det.tab$`det in StromalNormal` <- rowSums(UMI.Str.N != 0)
row.names(det.tab) <- row.names(UMI)
write.csv(det.tab, 'detaction in each category.csv')

# try to plot some genes
par(mfrow=c(2,3))
plotCounts(dds.all, gene = "CXCL9_149250", intgroup = "group")
plotCounts(dds.all, gene = "CCR7_133237", intgroup = "group")
plotCounts(dds.all, gene = "CCL21_432539", intgroup = "group")
plotCounts(dds.all, gene = "CCL22_79185", intgroup = "group")
plotCounts(dds.all, gene = "S100A8_239334", intgroup = "group")
plotCounts(dds.all, gene = "S100A9_174280", intgroup = "group")

# sample distance, did not look good
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Patient, vsd$Age, vsd$Tissue, vsd$Condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

####################################################
# use log2 transfromed counts for expression plot ##
####################################################
degene <- read.csv("DEseq2_significant.csv", stringsAsFactors = F)

#upregulated genes in epi
merged.epi.up <- merge(colData(dds.all), t(lognorm.all.counts[degene$Target.Epi.up,]), by="row.names")
gathered.epi.up <- gather(merged.epi.up, gene, expression, 9:ncol(merged.epi.up))

#rename the genes and sort by fold change
gathered.epi.up$gene <- mapvalues(gathered.epi.up$gene, from = degene$Target.Epi.up, to = degene$Gene.Epi.up)
gathered.epi.up$FC <- mapvalues(gathered.epi.up$gene, from = degene$Gene.Epi.up, to = degene$FC.Epi.up)
gathered.epi.up <- gathered.epi.up[order(gathered.epi.up$FC, decreasing = T),]

ggplot(gathered.epi.up, aes(Tissue, expression, fill=Condition)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  labs(x="Tissue", 
       y="Expression (log normalized counts)", 
       fill="Condition", 
       title="Top Results")
# plot manually saved
# save expression matrix
write.csv(merged.epi.up, "epi.upregulated.expression_matrix.lognorm.csv")

############################
#downregulated genes in epi#
############################
merged.epi.down <- merge(colData(dds.all), t(lognorm.all.counts[degene$Target.Epi.down[degene$Target.Epi.down != ""],]), by="row.names")
gathered.epi.down <- gather(merged.epi.down, gene, expression, 9:ncol(merged.epi.down))

#rename the genes and sort by fold change
gathered.epi.down$gene <- mapvalues(gathered.epi.down$gene, from = degene$Target.Epi.down, to = degene$Gene.Epi.down)
gathered.epi.down$FC <- mapvalues(gathered.epi.down$gene, from = degene$Gene.Epi.down, to = degene$FC.Epi.down)
gathered.epi.down <- gathered.epi.down[order(gathered.epi.down$FC, decreasing = T),]

ggplot(gathered.epi.down, aes(Tissue, expression, fill=Condition)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  labs(x="Tissue", 
       y="Expression (log normalized counts)", 
       fill="Condition", 
       title="Top Results")
# plot manually saved
# save expression matrix
write.csv(merged.epi.down, "epi.downregulated.expression_matrix.lognorm.csv")

##############################
#upregulated genes in stromal#
##############################
merged.stro.up <- merge(colData(dds.all), t(lognorm.all.counts[degene$Target.Stro.up[degene$Target.Stro.up != ""],]), by="row.names")
gathered.stro.up <- gather(merged.stro.up, gene, expression, 9:ncol(merged.stro.up))

#rename the genes and sort by fold change
gathered.stro.up$gene <- mapvalues(gathered.stro.up$gene, from = degene$Target.Stro.up, to = degene$Gene.Stro.up)
gathered.stro.up$FC <- mapvalues(gathered.stro.up$gene, from = degene$Gene.Stro.up, to = degene$FC.Stro.up)
gathered.stro.up <- gathered.stro.up[order(gathered.stro.up$FC, decreasing = T),]

ggplot(gathered.stro.up, aes(Tissue, expression, fill=Condition)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  labs(x="Tissue", 
       y="Expression (log normalized counts)", 
       fill="Condition", 
       title="Top Results")
# plot manually saved
# save expression matrix
write.csv(merged.stro.up, "stro.upregulated.expression_matrix.lognorm.csv")


################################
#downregulated genes in stromal#
################################
merged.stro.down <- merge(colData(dds.all), t(lognorm.all.counts[degene$Target.Stro.down[degene$Target.Stro.down != ""],]), by="row.names")
gathered.stro.down <- gather(merged.stro.down, gene, expression, 9:ncol(merged.stro.down))

#rename the genes and sort by fold change
gathered.stro.down$gene <- mapvalues(gathered.stro.down$gene, from = degene$Target.Stro.down, to = degene$Gene.Stro.down)
gathered.stro.down$FC <- mapvalues(gathered.stro.down$gene, from = degene$Gene.Stro.down, to = degene$FC.Stro.down)
gathered.stro.down <- gathered.stro.down[order(gathered.stro.down$FC, decreasing = T),]

ggplot(gathered.stro.down, aes(Tissue, expression, fill=Condition)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  labs(x="Tissue", 
       y="Expression (log normalized counts)", 
       fill="Condition", 
       title="Top Results")
# plot manually saved
# save expression matrix
write.csv(merged.stro.down, "stro.downregulated.expression_matrix.lognorm.csv")

#########################################
# perform PCA on only epithelium samples#
#########################################
coldata.epi <- coldata[coldata$Tissue == "Epithelial",]
UMI.epi <- UMI[,row.names(coldata.epi)]
dds.epi <- DESeqDataSetFromMatrix(
  countData = UMI.epi,
  colData = coldata.epi,
  design = ~ Patient + group
)

# get normalized counts
vsd.epi <- varianceStabilizingTransformation(dds.epi)
norm.epi <- assay(vsd.epi)
write.csv(norm.epi, "VST_normalized_count.epi.csv")

# PCA with all genes
pcadata.allgenes.epi <- plotPCA(vsd.epi, intgroup=c("Patient", "Condition", "Tissue"))
pca_output.allgenes.epi <- pcadata.allgenes.epi$data
write.csv(pca_output.allgenes.epi, "PCA_output.epi.all_genes.csv")

# PCA with top 50 genes
pcadata.topgenes.epi <- plotPCA(vsd.epi, intgroup=c("Patient", "Condition", "Tissue"), ntop = 50)
pca_output.topgenes.epi <- pcadata.topgenes.epi$data
write.csv(pca_output.topgenes.epi, "PCA_output.epi.50_genes.csv")


######################################
# perform PCA on only stromal samples#
######################################
coldata.str <- coldata[coldata$Tissue == "Stromal",]
UMI.str <- UMI[,row.names(coldata.str)]
dds.str <- DESeqDataSetFromMatrix(
  countData = UMI.str,
  colData = coldata.str,
  design = ~ Patient + group
)

# get normalized counts
vsd.str <- varianceStabilizingTransformation(dds.str)
norm.str <- assay(vsd.str)
write.csv(norm.str, "VST_normalized_count.str.csv")

# PCA with all genes
pcadata.allgenes.str <- plotPCA(vsd.str, intgroup=c("Patient", "Condition", "Tissue"))
pca_output.allgenes.str <- pcadata.allgenes.str$data
write.csv(pca_output.allgenes.str, "PCA_output.str.all_genes.csv")

# PCA with top 50 genes
pcadata.topgenes.str <- plotPCA(vsd.str, intgroup=c("Patient", "Condition", "Tissue"), ntop = 50)
pca_output.topgenes.str <- pcadata.topgenes.str$data
write.csv(pca_output.topgenes.str, "PCA_output.str.50_genes.csv")

# clean up expression matrix
colnames(merged.stro.down) <- mapvalues(colnames(merged.stro.down), from = genes$Target, to = genes$Gene, warn_missing = F)
colnames(merged.stro.up) <- mapvalues(colnames(merged.stro.up), from = genes$Target, to = genes$Gene, warn_missing = F)
colnames(merged.epi.down) <- mapvalues(colnames(merged.epi.down), from = genes$Target, to = genes$Gene, warn_missing = F)
colnames(merged.epi.up) <- mapvalues(colnames(merged.epi.up), from = genes$Target, to = genes$Gene, warn_missing = F)
write.csv(merged.stro.down, "stro.downregulated.expression_matrix.lognorm.cleaned.csv")
write.csv(merged.stro.up, "stro.upregulated.expression_matrix.lognorm.cleaned.csv")
write.csv(merged.epi.down, "epi.downregulated.expression_matrix.lognorm.cleaned.csv")
write.csv(merged.epi.up, "epi.upregulated.expression_matrix.lognorm.cleaned.csv")
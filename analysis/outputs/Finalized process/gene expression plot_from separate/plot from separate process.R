# expression plot revised

# gene expression over age
library(ggplot2)

# all data stored at merged.epi.up
# counts were normalized in pooled expression matrix (epi + stro), extracted by function counts()
# then log2 transformation with 1 pseudo count.
#################################
# all gene in epithelial#
#################################
setwd("C:/Users/Administrator/Desktop/ACF paper/Finalized process/gene expression plot_from separate")

counts.epi <- read.csv("merged_counts.epi.lognormalized.csv", stringsAsFactors = F, check.names = F)

colnames(counts.epi) <- mapvalues(colnames(counts.epi), from = genes$Target, to = genes$Gene, warn_missing = F)
setwd("C:/Users/Administrator/Desktop/ACF paper/Finalized process/gene expression plot_from separate/epi")
counts.epi[,1] <- NULL

for (i in 9:ncol(counts.epi)) {
  df <- cbind(counts.epi[,1:8], counts.epi[,i])
  
  # split 2 replicates
  df.1 <- df[seq(1, nrow(df), 2),]
  df.2 <- df[seq(2, nrow(df), 2),]
  
  # compute mean and stdev
  df.mean <- cbind(df.1[,2:6], (df.1[,9] + df.2[,9])/2)
  temp <- cbind(df.1[,9], df.2[,9])
  df.mean$stdev <- apply(temp, 1, sd)
  colnames(df.mean)[6] <- colnames(counts.epi)[i]
  
  # sort by age
  df.mean <- df.mean[order(df.mean$Age, df.mean$Patient),]
  
  # plot 
  Normalized_expression <- df.mean[,6]
  stdev <- df.mean[,7]
  
  df.mean$Patient <- as.factor(df.mean$Patient)
  ggplot(df.mean, aes(fill=Condition, y=Normalized_expression, x=Patient)) + 
    geom_bar(position="dodge", stat="identity", width=0.8) +
    ggtitle(colnames(counts.epi)[i]) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_errorbar(aes(ymin=Normalized_expression-stdev, ymax=Normalized_expression+stdev), width=.2,
                  position=position_dodge(.9))
  
  ggsave(paste(colnames(counts.epi)[i],"png", sep = "."), width = 5, height = 5)
  
}



#################################
# all gene in stromal#
#################################
setwd("C:/Users/Administrator/Desktop/ACF paper/Finalized process/gene expression plot_from separate")
counts.stro <- read.csv("merged_counts.stro.lognormalized.csv", stringsAsFactors = F, check.names = F)

colnames(counts.stro) <- mapvalues(colnames(counts.stro), from = genes$Target, to = genes$Gene, warn_missing = F)
setwd("C:/Users/Administrator/Desktop/ACF paper/Finalized process/gene expression plot_from separate/stro")
counts.stro[,1] <- NULL

for (i in 9:ncol(counts.stro)) {
  df <- cbind(counts.stro[,1:8], counts.stro[,i])
  
  # split 2 replicates
  df.1 <- df[seq(1, nrow(df), 2),]
  df.2 <- df[seq(2, nrow(df), 2),]
  
  # compute mean and stdev
  df.mean <- cbind(df.1[,2:6], (df.1[,9] + df.2[,9])/2)
  temp <- cbind(df.1[,9], df.2[,9])
  df.mean$stdev <- apply(temp, 1, sd)
  colnames(df.mean)[6] <- colnames(counts.stro)[i]
  
  # sort by age
  df.mean <- df.mean[order(df.mean$Age, df.mean$Patient),]
  
  # plot 
  Normalized_expression <- df.mean[,6]
  stdev <- df.mean[,7]
  
  df.mean$Patient <- as.factor(df.mean$Patient)
  ggplot(df.mean, aes(fill=Condition, y=Normalized_expression, x=Patient)) + 
    geom_bar(position="dodge", stat="identity", width=0.8) +
    ggtitle(colnames(counts.stro)[i]) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_errorbar(aes(ymin=Normalized_expression-stdev, ymax=Normalized_expression+stdev), width=.2,
                  position=position_dodge(.9))
  
  ggsave(paste(colnames(counts.stro)[i],"png", sep = "."), width = 5, height = 5)
  
}




########################################
# plot only significant genes          #
########################################




setwd("C:/Users/Administrator/Desktop/ACF paper/Finalized process/gene expression plot_from separate")
gene.meta <- read.csv("gene_metadata.csv", stringsAsFactors = F)
degene <- read.csv("DE genes.csv", stringsAsFactors = F)
degene <- apply(degene, 2, function(x) mapvalues(x, from = gene.meta$Gene, to = gene.meta$Target, warn_missing = F)) # convert to target name
degene <- as.data.frame(degene)

lognorm.epi <- read.csv("counts.epi.lognormalized.csv", stringsAsFactors = F, row.names = 1, check.names = F)
lognorm.stro <- read.csv("counts.stro.lognormalized.csv", stringsAsFactors = F, row.names = 1, check.names = F)
lognorm.all <- cbind(lognorm.epi, lognorm.stro)
sample.meta <- read.csv("sample_metadata.csv", stringsAsFactors = F)
# make sample metadata with 2 replicates
sample.meta.c1 <- sample.meta
sample.meta.c2 <- sample.meta
sample.meta.c1$Replication <- 1
sample.meta.c2$Replication <- 2
sample.meta.c1$Barcode <- paste(sample.meta.c1$Sample,"1",sep = "_")
sample.meta.c2$Barcode <- paste(sample.meta.c2$Sample,"2",sep = "_")
sample.meta.combine <- rbind(sample.meta.c1, sample.meta.c2)
row.names(sample.meta.combine) <- sample.meta.combine$Barcode


# epi upregulated
epi.up.genes <- as.character(degene$Epi.up[degene$Epi.up != ""])
merged.epi.up <- merge(sample.meta.combine, t(lognorm.all[epi.up.genes,]), by="row.names")
gathered.epi.up <- gather(merged.epi.up, gene, expression, 9:ncol(merged.epi.up))

#rename the genes 
gathered.epi.up$gene <- mapvalues(gathered.epi.up$gene, from = gene.meta$Target, to = gene.meta$Gene, warn_missing = F)

library(ggplot2)
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

# epi downregulated
epi.down.genes <- as.character(degene$Epi.down[degene$Epi.down != ""])
merged.epi.down <- merge(sample.meta.combine, t(lognorm.all[epi.down.genes,]), by="row.names")
gathered.epi.down <- gather(merged.epi.down, gene, expression, 9:ncol(merged.epi.down))

#rename the genes 
gathered.epi.down$gene <- mapvalues(gathered.epi.down$gene, from = gene.meta$Target, to = gene.meta$Gene, warn_missing = F)

library(ggplot2)
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



# stro upregulated
stro.up.genes <- as.character(degene$Stro.up[degene$Stro.up != ""])
merged.stro.up <- merge(sample.meta.combine, t(lognorm.all[stro.up.genes,]), by="row.names")
gathered.stro.up <- gather(merged.stro.up, gene, expression, 9:ncol(merged.stro.up))

#rename the genes 
gathered.stro.up$gene <- mapvalues(gathered.stro.up$gene, from = gene.meta$Target, to = gene.meta$Gene, warn_missing = F)

library(ggplot2)
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

# epi downregulated
stro.down.genes <- as.character(degene$Stro.down[degene$Stro.down != ""])
merged.stro.down <- merge(sample.meta.combine, t(lognorm.all[stro.down.genes,]), by="row.names")
gathered.stro.down <- gather(merged.stro.down, gene, expression, 9:ncol(merged.stro.down))

#rename the genes 
gathered.stro.down$gene <- mapvalues(gathered.stro.down$gene, from = gene.meta$Target, to = gene.meta$Gene, warn_missing = F)

library(ggplot2)
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

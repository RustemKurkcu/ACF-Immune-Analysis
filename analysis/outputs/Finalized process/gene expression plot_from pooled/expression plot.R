# gene expression over age
library(ggplot2)

# all data stored at merged.epi.up
# counts were normalized in pooled expression matrix (epi + stro), extracted by function counts()
# then log2 transformation with 1 pseudo count.
#################################
# upregulated gene in epithelial#
#################################
epi.ACF.up <- merged.epi.up[which(merged.epi.up$Tissue == "Epithelial"),]
colnames(epi.ACF.up) <- mapvalues(colnames(epi.ACF.up), from = genes$Target, to = genes$Gene, warn_missing = F)
setwd("C:/Users/Administrator/Desktop/ACF paper/Finalized process/gene expression plot")

for (i in 9:ncol(epi.ACF.up)) {
  df <- cbind(epi.ACF.up[,1:8], epi.ACF.up[,i])
  
  # split 2 replicates
  df.1 <- df[seq(1, nrow(df), 2),]
  df.2 <- df[seq(2, nrow(df), 2),]
  
  # compute mean and stdev
  df.mean <- cbind(df.1[,2:6], (df.1[,9] + df.2[,9])/2)
  temp <- cbind(df.1[,9], df.2[,9])
  df.mean$stdev <- apply(temp, 1, sd)
  colnames(df.mean)[6] <- colnames(epi.ACF.up)[i]
  
  # sort by age
  df.mean <- df.mean[order(df.mean$Age, df.mean$Patient),]

  # plot 
  Normalized_expression <- df.mean[,6]
  stdev <- df.mean[,7]
  
  ggplot(df.mean, aes(fill=Condition, y=Normalized_expression, x=Patient)) + 
    geom_bar(position="dodge", stat="identity", width=0.8) +
    ggtitle(colnames(epi.ACF.up)[i]) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_errorbar(aes(ymin=Normalized_expression-stdev, ymax=Normalized_expression+stdev), width=.2,
                  position=position_dodge(.9))
  
  ggsave(paste(colnames(epi.ACF.up)[i],"png", sep = "."), width = 5, height = 5)
  
}



#################################
# upregulated gene in stromal#
#################################
stro.ACF.up <- merged.stro.up[which(merged.stro.up$Tissue == "Stromal"),]
colnames(stro.ACF.up) <- mapvalues(colnames(stro.ACF.up), from = genes$Target, to = genes$Gene, warn_missing = F)
setwd("C:/Users/Administrator/Desktop/ACF paper/Finalized process/gene expression plot")

for (i in 9:ncol(stro.ACF.up)) {
  df <- cbind(stro.ACF.up[,1:8], stro.ACF.up[,i])
  
  # split 2 replicates
  df.1 <- df[seq(1, nrow(df), 2),]
  df.2 <- df[seq(2, nrow(df), 2),]
  
  # compute mean and stdev
  df.mean <- cbind(df.1[,2:6], (df.1[,9] + df.2[,9])/2)
  temp <- cbind(df.1[,9], df.2[,9])
  df.mean$stdev <- apply(temp, 1, sd)
  colnames(df.mean)[6] <- colnames(stro.ACF.up)[i]
  
  # sort by age
  df.mean <- df.mean[order(df.mean$Age, df.mean$Patient),]
  
  # plot 
  Normalized_expression <- df.mean[,6]
  stdev <- df.mean[,7]
  
  ggplot(df.mean, aes(fill=Condition, y=Normalized_expression, x=Patient)) + 
    geom_bar(position="dodge", stat="identity", width=0.8) +
    ggtitle(colnames(stro.ACF.up)[i]) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_errorbar(aes(ymin=Normalized_expression-stdev, ymax=Normalized_expression+stdev), width=.2,
                  position=position_dodge(.9))
  
  ggsave(paste(colnames(stro.ACF.up)[i],"png", sep = "."), width = 5, height = 5)
  
}



###################################
# downregulated gene in epithelial#
###################################
epi.ACF.down <- merged.epi.down[which(merged.epi.down$Tissue == "Epithelial"),]
colnames(epi.ACF.down) <- mapvalues(colnames(epi.ACF.down), from = genes$Target, to = genes$Gene, warn_missing = F)
setwd("C:/Users/Administrator/Desktop/ACF paper/Finalized process/gene expression plot")

for (i in 9:ncol(epi.ACF.down)) {
  df <- cbind(epi.ACF.down[,1:8], epi.ACF.down[,i])
  
  # split 2 replicates
  df.1 <- df[seq(1, nrow(df), 2),]
  df.2 <- df[seq(2, nrow(df), 2),]
  
  # compute mean and stdev
  df.mean <- cbind(df.1[,2:6], (df.1[,9] + df.2[,9])/2)
  temp <- cbind(df.1[,9], df.2[,9])
  df.mean$stdev <- apply(temp, 1, sd)
  colnames(df.mean)[6] <- colnames(epi.ACF.down)[i]
  
  # sort by age
  df.mean <- df.mean[order(df.mean$Age, df.mean$Patient),]
  
  # plot 
  Normalized_expression <- df.mean[,6]
  stdev <- df.mean[,7]
  
  ggplot(df.mean, aes(fill=Condition, y=Normalized_expression, x=Patient)) + 
    geom_bar(position="dodge", stat="identity", width=0.8) +
    ggtitle(colnames(epi.ACF.down)[i]) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_errorbar(aes(ymin=Normalized_expression-stdev, ymax=Normalized_expression+stdev), width=.2,
                  position=position_dodge(.9))
  
  ggsave(paste(colnames(epi.ACF.down)[i],"png", sep = "."), width = 5, height = 5)
  
}



#################################
# downregulated gene in stromal #
#################################
stro.ACF.down <- merged.stro.down[which(merged.stro.down$Tissue == "Stromal"),]
colnames(stro.ACF.down) <- mapvalues(colnames(stro.ACF.down), from = genes$Target, to = genes$Gene, warn_missing = F)
setwd("C:/Users/Administrator/Desktop/ACF paper/Finalized process/gene expression plot")

for (i in 9:ncol(stro.ACF.down)) {
  df <- cbind(stro.ACF.down[,1:8], stro.ACF.down[,i])
  
  # split 2 replicates
  df.1 <- df[seq(1, nrow(df), 2),]
  df.2 <- df[seq(2, nrow(df), 2),]
  
  # compute mean and stdev
  df.mean <- cbind(df.1[,2:6], (df.1[,9] + df.2[,9])/2)
  temp <- cbind(df.1[,9], df.2[,9])
  df.mean$stdev <- apply(temp, 1, sd)
  colnames(df.mean)[6] <- colnames(stro.ACF.down)[i]
  
  # sort by age
  df.mean <- df.mean[order(df.mean$Age, df.mean$Patient),]
  
  # plot 
  Normalized_expression <- df.mean[,6]
  stdev <- df.mean[,7]
  
  ggplot(df.mean, aes(fill=Condition, y=Normalized_expression, x=Patient)) + 
    geom_bar(position="dodge", stat="identity", width=0.8) +
    ggtitle(colnames(stro.ACF.down)[i]) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_errorbar(aes(ymin=Normalized_expression-stdev, ymax=Normalized_expression+stdev), width=.2,
                  position=position_dodge(.9))
  
  ggsave(paste(colnames(stro.ACF.down)[i],"png", sep = "."), width = 5, height = 5)
  
}

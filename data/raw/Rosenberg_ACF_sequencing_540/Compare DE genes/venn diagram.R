# compare DE genes from two approaches
library(VennDiagram)
library(RColorBrewer)
setwd("C:/Users/Administrator/Desktop/ACF paper/output_summary")
# read in DE genes from norm couts
DGE <- read.csv("collapsed DGE.csv", stringsAsFactors = F)

# get gene lists unique(x[x != ""])
up.epi.edger <- unique(DGE$Upregulated.epi.edger[DGE$Upregulated.epi.edger != ""]) 
down.epi.edger <- unique(DGE$Downregulated.epi.edger[DGE$Downregulated.epi.edger != ""])  
up.epi.deseq <- unique(DGE$Upregulated.epi.deseq2[DGE$Upregulated.epi.deseq2 != ""]) 
down.epi.deseq <- unique(DGE$Downregulated.epi.deseq2[DGE$Downregulated.epi.deseq2 != ""]) 
up.s.edger <- unique(DGE$Upregulated.s.edger[DGE$Upregulated.s.edger != ""]) 
down.s.edger <- unique(DGE$Downregulated.s.edger[DGE$Downregulated.s.edger != ""])  
up.s.deseq <- unique(DGE$Upregulated.s.deseq2[DGE$Upregulated.s.deseq2 != ""]) 
down.s.deseq <- unique(DGE$Downregulated.s.deseq2[DGE$Downregulated.s.deseq2 != ""])


myCol <- brewer.pal(4, "Pastel2")
# compare upregulated genes
venn.diagram(
  x = list(up.epi.edger,up.epi.deseq),
  category.names = c("EdgeR", "DEseq2"),
  filename = "Compare ACF Epi up_genes.png",
  output = T,
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c('pink','green'),
  # Numbers
  cex = .5,
  fontface = "bold",
  fontfamily = "sans",

)


# compare downregulated genes
venn.diagram(
  x = list(down.epi.edger,down.epi.deseq),
  category.names = c("EdgeR", "DEseq2"),
  filename = "Compare ACF Epi down_genes.png",
  output = T,
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c('pink','green'),
  # Numbers
  cex = .5,
  fontface = "bold",
  fontfamily = "sans",
  
)

# 4-way on epi
venn.diagram(
  x = list(up.epi.edger,up.epi.deseq,down.epi.edger,down.epi.deseq),
  category.names = c("Up.EdgeR", "Up.DESeq2", "Down.EdgeR", "Down.DESeq2"),
  filename = "Compare ACF Epi all_genes.png",
  output = T,
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .5,
  fontface = "bold",
  fontfamily = "sans",
  # Label font
  cat.cex = 0.3,
)


# 4-way on stromal
venn.diagram(
  x = list(up.s.edger,up.s.deseq,down.s.edger,down.s.deseq),
  category.names = c("Up.EdgeR", "Up.DESeq2", "Down.EdgeR", "Down.DESeq2"),
  filename = "Compare ACF Stromal all_genes.png",
  output = T,
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .5,
  fontface = "bold",
  fontfamily = "sans",
  # Label font
  cat.cex = 0.3,
)

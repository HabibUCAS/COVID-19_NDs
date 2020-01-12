# 1. Set your work folder
setwd("/home/habib/T2D/Server")

# 2. Load libraries for script two
library(GOSemSim)
library(readtext)
library(stringr)
library(factoextra)
library(dendextend)
library(corrplot)
library(RColorBrewer)


# 3. List the available correspondence files 

path <- c("T2D1_GSE20966_correspondence.txt","T2D2_GSE23343_correspondence.txt")
files <- readtext(path)


gene_grasp <- function(text_gene){
  aux_1 <- str_replace_all(text_gene,"GO:.{7} |\n|genes:", "")
  aux_2 <- strsplit(aux_1, " ")
  aux_3 <- strsplit(aux_2[[1]][2:6],",")
  aux_4 <- unique(unlist(aux_3))
  return(aux_4)
}

library(clusterProfiler)

corrisp <- lapply(gene_term_sub,function(x) bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))
corrisp <- lapply(corrisp,"[","ENTREZID")
corrisp <- lapply(corrisp, unlist)
names(corrisp) = substr(names(corrisp),1,4)

ccluster <- compareCluster(geneCluster = corrisp, fun = "enrichKEGG")
pdf(file = "enrich_KEGG.pdf", width = 15,height = 10)
dotplot(ccluster, font.size = 6)
dev.off()


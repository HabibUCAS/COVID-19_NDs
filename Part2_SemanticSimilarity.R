# 1. Set your work folder
setwd("/home/habib/T2D/Server")

library(clusterProfiler)

corrisp <- lapply(gene_term_sub,function(x) bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))
corrisp <- lapply(corrisp,"[","ENTREZID")
corrisp <- lapply(corrisp, unlist)
names(corrisp) = substr(names(corrisp),1,4)

ccluster <- compareCluster(geneCluster = corrisp, fun = "enrichKEGG")
pdf(file = "enrich_KEGG.pdf", width = 15,height = 10)
dotplot(ccluster, font.size = 6)
dev.off()

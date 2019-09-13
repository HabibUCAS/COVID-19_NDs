setwd("/home/habib/T2D/DB")

# 2. Load libraries for script one
library(RCurl)
library(GEOquery)
library(limma)
library(topGO)
library(genefilter)

#url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE11nnn/GSE11501/matrix/GSE84729_series_matrix.txt.gz" -P /home/habib/T2D/
#wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE11nnn/GSE11501/matrix/GSE84729_series_matrix.txt.gz' -P /home/habib/T2D/
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29231/matrix/"
#url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE11nnn/GSE11501/matrix/"
#url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE11nnn/GSE11501/matrix/"
#url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE11nnn/GSE11501/matrix/"
#dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
#print(dataset)

#for (ds in dataset){print(paste(url))
# }
gse <- getGEO(filename = "GSE29231_series_matrix.txt.gz",destdir=getwd())
d <- factor(c(rep('CTRL', 12),rep('T2D',12)))


#d <- factor(c('CTRL','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','CTRL','T2D','CTRL','T2D','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','CTRL','T2D','T2D'))
#print(d)
# [1] CVID   CVCTRL CVID   CVCTRL CVID   CVCTRL CVID   CVID   CVCTRL CVID   CVID   CVID   CVCTRL CVID   CVID   CVID   CVCTRL CVCTRL
# [19] CVCTRL CVID   CD     CD     CD     CD     CD     CD     CD     CD     CD     CD     CTRL   CTRL   CTRL   CTRL   CTRL   CTRL  
# [37] CTRL   CTRL   CTRL   CTRL   CTRL   CTRL   CTRL   CTRL   CTRL   CTRL   CTRL  
# Levels: CD CTRL CVCTRL CVID



# ok d <- factor(c(rep('CVID',20),rep('CD',10),rep('CTRL',17)))

mod <- model.matrix(~0+d)
# print(mod)
# dCD dCTRL dCVCTRL dCVID
# 1    0     0       0     1
# 2    0     0       1     0

fit_1 <- lmFit(gse, mod)
# #print(fit_1)
contr <- makeContrasts(dCTRL-dT2D,levels = mod)
#contr <- makeContrasts(dCD-dCTRL, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)
# 
# #> colnames(fit_3)
# #[1] "dT2D - dCTRL"   "dTREAT - dT2D"  "dTREAT - dCTRL"
# #print(fit_3)
# #show(fit_3)
# #n=length(fit_3)
# #print(n)
# #topTable(fit_3, coef=1, adjust="BH") #coef=1 means dT2D-dCTRL
# #topTable(fit_3, coef=2, adjust="BH") #coef=1 means dTRAET-dT2D
# 
#table_result <- topTable(fit_3, coef=1,p.value=1e-2, sort.by = "logFC")
#table_result <- topTable(fit_3, coef=1,p.value=1e-2, sort.by = "logFC")
#table_result <- topTable(fit_3, coef=1, n=Inf,p.value=5e-2, sort.by = "logFC")
#dim(table_result)
#[1] 156  32
table_result <- topTable(fit_3, coef=1,n=Inf,adjust="BH", sort.by = "logFC")
#> dim(table_result)
#[1] 18981    32
#[1] 156  32
#dim(table_result)
# # #[1]  9 11
subtable_result <- subset(table_result, select=c("ID","Symbol","adj.P.Val","P.Value","t","B","logFC"))

subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
# dim(subtable_result)
#[1] 16930     7
geneList <- subtable_result$logFC
#n=length(geneList)
#print(n)
#[1] 16930
names(geneList) <- subtable_result$Symbol
# # > print(names(geneList))
# names of 16930 genes
write.csv2(subtable_result,"GSE29231_table.csv")
nrow(subtable_result[subtable_result$P.Value<0.05,])
#[1] 1


topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}
# 
T2D_GOdata <- new("topGOdata",
                  description = "t2d study",
                  ontology = "BP", 
                  allGenes = geneList,
                  geneSel = topDiffGenes,
                  annot = annFUN.org, 
                  ID = "Symbol", 
                  mapping = "org.Hs.eg.db",
                  nodeSize = 10)
# Building most specific GOs .....
# ( 8180 GO terms found. )
# 
# Build GO DAG topology ..........
# ( 12438 GO terms and 29096 relations. )
# 
# Annotating nodes ...............
# ( 4905 genes annotated to the GO terms. )
# 
#print(CD_GOdata)
# ------------------------- topGOdata object -------------------------
#   
#   Description:
#   -  cd study 
# 
# Ontology:
#   -  BP 
# 
# 16930 available genes (all genes from the array):
#   - symbol:  HLA-DRB1 HLA-DRB5 CHURC1 CSF3R GZMH  ...
# - score :  -1.39492 -0.830915 0.4445 -0.424 0.4  ...
# - 1  significant genes. 
# 
# 12829 feasible genes (genes that can be used in the analysis):
#   - symbol:  HLA-DRB1 HLA-DRB5 CHURC1 CSF3R GZMH  ...
# - score :  -1.39492 -0.830915 0.4445 -0.424 0.4  ...
# - 1  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 5831 
# - number of edges = 13072 
# 
# ------------------------- topGOdata object -------------------------
#   
n_sg <- sum(topDiffGenes(geneList))
# print(n_sg)
# [1] 1
sg <- sigGenes(T2D_GOdata)
# print(sg)
#[1] "HLA-DRB1"
ug <- usedGO(T2D_GOdata)
#k=length(ug)
# print(k)
#[1] 5831 GO names
resultFisher <- runTest(T2D_GOdata, algorithm = "classic", statistic = "fisher")
# -- Classic Algorithm -- 
#   
#   the algorithm is scoring 54 nontrivial nodes
# parameters: 
#   test statistic: fisher

#print(resultFisher)
# Description: cd study 
# Ontology: BP 
# 'classic' algorithm with the 'fisher' test
# 5831 GO terms scored: 10 terms with p < 0.01
# Annotation data:
#   Annotated genes: 12829 
# Significant genes: 1 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 54 


resultKS <- runTest(T2D_GOdata, algorithm = "classic", statistic = "ks")

# -- Classic Algorithm -- 
#   
#   the algorithm is scoring 5831 nontrivial nodes
# parameters: 
#   test statistic: ks
# score order: increasing

#print(resultKS)
# Description: cd study 
# Ontology: BP 
# 'classic' algorithm with the 'ks' test
# 5831 GO terms scored: 1630 terms with p < 0.01
# Annotation data:
#   Annotated genes: 12829 
# Significant genes: 1 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 5831 


allRes <- GenTable(T2D_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 30)

# > print(allRes)
# GO.ID                                        Term Annotated Significant Expected Rank in KS classic      KS
# 1  GO:0060333 interferon-gamma-mediated signaling path...        59           1     0.00        884  0.0046 0.00056
# 2  GO:0019886 antigen processing and presentation of e...        68           1     0.01       1949  0.0053 0.02001
# 3  GO:0002495 antigen processing and presentation of p...        70           1     0.01       1789  0.0055 0.01456
# 4  GO:0002504 antigen processing and presentation of p...        71           1     0.01       1568  0.0055 0.00865
# 5  GO:0002478 antigen processing and presentation of e...        81           1     0.01       2073  0.0063 0.02611
# 6  GO:0019884 antigen processing and presentation of e...        87           1     0.01       1941  0.0068 0.01990
# 7  GO:0048002 antigen processing and presentation of p...        88           1     0.01       1727  0.0069 0.01260
# 8  GO:0050852           T cell receptor signaling pathway       102           1     0.01        492  0.0080 5.0e-06
# 9  GO:0019882         antigen processing and presentation       114           1     0.01        877  0.0089 0.00054
# 10 GO:0071346       cellular response to interferon-gamma       121           1     0.01        950  0.0094 0.00086
# 11 GO:0050851 antigen receptor-mediated signaling path...       135           1     0.01        409  0.0105 2.6e-07
# 12 GO:0034341                response to interferon-gamma       136           1     0.01        963  0.0106 0.00093
# 13 GO:0002429 immune response-activating cell surface ...       217           1     0.02        266  0.0169 3.0e-10
# 14 GO:0002768 immune response-regulating cell surface ...       242           1     0.02        276  0.0189 5.9e-10
# 15 GO:0002757 immune response-activating signal transd...       322           1     0.03        259  0.0251 2.5e-10
# 16 GO:0002764 immune response-regulating signaling pat...       348           1     0.03        251  0.0271 1.8e-10
# 17 GO:0002253               activation of immune response       376           1     0.03        302  0.0293 2.7e-09
# 18 GO:0050778      positive regulation of immune response       493           1     0.04        298  0.0384 2.5e-09
# 19 GO:0019221         cytokine-mediated signaling pathway       511           1     0.04        284  0.0398 1.0e-09
# 20 GO:0045087                      innate immune response       551           1     0.04        517  0.0429 7.1e-06
# 21 GO:0050776               regulation of immune response       643           1     0.05        301  0.0501 2.6e-09
# 22 GO:0034622 cellular protein-containing complex asse...       644           1     0.05        487  0.0502 4.1e-06
# 23 GO:0002684 positive regulation of immune system pro...       693           1     0.05        161  0.0540 7.8e-13
# 24 GO:0071345      cellular response to cytokine stimulus       732           1     0.06        212  0.0571 4.0e-11
# 25 GO:0034097                        response to cytokine       799           1     0.06        168  0.0623 1.2e-12
# 26 GO:0002682         regulation of immune system process      1001           1     0.08        133  0.0780 6.0e-15
# 27 GO:0006952                            defense response      1071           1     0.08        190  0.0835 8.4e-12
# 28 GO:0065003         protein-containing complex assembly      1209           1     0.09        349  0.0942 4.6e-08
# 29 GO:0006955                             immune response      1359           1     0.11         38  0.1059 2.9e-23
# 30 GO:0043933 protein-containing complex subunit organ...      1427           1     0.11        417  0.1112 3.5e-07

# showSigOfNodes(CD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
# $dag
# A graphNEL graph with directed edges
# Number of Nodes = 33 
# Number of Edges = 51 
# 
# $complete.dag
# [1] "A graph with 33 nodes."

printGraph(T2D_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "T2D1_GSE29231", useInfo = "all", pdfSW = TRUE)
terms <- allRes$GO.ID
# print(terms)
genes <- genesInTerm(T2D_GOdata,terms)
# print(genes)
for (i in 1:length(terms))
{
  term <- terms[i]
  genes_term <- genes[term][[1]]
  # find the genes that are in the list of genes of interest
  fact <- genes_term %in% sg
  genes_term_2 <- genes_term[fact == TRUE]
  genes_term_2 <- paste(genes_term_2, collapse=',')
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "T2D1_GSE29231_correspondence.txt" )
}

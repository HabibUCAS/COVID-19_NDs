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


subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]

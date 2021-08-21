library(gplots)
library(readxl)
library(plyr)

Heatmap <- read.csv("Heatmap.csv",row.names = 1)
df<-data.matrix(Heatmap)

heatmap.2(df, scale = "none", col = bluered(200), key.xlab = "logFC",
          trace = "none", density.info = "none",cexCol=1 #,cexRow = 1
          )


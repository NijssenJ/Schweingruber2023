library(tidyverse)
library(Vennerable)

load("Sig_MNs.Rdata")
MNsig <- Sig 

load("Sig_v2a.Rdata")
v2asig <- Sig

load("Sig_nonv2a.Rdata")
intersig <- Sig

p <- 0.05

inter <- intersig$TDP43M337V
list1 <- rownames(inter[which(inter$padj < p & inter$log2FoldChange > 0),])
list2 <- rownames(inter[which(inter$padj < p & inter$log2FoldChange < 0),])

motor <- MNsig$TDP43M337V
list3 <- rownames(motor[which(motor$padj < p & motor$log2FoldChange > 0),])
list4 <- rownames(motor[which(motor$padj < p & motor$log2FoldChange < 0),])

v2a <- v2asig$TDP43M337V
list5 <- rownames(v2a[which(v2a$padj < p & v2a$log2FoldChange > 0),])
list6 <- rownames(v2a[which(v2a$padj < p & v2a$log2FoldChange < 0),])


list.up <- list( "Motorneurons" = list3, "V2a interneurons" = list5, "Non-V2a interneurons"= list1)
list.down <- list("Motorneurons" = list4, "V2a interneurons" = list6, "Non-V2a interneurons"= list2)

venn.up <- Venn(list.up)
venn.down <- Venn(list.down)


plot(venn.up)


## Exporting excel data on Venn intersection lists
vennobj <- venn.up

maxintersect <- max(sapply(attr(vennobj, "IntersectionSets"), length))
vennint <- data.frame(row.names = 1:maxintersect)
for(i in 2:length(attr(vennobj, "IntersectionSets"))){
  vennint <- cbind(vennint, c(attr(vennobj, "IntersectionSets")[[i]], rep(NA, maxintersect - length(attr(vennobj, "IntersectionSets")[[i]]))))
}

topdf <- t(as.data.frame(vennobj@IndicatorWeight))[,-1]
rownames(topdf)[4] <- "counts"
colnames(vennint) <- colnames(topdf)

vennint <- rbind(topdf, vennint)

down <- vennint


WriteXLS::WriteXLS(list(up,down), ExcelFileName = "TARDBP_M337V_celltypes.xlsx", col.names = F, row.names = T, SheetNames = c("upregulated", "downregulated"))


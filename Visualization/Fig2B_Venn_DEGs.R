library(tidyverse)
library(Vennerable)


load("Sig_MNs.Rdata")

p <- 0.05

FUSKO <- Sig$FUSKO
list1 <- rownames(FUSKO[which(FUSKO$padj < p & FUSKO$log2FoldChange > 0),])
list2 <- rownames(FUSKO[which(FUSKO$padj < p & FUSKO$log2FoldChange < 0),])

FUSP525Lho <- Sig$FUSP525Lho
list3 <- rownames(FUSP525Lho[which(FUSP525Lho$padj < p & FUSP525Lho$log2FoldChange > 0),])
list4 <- rownames(FUSP525Lho[which(FUSP525Lho$padj < p & FUSP525Lho$log2FoldChange < 0),])

FUSR495X <- Sig$FUSR495X
list5 <- rownames(FUSR495X[which(FUSR495X$padj < p & FUSR495X$log2FoldChange > 0),])
list6 <- rownames(FUSR495X[which(FUSR495X$padj < p & FUSR495X$log2FoldChange < 0),])

TDP43M337V <- Sig$TDP43M337V
list7 <- rownames(TDP43M337V[which(TDP43M337V$padj < p & TDP43M337V$log2FoldChange > 0),])
list8 <- rownames(TDP43M337V[which(TDP43M337V$padj < p & TDP43M337V$log2FoldChange < 0),])

FUSP525Lhe <- Sig$FUSP525Lhe
list9 <- rownames(FUSP525Lhe[which(FUSP525Lhe$padj < p & FUSP525Lhe$log2FoldChange > 0),])
list10 <- rownames(FUSP525Lhe[which(FUSP525Lhe$padj < p & FUSP525Lhe$log2FoldChange < 0),])


# list.total <- list("TDP43 M337V" = list4, "FUS P525L"= list2, "FUS KO" = list1,  "FUS R495X" = list3)
list.fus.up <- list("FUS P525Lho"= list3, "FUS P525Lhe" = list9, "FUS KO" = list1,  "FUS R495X" = list5) 
list.fus.down <- list("FUS P525Lho"= list4, "FUS P525Lhe" = list10, "FUS KO" = list2,  "FUS R495X" = list6) 
list.all.up <- list("FUS P525Lho"= list3, "FUS P525Lhe" = list9, "FUS KO" = list1,  "FUS R495X" = list5, "TDP43 M337V" = list7) 
list.all.down <- list("FUS P525Lho"= list4, "FUS P525Lhe" = list10, "FUS KO" = list2,  "FUS R495X" = list6, "TDP43 M337V" = list8) 


venn.up <- Venn(list.fus.up)
venn.down <- Venn(list.fus.down)

cvenn <- compute.Venn(venn)
gp <- VennThemes(cvenn)
for(f in names(gp$Face[-1])){
  if(length(venn@IntersectionSets[[f]]) == 0){
    break
  }
  intersects <- sum(as.numeric(strsplit(f, "")[[1]]))
  if(intersects == 1){
    gp$Face[[f]]$fill <- "aliceblue"
  }
  if(intersects == 2){
    gp$Face[[f]]$fill <- "cyan2"
  }
  if(intersects == 3){
    gp$Face[[f]]$fill <- "darkseagreen2"
  }
  if(intersects == 4){
    gp$Face[[f]]$fill <- "gold2"
  }
  # if(intersects == 5){
  #   gp$Face[[f]]$fill <- "red"
  # }
}


plot(venn.down, type = "ellipses")
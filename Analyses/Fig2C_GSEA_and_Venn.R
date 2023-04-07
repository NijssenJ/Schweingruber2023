library(fgsea)
library(tidyverse)
library(Vennerable)

## Load desired DESeq2 output file (including all genes)
load("Res_MNs.Rdata")

## Load pathways from GO biological process
# downloaded from Broad Institute website: http://software.broadinstitute.org/gsea/msigdb/collections.jsp
pathways.GObp <- gmtPathways("c5.go.bp.v7.2.symbols.gmt")

# Generate empty list to store data into
Fgsea <- list()






#### FUS P525L ####

res <- Res$FUSP525Lho

ranklist <- res

# Remove genes with NA fold change
if(length(which(is.na(res$log2FoldChange))) != 0){
  ranklist <- res[-which(is.na(res$log2FoldChange)),]
}

# Order ranklist by test statistic
ranklist <- ranklist[order(ranklist$stat, decreasing=T),]

#### Make a vector of the fold-changes, and add the names of the genes ####
stats <- as.vector(ranklist$stat)
names(stats) <- rownames(ranklist)
  

#### Running the GSEA ####
fgsea <- fgsea(pathways = pathways, stats = stats, nperm = 10000, minSize = 10, maxSize = 500, nproc = 2, gseaParam = 0.25)

Fgsea$FUSP525Lho <- fgsea



#### TDP43 M337V ####

res <- Res$TDP43M337V

ranklist <- res

# Remove genes with NA fold change
if(length(which(is.na(res$log2FoldChange))) != 0){
  ranklist <- res[-which(is.na(res$log2FoldChange)),]
}

# Order ranklist by test statistic
ranklist <- ranklist[order(ranklist$stat, decreasing=T),]

#### Make a vector of the fold-changes, and add the names of the genes ####
stats <- as.vector(ranklist$stat)
names(stats) <- rownames(ranklist)


#### Running the GSEA ####
fgsea <- fgsea(pathways = pathways, stats = stats, nperm = 10000, minSize = 10, maxSize = 500, nproc = 2, gseaParam = 0.25)

Fgsea$TDP43M337V <- fgsea


#### FUS KO ####

res <- Res$FUSKO

ranklist <- res

# Remove genes with NA fold change
if(length(which(is.na(res$log2FoldChange))) != 0){
  ranklist <- res[-which(is.na(res$log2FoldChange)),]
}

# Order ranklist by test statistic
ranklist <- ranklist[order(ranklist$stat, decreasing=T),]

#### Make a vector of the fold-changes, and add the names of the genes ####
stats <- as.vector(ranklist$stat)
names(stats) <- rownames(ranklist)


#### Running the GSEA ####
fgsea <- fgsea(pathways = pathways, stats = stats, nperm = 10000, minSize = 10, maxSize = 500, nproc = 2, gseaParam = 0.25)

Fgsea$FUSKO <- fgsea




#### FUS R495X ####

res <- Res$FUSR495X

ranklist <- res

# Remove genes with NA fold change
if(length(which(is.na(res$log2FoldChange))) != 0){
  ranklist <- res[-which(is.na(res$log2FoldChange)),]
}

# Order ranklist by test statistic
ranklist <- ranklist[order(ranklist$stat, decreasing=T),]

#### Make a vector of the fold-changes, and add the names of the genes ####
stats <- as.vector(ranklist$stat)
names(stats) <- rownames(ranklist)


#### Running the GSEA ####
fgsea <- fgsea(pathways = pathways, stats = stats, nperm = 10000, minSize = 10, maxSize = 500, nproc = 2, gseaParam = 0.25)

Fgsea$FUSR495X <- fgsea




#### FUS P525Lhe ####

res <- Res$FUSP525Lhe

ranklist <- res

# Remove genes with NA fold change
if(length(which(is.na(res$log2FoldChange))) != 0){
  ranklist <- res[-which(is.na(res$log2FoldChange)),]
}

# Order ranklist by test statistic
ranklist <- ranklist[order(ranklist$stat, decreasing=T),]

#### Make a vector of the fold-changes, and add the names of the genes ####
stats <- as.vector(ranklist$stat)
names(stats) <- rownames(ranklist)

#### Running the GSEA ####
fgsea <- fgsea(pathways = pathways, stats = stats, nperm = 10000, minSize = 10, maxSize = 500, nproc = 2, gseaParam = 0.25)

Fgsea$FUSP525Lhe <- fgsea

save(Fgsea, file = "Fgsea_MN.Rdata")








###### Comparative Venn diagrams ####
p <- 0.05

FUSKO <- Fgsea$FUSKO
list1 <- FUSKO$pathway[FUSKO$padj < p & FUSKO$NES > 0]
list2 <- FUSKO$pathway[FUSKO$padj < p & FUSKO$NES < 0]

FUSP525Lho <- Fgsea$FUSP525Lho
list3 <- FUSP525Lho$pathway[FUSP525Lho$padj < p & FUSP525Lho$NES > 0]
list4 <- FUSP525Lho$pathway[FUSP525Lho$padj < p & FUSP525Lho$NES < 0]

FUSR495X <- Fgsea$FUSR495X
list5 <- FUSR495X$pathway[FUSR495X$padj < p & FUSR495X$NES > 0]
list6 <- FUSR495X$pathway[FUSR495X$padj < p & FUSR495X$NES < 0]

TDP43M337V <- Fgsea$TDP43M337V
list7 <- TDP43M337V$pathway[TDP43M337V$padj < p & TDP43M337V$NES > 0]
list8 <- TDP43M337V$pathway[TDP43M337V$padj < p & TDP43M337V$NES < 0]

FUSP525Lhe <- Fgsea$FUSP525Lhe
list9 <- FUSP525Lhe$pathway[FUSP525Lhe$padj < p & FUSP525Lhe$NES > 0]
list10 <- FUSP525Lhe$pathway[FUSP525Lhe$padj < p & FUSP525Lhe$NES < 0]

# Add desired list for Venn diagram creation like Fig1/2
list.fus.up <- list("FUS P525Lho"= list3, "FUS P525Lhe" = list9, "FUS KO" = list1,  "FUS R495X" = list5) 
list.fus.down <- list("FUS P525Lho"= list4, "FUS P525Lhe" = list10, "FUS KO" = list2,  "FUS R495X" = list6) 
list.all.up <- list("FUS P525Lho"= list3, "FUS P525Lhe" = list9, "FUS KO" = list1,  "FUS R495X" = list5, "TDP43 M337V" = list7) 
list.all.down <- list("FUS P525Lho"= list4, "FUS P525Lhe" = list10, "FUS KO" = list2,  "FUS R495X" = list6, "TDP43 M337V" = list8) 


venn <- Venn(list.fus.down)


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
  if(intersects == 5){
    gp$Face[[f]]$fill <- "red"
  }
}


plot(venn, type = "ellipses", gp = gp)


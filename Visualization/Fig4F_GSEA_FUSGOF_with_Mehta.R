library(tidyverse)
library(venn)

load("Fgsea_MN.Rdata")


###### COMPARISON ANALYSIS ####
p <- 0.05

FUSKO <- readxl::
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


# Reading in Mehta data and subsetting to overlap between 2x Mehta lines first
load("Fgsea_Mehta.Rdata")
mehta1 <- Fgsea$mehta1
mehta2 <- Fgsea$mehta2

list11 <- unlist(mehta1[which(mehta1$padj < 0.1 & mehta1$NES > 0),1])
list12 <- unlist(mehta1[which(mehta1$padj < 0.1 & mehta1$NES < 0),1])

list13 <- unlist(mehta2[which(mehta2$padj < 0.1 & mehta2$NES > 0),1])
list14 <- unlist(mehta2[which(mehta2$padj < 0.1 & mehta2$NES < 0),1])

list.mehta.up <- list(list11, list13)
list.mehta.down <- list(list12, list14)

vennmehta.up <- Venn(list.mehta.up) ; vennmehta.up <- vennmehta.up@IntersectionSets$`11`
vennmehta.down <- Venn(list.mehta.down) ; vennmehta.down <- vennmehta.down@IntersectionSets$`11`


list.all.up <- list("FUS P525Lho"= list3,  "Mehta" = vennmehta.up,  "FUS R495X" = list5, "TARDBP M337V" = list7, "FUS P525Lhe" = list9) 
list.all.down <- list("FUS P525Lho"= list4, "Mehta" = vennmehta.down,  "FUS R495X" = list6, "TARDBP M337V" = list8, "FUS P525Lhe" = list10) 

venn.up <- venn(list.all.up, ilab = T, zcolor = "royalblue3",
     opacity = 0.1, ilcs = 1.5, sncs = 1.5, borders = F, box = F)

venn.down <- venn(list.all.down, ilab = T, zcolor = "royalblue3",
                opacity = 0.1, ilcs = 1.5, sncs = 1.5, borders = F, box = F)



# Note that fig 4E is a minimized version of this Venn (4F)








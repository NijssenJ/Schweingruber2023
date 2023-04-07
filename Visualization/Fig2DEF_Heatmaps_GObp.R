library(pheatmap)
library(tidyverse)
source("Heatmap_colors.R")

### Loading GObp files for MNs
load("Fgsea_MN.Rdata")
p525lho <- Fgsea$FUSP525Lho
p525lhe <- Fgsea$FUSP525Lhe
ko <- Fgsea$FUSKO
r495x <- Fgsea$FUSR495X

MN <- left_join(p525lho, p525lhe, by = "pathway")
MN <- left_join(MN, ko, by = "pathway")
MN <- left_join(MN, r495x, by = "pathway")

MN <- as.data.frame(MN[,c(1,3,4,5,10,11,12, 17,18,19, 24,25,26)])
colnames(MN) <- c("pathway", "p525lho_padj", "p525lho_ES", "p525lho_NES", "p525lhe_padj", "p525lhe_ES", "p525lhe_NES", 
                    "ko_padj", "ko_ES", "ko_NES", "r495x_padj", "r495x_ES", "r495x_NES")

### Loading GObp files for V2as
load("Fgsea_V2a.Rdata")
p525lho <- Fgsea$FUSP525Lho
p525lhe <- Fgsea$FUSP525Lhe
ko <- Fgsea$FUSKO
r495x <- Fgsea$FUSR495X

v2a <- left_join(p525lho, p525lhe, by = "pathway")
v2a <- left_join(v2a, ko, by = "pathway")
v2a <- left_join(v2a, r495x, by = "pathway")
v2a <- as.data.frame(v2a[,c(1,3,4,5,10,11,12, 17,18,19, 24,25,26)])
colnames(v2a) <- c("pathway", "V2A_p525lho_padj", "V2A_p525lho_ES", "V2A_p525lho_NES", "V2A_p525lhe_padj", "V2A_p525lhe_ES", 
                    "V2A_p525lhe_NES", "V2A_ko_padj", "V2A_ko_ES", "V2A_ko_NES", "V2A_r495x_padj", "V2A_r495x_ES", "V2A_r495x_NES")

### Loading GObp files for Non-V2a inters
load("Fgsea_nonV2a.Rdata")
p525lho <- Fgsea$FUSP525Lho
p525lhe <- Fgsea$FUSP525Lhe
ko <- Fgsea$FUSKO
r495x <- Fgsea$FUSR495X

nonv2a <- left_join(p525lho, p525lhe, by = "pathway")
nonv2a <- left_join(nonv2a, ko, by = "pathway")
nonv2a <- left_join(nonv2a, r495x, by = "pathway")
nonv2a <- as.data.frame(nonv2a[,c(1,3,4,5,10,11,12, 17,18,19, 24,25,26)])
colnames(nonv2a) <- c("pathway", "Non-V2A_p525lho_padj", "Non-V2A_p525lho_ES", "Non-V2A_p525lho_NES", "Non-V2A_p525lhe_padj", "Non-V2A_p525lhe_ES", 
                   "Non-V2A_p525lhe_NES", "Non-V2A_ko_padj", "Non-V2A_ko_ES", "Non-V2A_ko_NES", "Non-V2A_r495x_padj", "Non-V2A_r495x_ES", "Non-V2A_r495x_NES")


combined <- left_join(MN, v2a, "pathway")
combined <- left_join(combined, nonv2a, "pathway")

# Clean up objects
rm(MN, v2a, nonv2a, p525lhe, p525lho, r495x, ko, Fgsea)
gc()



#  Subset to the desired pathways (example here is all common pathways)
for(col in c(2,5,8,11)){
  combined <- combined[which(combined[,col] < 0.05),]
}


# Subsetting p-value data
plotdata <- combined
rownames(plotdata) <- plotdata[,1]
plotdata <- plotdata[,grep("NES", colnames(plotdata))]

plotdata[is.na(plotdata)] <- 0

annotation <- data.frame(row.names = colnames(plotdata), "Line" = rep(c("FUS_P525Lho", "FUS_P525Lhe", "FUS_KO","FUS_R495X"), 3), "Celltype" = c(rep("Motor neuron", 4), rep("V2a interneuron", 4), rep("Other interneuron", 4)))
annotation <- annotation[order(annotation$Celltype, annotation$Line),]
annotation <- annotation[c(1:4, 9:12, 5:8),] # Custom ordering of V2a before other interneurons
annotation <- annotation[c(1,4,2,3,5,8,6,7,9,12,10,11),]
plotdata <- plotdata[,match(rownames(annotation), colnames(plotdata))]

# Changing GO-term names (Removing "GO_" and underscores throughout. Changing to only first letter caps)
firstup <- function(x) {
  x <- gsub("GO_", "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
rownames(plotdata) <- firstup(rownames(plotdata))

## Heatmap plotting
col_nonscaled(plotdata, cellheight = 3, aspect.ratio = 1, low = -3, high = 3, cent = 0, by = 0.05)
# col_scaled(plotdata, cellheight = 15, aspect.ratio = 1.5, cutoff = 0.8)

Line = c("darkorange2", "plum3", "skyblue", "royalblue1")
names(Line) = c("FUS_KO", "FUS_R495X", "FUS_P525Lhe", "FUS_P525Lho")
Celltype = c("tomato3", "blue", "gold2")
names(Celltype) = c("Motor neuron", "V2a interneuron", "Other interneuron")

ann_colors = list(Line = Line, Celltype = Celltype)


pheatmap(plotdata, scale="none", cellheight = cellheight, color = color, breaks = breaks,
         annotation_colors = ann_colors, cellwidth=cellwidth, treeheight_col = 15, 
         clustering_distance_cols = "euclidean", show_colnames = F, cluster_rows = T, 
         cluster_cols=F, fontsize_col=cellwidth*0.9, fontsize_row= cellheight*0.8, drop_levels = T,
         annotation_col = annotation, border_color=NA, annotation_names_col = T)


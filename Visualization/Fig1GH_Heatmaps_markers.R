library(pheatmap)
library(tidyverse)
source("Heatmap_colors.R")

### Loading RPKMs
load("rpkms_1415cells.Rdata")

### Loading annotations
load("annotation_1415cells.Rdata")
annotation <- annotation[-grep("Other", annotation$Celltype_genes),]

annotation$Line <- factor(annotation$Line, levels = c("CTRL", "FUS_KO", "FUS_R495X", "FUS_P525Lhe", "FUS_P525Lho"))
annotation <- annotation[order(annotation$Celltype_genes_v2a, annotation$Line),]


# Genes for Fig1 heatmap #1
genes <- c("ISL1", "ISL2", "MNX1", "UTS2", "CHAT", "SLC18A3", "SLC5A7", "PRPH", "MAPT", "NEFM", "SNAP25", "NGFR",
          "PAX6", "SOX2", "MKI67", "CDK1", "GAPDH", "ACTB")

# Genes for Fig1 heatmap #2
# genes <- c("ISL1", "SLC18A3", "SLC5A7", "SLC17A6", "SLC6A5",
#           "SLC32A1", "GAD2", "SOX14", "VSX2", "SOX21")

# Subsetting rpkms
filtered.rpkms <- rpkms[, match(rownames(annotation), colnames(rpkms))]
filtered.rpkms <- filtered.rpkms[match(genes, rownames(filtered.rpkms), nomatch = 0),]

# Log-transform
log.rpkms <- log2(filtered.rpkms + 1)

#######HEATMAP

col_nonscaled(log.rpkms, cellheight = 15, low = min(log.rpkms), high = max(log.rpkms), aspect.ratio = 1)
# col_scaled(log.rpkms, cellheight = 20, aspect.ratio = 1, cutoff = 0.25) # Only use if 'scale = row' in the pheatmap call

Line = c("Black", "tomato3", "royalblue1", "plum3", "darkorange2", "skyblue")
names(Line) = c("CTRL","TDP43_M337V", "FUS_P525Lho", "FUS_R495X", "FUS_KO", "FUS_P525Lhe")
Diff = c("lightblue", "orange", "green", "purple")
names(Diff) = c("D1", "D2", "D3", "D4")
Celltype_genes <- c("Grey70", "gold2", "tomato3")
names(Celltype_genes) <- c("Other", "Neuron", "Motorneuron")
Celltype_genes_v2a <- c("Grey70", "blue", "gold2", "tomato3")
names(Celltype_genes_v2a) <- c("Other", "v2a", "Neuron", "Motorneuron")
Protocol <- c("grey75", "black")
names(Protocol) <- c("Mixed", "Pure")
ramp = colorRampPalette(colors = c("green3", "lightyellow", "firebrick4"))(40)


ann_colors = list(Line = Line, Diff = Diff,  Celltype_genes = Celltype_genes,Celltype_genes_v2a = Celltype_genes_v2a,
                  hg38mapping = ramp, mapped_reads = ramp, detected_genes = ramp, Protocol = Protocol)


pheatmap(log.rpkms, scale="none", cellheight = cellheight, color=color, breaks=breaks, 
         annotation_colors = ann_colors, cellwidth=cellwidth, treeheight_col = 15, 
         clustering_distance_cols = "euclidean", show_colnames = F, cluster_rows = F, 
         cluster_cols=F, fontsize_col=cellwidth*0.9, fontsize_row= cellheight*0.8, drop_levels = T,
         annotation_col = annotation[,c(1,7),drop=F], border_color=NA, annotation_names_col = T)


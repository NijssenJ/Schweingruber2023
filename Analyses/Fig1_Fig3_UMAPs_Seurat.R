library(tidyverse)
library(readxl)
library(Seurat)


## Load in data
load("counts_1415cells.Rdata")
load("annotation_1415cells.Rdata")

# Subset metadata and counts to desired cells
annotation <- annotation[-which(annotation$Celltype_genes == "Other"),] # Other cells are excluded or included for respective analyses
counts <- counts[,match(rownames(annotation), colnames(counts))]


## Create Seurate Object
facs <- CreateSeuratObject(counts = counts, project = "FACS")

## Add meta data
all(colnames(counts) == rownames(annotation))
facs <- AddMetaData(object = facs, metadata = annotation$Line, col.name = "Cell_line")
facs <- AddMetaData(object = facs, metadata = annotation$Seqround, col.name = "NGI_seq")
facs <- AddMetaData(object = facs, metadata = annotation$Diff, col.name = "Diff")
facs <- AddMetaData(object = facs, metadata = annotation$detected_genes, col.name = "Detected_genes")
facs <- AddMetaData(object = facs, metadata = annotation$mapped_reads, col.name = "Mapped_reads")
facs <- AddMetaData(object = facs, metadata = annotation$Celltype_genes, col.name = "Celltype")


## QC plots
# facs <- SetIdent(facs, value = facs@meta.data$Cell_line)
# plot1 <- VlnPlot(facs, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
# plot2 <- FeatureScatter(facs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot3 <- FeatureScatter(facs, feature1 = "percent_mito", feature2 = "nFeature_RNA")

# plot1
# plot2
# plot3


facs1 <- subset(facs, subset = nFeature_RNA >= 5000 & nCount_RNA >= 100000)


## Batch correction by Biological Differentiation using CCA

## Split up dataset by differentiation
facs1_list <- SplitObject(facs1, split.by = "Diff")

## Normalize each subset using SCTransform (takes time)
facs1_list <- lapply(X = facs1_list, FUN = SCTransform)

## Select features that are repeatedly variable across the differentiation for integration
integration_features <- SelectIntegrationFeatures(object.list = facs1_list, nfeatures = 2000)
facs1_list <- PrepSCTIntegration(object.list = facs1_list, anchor.features = integration_features)
facs1_anchors <- FindIntegrationAnchors(object.list = facs1_list, normalization.method = "SCT", anchor.features = integration_features)

## Integrate Dataset
facs2 <- IntegrateData(anchorset = facs1_anchors, normalization.method = "SCT")

## Dim reduction after integration
facs2 <- RunPCA(facs2, npcs = 50, verbose = FALSE)
elbow <- ElbowPlot(facs2, ndims = 50)
elbow
pc_cutoff <- length(which(diff(elbow$data$stdev) < -0.05))


facs2 <- RunUMAP(facs2, reduction = "pca", dims = 1:pc_cutoff, verbose = FALSE, n.neighbors = 12)


## Export x,y coordinates of UMAP projections for plotting
coordinates <- as.data.frame(facs2@reductions$umap@cell.embeddings)
annotation <- cbind(annotation, coordinates)

save(annotation, file = "UMAP_neuronsOnly_annotation.Rdata")





##### Plot from annotation #######
load("UMAP_neuronsOnly_annotation.Rdata") # Neurons only (fig 1 and fig 3)
load("UMAP_allCells_annotation.Rdata") # All cells (fig 1 and fig 3)

annotation$UMAP_1 <- annotation$UMAP_1 * -1 # Custom flip of one axis for better fit within figure panels

# annotation <- annotation[-grep("TDP", annotation$Line),] # TDP43 M337V cells are excluded from fig1 and included in fig 3

# Plotting of UMAPs
p <- ggplot(annotation, aes(x = UMAP_1, y = UMAP_2))

# Plot per celltype
p + geom_point(aes(color = Celltype_genes), size = 1.5) +
  scale_color_manual(values = c("tomato3", "gold2", "darkblue")) +
  theme_light()


# Plot per line
p + geom_point(aes(color = annotation$Line), size = 1.5) +
  scale_color_manual(values = c("Black", "darkorange2", "skyblue", "royalblue1", "plum3", "tomato3")) +
  theme_light()

# Plot per line (TDP43 highlight)
p + geom_point(aes(color = annotation$Line), size = 1.5) +
  scale_color_manual(values = c("grey70", "grey70", "grey70", "grey70", "grey70", "tomato3")) +
  theme_light()




# Plot gene expression onto UMAP projection
library(gridExtra)
library(grid)

genes <- c("gene_name")

load("rpkms_1415cells.Rdata")
rpkms <- rpkms[,match(rownames(annotation), colnames(rpkms))]

plots <- list()

# Plot per celltype
plots[[1]] <- p + geom_point(aes(color = Celltype_genes), size = 1.5) +
  scale_color_manual(values = c("tomato3", "gold2", "darkblue")) +
  theme_light()


# Plot per line
plots[[2]] <- p + geom_point(aes(color = annotation$Line), size = 1.5) +
  scale_color_manual(values = c("Black", "darkorange2", "skyblue", "royalblue1", "plum3", "tomato3")) +
  theme_light()


counter <- 3
for(g in genelist){

  p2 <- p + geom_point(aes(color = !!log10(as.numeric(rpkms[g,])+0.1)), size = 1.5) +
    scale_color_gradientn(colors = c("grey85", "royalblue1", "firebrick2")) +
    theme_light() + ggtitle(g)
  
  plots[[counter]] <- p2
  counter <- counter + 1
  
}

plots <- lapply(plots, ggplotGrob)

ggsave(filename = "gene_name.pdf", device = "pdf",  
       width = unit(20, "cm"), height = unit(12, "cm"),
       plot = marrangeGrob(grobs = plots, layout_matrix = matrix(1:4, 2, 2, TRUE), top = NULL))
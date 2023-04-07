library(DESeq2)
library(BiocParallel)

#### Import annotation table ####
load("annotation_1415cells.Rdata")
annotation <- annotation[annotation$Line == "FUS_KO",] # Change to other cell types for respective analyses

load("Protein_coding_genes.Rdata")

#### Counts table ####
load("counts_1415cells.Rdata")
counts <- counts[,match(rownames(annotation), colnames(counts))]
counts <- counts[which(rownames(counts) %in% coding_genes),]
counts <- counts[which(rowSums(counts > 1) > 5),]


Res <- list()
Sig <- list()


##########################################

# Subset metadata and counts files
subanno <- annotation[grep("Motorneuron|Neuron", annotation$Celltype_genes_v2a),]
subcounts <- counts[,match(rownames(subanno), colnames(counts))]

#### generate DESeq file ####
deseq <- DESeqDataSetFromMatrix(subcounts, subanno, ~ Celltype_genes_v2a, tidy = F)

#### Running differential gene expression ####
deseq <- DESeq(deseq, parallel = T, BPPARAM = MulticoreParam(workers=2), fitType = "local")
res <- as.data.frame(results(deseq, contrast = c("Celltype_genes_v2a", "Motorneuron", "Neuron")))

#### Isolating significant genes ####
sig <- as.data.frame(subset(res, padj < 0.05))
sig <- sig[order(sig[,6], decreasing=F),]

Res$MNvsNeuron <- res
Sig$MNvsNeuron <- sig

##########################################

# Subset metadata and counts files
subanno <- annotation[grep("Motorneuron|v2a", annotation$Celltype_genes_v2a),]
subcounts <- counts[,match(rownames(subanno), colnames(counts))]

#### generate DESeq file ####
deseq <- DESeqDataSetFromMatrix(subcounts, subanno, ~ Celltype_genes_v2a, tidy = F)

#### Running differential gene expression ####
deseq <- DESeq(deseq, parallel = T, BPPARAM = MulticoreParam(workers=2), fitType = "local")
res <- as.data.frame(results(deseq), contrast = c("Celltype_genes_v2a", "Motorneuron", "v2a"))

#### Isolating significant genes ####
sig <- as.data.frame(subset(res, padj < 0.05))
sig <- sig[order(sig[,6], decreasing=F),]

Res$MNvsV2a <- res
Sig$MNvsV2a <- sig


##########################################

# Subset metadata and counts files
subanno <- annotation[grep("Neuron|v2a", annotation$Celltype_genes_v2a),]
subcounts <- counts[,match(rownames(subanno), colnames(counts))]

#### generate DESeq file ####
deseq <- DESeqDataSetFromMatrix(subcounts, subanno, ~ Celltype_genes_v2a, tidy = F)

#### Running differential gene expression ####
deseq <- DESeq(deseq, parallel = T, BPPARAM = MulticoreParam(workers=2), fitType = "local")
res <- as.data.frame(results(deseq), contrast = c("Celltype_genes_v2a", "v2a", "Neuron"))

#### Isolating significant genes ####
sig <- as.data.frame(subset(res, padj < 0.05))
sig <- sig[order(sig[,6], decreasing=F),]

Res$V2avsNeuron <- res
Sig$V2avsNeuron <- sig



save(Res, file = "../DEG_results/Between_celltypes_per_line/Res_celltypes_FUS_KO.Rdata")
save(Sig, file = "../DEG_results/Between_celltypes_per_line/Sig_celltypes_FUS_KO.Rdata")
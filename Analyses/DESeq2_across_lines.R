library(DESeq2)


## Import metadata and list of protein-coding genes
load("Protein_coding_genes.Rdata")
load("annotation_1415cells.Rdata")

## Subset metadata to desired cell type
annotation <- annotation[annotation$Celltype_genes == "Motorneuron",]

## Import counts table and subset
load("counts_1415cells.Rdata")
counts <- counts[,match(rownames(annotation), colnames(counts))]
counts <- counts[which(rownames(counts) %in% coding_genes),]
counts <- counts[which(rowSums(counts > 1) > 5),]


Res <- list()
Sig <- list()


##########################################

# Subset metadata and counts table
subanno <- annotation[grep("CTRL|KO", annotation$Line),]
subcounts <- counts[,match(rownames(subanno), colnames(counts))]

#### generate DESeq file ####
deseq <- DESeqDataSetFromMatrix(subcounts, subanno, ~ Diff + Line, tidy = F)

#### Running differential gene expression ####
deseq <- DESeq(deseq, parallel = T, BPPARAM = MulticoreParam(workers=2), fitType = "local")
res <- as.data.frame(results(deseq, contrast = c("Line", "FUS_KO", "CTRL")))

#### Isolating significant genes ####
sig <- as.data.frame(subset(res, padj < 0.05))
sig <- sig[order(sig[,6], decreasing=F),]

Res$FUSKO <- res
Sig$FUSKO <- sig

##########################################

# Subset metadata and counts table
subanno <- annotation[grep("CTRL|FUS_P525Lho", annotation$Line),]
subcounts <- counts[,match(rownames(subanno), colnames(counts))]

#### generate DESeq file ####
deseq <- DESeqDataSetFromMatrix(subcounts, subanno, ~ Diff + Line, tidy = F)

#### Running differential gene expression ####
deseq <- DESeq(deseq, BPPARAM = MulticoreParam(workers=2), fitType = "local")
res <- as.data.frame(results(deseq, contrast = c("Line", "FUS_P525Lho", "CTRL")))

#### Isolating significant genes ####
sig <- as.data.frame(subset(res, padj < 0.05))
sig <- sig[order(sig[,6], decreasing=F),]

Res$FUSP525Lho <- res
Sig$FUSP525Lho <- sig

##########################################

# Subset metadata and counts table
subanno <- annotation[grep("CTRL|M337V", annotation$Line),]
subcounts <- counts[,match(rownames(subanno), colnames(counts))]

#### generate DESeq file ####
deseq <- DESeqDataSetFromMatrix(subcounts, subanno, ~ Diff + Line, tidy = F)

#### Running differential gene expression ####
deseq <- DESeq(deseq, BPPARAM = MulticoreParam(workers=2), fitType = "local")
res <- as.data.frame(results(deseq, contrast = c("Line", "TDP43_M337V", "CTRL")))

#### Isolating significant genes ####
sig <- as.data.frame(subset(res, padj < 0.05))
sig <- sig[order(sig[,6], decreasing=F),]

Res$TDP43M337V <- res
Sig$TDP43M337V <- sig

##########################################

# Subset metadata and counts table
subanno <- annotation[grep("CTRL|FUS_R495X", annotation$Line),]
subcounts <- counts[,match(rownames(subanno), colnames(counts))]

#### generate DESeq file ####
deseq <- DESeqDataSetFromMatrix(subcounts, subanno, ~ Diff + Line, tidy = F)

#### Running differential gene expression ####
deseq <- DESeq(deseq, BPPARAM = MulticoreParam(workers=2), fitType = "local")
res <- as.data.frame(results(deseq, contrast = c("Line", "FUS_R495X", "CTRL")))

#### Isolating significant genes ####
sig <- as.data.frame(subset(res, padj < 0.05))
sig <- sig[order(sig[,6], decreasing=F),]

Res$FUSR495X <- res
Sig$FUSR495X <- sig


##########################################

# Subset metadata and counts table
subanno <- annotation[grep("CTRL|FUS_P525Lhe", annotation$Line),]
subcounts <- counts[,match(rownames(subanno), colnames(counts))]

#### generate DESeq file ####
deseq <- DESeqDataSetFromMatrix(subcounts, subanno, ~ Diff + Line, tidy = F)

#### Running differential gene expression ####
deseq <- DESeq(deseq, BPPARAM = MulticoreParam(workers=2), fitType = "local")
res <- as.data.frame(results(deseq, contrast = c("Line", "FUS_P525Lhe", "CTRL")))

#### Isolating significant genes ####
sig <- as.data.frame(subset(res, padj < 0.05))
sig <- sig[order(sig[,6], decreasing=F),]

Res$FUSP525Lhe <- res
Sig$FUSP525Lhe <- sig



# Export combined files including all analyzed genes (Res) or significant genes only (Sig)
save(Res, file = "Res_MNs.Rdata")
save(Sig, file = "Sig_MNs.Rdata")
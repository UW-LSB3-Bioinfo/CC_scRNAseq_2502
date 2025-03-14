## Single-cell RNA-seq data analysis

## set up the working direcoty 
setwd("~/Projects/Clemens_sc_202502/")
mainDir <- getwd()

## Step 1: Import data into R
## Use the read10×Counts function from DropletUtils package to import raw Cell Ranger output data into R
library(DropletUtils)
library(tidyverse)

## read the raw data
dir.name <- "/home/rstudio/Projects/Clemens_sc_202502/InputCounts/mud/"
list.files((dir.name))

## create an sce object (s4 class)
sce_mud <- read10xCounts(dir.name)
sce_mud
counts_mud <- counts(sce_mud)
class(counts_mud)


## Step 2: Filter empty droplets, using the R package DropletUtils
dim(sce_mud)

## use the barcodeRanks function to compute the ranks for all barcodes
bcrank_mud <- barcodeRanks(counts_mud)

## Only show unique points for plotting; identify knee and inflection
uniq_mud <- !duplicated(bcrank_mud$rank) & bcrank_mud$total > 0
bcrank_mud_uniq <- data.frame(bcrank_mud$rank[uniq_mud], bcrank_mud$total[uniq_mud])
colnames(bcrank_mud_uniq) <- c("rank", "total")
ggplot(bcrank_mud_uniq, aes(rank, total)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_hline(aes(yintercept = metadata(bcrank_mud)$inflection), col = "darkgreen", lty = 2, show.legend = TRUE) +
    geom_hline(aes(yintercept = metadata(bcrank_mud)$knee), col = "dodgerblue", lty = 2, show.legend = TRUE)  +
    geom_hline(aes(yintercept = 2000), col = "red", lty = 2)  +
    geom_vline(aes(xintercept = 10000), col = "orange", lty = 2)  +
    scale_fill_discrete(limits = c("knee", "inflection")) +
    labs(x = "log10(Rank)", y = "log10(Total UMI count)", title = "Drosophila _sens2l cells barcode ranks") +
    theme_classic()+
    theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 18),
        legend.text = element_text(size = 12)
    )

set.seed(42)
e.out_mud <- emptyDrops(counts_mud)
summary(e.out_mud$FDR <= 0.001)
table(Sig = e.out_mud$FDR <= 0.001, Limited = e.out_mud$Limited)

set.seed(42)
limit <- 100
all.out_mud <- emptyDrops(counts_mud, lower = limit, test.ambient = TRUE)
hist(all.out_mud$PValue[all.out_mud$Total <= limit & all.out_mud$Total > 0],
     xlab = "P-value", main = "", col = "grey80")

sce_mud <- sce_mud[, which(e.out_mud$FDR <= 0.001)]

## dimension of sce object after filtering empty droplets
dim(sce_mud)

## Step 3: Remove low quality cells

## use the scater package for removing low quality cells
library(ggplot2)
library(scater)
library(scuttle)

## definition of mitochondrial and ribosomal groups of genes as features of
is.mito <- grep("^mt:", rowData(sce_mud)$Symbol)
is.ribo <- grep("^Rp[SL]", rowData(sce_mud)$Symbol)
is.heatshock <- grep("^Hsp", rowData(sce_mud)$Symbol)
length(is.mito)
length(is.ribo)
length(is.heatshock)

## addPerCellQCmetrics to evaluate % of mitochondrial or ribosomal genes for each cell
sce_mud_QC <- addPerCellQC(sce_mud, subsets = list(Mt = is.mito, Hs = is.heatshock, Ribo = is.ribo))
colnames(colData(sce_mud_QC))

## plot the distribution of number of total counts
qplot(sce_mud_QC$sum, geom = "histogram", binwidth = 500,
      xlab = "number of counts", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22)

## plot the distribution of the log10-transformed number of total counts
qplot(log10(sce_mud_QC$sum), geom = "histogram", binwidth = 0.05,
      xlab = "log10 (number of counts)", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    geom_vline(xintercept = log10(1000), col = "red", size = 1) +
    theme_light(base_size = 22)

## plot the distribution of the total genes
qplot(sce_mud_QC$detected, geom = "histogram", binwidth = 75,
      xlab = "number of genes", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22)

## plot the distribution of the log10-transformed number of genes
qplot(log10(sce_mud_QC$detected), geom = "histogram", binwidth = 0.02,
      xlab = "log10 (number of genes)", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22) +
    geom_vline(xintercept = log10(500), color = "red", size = 1)


## plot the distribution of the log10-transformed total number of genes per UMI
sce_mud_QC$log10GenesPerUMI <- log10(sce_mud_QC$detected)/log10(sce_mud_QC$total)
qplot(sce_mud_QC$log10GenesPerUMI, geom = "density",
      xlab = "log10(number of genes ) / log10(number of counts)", ylab = "density",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22) +
    geom_vline(xintercept = 0.8, color = "red", size = 1)


## remove barcodes/cells with total log10(Genes/counts) < 0.8
log10GenesPerUMI.keep <- sce_mud_QC$log10GenesPerUMI >= 0.8
log10GenesPerUMI.remove <- sce_mud_QC$log10GenesPerUMI < 0.8
sce_mud_QC_filter <- sce_mud_QC[, log10GenesPerUMI.keep]
dim(sce_mud_QC_filter)


## remove barcodes/cells with total number of genes less than 500, and total UMIs greater than 1000
genes.keep <- sce_mud_QC$detected >= 500 & sce_mud_QC$total >= 1000 
sce_mud_QC_filter <- sce_mud_QC[, genes.keep]
dim(sce_mud_QC_filter)



## plot the distribution of the mitochondrial proportion
qplot(sce_mud_QC_filter$subsets_Mt_percent, geom = "histogram", binwidth = 0.5,
      xlab = "Mitochondrial prop. %", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22)

## plot the distribution of the ribosomal proportion
qplot(sce_mud_QC_filter$subsets_Ribo_percent, geom = "histogram", binwidth = 0.5,
      xlab = "Ribosome prop. %", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22) 

## plot the ribosomal proportion vs mitochodrial proportion
smoothScatter(sce_mud_QC_filter$subsets_Ribo_percent, sce_mud_QC_filter$subsets_Mt_percent,
              xlab = "Ribosome prop. (%)", ylab = "Mitochondrial prop. (%)",
              nrpoints = 1000, cex.axis = 1.5, cex.lab = 1.8)
abline(h = 18, lty = 1, col = "red")
abline(v = c(5, 40), lty = 1, col = "red")

## remove barcodes/cells with a mitochondrial percentage greater than 18%
mito.keep <- sce_mud_QC_filter$subsets_Mt_percent <= 18
mito.remove <- sce_mud_QC_filter$subsets_Mt_percent > 18
sce_mud_QC_filter_2 <- sce_mud_QC_filter[, mito.keep]
dim(sce_mud_QC_filter_2)

## remove barcodes/cells with a ribosomal percentage less than 5%, or greater than 40%
ribo.keep <- (sce_mud_QC_filter_2$subsets_Ribo_percent >= 5) & (sce_mud_QC_filter_2$subsets_Ribo_percent < 40)
sce_mud_QC_filter_2 <- sce_mud_QC_filter_2[, ribo.keep]
dim(sce_mud_QC_filter_2)

## plot the distribution of the heatshock proteins proportion
qplot(sce_mud_QC_filter_2$subsets_Hs_percent, geom = "histogram", binwidth = 0.3,
      xlab = "Heatshock proteins prop. %", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22) +
    geom_vline(xintercept = 5, color = "red", size = 1)


## remove barcodes/cells with a heatshock protein proportion greater than 5%
hs.keep <- sce_mud_QC_filter_2$subsets_Hs_percent <= 5
sce_mud_QC_filter_2 <-sce_mud_QC_filter_2[, hs.keep]
dim(sce_mud_QC_filter_2)


## Step 4: Identify highly expressed genes

## summarize gene_level information
summary(rowSums(assay(sce_mud_QC_filter_2)))
fea_greaterthan0 <- rowSums(assay(sce_mud_QC_filter_2)) > 0
summary(fea_greaterthan0)
summary(rowSums(assay(sce_mud_QC_filter_2)[fea_greaterthan0, ]))
fea_greaterthan2 <- rowSums(assay(sce_mud_QC_filter_2)) > 2
summary(fea_greaterthan2)
summary(rowSums(assay(sce_mud_QC_filter_2)[fea_greaterthan2, ]))
fea_greaterthan5 <- rowSums(assay(sce_mud_QC_filter_2)) > 5
summary(fea_greaterthan5)
summary(rowSums(assay(sce_mud_QC_filter_2)[fea_greaterthan5, ]))


## only keep genes have counts more than 2
sce_mud_QC_filter_2_gene <- sce_mud_QC_filter_2[fea_greaterthan2, ]
dim(sce_mud_QC_filter_2_gene)

## only keep genes express in at least 2 cells
keep_feature <- nexprs(sce_mud_QC_filter_2_gene, byrow = TRUE) >= 2
summary(keep_feature)
sce_mud_QC_filter_2_gene <- sce_mud_QC_filter_2_gene[keep_feature, ]
dim(sce_mud_QC_filter_2_gene)

## identify highly expressed genes
plotHighestExprs(sce_mud_QC_filter_2_gene, exprs_values = "counts", 
                 feature_names_to_plot = "Symbol", colour_cells_by = "detected")

summary(colData(sce_mud_QC_filter_2_gene)$detected)
summary(colData(sce_mud_QC_filter_2_gene)$total)


## Step 4: Normalization

## normalization by deconvolution, using the scran package
set.seed(42)
library(scran)
clust.dro__sens2 <- quickCluster(sce_mud_QC_filter_2_gene)
sce_mud_QC_filter_2_gene <- computeSumFactors(sce_mud_QC_filter_2_gene, 
                                                  cluster = clust.dro__sens2, 
                                                  min.mean = 0.1)
summary(sizeFactors(sce_mud_QC_filter_2_gene))

par(mar = c(6,6,6,6))
plot(sce_mud_QC_filter_2_gene$sum, 
     sizeFactors(sce_mud_QC_filter_2_gene),
     log = "xy", 
     xlab = "total counts",
     ylab = "size factors",
     cex.axis = 1.2,
     cex.lab = 1.2,
     cex =0.8,
     pch = 20,
     col = rgb(0.1, 0.2, 0.7, 0.3))

sce_mud_QC_filter_2_gene <- logNormCounts(sce_mud_QC_filter_2_gene)
assayNames(sce_mud_QC_filter_2_gene)

## assignment gene symbol from rowname
rownames(sce_mud_QC_filter_2_gene) <- rowData(sce_mud_QC_filter_2_gene)$Symbol
head(rownames(sce_mud_QC_filter_2_gene))

## save file
save.image("/home/rstudio/Projects/sens2_sc/Outputs/sce_mud_normalized.RData")


## Step 5: Conversion from sce objects to Seurat objects
library(dplyr)
library(Seurat)
seurat_mud <- as.Seurat(sce_mud_QC_filter_2_gene, counts = "counts", data = "logcounts")


## Step 6: Identify the most variable genes (feature selection)

## use the FindVariableFeatures function from the Seurat package to identify highly variable features
seurat_mud <- FindVariableFeatures(seurat_mud, selection.method = "vst", nfeatures = 2000)
top15 <- head(VariableFeatures(seurat_mud), 15)
top30 <- head(VariableFeatures(seurat_mud), 30)
plot1 <- VariableFeaturePlot(seurat_mud)
LabelPoints(plot = plot1, points = top30, repel = TRUE, xnudge = 0, ynudge = 0)


## Step 7: Scaling the data
all.genes <- rownames(seurat_mud)
seurat_mud <- ScaleData(seurat_mud, features = all.genes)


## Step 8: Principal component analysis (PCA)

## get the first 40 PCs
seurat_mud <- RunPCA(seurat_mud, npcs = 40, verbose = FALSE)
ElbowPlot(seurat_mud, ndims = 40)

VizDimLoadings(seurat_mud, dims = 1:4, reduction = "pca")
VizDimLoadings(seurat_mud, dims = 5:8, reduction = "pca")
VizDimLoadings(seurat_mud, dims = 9:12, reduction = "pca")
VizDimLoadings(seurat_mud, dims = 13:16, reduction = "pca")
VizDimLoadings(seurat_mud, dims = 17:20, reduction = "pca")
VizDimLoadings(seurat_mud, dims = 21:24, reduction = "pca")
VizDimLoadings(seurat_mud, dims = 25:28, reduction = "pca")

DimPlot(seurat_mud, reduction = "pca")

DimHeatmap(seurat_mud, dims = 1:9, cells = 300, balanced = TRUE)
DimHeatmap(seurat_mud, dims = 10:18, cells = 300, balanced = TRUE)
DimHeatmap(seurat_mud, dims = 19:27, cells = 300, balanced = TRUE)

## determine the 'dimensionality' of the dataset
seurat_mud <- JackStraw(seurat_mud, num.replicate = 100)
seurat_mud <- ScoreJackStraw(seurat_mud, dims = 1:20)
JackStrawPlot(seurat_mud, dims = 1:20)

## Step 9: Non-linear demensional reduction (UMAP/t-SNE)
## use the first 30 PCs
seurat_mud <- RunTSNE(seurat_mud, reduction = "pca", dims = 1:30)
DimPlot(seurat_mud, reduction = "tsne", label = FALSE)
seurat_mud <-RunUMAP(seurat_mud, reduction = "pca", dims = 1:30)
DimPlot(seurat_mud, reduction = "umap", label = FALSE)


## Step 10: cluster the cells, set up the resolution as 0.3 (it might need to be modified in the future)
seurat_mud <- FindNeighbors(seurat_mud, reduction = "pca", dims = 1:30)
seurat_mud <- FindClusters(seurat_mud, resolution = 0.3)
DimPlot(seurat_mud, reduction = "umap", label = TRUE)
DimPlot(seurat_mud, reduction = "tsne", label = TRUE)


## check the total counts and detected genes number for each cluster
counts_genes_cluster_sens2 <-  seurat_mud@meta.data %>%
    group_by(seurat_clusters) %>%
    dplyr::summarize(mean_total = mean(total),
                     mean_gene = mean(detected),
                     mean_size = mean(sizeFactor))

counts_genes_cluster_sens2

FeaturePlot(seurat_mud, features = "nrv3", label = TRUE)
FeaturePlot(seurat_mud, features = "brat", label = TRUE)
FeaturePlot(seurat_mud, features = "Dl", label = TRUE)
FeaturePlot(seurat_mud, features = "nrv2", label = TRUE)
FeaturePlot(seurat_mud, features = "pnt", label = TRUE)
FeaturePlot(seurat_mud, features = "jdp", label = TRUE)
FeaturePlot(seurat_mud, features = "repo", label = TRUE)
FeaturePlot(seurat_mud, features = "esg", label = TRUE)
FeaturePlot(seurat_mud, features = "stg", label = TRUE)

VlnPlot(seurat_mud, features = "detected")

table(seurat_mud@meta.data$seurat_clusters)

## plot clusters as a tree
library(ape)
seurat_mud <- BuildClusterTree(object = seurat_mud)
PlotClusterTree(object = seurat_mud)

## Step 11: predict doublets
suppressMessages(require(DoubletFinder))
seurat_mud <- RenameAssays(seurat_mud, originalexp = "RNA")
seurat_mud <- NormalizeData(seurat_mud)
seurat_mud <- FindVariableFeatures(seurat_mud, selection.method = "vst", nfeatures = 2000)
seurat_mud <- ScaleData(seurat_mud, features = all.genes)
seurat_mud <- RunPCA(seurat_mud, npcs = 30, verbose = FALSE)

## run parameter optimization with paramSweep
sweep.res <- paramSweep(seurat_mud, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    geom_point(col = "red") +
    geom_line(linetype = "dashed", col = "blue") +
    labs(title = "The BCmvn distributions")

## define the expected number of doublet cells, expect 5% doublets
nExp <- round(ncol(seurat_mud) * 0.05)
seurat_mud <- doubletFinder(seurat_mud, pN = 0.25, pK = 0.22, nExp = nExp, PCs = 1:30)

## extract correct column name of doublet finder prediction
DF.name <- colnames(seurat_mud@meta.data)[grepl("DF.classification", colnames(seurat_mud@meta.data))]
cowplot::plot_grid(
    ncol = 2,
    DimPlot(seurat_mud, reduction = "umap", label = TRUE) + NoAxes(),
    DimPlot(seurat_mud, group.by = DF.name) + NoAxes()
)

## check whether doublets have more detected genes
VlnPlot(seurat_mud, features = "total", group.by = DF.name, pt.size = 0.1)
summary(seurat_mud@meta.data$sizeFactor)

## remove all predicted doublets from our data
seurat_mud_singlet <- seurat_mud[, seurat_mud@meta.data[, DF.name] == "Singlet"]
dim(seurat_mud_singlet)

## 0205-2025 stop here
## Step 12: Repeat previous steps, re-perform clustering

## Step 6_2: Identify the most variable genes (feature selection)

## use the FindVariableFeatures function from the Seurat package to identify highly variable features
seurat_mud_singlet <- FindVariableFeatures(seurat_mud_singlet, selection.method = "vst", nfeatures = 2000)
top15 <- head(VariableFeatures(seurat_mud_singlet), 15)
top30 <- head(VariableFeatures(seurat_mud_singlet), 30)
plot1 <- VariableFeaturePlot(seurat_mud_singlet)
LabelPoints(plot = plot1, points = top30, repel = TRUE, xnudge = 0, ynudge = 0)


## Step 7_2: Scaling the data
all.genes <- rownames(seurat_mud_singlet)
seurat_mud_singlet <- ScaleData(seurat_mud_singlet, features = all.genes)


## Step 8_2: Principal component analysis (PCA)

## get the first 40 PCs
seurat_mud_singlet <- RunPCA(seurat_mud_singlet, npcs = 40, verbose = FALSE)
ElbowPlot(seurat_mud_singlet, ndims = 40)

VizDimLoadings(seurat_mud_singlet, dims = 1:4, reduction = "pca")
VizDimLoadings(seurat_mud_singlet, dims = 5:8, reduction = "pca")
VizDimLoadings(seurat_mud_singlet, dims = 9:12, reduction = "pca")
VizDimLoadings(seurat_mud_singlet, dims = 13:16, reduction = "pca")
VizDimLoadings(seurat_mud_singlet, dims = 17:20, reduction = "pca")
VizDimLoadings(seurat_mud_singlet, dims = 21:24, reduction = "pca")
VizDimLoadings(seurat_mud_singlet, dims = 25:28, reduction = "pca")

DimPlot(seurat_mud_singlet, reduction = "pca")

DimHeatmap(seurat_mud_singlet, dims = 1:9, cells = 300, balanced = TRUE)
DimHeatmap(seurat_mud_singlet, dims = 10:18, cells = 300, balanced = TRUE)


## determine the 'dimensionality' of the dataset
seurat_mud_singlet <- JackStraw(seurat_mud_singlet, num.replicate = 100)
seurat_mud_singlet <- ScoreJackStraw(seurat_mud_singlet, dims = 1:20)
JackStrawPlot(seurat_mud_singlet, dims = 1:20)



## Step 9_2: Non-linear demensional reduction (UMAP/t-SNE)
## use the first 30 PCs
seurat_mud_singlet <- RunTSNE(seurat_mud_singlet, reduction = "pca", dims = 1:30)
DimPlot(seurat_mud_singlet, reduction = "tsne", label = FALSE)
seurat_mud_singlet <-RunUMAP(seurat_mud_singlet, reduction = "pca", dims = 1:30)
DimPlot(seurat_mud_singlet, reduction = "umap", label = FALSE)


## Step 10_2: cluster the cells, set up the resolution as 0.3 (it might need to be modified in the future)
seurat_mud_singlet <- FindNeighbors(seurat_mud_singlet, reduction = "pca", dims = 1:30)
seurat_mud_singlet <- FindClusters(seurat_mud_singlet, resolution = 1.5)
DimPlot(seurat_mud_singlet, reduction = "umap", label = TRUE)
DimPlot(seurat_mud_singlet, reduction = "tsne", label = TRUE)

## Take a look at cluster sizes
table(seurat_mud_singlet@meta.data$seurat_clusters)

## check the total counts and detected genes number for each cluster
counts_genes_cluster_sens2_singlet <-  seurat_mud_singlet@meta.data %>%
    group_by(seurat_clusters) %>%
    dplyr::summarise(mean_total = mean(total),
                     mean_gene = mean(detected),
                     mean_size = mean(sizeFactor))

counts_genes_cluster_sens2_singlet

FeaturePlot(seurat_mud_singlet, features = "nrv3", label = TRUE)
FeaturePlot(seurat_mud_singlet, features = "nSyb", label = TRUE)
FeaturePlot(seurat_mud_singlet, features = "elav", label = TRUE)
FeaturePlot(seurat_mud_singlet, features = "brp", label = TRUE)

FeaturePlot(seurat_mud_singlet, features = "GFP", label = TRUE)
FeaturePlot(seurat_mud_singlet, features = "mCherry", label = TRUE)

FeaturePlot(seurat_mud_singlet, features = "sizeFactor", label = TRUE)

## plot clusters as a tree
library(ape)
seurat_mud_singlet <- BuildClusterTree(object = seurat_mud_singlet)
PlotClusterTree(object = seurat_mud_singlet)


## find markers for every cluster compared to all remaining cells, report only the positive ones
mud.singlet.markers <- FindAllMarkers(seurat_mud_singlet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
nrow(mud.singlet.markers)

## save the Drosophila peripheral neurons singlet seurat object
saveRDS(seurat_mud_singlet, file = "/home/rstudio/Projects/Clemens_sc_202502/Outputs/seurat_mud_singlet.rds")

## save the markers for all clusters
write.table(mud.singlet.markers, sep="\t", file.path(mainDir, paste0("Outputs/mud_singlet_markers_res1.5",".txt")), quote=FALSE, col.names=TRUE, row.names=FALSE)



## generate expression heatmap for all clusters
top5_markers_mud <- mud.singlet.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
View(top5_markers_mud)
DoHeatmap(seurat_mud_singlet, features = top5_markers_mud$gene) + NoLegend()

## save file
save.image("/home/rstudio/Projects/Clemens_sc_202502/Outputs/sce_mud_singletClusters.RData")

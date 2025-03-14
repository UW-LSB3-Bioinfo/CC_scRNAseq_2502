## Single-cell RNA-seq data analysis

## set up the working direcoty 
setwd("~/Projects/Clemens_sc_202502/")
mainDir <- getwd()

## Step 1: Import data into R
## Use the read10×Counts function from DropletUtils package to import raw Cell Ranger output data into R
library(DropletUtils)
library(tidyverse)

## read the raw data
dir.name <- "/home/rstudio/Projects/Clemens_sc_202502/InputCounts/nanobody/"
list.files((dir.name))

## create an sce object (s4 class)
sce_nanobody <- read10xCounts(dir.name)
sce_nanobody
counts_nanobody <- counts(sce_nanobody)
class(counts_nanobody)


## Step 2: Filter empty droplets, using the R package DropletUtils
dim(sce_nanobody)

## use the barcodeRanks function to compute the ranks for all barcodes
bcrank_nanobody <- barcodeRanks(counts_nanobody)

## Only show unique points for plotting; identify knee and inflection
uniq_nanobody <- !duplicated(bcrank_nanobody$rank) & bcrank_nanobody$total > 0
bcrank_nanobody_uniq <- data.frame(bcrank_nanobody$rank[uniq_nanobody], bcrank_nanobody$total[uniq_nanobody])
colnames(bcrank_nanobody_uniq) <- c("rank", "total")
ggplot(bcrank_nanobody_uniq, aes(rank, total)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_hline(aes(yintercept = metadata(bcrank_nanobody)$inflection), col = "darkgreen", lty = 2, show.legend = TRUE) +
    geom_hline(aes(yintercept = metadata(bcrank_nanobody)$knee), col = "dodgerblue", lty = 2, show.legend = TRUE)  +
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
e.out_nanobody <- emptyDrops(counts_nanobody)
summary(e.out_nanobody$FDR <= 0.001)
table(Sig = e.out_nanobody$FDR <= 0.001, Limited = e.out_nanobody$Limited)

set.seed(42)
limit <- 100
all.out_nanobody <- emptyDrops(counts_nanobody, lower = limit, test.ambient = TRUE)
hist(all.out_nanobody$PValue[all.out_nanobody$Total <= limit & all.out_nanobody$Total > 0],
     xlab = "P-value", main = "", col = "grey80")

sce_nanobody <- sce_nanobody[, which(e.out_nanobody$FDR <= 0.001)]

## dimension of sce object after filtering empty droplets
dim(sce_nanobody)

## Step 3: Remove low quality cells

## use the scater package for removing low quality cells
library(ggplot2)
library(scater)
library(scuttle)

## definition of mitochondrial and ribosomal groups of genes as features of
is.mito <- grep("^mt:", rowData(sce_nanobody)$Symbol)
is.ribo <- grep("^Rp[SL]", rowData(sce_nanobody)$Symbol)
is.heatshock <- grep("^Hsp", rowData(sce_nanobody)$Symbol)
length(is.mito)
length(is.ribo)
length(is.heatshock)

## addPerCellQCmetrics to evaluate % of mitochondrial or ribosomal genes for each cell
sce_nanobody_QC <- addPerCellQC(sce_nanobody, subsets = list(Mt = is.mito, Hs = is.heatshock, Ribo = is.ribo))
colnames(colData(sce_nanobody_QC))

## plot the distribution of number of total counts
qplot(sce_nanobody_QC$sum, geom = "histogram", binwidth = 500,
      xlab = "number of counts", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22)

## plot the distribution of the log10-transformed number of total counts
qplot(log10(sce_nanobody_QC$sum), geom = "histogram", binwidth = 0.05,
      xlab = "log10 (number of counts)", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    geom_vline(xintercept = log10(1000), col = "red", size = 1) +
    theme_light(base_size = 22)

## plot the distribution of the total genes
qplot(sce_nanobody_QC$detected, geom = "histogram", binwidth = 75,
      xlab = "number of genes", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22)

## plot the distribution of the log10-transformed number of genes
qplot(log10(sce_nanobody_QC$detected), geom = "histogram", binwidth = 0.02,
      xlab = "log10 (number of genes)", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22) +
    geom_vline(xintercept = log10(500), color = "red", size = 1)


## plot the distribution of the log10-transformed total number of genes per UMI
sce_nanobody_QC$log10GenesPerUMI <- log10(sce_nanobody_QC$detected)/log10(sce_nanobody_QC$total)
qplot(sce_nanobody_QC$log10GenesPerUMI, geom = "density",
      xlab = "log10(number of genes ) / log10(number of counts)", ylab = "density",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22) +
    geom_vline(xintercept = 0.8, color = "red", size = 1)


## remove barcodes/cells with total log10(Genes/counts) < 0.8
log10GenesPerUMI.keep <- sce_nanobody_QC$log10GenesPerUMI >= 0.8
log10GenesPerUMI.remove <- sce_nanobody_QC$log10GenesPerUMI < 0.8
sce_nanobody_QC_filter <- sce_nanobody_QC[, log10GenesPerUMI.keep]
dim(sce_nanobody_QC_filter)


## remove barcodes/cells with total number of genes less than 500, and total UMIs greater than 1000
genes.keep <- sce_nanobody_QC$detected >= 500 & sce_nanobody_QC$total >= 1000 
sce_nanobody_QC_filter <- sce_nanobody_QC[, genes.keep]
dim(sce_nanobody_QC_filter)



## plot the distribution of the mitochondrial proportion
qplot(sce_nanobody_QC_filter$subsets_Mt_percent, geom = "histogram", binwidth = 0.5,
      xlab = "Mitochondrial prop. %", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22)

## plot the distribution of the ribosomal proportion
qplot(sce_nanobody_QC_filter$subsets_Ribo_percent, geom = "histogram", binwidth = 0.5,
      xlab = "Ribosome prop. %", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22) 

## plot the ribosomal proportion vs mitochodrial proportion
smoothScatter(sce_nanobody_QC_filter$subsets_Ribo_percent, sce_nanobody_QC_filter$subsets_Mt_percent,
              xlab = "Ribosome prop. (%)", ylab = "Mitochondrial prop. (%)",
              nrpoints = 1000, cex.axis = 1.5, cex.lab = 1.8)
abline(h = 18, lty = 1, col = "red")
abline(v = c(5, 40), lty = 1, col = "red")

## remove barcodes/cells with a mitochondrial percentage greater than 18%
mito.keep <- sce_nanobody_QC_filter$subsets_Mt_percent <= 18
mito.remove <- sce_nanobody_QC_filter$subsets_Mt_percent > 18
sce_nanobody_QC_filter_2 <- sce_nanobody_QC_filter[, mito.keep]
dim(sce_nanobody_QC_filter_2)

## remove barcodes/cells with a ribosomal percentage less than 5%, or greater than 40%
ribo.keep <- (sce_nanobody_QC_filter_2$subsets_Ribo_percent >= 5) & (sce_nanobody_QC_filter_2$subsets_Ribo_percent < 40)
sce_nanobody_QC_filter_2 <- sce_nanobody_QC_filter_2[, ribo.keep]
dim(sce_nanobody_QC_filter_2)

## plot the distribution of the heatshock proteins proportion
qplot(sce_nanobody_QC_filter_2$subsets_Hs_percent, geom = "histogram", binwidth = 0.3,
      xlab = "Heatshock proteins prop. %", ylab = "number of cells",
      fill = I("#00ABFD"), col = I("black")) +
    theme_light(base_size = 22) +
    geom_vline(xintercept = 5, color = "red", size = 1)


## remove barcodes/cells with a heatshock protein proportion greater than 5%
hs.keep <- sce_nanobody_QC_filter_2$subsets_Hs_percent <= 5
sce_nanobody_QC_filter_2 <-sce_nanobody_QC_filter_2[, hs.keep]
dim(sce_nanobody_QC_filter_2)


## Step 4: Identify highly expressed genes

## summarize gene_level information
summary(rowSums(assay(sce_nanobody_QC_filter_2)))
fea_greaterthan0 <- rowSums(assay(sce_nanobody_QC_filter_2)) > 0
summary(fea_greaterthan0)
summary(rowSums(assay(sce_nanobody_QC_filter_2)[fea_greaterthan0, ]))
fea_greaterthan2 <- rowSums(assay(sce_nanobody_QC_filter_2)) > 2
summary(fea_greaterthan2)
summary(rowSums(assay(sce_nanobody_QC_filter_2)[fea_greaterthan2, ]))
fea_greaterthan5 <- rowSums(assay(sce_nanobody_QC_filter_2)) > 5
summary(fea_greaterthan5)
summary(rowSums(assay(sce_nanobody_QC_filter_2)[fea_greaterthan5, ]))


## only keep genes have counts more than 2
sce_nanobody_QC_filter_2_gene <- sce_nanobody_QC_filter_2[fea_greaterthan2, ]
dim(sce_nanobody_QC_filter_2_gene)

## only keep genes express in at least 2 cells
keep_feature <- nexprs(sce_nanobody_QC_filter_2_gene, byrow = TRUE) >= 2
summary(keep_feature)
sce_nanobody_QC_filter_2_gene <- sce_nanobody_QC_filter_2_gene[keep_feature, ]
dim(sce_nanobody_QC_filter_2_gene)

## identify highly expressed genes
plotHighestExprs(sce_nanobody_QC_filter_2_gene, exprs_values = "counts", 
                 feature_names_to_plot = "Symbol", colour_cells_by = "detected")

summary(colData(sce_nanobody_QC_filter_2_gene)$detected)
summary(colData(sce_nanobody_QC_filter_2_gene)$total)


## Step 4: Normalization

## normalization by deconvolution, using the scran package
set.seed(42)
library(scran)
clust.dro__sens2 <- quickCluster(sce_nanobody_QC_filter_2_gene)
sce_nanobody_QC_filter_2_gene <- computeSumFactors(sce_nanobody_QC_filter_2_gene, 
                                              cluster = clust.dro__sens2, 
                                              min.mean = 0.1)
summary(sizeFactors(sce_nanobody_QC_filter_2_gene))

par(mar = c(6,6,6,6))
plot(sce_nanobody_QC_filter_2_gene$sum, 
     sizeFactors(sce_nanobody_QC_filter_2_gene),
     log = "xy", 
     xlab = "total counts",
     ylab = "size factors",
     cex.axis = 1.2,
     cex.lab = 1.2,
     cex =0.8,
     pch = 20,
     col = rgb(0.1, 0.2, 0.7, 0.3))

sce_nanobody_QC_filter_2_gene <- logNormCounts(sce_nanobody_QC_filter_2_gene)
assayNames(sce_nanobody_QC_filter_2_gene)

## assignment gene symbol from rowname
rownames(sce_nanobody_QC_filter_2_gene) <- rowData(sce_nanobody_QC_filter_2_gene)$Symbol
head(rownames(sce_nanobody_QC_filter_2_gene))

## save file
save.image("/home/rstudio/Projects/sens2_sc/Outputs/sce_nanobody_normalized.RData")


## Step 5: Conversion from sce objects to Seurat objects
library(dplyr)
library(Seurat)
seurat_nanobody <- as.Seurat(sce_nanobody_QC_filter_2_gene, counts = "counts", data = "logcounts")


## Step 6: Identify the most variable genes (feature selection)

## use the FindVariableFeatures function from the Seurat package to identify highly variable features
seurat_nanobody <- FindVariableFeatures(seurat_nanobody, selection.method = "vst", nfeatures = 2000)
top15 <- head(VariableFeatures(seurat_nanobody), 15)
top30 <- head(VariableFeatures(seurat_nanobody), 30)
plot1 <- VariableFeaturePlot(seurat_nanobody)
LabelPoints(plot = plot1, points = top30, repel = TRUE, xnudge = 0, ynudge = 0)


## Step 7: Scaling the data
all.genes <- rownames(seurat_nanobody)
seurat_nanobody <- ScaleData(seurat_nanobody, features = all.genes)


## Step 8: Principal component analysis (PCA)

## get the first 40 PCs
seurat_nanobody <- RunPCA(seurat_nanobody, npcs = 40, verbose = FALSE)
ElbowPlot(seurat_nanobody, ndims = 40)

VizDimLoadings(seurat_nanobody, dims = 1:4, reduction = "pca")
VizDimLoadings(seurat_nanobody, dims = 5:8, reduction = "pca")
VizDimLoadings(seurat_nanobody, dims = 9:12, reduction = "pca")
VizDimLoadings(seurat_nanobody, dims = 13:16, reduction = "pca")
VizDimLoadings(seurat_nanobody, dims = 17:20, reduction = "pca")
VizDimLoadings(seurat_nanobody, dims = 21:24, reduction = "pca")
VizDimLoadings(seurat_nanobody, dims = 25:28, reduction = "pca")

DimPlot(seurat_nanobody, reduction = "pca")

DimHeatmap(seurat_nanobody, dims = 1:9, cells = 300, balanced = TRUE)
DimHeatmap(seurat_nanobody, dims = 10:18, cells = 300, balanced = TRUE)
DimHeatmap(seurat_nanobody, dims = 19:27, cells = 300, balanced = TRUE)

## determine the 'dimensionality' of the dataset
seurat_nanobody <- JackStraw(seurat_nanobody, num.replicate = 100)
seurat_nanobody <- ScoreJackStraw(seurat_nanobody, dims = 1:20)
JackStrawPlot(seurat_nanobody, dims = 1:20)



## Step 9: Non-linear demensional reduction (UMAP/t-SNE)
## use the first 30 PCs
seurat_nanobody <- RunTSNE(seurat_nanobody, reduction = "pca", dims = 1:30)
DimPlot(seurat_nanobody, reduction = "tsne", label = FALSE)
seurat_nanobody <-RunUMAP(seurat_nanobody, reduction = "pca", dims = 1:30)
DimPlot(seurat_nanobody, reduction = "umap", label = FALSE)


## Step 10: cluster the cells, set up the resolution as 0.3 (it might need to be modified in the future)
seurat_nanobody <- FindNeighbors(seurat_nanobody, reduction = "pca", dims = 1:30)
seurat_nanobody <- FindClusters(seurat_nanobody, resolution = 0.3)
DimPlot(seurat_nanobody, reduction = "umap", label = TRUE)
DimPlot(seurat_nanobody, reduction = "tsne", label = TRUE)


## check the total counts and detected genes number for each cluster
counts_genes_cluster_sens2 <-  seurat_nanobody@meta.data %>%
    group_by(seurat_clusters) %>%
    dplyr::summarize(mean_total = mean(total),
                     mean_gene = mean(detected),
                     mean_size = mean(sizeFactor))

counts_genes_cluster_sens2

FeaturePlot(seurat_nanobody, features = "nrv3", label = TRUE)
FeaturePlot(seurat_nanobody, features = "brat", label = TRUE)
FeaturePlot(seurat_nanobody, features = "Dl", label = TRUE)
FeaturePlot(seurat_nanobody, features = "nrv2", label = TRUE)
FeaturePlot(seurat_nanobody, features = "pnt", label = TRUE)
FeaturePlot(seurat_nanobody, features = "jdp", label = TRUE)
FeaturePlot(seurat_nanobody, features = "repo", label = TRUE)
FeaturePlot(seurat_nanobody, features = "esg", label = TRUE)
FeaturePlot(seurat_nanobody, features = "stg", label = TRUE)

VlnPlot(seurat_nanobody, features = "detected")

table(seurat_nanobody@meta.data$seurat_clusters)

## plot clusters as a tree
library(ape)
seurat_nanobody <- BuildClusterTree(object = seurat_nanobody)
PlotClusterTree(object = seurat_nanobody)

## Step 11: predict doublets
suppressMessages(require(DoubletFinder))
seurat_nanobody <- RenameAssays(seurat_nanobody, originalexp = "RNA")
seurat_nanobody <- NormalizeData(seurat_nanobody)
seurat_nanobody <- FindVariableFeatures(seurat_nanobody, selection.method = "vst", nfeatures = 2000)
seurat_nanobody <- ScaleData(seurat_nanobody, features = all.genes)
seurat_nanobody <- RunPCA(seurat_nanobody, npcs = 30, verbose = FALSE)

## run parameter optimization with paramSweep
sweep.res <- paramSweep(seurat_nanobody, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    geom_point(col = "red") +
    geom_line(linetype = "dashed", col = "blue") +
    labs(title = "The BCmvn distributions")

## define the expected number of doublet cells, expect 5% doublets
nExp <- round(ncol(seurat_nanobody) * 0.05)
seurat_nanobody <- doubletFinder(seurat_nanobody, pN = 0.25, pK = 0.20, nExp = nExp, PCs = 1:30)

## extract correct column name of doublet finder prediction
DF.name <- colnames(seurat_nanobody@meta.data)[grepl("DF.classification", colnames(seurat_nanobody@meta.data))]
cowplot::plot_grid(
    ncol = 2,
    DimPlot(seurat_nanobody, reduction = "umap", label = TRUE) + NoAxes(),
    DimPlot(seurat_nanobody, group.by = DF.name) + NoAxes()
)

## check whether doublets have more detected genes
VlnPlot(seurat_nanobody, features = "total", group.by = DF.name, pt.size = 0.1)
summary(seurat_nanobody@meta.data$sizeFactor)

## remove all predicted doublets from our data
seurat_nanobody_singlet <- seurat_nanobody[, seurat_nanobody@meta.data[, DF.name] == "Singlet"]
dim(seurat_nanobody_singlet)

## 0205-2025 stop here
## Step 12: Repeat previous steps, re-perform clustering

## Step 6_2: Identify the most variable genes (feature selection)

## use the FindVariableFeatures function from the Seurat package to identify highly variable features
seurat_nanobody_singlet <- FindVariableFeatures(seurat_nanobody_singlet, selection.method = "vst", nfeatures = 2000)
top15 <- head(VariableFeatures(seurat_nanobody_singlet), 15)
top30 <- head(VariableFeatures(seurat_nanobody_singlet), 30)
plot1 <- VariableFeaturePlot(seurat_nanobody_singlet)
LabelPoints(plot = plot1, points = top30, repel = TRUE, xnudge = 0, ynudge = 0)


## Step 7_2: Scaling the data
all.genes <- rownames(seurat_nanobody_singlet)
seurat_nanobody_singlet <- ScaleData(seurat_nanobody_singlet, features = all.genes)


## Step 8_2: Principal component analysis (PCA)

## get the first 40 PCs
seurat_nanobody_singlet <- RunPCA(seurat_nanobody_singlet, npcs = 40, verbose = FALSE)
ElbowPlot(seurat_nanobody_singlet, ndims = 40)

VizDimLoadings(seurat_nanobody_singlet, dims = 1:4, reduction = "pca")
VizDimLoadings(seurat_nanobody_singlet, dims = 5:8, reduction = "pca")
VizDimLoadings(seurat_nanobody_singlet, dims = 9:12, reduction = "pca")
VizDimLoadings(seurat_nanobody_singlet, dims = 13:16, reduction = "pca")
VizDimLoadings(seurat_nanobody_singlet, dims = 17:20, reduction = "pca")
VizDimLoadings(seurat_nanobody_singlet, dims = 21:24, reduction = "pca")
VizDimLoadings(seurat_nanobody_singlet, dims = 25:28, reduction = "pca")

DimPlot(seurat_nanobody_singlet, reduction = "pca")

DimHeatmap(seurat_nanobody_singlet, dims = 1:9, cells = 300, balanced = TRUE)
DimHeatmap(seurat_nanobody_singlet, dims = 10:18, cells = 300, balanced = TRUE)


## determine the 'dimensionality' of the dataset
seurat_nanobody_singlet <- JackStraw(seurat_nanobody_singlet, num.replicate = 100)
seurat_nanobody_singlet <- ScoreJackStraw(seurat_nanobody_singlet, dims = 1:20)
JackStrawPlot(seurat_nanobody_singlet, dims = 1:20)



## Step 9_2: Non-linear demensional reduction (UMAP/t-SNE)
## use the first 30 PCs
seurat_nanobody_singlet <- RunTSNE(seurat_nanobody_singlet, reduction = "pca", dims = 1:30)
DimPlot(seurat_nanobody_singlet, reduction = "tsne", label = FALSE)
seurat_nanobody_singlet <-RunUMAP(seurat_nanobody_singlet, reduction = "pca", dims = 1:30)
DimPlot(seurat_nanobody_singlet, reduction = "umap", label = FALSE)


## Step 10_2: cluster the cells, set up the resolution as 0.3 (it might need to be modified in the future)
seurat_nanobody_singlet <- FindNeighbors(seurat_nanobody_singlet, reduction = "pca", dims = 1:30)
seurat_nanobody_singlet <- FindClusters(seurat_nanobody_singlet, resolution = 1.5)
DimPlot(seurat_nanobody_singlet, reduction = "umap", label = TRUE)
DimPlot(seurat_nanobody_singlet, reduction = "tsne", label = TRUE)

## Take a look at cluster sizes
table(seurat_nanobody_singlet@meta.data$seurat_clusters)

## check the total counts and detected genes number for each cluster
counts_genes_cluster_sens2_singlet <-  seurat_nanobody_singlet@meta.data %>%
    group_by(seurat_clusters) %>%
    dplyr::summarise(mean_total = mean(total),
                     mean_gene = mean(detected),
                     mean_size = mean(sizeFactor))

counts_genes_cluster_sens2_singlet

FeaturePlot(seurat_nanobody_singlet, features = "nrv3", label = TRUE)
FeaturePlot(seurat_nanobody_singlet, features = "nSyb", label = TRUE)
FeaturePlot(seurat_nanobody_singlet, features = "elav", label = TRUE)
FeaturePlot(seurat_nanobody_singlet, features = "brp", label = TRUE)

FeaturePlot(seurat_nanobody_singlet, features = "GFP", label = TRUE)
FeaturePlot(seurat_nanobody_singlet, features = "mCherry", label = TRUE)

FeaturePlot(seurat_nanobody_singlet, features = "sizeFactor", label = TRUE)

## plot clusters as a tree
library(ape)
seurat_nanobody_singlet <- BuildClusterTree(object = seurat_nanobody_singlet)
PlotClusterTree(object = seurat_nanobody_singlet)


## find markers for every cluster compared to all remaining cells, report only the positive ones
nanobody.singlet.markers <- FindAllMarkers(seurat_nanobody_singlet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
nrow(nanobody.singlet.markers)

## save the Drosophila peripheral neurons singlet seurat object
saveRDS(seurat_nanobody_singlet, file = "/home/rstudio/Projects/Clemens_sc_202502/Outputs/seurat_nanobody_singlet.rds")

## save the markers for all clusters
write.table(nanobody.singlet.markers, sep="\t", file.path(mainDir, paste0("Outputs/nanobody_singlet_markers_res1.5",".txt")), quote=FALSE, col.names=TRUE, row.names=FALSE)



## generate expression heatmap for all clusters
top5_markers_nanobody <- nanobody.singlet.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
View(top5_markers_nanobody)
DoHeatmap(seurat_nanobody_singlet, features = top5_markers_nanobody$gene) + NoLegend()

## save file
save.image("/home/rstudio/Projects/Clemens_sc_202502/Outputs/sce_nanobody_singletClusters.RData")

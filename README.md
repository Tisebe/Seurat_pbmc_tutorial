# Seurat_pbmc_tutorial
Single_cell RNA tutorial


library(dplyr)
library(Seurat)
 ## load data

pbmc.data <- Read10X(data.dir = 'C:/Users/Tony Isebe/Desktop/R tutorials/filtered_gene_bc_matrices/hg19')

### initialize Seurat object

pbmc <-  CreateSeuratObject(counts = pbmc.data, project = 'pbmc3k', min.cells = 3, min.features = 200)
pbmc


## Pre-processing workflow based on QC metrics, data normalization and scaling, detection of highly variable features

pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc, pattern = '^MT-') ## Mitochondrial QC metrics

## visualize QC metrics in violin plot

VlnPlot(pbmc, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

## visualizing feature-feature relationships using FeatureScatter

plot1 <- FeatureScatter(pbmc, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
plot2 <- FeatureScatter(pbmc, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
CombinePlots(plots = list(plot1, plot2))

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 &  nFeature_RNA <2500 &  percent.mt < 5)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#### Normalize data

pbmc <- NormalizeData(pbmc)


## Identify highly variable features

pbmc <- FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)

## identify top 10 most variable features


top10 <- head(VariableFeatures(pbmc),10)

View(top10)


### plot variable features with and withou labels

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

CombinePlots(plots = list(plot1,plot2))

### scaling the data

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- ScaleData(pbmc)
#### perform linera dimensional reduction using defined features

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


## visualize PCA results in different ways

print(pbmc[['pca']], dims=1:5, nfeatures=5)

plot3 <- VizDimLoadings(pbmc, dims=1:2, reduction = 'pca')
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)


DimPlot(pbmc, reduction = 'pca')

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)


pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc),5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = 'umap')


saveRDS(pbmc, file = '../output/pbmc_tutorial.rds')
## find all markers of cluster1

cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n=5)
####assign IDs to each cluster
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

### intgration of a new data set

library(Seurat)
library(SeuratData)
InstallData('panc8')

## constructing a new dataset by identifying new anchors

data('panc8')

pancreas.list <- SplitObject(panc8, split.by = 'tech')
pancreas.list <- pancreas.list[c('celseq','celseq2','fluidigmc1','smartseq2')]

## perform standard pre-processing mainly log-normalization and identification of variable features

for(i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]],selection.method = 'vst',
nfeatures=2000, verbose = FALSE)
}

for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

## integration of 3 datasets

reference.list <- pancreas.list[c('celseq','celseq2','smartseq2')]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)


## pass anchors to integrate data function

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

## use integrated data for downstream analysis

library(ggplot2)
library(cowplot)

## set to integrated data

DefaultAssay(pancreas.integrated) <- 'integrated'


## Running standard workflow

pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = 'pca',dims = 1:30)


p1 <- DimPlot(pancreas.integrated, reduction = 'umap',group.by = 'tech')
p2 <- DimPlot(pancreas.integrated, reduction = 'umap', group.by = 'celltype',label = TRUE, repel = TRUE)+ NoLegend()

plot_grid(p1,p2)






library(dplyr)
library(Seurat)
library(tidydr)
library(patchwork)
library(ggplot2)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/Users/xiaofei/Desktop/SeuratApp/app/data/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


pbmc <- FindNeighbors(pbmc, dims = 1:35)
pbmc <- FindClusters(pbmc, resolution = 0.5)


pbmc.umap <- RunUMAP(pbmc, dims = 1:15)
pbmc.tsne <- RunTSNE(pbmc.umap, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap", label = TRUE) + theme_dr()

saveRDS(pbmc, file = "/Users/xiaofei/Desktop/SeuratApp/app/data/pbmc_tutorial.rds")

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
topN <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

library(scCustomize)
library(ComplexHeatmap)
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', 
               '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', 
               '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', 
               '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', 
               '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', 
               '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

ht <- Clustered_DotPlot(seurat_object = pbmc, features = unique(topN$gene), 
                  colors_use_exp = c('#330066','#336699','#66CC66','#FFCC33'),
                  colors_use_idents = my36colors,
                  flip = TRUE, cluster_ident = FALSE, plot_km_elbow = FALSE,
                  x_lab_rotate = 90, row_label_size = 16)

draw(ht, merge_legend=TRUE)

Cluster_Highlight_Plot(seurat_object = pbmc, cluster_name = "1")

DoHeatmap(object = pbmc, features = unique(topN$gene)) +
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))





















# https://satijalab.org/seurat/articles/project3k_tutorial.html

path = "/Users/miyang/Documents/STANFORD/TUTORIALS/Seurat/"

# create a seurat object

library(dplyr)
library(Seurat)
library(patchwork)

# Load the project dataset
project.data <- Read10X(data.dir = paste0(path,"data/project3k/filtered_gene_bc_matrices/hg19/"))
# Initialize the Seurat object with the raw (non-normalized data).
project <- CreateSeuratObject(counts = project.data, project = "project3k", min.cells = 3, min.features = 200)
# counts: Either a matrix-like object with unnormalized data with cells as columns and features as rows or an Assay-derived object
# min.cells: Include features detected in at least this many cells. Will subset the counts matrix as well. 
#            To reintroduce excluded features, create a new object with a lower cutoff.
# min.features: Include cells where at least this many features are detected.

project

# Standard pre-processing workflow

######################## QC and selecting cells for further analysis
#The percentage of reads that map to the mitochondrial genome
#Low-quality / dying cells often exhibit extensive mitochondrial contamination
#We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
#We use the set of all genes starting with MT- as a set of mitochondrial genes
project[["percent.mt"]] <- PercentageFeatureSet(project, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(project, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(project, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(project, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

project <- subset(project, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


######################## Normalizing the data
# After removing unwanted cells from the dataset, the next step is to normalize the data. 
# By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression 
# measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
# Normalized values are stored in project[["RNA"]]@data.
project <- NormalizeData(project, normalization.method = "LogNormalize", scale.factor = 10000)


######################## Identification of highly variable features (feature selection)

project <- FindVariableFeatures(project, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(project), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(project)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

######################## Scaling the data

all.genes <- rownames(project)
project <- ScaleData(project, features = all.genes)


######################## Perform linear dimensional reduction

project <- RunPCA(project, features = VariableFeatures(object = project))

# Examine and visualize PCA results a few different ways
print(project[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(project, dims = 1:2, reduction = "pca")
DimPlot(project, reduction = "pca")
DimHeatmap(project, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(project, dims = 1:15, cells = 500, balanced = TRUE)


######################## Determine the ‘dimensionality’ of the dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
project <- JackStraw(project, num.replicate = 100)
project <- ScoreJackStraw(project, dims = 1:20)
JackStrawPlot(project, dims = 1:15)

ElbowPlot(project)

######################## Cluster the cells

project <- FindNeighbors(project, dims = 1:10)
project <- FindClusters(project, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(project), 5)

######################## Run non-linear dimensional reduction (UMAP/tSNE)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
project <- RunUMAP(project, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(project, reduction = "umap")
saveRDS(project, file = paste0(path,"output/project_tutorial.rds"))
 

######################## Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(project, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(project, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
project.markers <- FindAllMarkers(project, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
project.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(project, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(project, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(project, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(project, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A" ))

project.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(project, features = top10$gene) + NoLegend()


######################## Assigning cell type identity to clusters

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(project)
project <- RenameIdents(project, new.cluster.ids)
DimPlot(project, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(project, file = paste0(path,"output/project3k_final.rds"))

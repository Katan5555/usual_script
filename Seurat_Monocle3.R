
# https://satijalab.org/signac/articles/monocle.html

library(Signac) 
library(Seurat)
library(SeuratWrappers) # remotes::install_github('satijalab/seurat-wrappers')
library(monocle3) # devtools::install_github('cole-trapnell-lab/monocle3')  devtools::install_github('cole-trapnell-lab/leidenbase')
library(Matrix) 
library(ggplot2)
library(patchwork)
set.seed(1234)

path <- "/Users/miyang/Documents/STANFORD/TUTORIALS/Seurat/GSE129785/"
filepath <- "/Users/miyang/Documents/STANFORD/TUTORIALS/Seurat/GSE129785/GSE129785_scATAC-Hematopoiesis-CD34"
setwd(path)

peaks <- read.table(paste0(filepath, ".peaks.txt.gz"), header = TRUE)
cells <- read.table(paste0(filepath, ".cell_barcodes.txt.gz"), header = TRUE, stringsAsFactors = FALSE)
rownames(cells) <- make.unique(cells$Barcodes)

mtx <- readMM(file = paste0(filepath, ".mtx.gz"))
mtx <- as(object = mtx, Class = "dgCMatrix")
colnames(mtx) <- rownames(cells)
rownames(mtx) <- peaks$Feature



# indexing with tabix https://vcf.iobio.io/help.html
# cd /Users/miyang/Documents/STANFORD/
# PACKAGES/TabixBinary/tabix -p bed TUTORIALS/Seurat/GSM3722029_CD34_Progenitors_Rep1_fragments.tsv.gz

setwd(path)
# folder containing both: GSM3722029_CD34_Progenitors_Rep1_fragments.tsv.gz 
#                         GSM3722029_CD34_Progenitors_Rep1_fragments.tsv.gz.tbi
bone_assay <- CreateChromatinAssay(
  counts = mtx,
  min.cells = 5,
  fragments =  "GSM3722029_CD34_Progenitors_Rep1_fragments.tsv.gz",
  sep = c("_", "_"),
  genome = "hg19"
)
bone <- CreateSeuratObject(
  counts = bone_assay,
  meta.data = cells,
  assay = "ATAC"
)

# The dataset contains multiple cell types
# We can subset to include just one replicate of CD34+ progenitor cells
bone <- bone[, bone$Group_Barcode == "CD34_Progenitors_Rep1"]

# add cell type annotations from the original paper
cluster_names <- c("HSC",   "MEP",  "CMP-BMP",  "LMPP", "CLP",  "Pro-B",    "Pre-B",    "GMP",
                   "MDP",    "pDC",  "cDC",  "Monocyte-1",   "Monocyte-2",   "Naive-B",  "Memory-B",
                   "Plasma-cell",    "Basophil", "Immature-NK",  "Mature-NK1",   "Mature-NK2",   "Naive-CD4-T1",
                   "Naive-CD4-T2",   "Naive-Treg",   "Memory-CD4-T", "Treg", "Naive-CD8-T1", "Naive-CD8-T2",
                   "Naive-CD8-T3",   "Central-memory-CD8-T", "Effector-memory-CD8-T",    "Gamma delta T")
num.labels <- length(cluster_names)
names(cluster_names) <- paste0( rep("Cluster", num.labels), seq(num.labels) )
bone$celltype <- cluster_names[as.character(bone$Clusters)]

bone[["ATAC"]]


# add gene annotations for the hg19 genome to the object. 
# This will be useful for computing quality control metrics (TSS enrichment score) and plotting
library(EnsDb.Hsapiens.v75)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(bone) <- annotations


###################################### Quality control

# compute TSS enrichment, nucleosome signal score, and the percentage of counts in genomic blacklist regions for each cell, 
# and use these metrics to help remove low quality cells from the datasets.

bone <- TSSEnrichment(bone)
bone <- NucleosomeSignal(bone)
bone$blacklist_fraction <- FractionCountsInRegion(bone, regions = blacklist_hg19)

VlnPlot(
  object = bone,
  features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "blacklist_fraction"),
  pt.size = 0.1,
  ncol = 4
)


bone <- bone[, (bone$nCount_ATAC < 50000) &
               (bone$TSS.enrichment > 2) & 
               (bone$nucleosome_signal < 5)]



###################################### Dataset preprocessing

bone <- FindTopFeatures(bone, min.cells = 10)
bone <- RunTFIDF(bone)
bone <- RunSVD(bone, n = 100)
DepthCor(bone)

bone <- RunUMAP(
  bone,
  reduction = "lsi",
  dims = 2:50,
  reduction.name = "UMAP"
)

bone <- FindNeighbors(bone, dims = 2:50, reduction = "lsi")
bone <- FindClusters(bone, resolution = 0.8, algorithm = 3)

DimPlot(bone, label = TRUE) + NoLegend()

for(i in levels(bone)) {
  cells_to_reid <- WhichCells(bone, idents = i)
  newid <- names(sort(table(bone$celltype[cells_to_reid]),decreasing=TRUE))[1]
  Idents(bone, cells = cells_to_reid) <- newid
}
bone$assigned_celltype <- Idents(bone)

DimPlot(bone, label = TRUE)

# subset the different lineages and create a trajectory for each lineage. 
# Another way to build the trajectories is to use the whole dataset and build separate pseudotime trajectories 
# for the different cell partitions found by Monocle 3.

DefaultAssay(bone) <- "ATAC"

erythroid <- bone[,  bone$assigned_celltype %in% c("HSC", "MEP", "CMP-BMP")]
lymphoid <- bone[, bone$assigned_celltype %in% c("HSC", "LMPP", "GMP", "CLP", "Pro-B", "pDC", "MDP", "GMP")]


###################################### Building trajectories with Monocle 3
 
erythroid.cds <- as.cell_data_set(erythroid)
erythroid.cds <- cluster_cells(cds = erythroid.cds, reduction_method = "UMAP")
erythroid.cds <- learn_graph(erythroid.cds, use_partition = TRUE)

lymphoid.cds <- as.cell_data_set(lymphoid)
lymphoid.cds <- cluster_cells(cds = lymphoid.cds, reduction_method = "UMAP")
lymphoid.cds <- learn_graph(lymphoid.cds, use_partition = TRUE)

# load the pre-selected HSCs
# Decide what the start of each trajectory is. We know that the hematopoietic stem cells are the progenitors of other cell types 
hsc <- readLines("hsc_cells.txt")
 
# order cells
# in the trajectory, so we can set these cells as the root of the trajectory. 
erythroid.cds <- order_cells(erythroid.cds, reduction_method = "UMAP", root_cells = hsc)
lymphoid.cds  <- order_cells(lymphoid.cds, reduction_method = "UMAP", root_cells = hsc)

# plot trajectories colored by pseudotime
plot_cells(
  cds = erythroid.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)


plot_cells(
  cds = lymphoid.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

# Extract the pseudotime values and add to the Seurat object
bone <- AddMetaData(
  object = bone,
  metadata = erythroid.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Erythroid"
)

bone <- AddMetaData(
  object = bone,
  metadata = lymphoid.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Lymphoid"
)

FeaturePlot(bone, c("Erythroid", "Lymphoid"), pt.size = 0.1) & scale_color_viridis_c()









 
# read and write feature file
path_data <- "/seurat/"
features <- read.table(gzfile("/seurat/features.tsv.gz"))
features$V2 <- gsub("\\_","-",features$V2)
features$V5 <- vapply(strsplit( features$V2 , split="[-]"), "[", "", 2)

write.table(features, file=gzfile("/seurat/features.tsv.gz"),col.names = FALSE,row.names=FALSE, sep="\t")


# Load the dataset
project.data <- Read10X(data.dir = path_data )
# Initialize the Seurat object with the raw (non-normalized data).
project <- CreateSeuratObject(counts = project.data , min.cells = 3, min.features = 200)


# subset of certain genes
features <- read.table(gzfile("/seurat/features.tsv.gz"))
features <- features[!duplicated(features$V5), ]
project <- subset(x=project, features=features$V2  )

# change gene names
# https://github.com/satijalab/seurat/issues/1049
genes <- rownames(project)
features <- features[features$V2 %in% genes, ]
new_gene_names <- features$V5

project@assays$RNA@counts@Dimnames[[1]] <- new_gene_names
project@assays$RNA@data@Dimnames[[1]] <- new_gene_names
rownames(project@assays$RNA@meta.features) <- new_gene_names


# add cell type names 
sample_info <- read.csv("/Users/i0535027/Documents/Sanofi/DATABASE/SingleCell/MM/DATASET_Cho_Kim/DATA/GSE155795_metadata_1656537305789.tsv",check.names = F,sep="\t")
sample_info <- sample_info[grep("GSM4743044_Pt1",sample_info$Barcodes), ]
sample_info$Barcodes <- vapply(strsplit( sample_info$Barcodes , split="[_]"), "[", "", 3)
rownames(sample_info) <- sample_info$Barcodes  
sample_info <- sample_info[colnames(project), ]
sample_info$`Author's cell type`[is.na(sample_info$`Author's cell type`)==T] <- "other"

head( testdata@meta.data )
testdata@active.ident
head( project@meta.data ) 

project@meta.data$seurat_annotations <-  as.factor( sample_info$`Author's cell type` ) 
project@active.ident <- as.factor( sample_info$`Author's cell type` ) 
names(project@active.ident) <- colnames(project)

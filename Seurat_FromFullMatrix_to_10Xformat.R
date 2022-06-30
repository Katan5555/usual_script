

path_data <- "XXX"
mat <- read.table(gzfile(paste0(path_data,"matrix.txt.gz")), header = T,row.names = 1)
 
barcodes <- data.frame(colnames(mat))
write.table(barcodes, file=gzfile(paste0(path_data,"/seurat/barcodes.tsv.gz")),col.names = FALSE,row.names=FALSE, sep="\t")

features <- data.frame(rownames(mat))
write.table(features, file=gzfile(paste0(path_data,"/seurat/features.tsv.gz")),col.names = FALSE,row.names=FALSE, sep="\t")
 

library(Matrix)
mat <- data.matrix(mat)
sp_matrix <- as(mat, "sparseMatrix")
writeMM(obj=sp_matrix, file= paste0(path_data,"seurat/matrix.mtx.gz") ,col.names = FALSE,row.names=FALSE, sep="\t")



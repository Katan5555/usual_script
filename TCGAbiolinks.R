
# https://bioc.ism.ac.jp/packages/3.2/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html
 
library(TCGAbiolinks)
directory = "/Users/miyang/Documents/STANFORD/DATABASE/TCGA/TCGAbiolinks/"
 
tissue <- "TCGA-BRCA"

#### if you don't know the list of possible entry, just write a wrong one, it will give you the entire list of entries to choose from. ###
query <- GDCquery(project = tissue,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification")

# find mapping between download IDs and TCGA patient IDs. 
result <- getResults(query)

# download the data
GDCdownload(query, method = "api",directory=directory)

files <- list.dirs(paste0(directory,tissue,"/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/"))
files <- vapply(strsplit( files , split="[/]"), "[", "",  14)
files <- files[!is.na(files)]

# map IDs
common <- intersect(files, result$id) 
result <- result[result$id %in% common,c("id","cases.submitter_id","sample_type")]
result <- result[result$sample_type=="Primary Tumor", ]
result <- result[!duplicated(result$cases.submitter_id), ]
rownames(result) <- result$cases.submitter_id

# clinical data
clinical <- GDCquery_clinic(project = tissue, type = "clinical", save.csv = F)
 
common <- intersect(clinical$submitter_id, result$cases.submitter_id) 
clinical <- clinical[clinical$submitter_id %in% common, ]
result  <- result[clinical$submitter_id, ]

# assemble gene expression data
GEX <- c()
for(s in rownames(result)) {
  print(s)
  work_dir <- paste0(directory,tissue,"/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/",result[s,c("id")],"/")
  file <- list.files(work_dir)
  df <- read.delim(paste0(work_dir,file), comment.char="#")
  df <- df[df$gene_id %not in% c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous"), ]
  GEX <- cbind(GEX, df$unstranded)
}

rownames(GEX) <- df$gene_name
colnames(GEX) <- rownames(result)

# normalize by library size
GEX_library_size <- colSums(GEX)
GEX_norm <- sweep(GEX, 2, GEX_library_size, "/")
GEX_norm <- GEX_norm*mean(GEX_library_size)

# log2
GEX_prot <- convert_protein_product_df_row(GEX_norm)
GEX_prot_log2 <- log2(GEX_prot+1)

# save data
dir.create(file.path(paste0(directory,tissue,"/DATA/")), recursive = T, showWarnings = F)
GEX_list <- list(GEX_norm,GEX_prot,GEX_prot_log2)
names(GEX_list) <- c("GEX","GEX_prot","GEX_prot_log2")
save(GEX_list,file=paste0(directory,tissue,"/DATA/GEX_list.Rdata"))
write.csv(GEX_prot_log2,paste0(directory,tissue,"/DATA/GEX.csv"))

write.csv(clinical, paste0(directory,tissue,"/DATA/clinical_original.csv") )





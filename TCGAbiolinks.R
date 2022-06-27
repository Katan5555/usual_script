
# https://bioc.ism.ac.jp/packages/3.2/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html
 
library(TCGAbiolinks)
directory = "/Users/miyang/Documents/STANFORD/DATABASE/TCGA/TCGAbiolinks/"


#### if you don't know the list of possible entry, just write a wrong one, it will give you the entire list of entries to choose from. ###
query <- GDCquery(project = "TCGA-DLBC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification")

# find mapping between download IDs and TCGA patient IDs. 
result <- getResults(query)

# download the data
GDCdownload(query, method = "api",directory=directory)

files <- list.dirs("/Users/miyang/Documents/STANFORD/DATABASE/TCGA/TCGAbiolinks/TCGA-DLBC/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/")
files <- vapply(strsplit( files , split="[/]"), "[", "",  14)
files <- files[!is.na(files)]
intersect(files, result$id) 


# clinical data

clinical <- GDCquery_clinic(project = "TCGA-DLBC", type = "clinical", save.csv = F)
clinical_id <- clinical[ ,grep("_id", colnames(clinical))]

query <- GDCquery(
  project = "TCGA-DLBC", 
  data.category = "Clinical",
  data.type = "Clinical Supplement"
)
result <- getResults(query)
GDCdownload(query, method = "api",directory=directory)





load.lib<-c("Rdpack","TFisher","mutoss")

load.lib<-c("glmnet","doParallel","ggplot2","caret","Hmisc","randomForest","randomForestSRC","survminer",
            "hydroGOF","pROC","glmnet","rsample","dplyr","GSA","SuppDists","nnls","e1071","devtools","corrplot",
            "HiveR","pheatmap","beeswarm","plotflow","colorRamps","metap","HiClimR","Rtsne","Seurat","xgboost")

install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")
BiocManager::install("GSA")
BiocManager::install("survcomp")
BiocManager::install("preprocessCore")
BiocManager::install("sva")
BiocManager::install("multtest")
BiocManager::install("affyPLM")
BiocManager::install("qvalue")
BiocManager::install("aff")  # package ‘aff’ is not available (for R version 3.5.1)
BiocManager::install("biomaRt")


remove.packages("devtools")
install.packages("devtools")


library(devtools)
install_github("trinker/plotflow")
path <- "/home/users/miyang/STANFORD"  # Sherlock
path <- "/home/alizadehlab/miyang/STANFORD" # server lung
install.packages(paste0(path,"/PACKAGES/bapred_1.0.tgz"),repos=NULL,type="source")
install.packages(paste0(path,"/PACKAGES/bapred_0.1.tar.gz"),repos=NULL,type="source")

install.packages("lme4.0",repos=c("http://lme4.r-forge.r-project.org/repos",getOption("repos")[["CRAN"]]))
install.packages("lme4",repos=c("http://lme4.r-forge.r-project.org/repos",getOption("repos")[["CRAN"]]), dependencies = TRUE)


path <- "/home/alizadehlab/miyang/STANFORD" # server lung
install.packages(paste0(path,"/PACKAGES/lme4_1.1-23.tar.gz"),repos=NULL,type="source")

install.packages(paste0(path,"/PACKAGES/metap_1.4.tar.gz"),repos=NULL,type="source", dependencies = TRUE )
install.packages(paste0(path,"/PACKAGES/bapred_0.1.tar.gz"), repos = NULL, type = "source")
load.lib<-c("Rdpack","TFisher","mutoss","multtest","RcppArmadillo","conquer","quantreg","sn")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)


####################### install HiClimR on Sherlock #######################
# remove copy of ncdf4 since R is looking here - /home/users/miyang/R/x86_64-pc-linux-gnu-library/3.5/ncdf4/libs/ncdf4.so
# ml load R
# ml load gcc
# ml load devel netcdf
# R
# install.packages("HiClimR", repos = "https://cloud.r-project.org/", dependencies = TRUE)


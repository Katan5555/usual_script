# https://omnideconv.org/immunedeconv/articles/immunedeconv.html
# https://github.com/omnideconv/immunedeconv/blob/master/R/immune_deconvolution_methods.R
# https://github.com/omnideconv/immunedeconv

# QuantiSeq: https://icbi.i-med.ac.at/software/quantiseq/doc/index.html
# ABIS:https://github.com/giannimonaco/ABIS



# Input Data
# The input data is a gene × sample gene expression matrix. In general values should be
# 
# TPM-normalized, not log-transformed.
# For xCell and MCP-counter this is not so important. xCell works on the ranks of the gene expression only and MCP-counter sums up the gene expression values.
# 
# Rownames are expected to be HGNC gene symbols. Instead of a matrix, immunedeconv also supports ExpressionSets (see below).


# Immune cell fractions
library(immunedeconv)

df_immune <- deconvolute(df, "mcp_counter")
df_immune <- deconvolute(df, "quantiseq" )
df_immune <- deconvolute(df, "abis" )

#' Deconvolute using TIMER
#'
#' Unlike the other methods, TIMER needs the specification of the
#' cancer type for each sample.
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples
#' @param indications a n-vector giving and indication string (e.g. 'brca') for each sample.
#'     Accepted indications are 'kich', 'blca', 'brca', 'cesc', 'gbm', 'hnsc', 'kirp', 'lgg',
#'     'lihc', 'luad', 'lusc', 'prad', 'sarc', 'pcpg', 'paad', 'tgct',
#'     'ucec', 'ov', 'skcm', 'dlbc', 'kirc', 'acc', 'meso', 'thca',
#'     'uvm', 'ucs', 'thym', 'esca', 'stad', 'read', 'coad', 'chol'
#' @export


df_immune <- deconvolute(df, "timer", indications=indications )
df_immune <- deconvolute(df, "consensus_tme", indications=indications )










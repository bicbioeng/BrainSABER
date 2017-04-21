#' Get an example vector for specified genes
#'
#' This function returns a named example vector of gene expression values for
#' the specified genes, taken from the 1st row of dataSet, for use in
#' demonstrating getSimScores.
#'
#' @param genes a character vector of HGNC-compliant gene names
#' @param dataSet a Biobase ExpressionSet
#'
#' @return a named character vector of gene-expression values
#' @import Biobase
#' @export

getExampleVector <- function(genes, dataSet) {
  #get indices
  geneIdx <- which(fData(dataSet)$gene_symbol %in% genes, arr.ind = TRUE)
  #get relevantGenes of dataSet
  relevantGenes <- dataSet[geneIdx]
  #get 8pcw exprs
  v <- as.vector(exprs(relevantGenes[, 1]))
  #name v
  names(v) <- as.character(fData(relevantGenes)$gene_symbol)
  #return v
  return(v)
}

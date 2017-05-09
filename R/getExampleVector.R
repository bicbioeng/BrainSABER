#' Get an example vector for specified genes
#'
#' This function returns a named example vector of gene expression values for
#' the specified genes, taken from the 1st row of AIBSARNA, for use in
#' demonstrating getSimScores.
#'
#' @param genes a character vector of HGNC-compliant gene names
#'
#' @return a named character vector of gene-expression values
#' @import Biobase
#' @export
#' @examples
#' myGenes <- c("TNFRSF1A", "BCL3", "NEFH")
#' myExampleVector <- getExampleVector(myGenes)

getExampleVector <- function(genes) {
  #get relevant genes
  names(genes) <- genes
  relevantGenes <- getRelevantGenes(genes, gene_names = "HGNC")
  #get 8pcw exprs
  v <- as.vector(exprs(relevantGenes[, 1]))
  #name v
  names(v) <- as.character(fData(relevantGenes)$gene_symbol)
  #return v
  return(v)
}

#' Get a subset of AIBSARNA using a Gene Expression Vector
#'
#' This function returns a subset of the AIBSARNA dataset, containing only
#' the genes listed in the parameter vector, v.  The vector must consist of
#' numerical gene expression values named with HGNC-conformant gene names.
#'
#' @param v a vector of HGNC-conformant named gene expression values
#'
#' @return a Biobase ExpressionSet consisting of genes in v
#' @import Biobase
#' @export
#' @examples
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes)
#'
getRelevantGenes <- function(v){
  #relevantGenes <- subset of AIBSARNA containing only genes specified in V
  vInd <- which(fData(AIBSARNA::AIBSARNA)$gene_symbol %in% names(v),
                arr.ind = TRUE)
  relevantGenes <- AIBSARNA::AIBSARNA[vInd]
  return(relevantGenes)
}

#' Get a subset of AIBSARNA using a Gene Expression Vector
#'
#' This function returns a subset of the AIBSARNA dataset, containing only
#' the genes listed in the parameter vector, v.  The vector must consist of
#' numerical gene expression values named with HGNC-conformant gene names for
#' \code{gene_names = "HGNC"}, or named with RefSeq IDs for \code{gene_names =
#' "RefSeq"}.
#'
#' @param v a vector of HGNC-conformant or RefSeq ID named gene expression values
#' @param gene_names either "HGNC" or "RefSeq", defaults to "HGNC"
#'
#' @return a Biobase ExpressionSet consisting of genes in v
#' @import Biobase
#' @export
#' @examples
#'
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes, gene_names = "HGNC")
#'
getRelevantGenes <- function(v, gene_names = "HGNC"){
  #relevantGenes <- subset of AIBSARNA containing only genes specified in v
  if (gene_names == "RefSeq"){
    vInd <- which(fData(AIBSARNA::AIBSARNA)$refseq_ids %in% names(v),
                  arr.ind = TRUE)
  } else {
    vInd <- which(fData(AIBSARNA::AIBSARNA)$gene_symbol %in% names(v),
                  arr.ind = TRUE)
  }
  # subset feature data
  relfd <- fData(AIBSARNA::AIBSARNA)[vInd, ]
  # convert to Annotated Data Frame
  relfd <- new("AnnotatedDataFrame", data = relfd)
  # subset exprs
  relexprs <- exprs(AIBSARNA::AIBSARNA)[vInd, ]
  # put together ExpressionSet
  relevantGenes <- ExpressionSet(assayData = relexprs,
                                 phenoData = phenoData(AIBSARNA::AIBSARNA),
                                 featureData = relfd)
  return(relevantGenes)
}

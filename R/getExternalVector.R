#' getExternalVector
#'
#' Get a named vector of gene expression values from a single sample in an
#' outside Biobase ExpressionSet, for use in creating subsets of AIBSARNA with
#'  \code{getRelevantGenes} and comparison with that subset with
#'  \code{getSimScores}
#'
#' @param dataSet a Biobase ExpressionSet
#' @param index the integer index of the sample of dataSet to be used
#' @param dataSetId the name of the column of gene identifiers in fData(dataSet)
#'     to be used to compare dataSet to AIBSARNA.
#' @param AIBSARNAid the name of the column of fData(AIBSARNA) that is
#'     comparable to dataSetId.  One of "gene_id", "ensembl_gene_id",
#'     "gene_symbol", "entrez_id", "refseq_ids"
#'
#' @return a named vector of gene expression values
#' @export
#' @import Biobase
#'
#' @examples
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes, gene_names = "HGNC")
#' myGeneSampleVector <- getExternalVector(myGeneSet, index = 1,
#'      dataSetId = "gene_symbol", AIBSARNAid = "gene_symbol")
#'
getExternalVector <- function(dataSet, index = 1, dataSetId = "gene_symbol",
                              AIBSARNAid = c("gene_id", "ensembl_gene_id",
                                                  "gene_symbol", "entrez_id",
                                                  "refseq_ids")){
  v <- exprs(dataSet)[, index]
  names(v) <- as.character(fData(dataSet)[[dataSetId]])
  # get gene identifiers common to v and AIBSARNA
  vInAIBSARNA <- which(fData(AIBSARNA::AIBSARNA)[[AIBSARNAid]] %in%
                         names(v), arr.ind = TRUE)
  vInAIBSARNA <- as.character(
                      fData(AIBSARNA::AIBSARNA)[[AIBSARNAid]][vInAIBSARNA])
  # get indices of v that are present in vInAIBSARNA
  genesToKeep <- which(names(v) %in% vInAIBSARNA, arr.ind = TRUE)
  # update v to only include comparable genes (genes in vInAIBSARNA)
  v <- v[genesToKeep]
  # remove any duplicate genes
  genesToKeep <- unique(names(v))
  v <- v[genesToKeep]
  return(v)
}

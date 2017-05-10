#' Get a subset of AIBSARNA using a Gene Expression Vector
#'
#' This function returns a subset of the AIBSARNA dataset, containing only
#' the genes \code{data}, which may be a vector or a Biobase ExpressionSet, or a
#' derivative of ExpressionSet that supports the Biobase methods \code{exprs()}
#' and \code{fData()}.  If a vector is used, it must consist of numerical gene
#' expression values with names comparable to one columns of identifiers present
#' in AIBSARNA.
#'
#' @param data a vector of named gene expression values, or a compatible data set
#' @param dataSetId (Optional) If \code{data} is not a vector, the name of the
#'     column of gene identifiers in fData(dataSet) to be used to compare
#'     \code{data} to AIBSARNA.
#' @param AIBSARNAid the name of the column of fData(AIBSARNA) that is
#'     comparable to dataSetId.  One of "gene_id", "ensembl_gene_id",
#'     "gene_symbol", "entrez_id", "refseq_ids"
#'
#' @return a Biobase ExpressionSet consisting of genes in v
#' @import Biobase
#' @import AIBSARNA
#' @export
#' @examples
#'
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes, AIBSARNAid = "gene_symbol")
#'
getRelevantGenes <- function(data, dataSetId = NULL,
                             AIBSARNAid = c("gene_id", "ensembl_gene_id",
                                            "gene_symbol", "entrez_id",
                                            "refseq_ids")){
  # get a single sample vector if data is not a vector
  if(is.vector(data)){
    v <- data
  } else {
    #stop if data is not a vector and dataSetId == NULL
    stopifnot(!is.null(dataSetId))
    # get a trimmed vector faciliate things
    v <- getExternalVector(dataSet = data, index = 1, dataSetId = dataSetId,
                           AIBSARNAid = AIBSARNAid)
  }
  #relevantGenes <- subset of AIBSARNA containing only genes specified in v
  # get indices of genes common to v and AIBSARNA
  vInd <- which(fData(AIBSARNA::AIBSARNA)[[AIBSARNAid]] %in%
                    names(v), arr.ind = TRUE)
  # subset feature data
  relfd <- fData(AIBSARNA::AIBSARNA)[vInd, ]
  # remove any unused factor levels
  relfd <- as.data.frame(apply(relfd, 2, function(x) {x[drop = TRUE]}))
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

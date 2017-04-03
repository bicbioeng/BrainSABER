#' Get Cosine Similarity Scoring for a Gene Expression Vector
#'
#' This function computes the cosine similarity score of a gene expression
#' vector as compared to the AIBSARNA dataset, though other dataSets may be
#' supported in the future.  The vector must consist of numerical expression
#' values named with HGNC-conformant gene names.
#'
#' @param v a vector of HGNC-conformant named gene expression values
#' @param dataSet the AIBSARNA dataSet
#'
#' @return a Biobase ExpressionSet, with similarity scores stored as part of the
#'    sample data (phenoData)
#' @import Biobase
#' @export

getCosineSim <- function(v, dataSet){
  #relevantGenes <- subset of AIBSARNA containing only genes specified in V
  relevantGenes <- dataSet[featureData(dataSet)$gene_symbol == names(v), ]
  nrows <- length(v)
  ncols <- lengths(relevantGenes[1,]) #returns row length
  #initialize Similarity_Score vector
  cosine_similarity <- lsa::cosine(exprs(relevantGenes)[, 1], v)
  #get similarity score of v with each gene of relevantGenes
  for(i in 2:ncols){
      cosine_similarity <- c(cosine_similarity, lsa::cosine(exprs(
        relevantGenes)[, i], v))
  }
  #get similarity matrix
  #bind cosine_similarity to relevantGenes
  phenoData(relevantGenes)@data <- cbind(phenoData(relevantGenes)@data,
                                           cosine_similarity)
  return(relevantGenes)
}

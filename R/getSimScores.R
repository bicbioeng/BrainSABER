#' Get Similarity Scoring for a Gene Expression Vector
#'
#' This function computes the similarity score of a gene expression
#' vector as compared to an ExpressionSet of genes, obtained by
#' \code{getRelevantGenes}.  The vector must consist of numerical expression
#' values named with HGNC-conformant gene names, and should be the same vector
#' used to generate the ExpressionSet.
#'
#' @param v a vector of HGNC-conformant named gene expression values
#' @param relevantGenes an ExpressionSet returned by \code{getRelevantGenes}
#' @param similarity_method currently supported similarity methods are "cosine"
#'     and "euclidean", defaults to "cosine"
#'
#' @return relevantGenes, with similarity scores stored as part of the sample
#'     data (phenoData)
#' @import Biobase
#' @export
#' @examples
#' \dontrun{myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes)
#' myGeneSet <- getSimScores(myGenes, myGeneSet, similarity_method = "cosine")
#' myGeneSet <- getSimScores(myGenes, myGeneSet, similarity_method = "euclidean")}

getSimScores <- function(v, relevantGenes, similarity_method = "cosine"){
  #get row length
  ncols <- lengths(relevantGenes[1,])
  if (similarity_method == "euclidean"){
    #initialize Similarity_Score vector
    euclidean_similarity <- euclideanSim(exprs(relevantGenes)[, 1], v)
    #get similarity score of v with each gene of relevantGenes
    for(i in 2:ncols){
      euclidean_similarity <- c(euclidean_similarity, euclideanSim(exprs(
        relevantGenes)[, i], v))
    }
    #bind euclidean_similarity to relevantGenes
    phenoData(relevantGenes)@data <- cbind(phenoData(relevantGenes)@data,
                                           euclidean_similarity)
  } else {
    #initialize Similarity_Score vector
    cosine_similarity <- lsa::cosine(exprs(relevantGenes)[, 1], v)
    #get similarity score of v with each gene of relevantGenes
    for(i in 2:ncols){
        cosine_similarity <- c(cosine_similarity, lsa::cosine(exprs(
          relevantGenes)[, i], v))
    }
    #bind cosine_similarity to relevantGenes
    phenoData(relevantGenes)@data <- cbind(phenoData(relevantGenes)@data,
                                           cosine_similarity)
  }
  return(relevantGenes)
}

#helper function
euclideanSim <- function(x,y){
  tmp <- sqrt(sum((x-y)^2))
  return(1/(1+tmp))
}

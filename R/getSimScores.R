#' Get Similarity Scoring for a Gene Expression Vector
#'
#' This function computes the similarity score of a gene expression
#' vector returned by \code{getExternalVector} or a trimmed data set returned by
#' \code{getTrimmedExternalSet}, compared to a subset of AIBSARNA,
#' obtained by \code{getRelevantGenes}.
#'
#' @param data a named vector of gene expression values returned by
#'     \code{getExternalVector}, or an ExpressionSet returned by
#'     \code{getTrimmedExternalSet}
#' @param relevantGenes an ExpressionSet returned by \code{getRelevantGenes}
#' @param similarity_method currently supported similarity methods are "cosine"
#'     and "euclidean", defaults to "cosine"
#'
#' @return If \code{data} is a vector, returns a vector of similarity scores for
#'     each sample in relevantGenes. If \code{data} is an ExpressionSet, returns
#'     a data frame, with columns containing the similarity scores for and named
#'     after each sample in \code{data}, and rows named after each sample in
#'     \code{relevantGenes}
#' @import Biobase
#' @export
#' @examples
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes)
#' CosScores <- getSimScores(myGenes, myGeneSet, similarity_method = "cosine")
#' EucScores <- getSimScores(myGenes, myGeneSet, similarity_method = "euclidean")
#' myCosSimFrame <- getSimScores(myGeneSet, myGeneSet)

getSimScores <- function(data, relevantGenes, similarity_method = "cosine"){
  if (is.vector(data)){
    if (similarity_method == "euclidean"){
      euclidean_similarity <- apply(X = exprs(relevantGenes), MARGIN = 2,
                                    FUN = euclideanSim, y = data)
      return(euclidean_similarity)
    } else {
      cosine_similarity <- apply(X = exprs(relevantGenes), MARGIN = 2,
                                 FUN = lsa::cosine, y = data)
      return(cosine_similarity)
    }
  } else {
    simFrame <- as.data.frame(apply(X = exprs(data), MARGIN = 2,
                                    FUN = getSimScores,
                                    relevantGenes = relevantGenes,
                                    similarity_method = similarity_method))
    return(simFrame)
  }
}

#helper function
euclideanSim <- function(x,y){
  tmp <- sqrt(sum((x-y)^2))
  return(1/(1+tmp))
}

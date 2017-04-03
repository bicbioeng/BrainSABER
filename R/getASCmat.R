#create matrix with rows = age, cols = structure_acronym, value = cosine_similarity
#' Get Age, Structure Acronym, and Cosine Similarity Matrix
#'
#' This function takes a \code{Biobase ExpressionSet} as returned by
#' \code{getCosineSim} and returns a numeric matrix of cosine_similarity with
#' rows labeled by age, and columns labeled by structure_acronym.
#'
#' @param relevantGenes a Biobase ExpressionSet
#'
#' @return a numeric matrix of cosine similaritys
#' @import Biobase
#' @export

getASCmat <- function(relevantGenes){
  #get data vectors
  age <- phenoData(relevantGenes)$age
  structure_acronym <- phenoData(relevantGenes)$structure_acronym
  cosine_similarity <- phenoData(relevantGenes)$cosine_similarity
  #get the dimensions
  nrows <- nlevels(age)
  ncols <- nlevels(structure_acronym)
  #create and prefill the matrix, in case not all ages have all structures
  mat <- matrix(data = rep(0, nrows * ncols), nrow = nrows,
                ncol = ncols)
  rownames(mat) <- as.character(levels(age))
  colnames(mat) <- as.character(levels(structure_acronym))
  #populate the matrix with cosine_similarity
  for(i in 1:length(age)){
    a <- as.character(age[i])
    s <- as.character(structure_acronym[i])
    mat[a,s] <- cosine_similarity[i]
  }
  return(mat)
}

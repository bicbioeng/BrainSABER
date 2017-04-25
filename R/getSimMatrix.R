#' Get Age, Structure Acronym, and Similarity Score Matrix
#'
#' This function takes a \code{Biobase ExpressionSet} as returned by
#' \code{getSimScores} and returns a numeric matrix of similarity scores with
#' rows labeled by age, and columns labeled by structure_acronym.
#'
#' @param relevantGenes a Biobase ExpressionSet
#' @param similarity_method currently supported similarity methods are "cosine"
#'     and "euclidean", defaults to "cosine"
#'
#' @return a numeric matrix of similarity scores
#' @import Biobase
#' @export
#' @examples
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes)
#' myGeneSet <- getSimScores(myGenes, myGeneSet, similarity_method = "cosine")
#' myGeneSet <- getSimScores(myGenes, myGeneSet, similarity_method = "euclidean")
#' myCosineMatrix <- getSimMatrix(myGeneSet, similarity_method = "cosine")
#' myEuclideanMatrix <- getSimMatrix(myGeneSet, similarity_method = "euclidean")

getSimMatrix <- function(relevantGenes, similarity_method = "cosine"){
  #get data vectors
  age <- phenoData(relevantGenes)$age
  structure_acronym <- phenoData(relevantGenes)$structure_acronym
  if (similarity_method == "euclidean"){
    sim_score <- phenoData(relevantGenes)$euclidean_similarity
  } else {
    sim_score <- phenoData(relevantGenes)$cosine_similarity
  }
  #get the dimensions
  nrows <- nlevels(age)
  ncols <- nlevels(structure_acronym)
  #create and prefill the matrix, in case not all ages have all structures
  mat <- matrix(data = rep(NA, nrows * ncols), nrow = nrows,
                ncol = ncols)
  rownames(mat) <- as.character(levels(age))
  colnames(mat) <- as.character(levels(structure_acronym))
  checkedMat <- mat
  #populate the matrix with similarity scores
  for(i in 1:length(age)){
    a <- as.character(age[i])
    s <- as.character(structure_acronym[i])
    mat[a,s] <- sim_score[i]
  }
  # the following columns and rows are removed because they have 3 or less
  # values in them
  badcols <- c("CB", "CGE", "DTH", "LGE", "M1C-S1C", "MGE", "Ocx", "PCx", "TCx", "URL")
  badrows <- c("25 pcw", "26 pcw","35 pcw")
  inn <- which(colnames(mat) %in% badcols, arr.ind = TRUE)
  inr <- which(rownames(mat) %in% badrows, arr.ind = TRUE)
  mat <- mat[,-inn]
  mat <- mat[-inr,]
  # the following columns contain less than 5 different
  #keepCols <- NULL
  #for(i in 1:ncols){
  #  if (sum(mat[, i]) > 1){
  #    checkedMat[, i] <- mat[, i]
  #    keepCols <- c(keepCols, colnames(mat)[i])
  #  }
  #}

 # if (!(is.null(keepCols))){
  #  mat <- checkedMat[, keepCols]
  #}
  return(mat)
}

#' Get Age, Structure Acronym, and Cosine Similarity Data Frame
#'
#' This function takes a \code{Biobase ExpressionSet} as returned by
#' \code{getCosineSim} and constructs a \code{data.frame} with columns age,
#' structure_acroynym, and cosine_similarity, sorted by cosine_similarity in
#' decreasing order.
#'
#' @param relevantGenes a Biobase Expression set
#'
#' @return a three-column data.frame
#' @import Biobase
#' @export
#'

getASCdf <- function(relevantGenes){
  #get data vectors
  age <- phenoData(relevantGenes)$age
  structure_acronym <- phenoData(relevantGenes)$structure_acronym
  cosine_similarity <- phenoData(relevantGenes)$cosine_similarity
  #get column numbers to allow comparison to original dataset
  column_num <- phenoData(relevantGenes)$column_num
  #construct data frame, using column_num as rownames
  df <- data.frame(age, structure_acronym, cosine_similarity,
                   row.names = column_num)
  #sort df by cosine_similarity
  df <- df[order(df$cosine_similarity, decreasing = TRUE), ]
  return(df)
}

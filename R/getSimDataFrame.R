#' Get Age, Structure Acronym, and Similarity Scores Data Frame
#'
#' This function takes a \code{Biobase ExpressionSet} as returned by
#' \code{getSimScores} and constructs a \code{data.frame} with columns age,
#' structure_acroynym, and either cosine_similarity or euclidean_similarity,
#' sorted by similarity score in decreasing order.
#'
#' @param relevantGenes a Biobase Expression set
#' @param similarity_method currently supported similarity methods are "cosine"
#'     and "euclidean", defaults to "cosine"
#'
#' @return a three-column data.frame
#' @import Biobase
#' @export
#' @examples
#' \dontrun{myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes)
#' myGeneSet <- getSimScores(myGenes, myGeneSet, similarity_method = "cosine")
#' myGeneSet <- getSimScores(myGenes, myGeneSet, similarity_method = "euclidean")
#' myCosineDF <- getSimDataFrame(myGeneSet, similarity_method = "cosine")
#' myEuclideanDF <- getSimDataFrame(myGeneSet, similarity_method = "euclidean")}

getSimDataFrame <- function(relevantGenes, similarity_method = "cosine"){
  #get column numbers to allow comparison to original dataset
  column_num <- phenoData(relevantGenes)$column_num
  #get data vectors
  age <- phenoData(relevantGenes)$age
  structure_acronym <- phenoData(relevantGenes)$structure_acronym
  if (similarity_method == "euclidean"){
    euclidean_similarity <- phenoData(relevantGenes)$euclidean_similarity
    #construct data frame, using column_num as rownames
    df <- data.frame(age, structure_acronym, euclidean_similarity,
                     row.names = column_num)
    #sort df by euclidean_similarity
    df <- df[order(df$euclidean_similarity, decreasing = TRUE), ]
  } else {
  cosine_similarity <- phenoData(relevantGenes)$cosine_similarity
  #construct data frame, using column_num as rownames
  df <- data.frame(age, structure_acronym, cosine_similarity,
                   row.names = column_num)
  #sort df by cosine_similarity
  df <- df[order(df$cosine_similarity, decreasing = TRUE), ]
  }

  return(df)
}

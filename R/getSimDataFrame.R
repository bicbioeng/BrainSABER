#' Get Age, Structure Acronym, and Similarity Scores Data Frame
#'
#' This function takes a similarity vector or data frame as returned by
#' \code{getSimScores} and a subset of AIBSARNA as returned by
#' \code{getRelevantGenes} and constructs a \code{data.frame} with columns age,
#' structure_acroynym, and either cosine_similarity or euclidean_similarity,
#' sorted by similarity score in decreasing order.  In the case of a similarity
#' data frame, a list of data frames is returned.
#'
#' @param sim_score a vector or data frame of similarity scores
#' @param relevantGenes a Biobase Expression set
#' @param similarity_method currently supported similarity methods are "cosine"
#'     and "euclidean", defaults to "cosine"
#'
#' @return a three-column data.frame or list of data frames
#' @import Biobase
#' @importFrom methods is
#' @export
#' @examples
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes, AIBSARNAid = "gene_symbol")
#' myCosScore <- getSimScores(myGenes, myGeneSet, similarity_method = "cosine")
#' myEucScore <- getSimScores(myGenes, myGeneSet, similarity_method = "euclidean")
#' myCosineDF <- getSimDataFrame(myCosScore, myGeneSet,
#'   similarity_method = "cosine")
#' myEuclideanDF <- getSimDataFrame(myEucScore, myGeneSet,
#'   similarity_method = "euclidean")

getSimDataFrame <- function(sim_score, relevantGenes,
                            similarity_method = "cosine") {
    # sanity checks
    stopifnot(is.vector(sim_score) || is.data.frame(sim_score))
    stopifnot(is(relevantGenes, "ExpressionSet"))
    if (is.vector(sim_score)) {
        #get column numbers to allow comparison to original dataset
        column_num <- phenoData(relevantGenes)$column_num
        #get data vectors
        age <- phenoData(relevantGenes)$age
        structure_acronym <- phenoData(relevantGenes)$structure_acronym
        if (similarity_method == "euclidean") {
            euclidean_similarity <- sim_score
            #construct data frame, using column_num as rownames
            df <- data.frame(age, structure_acronym, euclidean_similarity,
                                row.names = column_num)
            #sort df by euclidean_similarity
            df <- df[order(df$euclidean_similarity, decreasing = TRUE),]
        } else {
            cosine_similarity <- sim_score
            #construct data frame, using column_num as rownames
            df <- data.frame(age, structure_acronym, cosine_similarity,
                                row.names = column_num)
            #sort df by cosine_similarity
            df <- df[order(df$cosine_similarity, decreasing = TRUE),]
        }
        return(df)
    } else {
        dfList <- as.list(apply(X = sim_score, MARGIN = 2,
                                FUN = getSimDataFrame,
                                relevantGenes = relevantGenes,
                                similarity_method = similarity_method))
        return(dfList)
    }
}

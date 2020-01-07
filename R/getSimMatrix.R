#' Get Age, Structure Acronym, and Similarity Score Matrix
#'
#' This function takes a similarity vector or data frame as returned by
#' \code{getSimScores} and returns a numeric matrix or list of matrices of
#' similarity scores with rows labeled by age, and columns labeled by
#' structure_acronym.
#'
#' @param sim_score a vector or data frame of similarity scores
#' @param relevantGenes a Biobase Expression set created using the 
#'     \code{getRelevantGenes()} function
#'
#' @return a numeric matrix of similarity scores, or a list of matrices
#' @import Biobase
#' @importFrom methods is
#' @export
#' @examples
#' AIBSARNA <- buildAIBSARNA(mini = TRUE)
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TSPAN6", "DPM1", "C1orf112")
#' myGeneSet <- getRelevantGenes(myGenes, AIBSARNA = AIBSARNA, 
#'     AIBSARNAid = "gene_symbol")
#' myCosScore <- getSimScores(myGenes, myGeneSet, similarity_method = "cosine")
#' myEucScore <- getSimScores(myGenes, myGeneSet, 
#'     similarity_method = "euclidean")
#' myCosineMatrix <- getSimMatrix(myCosScore, myGeneSet)
#' myEuclideanMatrix <- getSimMatrix(myEucScore, myGeneSet)

getSimMatrix <- function(sim_score, relevantGenes) {
    # sanity checks
    stopifnot(is.vector(sim_score) || is.data.frame(sim_score))
    stopifnot(is(relevantGenes, "ExpressionSet"))
    if (is.vector(sim_score)) {
        #get data vectors
        age <- phenoData(relevantGenes)$age
        structure_acronym <- phenoData(relevantGenes)$structure_acronym
        #get the dimensions
        nrows <- nlevels(age)
        ncols <- nlevels(structure_acronym)
        #create and prefill the matrix, in case not all ages have all structures
        mat <- matrix(data=rep(NA, nrows * ncols), nrow=nrows, ncol=ncols)
        rownames(mat) <- as.character(levels(age))
        colnames(mat) <- as.character(levels(structure_acronym))
        #populate the matrix with similarity scores
        for (i in seq_along(age)) {
            a <- as.character(age[i])
            s <- as.character(structure_acronym[i])
            mat[a, s] <- sim_score[i]
        }
        if(length(age) > 5){
            # the following columns and rows are removed because they have 
            # 3 or less values in them, crashing heatmaps later on
            badcols <- c("CB", "CGE", "DTH", "LGE", "M1C-S1C", "MGE", 
                        "Ocx", "PCx", "TCx", "URL")
            badrows <- c("25 pcw", "26 pcw", "35 pcw")
            inn <- which(colnames(mat) %in% badcols, arr.ind = TRUE)
            inr <- which(rownames(mat) %in% badrows, arr.ind = TRUE)
            mat <- mat[, -inn]
            mat <- mat[-inr, ]
        }
        return(mat)
    } else {
        nSamples <- length(sim_score[1,])
        # initialize matList
        matList <- vector("list", nSamples)
        for(x in seq_len(nSamples)){
            matList[[x]] <- getSimMatrix(sim_score[, x], relevantGenes)
        }
        return(matList)
    }
}

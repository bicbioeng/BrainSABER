#' Get Age, Structure Acronym, and Similarity Score Matrix
#'
#' This function takes a similarity vector or data frame as returned by
#' \code{getSimScores} and either a CellScabbard object or a subset of AIBSARNA
#' as returned by \code{getRelevantGenes}. Returns a numeric matrix or list
#' of matrices of similarity scores with rows labeled by age, and columns
#' labeled by structure_acronym.
#'
#' @param data a CellScabbard object with non-empty relevantGenes slot
#' @param sim_score a vector or data frame of similarity scores
#' @param relevantGenes a SummarizedExperiment created using the
#'     \code{getRelevantGenes()} function
#'
#' @return a numeric matrix of similarity scores, or a list of matrices
#' @import SummarizedExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom methods is
#' @export
#' @examples
#' AIBSARNA <- buildAIBSARNA(mini = TRUE)
#' # Example 1 - using CellScabbard class
#' # get a random sample of 3 genes
#' totalGenes <- nrow(AIBSARNA)
#' gene_idx <- sample.int(totalGenes, 3)
#' sample_idx <- c(1,3,5)
#' # Subset AIBSARNA
#' exprs <- assay(AIBSARNA)[gene_idx, sample_idx]
#' fd <- rowData(AIBSARNA)[gene_idx, ]
#' pd <- colData(AIBSARNA)[sample_idx, ]
#' # build a trimmed data set
#' myGenes <- CellScabbard(exprsData = exprs, phenoData = pd, featureData = fd,
#'                         AIBSARNA = AIBSARNA, autoTrim = TRUE)
#' similarityScores(myGenes) <- getSimScores(data = myGenes,
#'                                           similarity_method = "cosine")
#' similarityMatrices(myGenes) <- getSimMatrix(data = myGenes)
#' similarityMatrices(myGenes)
#' similarityScores(myGenes) <- getSimScores(data = myGenes,
#'                                           similarity_method = "euclidean")
#' similarityMatrices(myGenes) <- getSimMatrix(data = myGenes)
#' similarityMatrices(myGenes)
#'
#' # Example 2 - manual gene selection and relevant gene extraction
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TSPAN6", "DPM1", "C1orf112")
#' myGeneSet <- getRelevantGenes(myGenes, AIBSARNA = AIBSARNA,
#'     AIBSARNAid = "gene_symbol")
#' myCosScore <- getSimScores(myGenes, myGeneSet, similarity_method = "cosine")
#' myEucScore <- getSimScores(myGenes, myGeneSet,
#'                            similarity_method = "euclidean")
#' myCosineMatrix <- getSimMatrix(sim_score = myCosScore,
#'                                relevantGenes = myGeneSet)
#' myEuclideanMatrix <- getSimMatrix(sim_score = myEucScore,
#'                                   relevantGenes = myGeneSet)

getSimMatrix <- function(data = NULL, sim_score = NULL, relevantGenes = NULL) {
    # extract relevant genes from cellscabbard object
    if(is(data, "CellScabbard")){
      # extract sim scores and relevant genes from data
      sim_score <- similarityScores(data)
      relevantGenes <- relevantGenes(data)
    } else {
      # raise error if relevant genes not provided
      if(!is(relevantGenes,"SummarizedExperiment")){
        em <-
          "Relevant genes required and must be a SummarizedExperiment
          (or subclass of) if not using a CellScabbard data set."
        stop(em)
      }
    }
    # check for proper input
    stopifnot(is.vector(sim_score) || is.data.frame(sim_score))
    if (is.vector(sim_score)) {
        #get data vectors
        age <- colData(relevantGenes)$age
        structure_acronym <- factor(colData(relevantGenes)$structure_acronym)
        
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
            matList[[x]] <- getSimMatrix(sim_score = sim_score[, x],
                                         relevantGenes = relevantGenes)
        }
        matList <- SimpleList(matList)
        return(matList)
    }
}

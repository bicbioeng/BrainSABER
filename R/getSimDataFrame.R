#' Get Age, Structure Acronym, and Similarity Scores Data Frame
#'
#' This function takes in a CellScabbard object, or both a similarity vector
#' or data frame as returned by \code{getSimScores} and a subset of AIBSARNA
#' as returned by \code{getRelevantGenes}. Constructs a \code{data.frame}
#' with columns age, structure_acroynym, and either cosine_similarity or
#' euclidean_similarity, sorted by similarity score in decreasing order.
#' In the case of a similarity data frame, a list of data frames is returned.
#'
#' @param data a CellScabbard object with non-empty relevantGenes and
#' similarityScores slots, or a SummarizedExperiment created using the
#'     \code{getRelevantGenes()} function
#' @param sim_score a vector or data frame of similarity scores
#' @param relevantGenes a SummarizedExperiment object created using the
#'     \code{getRelevantGenes()} function
#' @param similarity_method currently supported similarity methods are "cosine"
#'     and "euclidean", defaults to "cosine"
#'
#' @return a three-column data.frame or list of data frames
#' @import SummarizedExperiment
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
#' # cosine similarity method
#' similarityScores(myGenes) <- getSimScores(data = myGenes,
#'                                           similarity_method = "cosine")
#' similarityDFs(myGenes) <- getSimDataFrame(data = myGenes,
#'                                           similarity_method = "cosine")
#' similarityDFs(myGenes)
#' # euclidean similarity method
#' similarityScores(myGenes) <- getSimScores(data = myGenes,
#'                                           similarity_method = "euclidean")
#' similarityDFs(myGenes) <- getSimDataFrame(data = myGenes,
#'                                           similarity_method = "euclidean")
#' similarityDFs(myGenes)
#'
#' # Example 2 - manual gene selection and relevant gene extraction
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TSPAN6", "DPM1", "C1orf112")
#' myGeneSet <- getRelevantGenes(myGenes, AIBSARNA = AIBSARNA,
#'     AIBSARNAid = "gene_symbol")
#' myCosScore <- getSimScores(myGenes, relevantGenes = myGeneSet,
#'                            similarity_method = "cosine")
#' myEucScore <- getSimScores(myGenes, relevantGenes = myGeneSet,
#'     similarity_method = "euclidean")
#' myCosineDF <- getSimDataFrame(sim_score = myCosScore,
#'                               relevantGenes = myGeneSet,
#'                               similarity_method = "cosine")
#' myEuclideanDF <- getSimDataFrame(sim_score = myEucScore,
#'                                  relevantGenes = myGeneSet,
#'                                  similarity_method = "euclidean")

getSimDataFrame <- function(data = NULL, sim_score = NULL, relevantGenes = NULL,
                            similarity_method = "cosine") {
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
        #get column numbers to allow comparison to original data
        column_num <- colData(relevantGenes)$column_num
        #get data vectors
        age <- colData(relevantGenes)$age
        structure_acronym <- colData(relevantGenes)$structure_acronym
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
                                FUN = function(x){ getSimDataFrame(
                                  sim_score = x,
                                  relevantGenes = relevantGenes,
                                  similarity_method = similarity_method)
                                }))
        return(dfList)
    }
}

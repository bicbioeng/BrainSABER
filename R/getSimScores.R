#' Get Similarity Scoring for a Gene Expression Vector
#'
#' This function computes the similarity score of a gene expression
#' vector returned by \code{getExternalVector} or a trimmed data set returned
#' by \code{getTrimmedExternalSet}, compared to a subset of AIBSARNA,
#' obtained by \code{getRelevantGenes}.
#'
#' @param data a named vector of gene expression values returned by
#'     \code{getExternalVector}, a SummarizedExperiment returned by
#'     \code{getTrimmedExternalSet}, or a CellScabbard object.
#' @param relevantGenes a SummarizedExperiment object created using the
#'     \code{getRelevantGenes()} function
#' @param similarity_method currently supported similarity methods are
#'     "cosine" and "euclidean", defaults to "cosine"
#'
#' @return If \code{data} is a vector, returns a vector of similarity scores
#'     for each sample in relevantGenes. If \code{data} is a
#'     SummarizedExperiment, returns a data frame, with columns containing  the
#'     similarity scores for and named after each sample in \code{data}, and
#'     rows named after each sample in \code{relevantGenes}. If \code{data} is
#'     a CellScabbard, the results will be stored in its similarityScores slot.
#'
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
#' similarityScores(myGenes) <- getSimScores(data = myGenes,
#'                                           similarity_method = "cosine")
#' similarityScores(myGenes)
#' similarityScores(myGenes) <- getSimScores(data = myGenes,
#'                                           similarity_method = "euclidean")
#' similarityScores(myGenes)
#'
#' # Example 2 - manual gene selection and relevant gene extraction
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TSPAN6", "DPM1", "C1orf112")
#' myGeneSet <- getRelevantGenes(myGenes, AIBSARNA = AIBSARNA,
#'     AIBSARNAid = "gene_symbol")
#' CosScores <- getSimScores(myGenes, myGeneSet,
#'     similarity_method = "cosine")
#' EucScores <- getSimScores(myGenes, myGeneSet,
#'     similarity_method = "euclidean")

getSimScores <-
    function(data, relevantGenes = NULL, similarity_method = "cosine") {
        # extract relevant genes from cellscabbard object
        if(is(data, "CellScabbard")){
          relevantGenes <- relevantGenes(data)
        } else {
          # raise error if relevant genes not provided
          if(!is(relevantGenes,"SummarizedExperiment")){
            em <-
              "Relevant genes required and must be a SummarizedExperiment
              (or subclass of) if not using a CellScabbard object."
            stop(em)
          }
        }
        if (is.vector(data)) {
            if (similarity_method == "euclidean") {
                #euclidean_similarity <- as.vector(apply(X=exprs(relevantGenes),
                euclidean_similarity <- as.vector(apply(X=assay(relevantGenes),
                    MARGIN = 2, FUN = euclideanDist, y = data))
                return(euclidean_similarity)
            } else {
                # cosine_similarity <- as.vector(apply(X = exprs(relevantGenes),
                cosine_similarity <- as.vector(apply(X = assay(relevantGenes),
                    MARGIN = 2, FUN = lsa::cosine,y = data))
                return(cosine_similarity)
            }
        } else {
            # simFrame <- as.data.frame(apply(X = exprs(data), MARGIN = 2,
            simFrame <- as.data.frame(apply(X = assay(data), MARGIN = 2,
                FUN = getSimScores, relevantGenes = relevantGenes,
                similarity_method = similarity_method))
            return(simFrame)
        }
    }

#helper function
euclideanDist <- function(x, y) {
    tmp <- sqrt(sum((x - y) ^ 2))
    return(1 / (1 + tmp))
}

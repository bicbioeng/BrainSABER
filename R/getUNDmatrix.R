#' getUNDmatrix
#'
#' This function returns a matrix showing whether gene expression values in
#'  \code{dataSet} are up-regulated, down-regulated, or normal.
#'  \code{method = "discrete"} will function on any CellScabbard object, while
#'  \code{method = "log2FC"} requires a trimmed data set as returned by
#'  \code{getTrimmedExternalSet} and a matching subset of AIBSARNA as
#'  returned by \code{getRelevantGenes}. Results are stored in the 'UNDmatrices'
#'  slot of the \code{dataSet} if it's a CellScabbard object.
#'
#' @param dataSet a CellScabbard or SummarizedExperiment object
#' @param relevantGenes (optional) a SummarizedExperiment and subset of AIBSARNA
#' @param method \code{"discrete"} applies thresholds directly to expression
#'     data. \code{"log2FC"} applies thresholds to the log2 fold-change
#'     between the expression data of each sample from \code{dataSet} and
#'     \code{relevantGenes}.
#' @param up_threshold a numerical value defining the lower bound
#'     (inclusive) by which to consider a gene up-regulated, defaults to 0.5
#' @param down_threshold a numerical value defining the upper bound
#'     (inclusive) by which to consider a gene down-regulated, defaults
#'     to -0.5
#' @param matrix_type either \code{"num"} for a numerical matrix with -1
#'     indicating down-regulation, 1 indicating up-regulation, and 0
#'     indicating normal, or \code{"char"} for a character matrix with "D"
#'     indicating down-regulation, "U" indicating up-regulation, and
#'     "N" indicating normal
#'
#' @return a list containing as many numerical or character matrices as samples
#'     in \code{dataSet}, with each matrix having as many rows as genes in
#'     \code{dataSet} and as many columns as samples in \code{relevantGenes}
#' @export
#' @import SummarizedExperiment
#' @importFrom methods is
#' @importFrom S4Vectors SimpleList
#'
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
#' UNDmatrices(myGenes) <- getUNDmatrix(myGenes, method = "discrete",
#'                                      up_threshold = 3,
#' down_threshold = 1, matrix_type = "char")
#' UNDmatrices(myGenes)
#' UNDmatrices(myGenes) <- getUNDmatrix(myGenes, method = "log2FC",
#'                                      up_threshold = 3,
#' down_threshold = 1, matrix_type = "num")
#' UNDmatrices(myGenes)
#'
#' # Example 2 - manual gene selection and relevant gene extraction
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TSPAN6", "DPM1", "C1orf112")
#' myGeneSet <- getRelevantGenes(myGenes, AIBSARNA = AIBSARNA,
#'     AIBSARNAid = "gene_symbol")
#' myTrimmedGeneSet <- getTrimmedExternalSet(myGeneSet,
#'     dataSetId = "gene_symbol", AIBSARNA, AIBSARNAid = "gene_symbol")
#' myUNDnumericalMatrix <- getUNDmatrix(myTrimmedGeneSet, method = "discrete",
#'     up_threshold = 3, down_threshold = 1, matrix_type = "num")
#' myUNDcharacterMatrix <- getUNDmatrix(myTrimmedGeneSet, myGeneSet,
#'                                      method = "log2FC",
#'     up_threshold = 3, down_threshold = 1, matrix_type = "char")

getUNDmatrix <- function(dataSet, relevantGenes = NULL,
        method=c("discrete", "log2FC"), up_threshold=0.5, down_threshold=-0.5,
        matrix_type = c("num", "char")) {
    # set up und vector for type of matrix
    if (matrix_type == "char") {
        und <- c("U", "N", "D")
    } else {
        und <- c(-1, 0, 1)
    }
    if (method == "discrete") {
        # get expression matrix from dataSet
        # exprs <- Biobase::exprs(dataSet)
        exprs <- assay(dataSet)
        # create empty matrix
        nrows <- length(exprs[, 1])
        ncols <- length(exprs[1,])
        mat <- matrix(data=rep(und[2], nrows * ncols), nrow=nrows, ncol=ncols)
        mat <-as.matrix(apply(X = exprs, MARGIN = c(1, 2), FUN = classifyValue,
                                und=und, u=up_threshold, d=down_threshold))
        rownames(mat) <- rownames(exprs)
        colnames(mat) <- colnames(exprs)
        matList <- SimpleList(mat)
        return(matList)
    } else {
        # extract relevant genes from cellscabbard object
        if(is(dataSet, "CellScabbard") ){
          relevantGenes <- relevantGenes(dataSet)
        } else {
          # raise error if relevant genes not provided
          if(!is(relevantGenes,"SummarizedExperiment")){
            em <-
              "Relevant genes required and must be a SummarizedExperiment
              (or subclass of) if not using a CellScabbard object."
            stop(em)
          }
        }
        dataExprs <- assay(dataSet)
        relExprs <- assay(relevantGenes)
        nSamples <- length(dataExprs[1,])
        # initialize matList
        matList <- lapply(seq_len(nSamples), FUN = function(z){
            ret <- as.matrix(apply(X=relExprs, MARGIN=2,
                                FUN=function(x, v){log2(v) - log2(x)},
                                v = dataExprs[, z]))
            rownames(ret) <- rownames(dataExprs)
            colnames(ret) <- colnames(relExprs)
            return(ret)
        })
        matList <- SimpleList(matList)
        return(matList)
    }
}

#helper function
classifyValue <- function(x, und, u, d) {
    if (is.finite(x) && x >= u) {
        return(und[3])
    } else if (is.finite(x) && x <= d) {
        return(und[1])
    } else {
        return(und[2])
    }
}

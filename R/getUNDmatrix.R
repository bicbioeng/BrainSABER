#' getUNDmatrix
#'
#' This function returns a matrix showing whether gene expression values in
#'  \code{dataSet} are up-regulated, down-regulated, or normal.
#'  \code{method = "discrete"} will function on any ExpressionSet, while
#'  \code{method = "log2FC"} requires a trimmed data set as returned by
#'  \code{getTrimmedExternalSet} and a matching subset of AIBSARNA as
#'  returned by \code{getRelevantGenes}.
#'
#' @param dataSet a Biobase ExpressionSet
#' @param relevantGenes (optional) an ExpressionSet that is a subset of AIBSARNA
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
#' @return for \code{method = "discrete"} either a numerical or character
#'     matrix with as many rows as genes in \code{dataSet} and as many
#'     columns as samples in \code{dataSet} or, for \code{method = "log2FC"},
#'     a list containing as many numerical or character matrices as samples
#'     in \code{dataSet}, with each matrix having as many rows as genes in
#'     \code{dataSet} and as many columns as samples in \code{relevantGenes}
#' @export
#'
#' @examples
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TNFRSF1A", "BCL3", "NEFH")
#' myGeneSet <- getRelevantGenes(myGenes, AIBSARNAid = "gene_symbol")
#' myUNDnumericalMatrix <- getUNDmatrix(myGeneSet, method = "discrete",
#'     up_threshold = 3, down_threshold = 1, matrix_type = "num")
#' myUNDcharacterMatrix <- getUNDmatrix(myGeneSet, method = "discrete",
#'     up_threshold = 3, down_threshold = 1, matrix_type = "char")

getUNDmatrix <- function(dataSet, relevantGenes = NULL,
                         method = c("discrete", "log2FC"),
                         up_threshold = 0.5,
                         down_threshold = -0.5,
                         matrix_type = c("num", "char")){
  # set up und vector for type of matrix
  if (matrix_type == "char"){
    und <- c("U", "N", "D")
  } else {
    und <- c(-1, 0, 1)
  }
  if (method == "discrete"){
    # get expression matrix from dataSet
    exprs <- Biobase::exprs(dataSet)
    # create empty matrix
    nrows <- length(exprs[ , 1])
    ncols <- length(exprs[1, ])
    mat <- matrix(data = rep(und[2], nrows * ncols), nrow = nrows,
                  ncol = ncols)
    mat <- as.matrix(apply(X = exprs, MARGIN = c(1,2), FUN = classifyValue,
                           und = und, u = up_threshold, d = down_threshold))
    rownames(mat) <- rownames(exprs)
    colnames(mat) <- colnames(exprs)
    return(mat)
  } else {
    dataExprs <- exprs(dataSet)
    relExprs <- exprs(relevantGenes)
    nSamples <- length(dataExprs[1, ])
    # initialize matList
    matList <- vector("list", nSamples)
    for(i in 1:nSamples){
      # populate matrix with log2 fold change against each sample
      # of relExprs
      matList[[i]] <- as.matrix(apply(X = relExprs, MARGIN = 2,
                                      FUN = function(x, v){
                                        log2(v) - log2(x)
                                      }, v = dataExprs[, i]))
      # convert log2 fold change matrix to UND
      matList[[i]] <- as.matrix(apply(X = matList[[i]], MARGIN = c(1,2),
                                      FUN = classifyValue, und = und,
                                      u = up_threshold, d = down_threshold))
      rownames(matList[[i]]) <- rownames(dataExprs)
      colnames(matList[[i]]) <- colnames(relExprs)
    }
    return(matList)
  }
}

#helper function
classifyValue <- function(x, und, u, d){
  if (is.finite(x) && x >= u){
    return(und[3])
  } else if (is.finite(x) && x <= d){
    return(und[1])
  } else {
    return(und[2])
  }
}

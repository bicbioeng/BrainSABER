#' getUNDmatrix
#'
#' This function returns a matrix showing whether gene expression values in
#'  \code{dataSet} are up-regulated, down-regulated, or normal.
#'  \code{method = "discrete"} will function on any ExpressionSet, while
#'  \code{method = "log2FC"} requires a trimmed data set as returned by
#'  \code{getTrimmedExternalSet} and a matching subset of AIBSARNA as returned by
#'  \code{getRelevantGenes}.
#'
#' @param dataSet a Biobase ExpressionSet
#' @param relevantGenes (optional) an ExpressionSet that is a subset of AIBSARNA
#' @param method "discrete" applies thresholds directly to expression data.
#'     "log2FC" applies thresholds to the log2 fold-change between the expression
#'     data from dataSet and AIBSARNA.
#' @param up_threshold a numerical value defining the lower bound (inclusive) by
#'     which to consider a gene up-regulated, defaults to 0.5
#' @param down_threshold a numerical value defining the lower bound (inclusive)
#'     by which to consider a gene down-regulated, defaults to -0.5
#' @param matrix_type either "num" for a numerical matrix with -1 indicating
#'     down-regulation, 1 indicating up-regulation, and 0 indicating normal, or
#'     "char" for a character matrix with "D" indicating down-regulation, "U"
#'     indicating up-regulation, and "N" indicating normal
#'
#' @return either a numerical or character matrix
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
                         method = c("discrete", "log2FC"), up_threshold = 0.5,
                         down_threshold = -0.5, matrix_type = c("num", "char")){
  # get expression matrix from dataSet
  exprs <- Biobase::exprs(dataSet)
  # set up und vector for type of matrix
  if (matrix_type == "char"){
    und <- c("U", "N", "D")
  } else {
    und <- c(-1, 0, 1)
  }
  # create empty matrix
  nrows <- length(exprs[ , 1])
  ncols <- length(exprs[1, ])
  mat <- matrix(data = rep(und[2], nrows * ncols), nrow = nrows,
                ncol = ncols)
  if (method == "discrete"){
    mat <- as.matrix(apply(X = exprs, MARGIN = c(1,2), FUN = classifyValue,
                           und = und, u = up_threshold, d = down_threshold))

    #for(r in 1:nrows){
     # for(c in 1:ncols){
      #  if (exprs[r, c] >= up_threshold){
       #   mat[r, c] <- und[3] # 1 or "U", depending on matrix type
      #  } else if (exprs[r, c] <= down_threshold)
       #   mat[r, c] <- und[1] # -1 or "D", depending on matrix type
      #}
    #}
  } else {
    mat <- apply(X = exprs, MARGIN = c(1,2),
                 FUN = function(x, und, u, d, rgexprs){
                   if (x >= u){
                     return(und[3])
                   } else if (x <= d){
                     return(und[1])
                   } else {
                     return(und[2])
                   }
                 }, und = und, u = up_threshold, d = down_threshold)
    for(r in 1:nrows){
      for(c in 1:ncols){
        val <- log2(exprs[r, c]) - log2(exprs(relevantGenes)[r, c])
        if (is.finite(val) && val >= up_threshold){
          mat[r, c] <- und[3] # 1 or "U", depending on matrix type
        } else if (is.finite(val) && val <= down_threshold)
          mat[r, c] <- und[1] # -1 or "D", depending on matrix type
      }
    }
  }
  return(mat)
}

#helper function
classifyValue <- function(x, und, u, d){
  if (x >= u){
    return(und[3])
  } else if (x <= d){
    return(und[1])
  } else {
    return(und[2])
  }
}

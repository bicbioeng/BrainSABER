#' getUNDmatrix
#'
#' This function returns a matrix showing whether gene expression values in
#'  \code{dataSet} are up-regulated, down-regulated, or normal. \code{method = "discrete"} requires \code{up_threshold} and \code{down_threshold} to be defined.
#'
#' @param dataSet a Biobase ExpressionSet
#' @param method currently only "discrete" is available
#' @param up_threshold a numerical value defining the lower bound (inclusive) by
#'     which to consider a gene up-regulated
#' @param down_threshold a numerical value defining the lower bound (inclusive)
#'     by which to consider a gene down-regulated
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
#' myGeneSet <- getRelevantGenes(myGenes, gene_names = "HGNC")
#' myUNDnumericalMatrix <- getUNDmatrix(myGeneSet, method = "discrete",
#'     up_threshold = 3, down_threshold = 1, matrix_type = "num")
#' myUNDcharacterMatrix <- getUNDmatrix(myGeneSet, method = "discrete",
#'     up_threshold = 3, down_threshold = 1, matrix_type = "char")

getUNDmatrix <- function(dataSet, method = "discrete", up_threshold,
                         down_threshold, matrix_type = c("num", "char")){
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
    for(r in 1:nrows){
      for(c in 1:ncols){
        if (exprs[r, c] >= up_threshold){
          mat[r, c] <- und[3] # 1 or "U", depending on matrix type
        } else if (exprs[r, c] <= down_threshold)
          mat[r, c] <- und[1] # -1 or "D", depending on matrix type
      }
    }
  }
  return(mat)
}

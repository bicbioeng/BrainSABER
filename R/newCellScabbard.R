#' Creates a new CellScabbard object.
#'
#' @param exprsData expression data matrix for an experiment
#' @param phenoData an AnnotatedDataFrame containing attributes of individual
#'     cells
#' @param featureData an AnnotatedDataFrame containing attributes of features
#'     (e.g. genes)
#' @return a new CellScabbard object
#' @importFrom Biobase annotatedDataFrameFrom assayDataNew
#' @export
#'
newCellScabbard <- function( exprsData,
                            phenoData = NULL,
                            featureData = NULL)
{
  if (class(exprsData) != "matrix"){
    stop("Error: argument exprsData must be a matrix")
  }

  if(is.null(phenoData))
    phenoData <- annotatedDataFrameFrom(exprsData, byrow=FALSE)
  if(is.null(featureData))
    featureData <- annotatedDataFrameFrom(exprsData, byrow=TRUE)

  cs <- new( "CellScabbard",
              assayData = assayDataNew("environment", exprs=exprsData),
              phenoData = phenoData,
              featureData = featureData)
  cs
}

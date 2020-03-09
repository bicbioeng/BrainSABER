#' Methods for the CellScabbard class
#'
#' These methods operate on CellScabbard objects.  They are used to access the
#' results of the \pkg{BrainSABER} workflow stored within a CellScabbard.
#'
#' @param cs A CellScabbard object
#' @param value data type, any of matrix, data.frame,list, or SimpleList
#' @return The contents of a slot of the CellScabbard object
#' @importFrom methods as is new
#' @import SummarizedExperiment
#' @examples
#' # construct example data set
#' AIBSARNA <- buildAIBSARNA(mini = TRUE)
#'
#' # get a random sample of 3 genes
#' totalGenes <- nrow(AIBSARNA)
#' gene_idx <- sample.int(totalGenes, 3)
#' sample_idx <- c(1,3,5)
#' 
#' # Subset AIBSARNA
#' exprs <- assay(AIBSARNA)[gene_idx, sample_idx]
#' fd <- rowData(AIBSARNA)[gene_idx, ]
#' pd <- colData(AIBSARNA)[sample_idx, ]
#' 
#' # construct a CellScabbard data set
#' myGenes <- CellScabbard(exprsData = exprs, phenoData = pd, featureData = fd, 
#'                         AIBSARNA = AIBSARNA, autoTrim = TRUE)
#' relevantGenes(myGenes)
#' 
#' # the following fields will be empty as output must be assigned to 
#' # them first
#' similarityScores(myGenes)
#' similarityMatrices(myGenes)
#' similarityDFs(myGenes)
#' UNDmatrices(myGenes)
#' @name CellScabbard-methods
NULL
#'
#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
dataSetId <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@dataSetId
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
`dataSetId<-` <- function(cs, value){
  stopifnot( is( cs, "CellScabbard" ) )
  stopifnot(is(value, "character") || length(value) == 1)
  cs@dataSetId <- value
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
AIBSARNAid <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@AIBSARNAid
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
`AIBSARNAid<-` <- function(cs, value){
  stopifnot( is( cs, "CellScabbard" ) )
  stopifnot(is(value, "character") || length(value) == 1)
  cs@AIBSARNAid <- value
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
relevantGenes <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@relevantGenes
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
`relevantGenes<-` <- function(cs, value){
  stopifnot( is( cs, "CellScabbard" ) )
  stopifnot(is(value, "SummarizedExperiment"))
  cs@relevantGenes <- value
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
similarityScores <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@similarityScores
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
`similarityScores<-` <- function(cs, value){
  stopifnot( is( cs, "CellScabbard" ) )
  stopifnot(is(value, "data.frame"))
  cs@similarityScores <- value
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
similarityDFs <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@similarityDFs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
`similarityDFs<-` <- function(cs, value){
  stopifnot( is( cs, "CellScabbard" ) )
  stopifnot(is(value, "list"))
  cs@similarityDFs <- value
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
similarityMatrices <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@similarityMatrices
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
`similarityMatrices<-` <- function(cs, value){
  stopifnot( is( cs, "CellScabbard" ) )
  stopifnot(is(value, "SimpleList"))
  cs@similarityMatrices <- value
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
UNDmatrices <- function( cs ) {
  stopifnot(is( cs, "CellScabbard" ) )
  cs@UNDmatrices
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
`UNDmatrices<-` <- function(cs, value){
  stopifnot(is( cs, "CellScabbard" ))
  stopifnot(is(value, "SimpleList"))
  cs@UNDmatrices <- value
  cs
}
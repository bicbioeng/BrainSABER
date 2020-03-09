#' Methods for the CellScabbard class
#'
<<<<<<< HEAD
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
=======
#' These methods operate on CellScabbard objects. Please note that
#' treeList<- and originalTrees<- are not intended to be
#' called directly.
#'
#' @param cs A CellScabbard object
#' @param value
#'    \code{cellTreeInfo(cs)<-}: A character vector of the name of the
#'       column of \code{phenoData()} to use for grouping data in the
#'       \pkg{cellTreeGenerator} workflow.
#'
#'    \code{monocleInfo(cs)<-}: A charater vector of parameters used by
#'       \code{generate_tree(treeType = "monocle")} in the
#'       \pkg{cellTreeGenerator} workflow.
#'
#'    \code{TSCANinfo(cs)<-}: A character vector of the row name of a
#'       single gene in \code{exprs()} to use for a single gene vs.
#'       pseudotime plot for \code{generate_tree(treeType = "TSCAN")} in the
#'       \pkg{cellTreeGenerator} workflow.
#'
#'    \code{sincellInfo(cs)<-}: The value to use as a named parameter for
#'       sincell, used by \code{generate_tree(treeType = "sincell")} in the
#'       \pkg{cellTreeGenerator} workflow
#'
#'    \code{oncoNEMdata(cs)<-}: A matrix containing the data to be used
#'       by \code{generate_tree(treeType = "oncoNEM")} in the
#'       \pkg{cellTreeGenerator} workflow.
#'
#'    \code{CanopyData(cs)<-}: A list containing the data to be used by
#'       \code{generate_tree(treeType = "Canopy")} in the
#'       \pkg{cellTreeGenerator} workflow.
#'
#'    \code{cellscapeData(cs, pt)<-}: The value to use as a named parameter for
#'       cellscape, used by \code{generate_tree(treeType = "cellscape")} in the
#'       \pkg{cellTreeGenerator} workflow
#'
#'    \code{timescapeData(cs, pt)<-}: The value to use as a named parameter for
#'       timescape used by \code{generate_tree(treeType = "timescape")} in the
#'       \pkg{cellTreeGenerator} workflow
#'
#'    \code{mapscapeData(cs)<-}: The value to use as a named parameter for
#'       mapscape, used by \code{generate_tree(treeType = "mapscape")} in the
#'       \pkg{cellTreeGenerator} workflow
#'
#' @param tt The type of tree being stored
#' @param pt The name of the \pkg{sincell}, \pkg{cellscape}, \pkg{timescape}, or
#'     \pkg{mapscape} parameter to store
#' @return An updated CellScabbard object, or the contents of a slot of the
#'      CellScabbard object
#' @importFrom methods as is new
#' @import Biobase
>>>>>>> refs/remotes/origin/master
#' @name CellScabbard-methods
NULL
#'
#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
<<<<<<< HEAD
dataSetId <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@dataSetId
=======
cellTreeInfo <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@cellTreeInfo
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
`cellTreeInfo<-` <- function(cs, value){
  stopifnot( is( cs, "CellScabbard" ) )
  if(length(value) != 1 || class(value) != "character" ||
     ! value %in% colnames(pData(cs))){
    stop("cellTreeInfo must be a character vector of length 1, and a column
         name in phenoData().")
  }
  cs@cellTreeInfo <- value
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
monocleInfo <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@monocleInfo
>>>>>>> refs/remotes/origin/master
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
<<<<<<< HEAD
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
=======
`monocleInfo<-` <- function(cs, value) {
  stopifnot( is( cs, "CellScabbard" ) )
  ex_type <- c("UMI", "TC", "FPKM", "TPM", "LTFPKM", "LTTPM")
  if(length(value) != 4 || class(value) != "character" ||
    !(value[1] %in% colnames(fData(cs))) ||
    !(value[2] %in% as.character(fData(cs)[[value[1]]])) ||
    !(value[3] %in% as.character(fData(cs)[[value[1]]])) ||
    !(value[4] %in% ex_type)){
      stop("monocleInfo must be a character vector of length 4.
           monocleInfo[1] must be a column name in featureData(),  and
           monocleInfo[2] and monocleInfo[3] must be gene identifiers found
           in that column. monocleInfo[4] must be one of 'UMI', 'TC',
           'FPKM', 'TPM', 'LTFPKM', 'LTTPM'")
  }
  cs@monocleInfo <- value
  cs
}


#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export

TSCANinfo <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@TSCANinfo
>>>>>>> refs/remotes/origin/master
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
<<<<<<< HEAD
`AIBSARNAid<-` <- function(cs, value){
  stopifnot( is( cs, "CellScabbard" ) )
  stopifnot(is(value, "character") || length(value) == 1)
  cs@AIBSARNAid <- value
=======
`TSCANinfo<-` <- function(cs, value) {
  stopifnot( is( cs, "CellScabbard" ) )
  if(length(value) != 1 || class(value) != "character" ||
     !(value %in% row.names(exprs(cs)))){
    stop("TSCANinfo must be a character vector of length 1, and a row name
         in exprs().")
  }
  cs@TSCANinfo <- value
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export

sincellInfo <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@sincellInfo
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @keywords internal
#' @export

`sincellInfo<-` <- function(cs, pt, value) {
  stopifnot(is(cs, "CellScabbard"))
  cs@sincellInfo[[pt]] <- value
>>>>>>> refs/remotes/origin/master
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
<<<<<<< HEAD
relevantGenes <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@relevantGenes
=======

oncoNEMdata <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@oncoNEMdata
>>>>>>> refs/remotes/origin/master
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
<<<<<<< HEAD
`relevantGenes<-` <- function(cs, value){
  stopifnot( is( cs, "CellScabbard" ) )
  stopifnot(is(value, "SummarizedExperiment"))
  cs@relevantGenes <- value
=======
#' @export
`oncoNEMdata<-` <- function(cs, value) {
  stopifnot( is( cs, "CellScabbard" ) )
  if(class(value) != "matrix" || class(value[1, ]) != "numeric" ||
     any(value > 2) || any(value < 0) ){
    stop("oncoNEMdata must be a numeric matrix, containing values from 0
         to 2.")
  }
  cs@oncoNEMdata <- value
>>>>>>> refs/remotes/origin/master
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
<<<<<<< HEAD
similarityScores <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@similarityScores
=======

CanopyData <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@CanopyData
>>>>>>> refs/remotes/origin/master
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
<<<<<<< HEAD
`similarityScores<-` <- function(cs, value){
  stopifnot( is( cs, "CellScabbard" ) )
  stopifnot(is(value, "data.frame"))
  cs@similarityScores <- value
=======
`CanopyData<-` <- function(cs, value) {
  stopifnot( is( cs, "CellScabbard" ) )
  if(class(value) != "list" || class(value$R) != "matrix" ||
     class(value$X) != "matrix"){
    stop("CanopyData must be a list containing two matrices and other data.
         See vignette for details.")
  }
  cs@CanopyData <- value
>>>>>>> refs/remotes/origin/master
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
<<<<<<< HEAD
similarityDFs <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@similarityDFs
=======

cellscapeData <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@cellscapeData
>>>>>>> refs/remotes/origin/master
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
<<<<<<< HEAD
#' @export
`similarityDFs<-` <- function(cs, value){
  stopifnot( is( cs, "CellScabbard" ) )
  stopifnot(is(value, "list"))
  cs@similarityDFs <- value
=======
#' @keywords internal
#' @export

`cellscapeData<-` <- function(cs, pt, value) {
  stopifnot(is(cs, "CellScabbard"))
  cs@cellscapeData[[pt]] <- value
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export

timescapeData <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@timescapeData
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @keywords internal
#' @export

`timescapeData<-` <- function(cs, pt, value) {
  stopifnot(is(cs, "CellScabbard"))
  cs@timescapeData[[pt]] <- value
>>>>>>> refs/remotes/origin/master
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
<<<<<<< HEAD
similarityMatrices <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@similarityMatrices
=======

mapscapeData <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@mapscapeData
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @keywords internal
#' @export

`mapscapeData<-` <- function(cs, pt, value) {
  stopifnot(is(cs, "CellScabbard"))
  cs@mapscapeData[[pt]] <- value
  cs
>>>>>>> refs/remotes/origin/master
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
<<<<<<< HEAD
`similarityMatrices<-` <- function(cs, value){
  stopifnot( is( cs, "CellScabbard" ) )
  stopifnot(is(value, "SimpleList"))
  cs@similarityMatrices <- value
=======

treeList <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@treeList
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @keywords internal
#' @export

`treeList<-` <- function(cs, tt, value) {
  stopifnot(is(cs, "CellScabbard"))
  cs@treeList[[tt]] <- value
>>>>>>> refs/remotes/origin/master
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
<<<<<<< HEAD
UNDmatrices <- function( cs ) {
  stopifnot(is( cs, "CellScabbard" ) )
  cs@UNDmatrices
=======

originalTrees <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@originalTrees
>>>>>>> refs/remotes/origin/master
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
<<<<<<< HEAD
#' @export
`UNDmatrices<-` <- function(cs, value){
  stopifnot(is( cs, "CellScabbard" ))
  stopifnot(is(value, "SimpleList"))
  cs@UNDmatrices <- value
  cs
}
=======
#' @keywords internal
#' @export

`originalTrees<-` <- function(cs, tt, value) {
  stopifnot(is(cs, "CellScabbard"))
  cs@originalTrees[[tt]] <- value
  cs
}
>>>>>>> refs/remotes/origin/master

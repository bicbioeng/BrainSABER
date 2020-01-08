#' Methods for the CellScabbard class
#'
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
#' @name CellScabbard-methods
NULL
#'
#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
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
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
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
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
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
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export

oncoNEMdata <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@oncoNEMdata
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
`oncoNEMdata<-` <- function(cs, value) {
  stopifnot( is( cs, "CellScabbard" ) )
  if(class(value) != "matrix" || class(value[1, ]) != "numeric" ||
     any(value > 2) || any(value < 0) ){
    stop("oncoNEMdata must be a numeric matrix, containing values from 0
         to 2.")
  }
  cs@oncoNEMdata <- value
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export

CanopyData <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@CanopyData
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export
`CanopyData<-` <- function(cs, value) {
  stopifnot( is( cs, "CellScabbard" ) )
  if(class(value) != "list" || class(value$R) != "matrix" ||
     class(value$X) != "matrix"){
    stop("CanopyData must be a list containing two matrices and other data.
         See vignette for details.")
  }
  cs@CanopyData <- value
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export

cellscapeData <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@cellscapeData
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
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
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export

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
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export

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
  cs
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @export

originalTrees <- function( cs ) {
  stopifnot( is( cs, "CellScabbard" ) )
  cs@originalTrees
}

#' @rdname CellScabbard-methods
#' @aliases CellScabbard,ANY,ANY-method
#' @keywords internal
#' @export

`originalTrees<-` <- function(cs, tt, value) {
  stopifnot(is(cs, "CellScabbard"))
  cs@originalTrees[[tt]] <- value
  cs
}

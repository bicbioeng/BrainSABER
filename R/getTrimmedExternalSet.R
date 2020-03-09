#' Get a trimmed version of an external data set
#'
<<<<<<< HEAD
#' Returns a SummarizedExperiment that is a subset of \code{dataSet} containing
#' only genes that are present in AIBSARNA, for use in \code{getSimScores} or
#' \code{getUNDmatrix}.
#'
#' @param dataSet a CellScabbard or SummarizedExperiment object
#' @param dataSetId the name of the column of gene identifiers in
#'     rowData(dataSet) to be used to compare dataSet to AIBSARNA.
#' @param AIBSARNA an instance of the AIBSARNA dataset, built using the
#'     \code{buildAIBSARNA()} function
#' @param AIBSARNAid the name of the column of rowData(AIBSARNA) that is
#'     comparable to dataSetId.  One of "gene_id", "ensembl_gene_id",
#'     "gene_symbol", "entrez_id", "refseq_ids"
#'
#' @return a SummarizedExperiment object containing trimmed data set
#' @export
#' @importFrom SummarizedExperiment rowData colData
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
#' # use the appropriate id's to extract the trimmed gene set from the data
#' dataSetId = dataSetId(myGenes)
#' AIBSARNAid = AIBSARNAid(myGenes)
#' myTrimmedGeneSet <- getTrimmedExternalSet(myGenes,
#'     dataSetId = dataSetId, AIBSARNA, AIBSARNAid = AIBSARNAid)
#'
#' # Example 2 - manual gene selection and relevant gene extraction
=======
#' Returns an ExpressionSet that is a subset of \code{dataSet} containing only
#' genes that are present in AIBSARNA, for use in \code{getSimScores} or
#' \code{getUNDmatrix}.
#'
#' @param dataSet a Biobase ExpressionSet compatible data set
#' @param dataSetId the name of the column of gene identifiers in 
#'     fData(dataSet) to be used to compare dataSet to AIBSARNA.
#' @param AIBSARNA an instance of the AIBSARNA dataset, built using the 
#'     \code{buildAIBSARNA()} function
#' @param AIBSARNAid the name of the column of fData(AIBSARNA) that is
#'     comparable to dataSetId.  One of "gene_id", "ensembl_gene_id",
#'     "gene_symbol", "entrez_id", "refseq_ids"
#'
#' @return a Biobase ExpressionSet
#' @export
#' @import Biobase
#'
#' @examples
#' AIBSARNA <- buildAIBSARNA(mini = TRUE)
>>>>>>> refs/remotes/origin/master
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TSPAN6", "DPM1", "C1orf112")
#' myGeneSet <- getRelevantGenes(myGenes, AIBSARNA = AIBSARNA,
#'     AIBSARNAid = "gene_symbol")
#' myTrimmedGeneSet <- getTrimmedExternalSet(myGeneSet,
#'     dataSetId = "gene_symbol", AIBSARNA, AIBSARNAid = "gene_symbol")

getTrimmedExternalSet <- function(dataSet, dataSetId = "gene_symbol",
<<<<<<< HEAD
    AIBSARNA = NULL, AIBSARNAid = c("gene_id", "ensembl_gene_id",
    "gene_symbol", "entrez_id", "refseq_ids")) {

    # get a trimmed vector to faciliate things
    v <- getExternalVector(dataSet = dataSet, index = 1, AIBSARNA,
                dataSetId = dataSetId, AIBSARNAid = AIBSARNAid)
    vInd <- which(match(rowData(dataSet)[[dataSetId]],
                        names(v), nomatch = 0, incomparables = c(NA, "")) > 0)
    # subset feature data
    relfd <- rowData(dataSet)[vInd,]
=======
    AIBSARNA = NULL, AIBSARNAid = c("gene_id", "ensembl_gene_id", 
    "gene_symbol", "entrez_id", "refseq_ids")) {
    # get a trimmed vector faciliate things
    v <- getExternalVector(dataSet = dataSet, index = 1, AIBSARNA, 
                dataSetId = dataSetId, AIBSARNAid = AIBSARNAid)
    vInd <- which(match(fData(dataSet)[[dataSetId]], 
                        names(v), nomatch = 0, incomparables = c(NA, "")) > 0)
    # subset feature data
    relfd <- fData(dataSet)[vInd,]
>>>>>>> refs/remotes/origin/master
    # get index(es) of any duplicates
    dup <- which(duplicated(relfd[[dataSetId]]), arr.ind = TRUE)
    # if there are duplicates, remove them from vInd and regenerate relfd
    if (length(dup) > 0) {
        vInd <- vInd[-dup]
<<<<<<< HEAD
        relfd <- rowData(dataSet)[vInd,]
=======
        relfd <- fData(dataSet)[vInd,]
>>>>>>> refs/remotes/origin/master
    }
    # remove any unused factor levels
    relfd <- as.data.frame(apply(relfd, 2, function(x) {x[drop = TRUE]}))
    # convert to Annotated Data Frame
<<<<<<< HEAD
    # subset exprs
    relexprs <- assay(dataSet)[vInd,]
    # put together SummarizedExperiment
    relevantGenes <- SummarizedExperiment(assays = SimpleList(relexprs),
                                          colData = colData(dataSet),
                                          rowData = relfd)
=======
    relfd <- AnnotatedDataFrame(data = relfd)
    # subset exprs
    relexprs <- exprs(dataSet)[vInd,]
    # put together ExpressionSet
    relevantGenes <- ExpressionSet(assayData = relexprs,
                        phenoData = phenoData(dataSet), featureData = relfd)
>>>>>>> refs/remotes/origin/master
    return(relevantGenes)
}

#' getExternalVector
#'
#' Get a named vector of gene expression values from a single sample in an
#' outside SummarizedExperiment, for use in creating subsets of AIBSARNA with
#'  \code{getRelevantGenes} and comparison with that subset with
#'  \code{getSimScores}
#'
#' @param dataSet a CellScabbard or SummarizedExperiment object
#' @param index the integer index of the sample of dataSet to be used
#' @param AIBSARNA an instance of the AIBSARNA dataset, built using the 
#'     \code{buildAIBSARNA()} function
#' @param dataSetId the name of the column of gene identifiers in 
#'     rowData(dataSet) to be used to compare dataSet to AIBSARNA.
#' @param AIBSARNAid the name of the column of rowData(AIBSARNA) that is
#'     comparable to dataSetId.  One of "gene_id", "ensembl_gene_id",
#'     "gene_symbol", "entrez_id", "refseq_ids"
#'
#' @return a named vector of gene expression values
#' @export
#' @import SummarizedExperiment
#'
#' @examples
#' miniAIBSARNA <- buildAIBSARNA(mini = TRUE)
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TSPAN6", "DPM1", "C1orf112")
#' myGeneSet <- getRelevantGenes(myGenes, "gene_symbol", miniAIBSARNA,
#'     AIBSARNAid = "gene_symbol")
#' myGeneSampleVector <- getExternalVector(myGeneSet, index = 1, miniAIBSARNA,
#'     dataSetId = "gene_symbol", AIBSARNAid = "gene_symbol")
#'
getExternalVector <- function(dataSet, index = 1, AIBSARNA = NULL, dataSetId,
            AIBSARNAid = c("gene_id",
                            "ensembl_gene_id",
                            "gene_symbol",
                            "entrez_id",
                            "refseq_ids")) {
    # check for proper data input
    if(!is(dataSet,"SummarizedExperiment")){
        stop("dataSet must be a CellScabbard or other SummarizedExperiment 
             object")
    }
    if(!is.character(dataSetId) | !is.character(AIBSARNAid)){
        stop("dataSetId and AIBSARNAid must be a character")
    }
    # extract id's from CellScabbard object
    if(is(dataSet, "CellScabbard")){ 
      dataSetId <- dataSetId(dataSet)
      AIBSARNAid <- AIBSARNAid(dataSet)
    }
    if(is.null(AIBSARNA)){
        em <-
            "AIBSARNA is required and must be built using buildAIBSARNA()."
        stop(em)
    }
    v <- assay(dataSet)[, index]
    names(v) <- as.character(rowData(dataSet)[[dataSetId]])
    
    genesToKeep = rowData(AIBSARNA)[[AIBSARNAid]][match(names(v), rowData(AIBSARNA)[[AIBSARNAid]])]
    uniqueGTK = unique(genesToKeep[!is.na(genesToKeep)])
    
    v <- v[uniqueGTK]
    
    return(v)
}

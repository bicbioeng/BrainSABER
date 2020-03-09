#' Get an example vector for specified genes
#'
#' This function returns a named example vector of gene expression values for
#' the specified genes, taken from the 1st row of AIBSARNA, for use in
#' demonstrating \code{getSimScores}.
#'
#' @param genes a character vector of HGNC-compliant gene names
#' @param AIBSARNA an instance of the AIBSARNA dataset, built using the 
#'     \code{buildAIBSARNA()} function
#'
#' @return a named character vector of gene-expression values
#' @import SummarizedExperiment
#' @export
#' @examples
#' AIBSARNA <- buildAIBSARNA(mini = TRUE)
#' myGenes <- c("TSPAN6", "DPM1", "C1orf112")
#' myExampleVector <- getExampleVector(myGenes, AIBSARNA)

getExampleVector <- function(genes, AIBSARNA = NULL) {
    if(is.null(AIBSARNA)){
        em <-
            "AIBSARNA is required and must be built using buildAIBSARNA()."
        stop(em)
    }
    #get relevant genes
    names(genes) <- genes
    relevantGenes <- getRelevantGenes(genes, AIBSARNA = AIBSARNA, 
                                        AIBSARNAid = "gene_symbol")
    #get 8pcw exprs
    v <- as.vector(assay(relevantGenes[, 1]))
    #name v
    names(v) <- as.character(rowData(relevantGenes)$gene_symbol)
    #return v
    return(v)
}

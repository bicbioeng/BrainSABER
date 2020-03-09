#' getExternalVector
#'
#' Get a named vector of gene expression values from a single sample in an
<<<<<<< HEAD
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
=======
#' outside Biobase ExpressionSet, for use in creating subsets of AIBSARNA with
#'  \code{getRelevantGenes} and comparison with that subset with
#'  \code{getSimScores}
#'
#' @param dataSet a Biobase ExpressionSet
#' @param index the integer index of the sample of dataSet to be used
#' @param AIBSARNA an instance of the AIBSARNA dataset, built using the 
#'     \code{buildAIBSARNA()} function
#' @param dataSetId the name of the column of gene identifiers in fData(dataSet)
#'     to be used to compare dataSet to AIBSARNA.
#' @param AIBSARNAid the name of the column of fData(AIBSARNA) that is
>>>>>>> refs/remotes/origin/master
#'     comparable to dataSetId.  One of "gene_id", "ensembl_gene_id",
#'     "gene_symbol", "entrez_id", "refseq_ids"
#'
#' @return a named vector of gene expression values
#' @export
#' @import Biobase
<<<<<<< HEAD
#' @import SummarizedExperiment
=======
>>>>>>> refs/remotes/origin/master
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
<<<<<<< HEAD
        if(is(dataSet, "CellScabbard")){ # extract id's from CellScabbard object
          dataSetId <- dataSetId(dataSet)
          AIBSARNAid <- AIBSARNAid(dataSet)
        }
=======
>>>>>>> refs/remotes/origin/master
        if(is.null(AIBSARNA)){
            em <-
                "AIBSARNA is required and must be built using buildAIBSARNA()."
            stop(em)
        }
<<<<<<< HEAD
        v <- assay(dataSet)[, index]
        names(v) <- as.character(rowData(dataSet)[[dataSetId]])
        # get gene identifiers common to v and AIBSARNA
        matchIdx <- which(match(rowData(AIBSARNA)[[AIBSARNAid]], 
                        names(v), nomatch = 0, incomparables = c(NA, "")) > 0)
        vInAIBSARNA <- 
            as.character(rowData(AIBSARNA)[[AIBSARNAid]][matchIdx])
=======
        v <- exprs(dataSet)[, index]
        names(v) <- as.character(fData(dataSet)[[dataSetId]])
        # get gene identifiers common to v and AIBSARNA
        matchIdx <- which(match(fData(AIBSARNA)[[AIBSARNAid]], 
                        names(v), nomatch = 0, incomparables = c(NA, "")) > 0)
        vInAIBSARNA <- 
            as.character(fData(AIBSARNA)[[AIBSARNAid]][matchIdx])
>>>>>>> refs/remotes/origin/master
        # get indices of v that are present in vInAIBSARNA
        genesToKeep <- which(match(names(v), vInAIBSARNA, nomatch = 0,
                             incomparables = c(NA, "")) > 0)
        # update v to only include comparable genes (genes in vInAIBSARNA)
        v <- v[genesToKeep]
        # remove any duplicate genes
        genesToKeep <- unique(names(v))
        v <- v[genesToKeep]
        return(v)
    }

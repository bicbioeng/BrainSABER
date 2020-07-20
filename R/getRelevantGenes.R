#' Get a subset of AIBSARNA using a Gene Expression Vector
#'
#' This function returns a subset of the AIBSARNA dataset, containing only
#' the genes in \code{data}, which may be a vector, a SummarizedExperiment
#' or derivative \code{assay()} and \code{rowData()}, or a CellScabbard.  
#' If a vector is used, it must consist of numerical gene expression values 
#' with names comparable to one column of identifiers present in AIBSARNA.
#' If \code{data} is a CellScabbard, results are stored in the relevantGenes
#' slot of the object.
#'
#' @param data a vector of named gene expression values, or a compatible
#'     data set
#' @param dataSetId (Optional) If \code{data} is not a vector, the name of
#'     the column of gene identifiers in rowData(dataSet) to be used to compare
#'     \code{data} to AIBSARNA.
#' @param AIBSARNA an instance of the AIBSARNA dataset, built using the 
#'     \code{buildAIBSARNA()} function
#' @param AIBSARNAid the name of the column of rowData(AIBSARNA) that is
#'     comparable to dataSetId.  One of "gene_id", "ensembl_gene_id",
#'     "gene_symbol", "entrez_id", "refseq_ids"
#'
#' @return a SummarizedExperiment consisting of genes in \code{data}, sorted
#'     to match the order of the genes in \code{data}
#' @import SummarizedExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom methods is
#' @export
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
#' relevantGenes(myGenes)
#' 
#' # Example 2 - manual gene selection and relevant gene extraction
#' myGenes <- c(4.484885, 0.121902, 0.510035)
#' names(myGenes) <- c("TSPAN6", "DPM1", "C1orf112")
#' myGeneSet <- getRelevantGenes(myGenes, AIBSARNA = AIBSARNA, 
#'     AIBSARNAid = "gene_symbol")
#'
getRelevantGenes <- function(data, dataSetId = NULL, AIBSARNA = NULL,
                    AIBSARNAid = c("gene_id", "ensembl_gene_id", "gene_symbol",
                                    "entrez_id", "refseq_ids")) {
    # check for proper input
    if(is.null(AIBSARNA)){
        em <-
            "AIBSARNA is required and must be built using buildAIBSARNA()."
        stop(em)
    }
    if(!is.character(AIBSARNAid)){
        stop("AIBSARNAid must be a character")
    }
    
    # get a single sample vector if data is not a vector
    if (is.vector(data)) {
        v <- data
    } else {
        #stop if data is not a vector and dataSetId == NULL
        stopifnot(!is.null(dataSetId))
        # get a trimmed vector faciliate things
        v <- getExternalVector(dataSet=data, index=1, AIBSARNA = AIBSARNA, 
                        dataSetId=dataSetId, AIBSARNAid = AIBSARNAid)
    }
    #relevantGenes <- subset of AIBSARNA containing only genes specified in v
    # get indices of genes common to v and AIBSARNA
    vInd <- which(match(rowData(AIBSARNA)[[AIBSARNAid]], 
                        names(v), nomatch = 0, incomparables = c(NA, "")) > 0)
    relfd <- rowData(AIBSARNA)[vInd,] # subset feature data
    # get indices of any duplicates, remove them from vInd and regenerate relfd
    dup <- which(duplicated(relfd[[AIBSARNAid]]), arr.ind = TRUE)
    if (length(dup) > 0) {
        vInd <- vInd[-dup]
        relfd <- rowData(AIBSARNA)[vInd,]
    }
    # remove any unused factor levels
    relfd <- as.data.frame(apply(relfd, 2, function(x) {x[drop = TRUE]}))
    relexprs <- assay(AIBSARNA)[vInd,] # subset exprs
    # sort relfd and relexprs so the genes are in the same order as dataSet
    relfd[[AIBSARNAid]] <- factor(relfd[[AIBSARNAid]], levels=unique(names(v)))
    relfd <- relfd[order(relfd[[AIBSARNAid]]), ]
    relfd$row_num <- factor(relfd$row_num, levels = unique(relfd$row_num))
    
    #order relexprs genes by dataset order
    rownames(relexprs) <- rowData(AIBSARNA)[[AIBSARNAid]][vInd]
    relexprs <- relexprs[as.character(relfd[[AIBSARNAid]]),]
    
    # convert relexprs to data.frame and add a column to do factor manipulation
    #relexprs <- as.data.frame(relexprs)
    #row_num <- row.names(relexprs)
    #relexprs <- cbind(relexprs, row_num)
    #relexprs$row_num <- as.character(relexprs$row_num)
    #relexprs$row_num <- factor(relexprs$row_num,
    #                            levels = trimws(as.character(relfd$row_num)))
    #relexprs <- relexprs[order(relexprs$row_num), ]
    # convert relexprs back to a matrix, remove the added column
    #lastcol <- length(relexprs[1,])
    #relexprs <- as.matrix(relexprs[,-lastcol])
    # put together SummarizedExperiment
    relevantGenes <- SummarizedExperiment(assays = SimpleList(relexprs),
        colData = colData(AIBSARNA), rowData = relfd)
    return(relevantGenes)
}

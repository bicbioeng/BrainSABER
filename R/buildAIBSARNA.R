#' Function to create a SummarizedExperiment containing BrainSpan Data
#'
#' This function is used to build the AIBSARNA SummarizedExperiment object, and
#' must be run  prior to running any other function in \pkg{BrainSABER}. This
#' function will download the data from http://brainspan.org and may take
#' several minutes, depending on internet connection speeds.
#' @param mini Default is FALSE.  If \code{mini=TRUE}, build a miniature
#'     version of AIBSARNA that does not require internet connectivity and is
#'     suitable for example purposes only
#' @return A SummarizedExperiment containing BrainSpan data, with the addition
#'     of RefSeq IDs via biomaRt
#' @import SummarizedExperiment
#' @importFrom utils download.file read.csv
#' @importFrom biomaRt useEnsembl getBM
#' @export
#' @examples
#' AIBSARNA <- buildAIBSARNA(mini = TRUE)
buildAIBSARNA <- function(mini = FALSE){
    # validate appropriate input
    if(!is.logical(mini)){
        stop("argument must be a logical")
    }
    if(mini){
        miniexprs <- matrix(
            data = c(36.447128, 24.251253, 19.330479, 27.668607, 19.998231,
                     0.044081, 0.067338, 0.000000, 0.145466, 0.185188,
                     34.373239, 20.765661, 18.734947, 22.366394, 19.228431,
                     4.379337, 4.227521, 2.551825, 3.603764, 2.948976,
                     3.957119, 3.520794, 2.037805, 3.487035, 2.177235),
            nrow = 5, ncol = 5, byrow = TRUE)
        rownames(miniexprs) <- c("1", "2", "3", "4", "5")
        colnames(miniexprs) <- c("1", "2", "3", "4", "5")
        minifd <- as.data.frame(
            rbind(c(1, 7062, "ENSG00000000003", "TSPAN6", 7105),
                c(2, 40735, "ENSG00000000005", "TNMD", 64102),
                c(3, 8736, "ENSG00000000419", "DPM1", 8813),
                c(4, 36423, "ENSG00000000457", "SCYL3", 57147),
                c(5, 35021, "ENSG00000000460", "C1orf112", 55732)),
             stringsAsFactors = FALSE)
        colnames(minifd) <- c("row_num", "gene_id", "ensembl_gene_id",
                                "gene_symbol", "entrez_id")
        rownames(minifd) <- c("1", "2", "3", "4", "5")
        minipd <- as.data.frame(
            rbind(c(1, 3058, "H376.IIA.51", "8 pcw", "M", 10268, "Ocx"),
                c(2, 13058, "H376.IIA.51", "8 pcw", "M", 10291, "M1C-S1C"),
                c(3, 13058, "H376.IIA.51", "8 pcw", "M", 10361, "AMY"),
                c(4, 13058, "H376.IIA.51", "8 pcw", "M", 10550, "MGE"),
                c(5, 13058, "H376.IIA.51", "8 pcw", "M", 10243, "STC")),
            stringsAsFactors = TRUE)
        colnames(minipd) <- c("column_num", "donor_id", "donor_name",
                        "age", "gender", "structure_id", "structure_acronym")
        rownames(minipd) <- c("1", "2", "3", "4", "5")
        miniAIB <- SummarizedExperiment(assays = SimpleList(miniexprs),
                                        colData = minipd,
                                        rowData = minifd)
        return(miniAIB)
    }
    # download AIBSARNA
    temp <- tempfile()
    download.file(
        "http://www.brainspan.org/api/v2/well_known_file_download/267666525",
        temp, mode = "wb")
    #read in csv files and clean them up
    exprs <- as.matrix(read.csv(
        file=unz(temp, "expression_matrix.csv"),
        header=FALSE, sep=","))
    exprs <- exprs[, 2:525]
    colnames(exprs) <- paste0("Sample_",as.character(c(seq(1,524))))
    rownames(exprs) <- paste0("Gene_",as.character(c(seq(1,length(exprs[,1])))))
    pd <- read.csv(file=unz(temp, "columns_metadata.csv"),
                   header=TRUE, sep=",")
    #reorder factor levels in pd$age so they are interpreted chronologically
    pd$age <- factor(pd$age, levels = unique(pd$age))
    rownames(pd) <- colnames(exprs)
    fd <- read.csv(file=unz(temp, "rows_metadata.csv"),
                   header=TRUE, sep=",", stringsAsFactors = FALSE)
    rownames(fd) <- rownames(exprs)
    #set up bioMaRt
    ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                                version = 67)
    #get ensembl.ids from AIBSARNA
    ensemblIds <- as.character(fd$ensembl_gene_id)
    #get map of applicable refseq_mrna ids
    geneMap <- biomaRt::getBM(attributes = c("ensembl_gene_id", "refseq_mrna"),
                              filters = "ensembl_gene_id",
                              values = ensemblIds,
                              mart = ensembl)
    idx <- match(ensemblIds, geneMap$ensembl_gene_id)
    refseq_ids <- geneMap$refseq_mrna[idx]
    #bind refseq_ids to AIBSARNA rowData
    fd <- cbind(fd, refseq_ids)
    #create SummarizedExperiment
    AIBSARNA <- SummarizedExperiment(assays = SimpleList(exprs),rowData=fd,
                                     colData=pd)
    #return AIBSARNA
    return(AIBSARNA)
}

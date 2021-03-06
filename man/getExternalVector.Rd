% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getExternalVector.R
\name{getExternalVector}
\alias{getExternalVector}
\title{getExternalVector}
\usage{
getExternalVector(
  dataSet,
  index = 1,
  AIBSARNA = NULL,
  dataSetId,
  AIBSARNAid = c("gene_id", "ensembl_gene_id", "gene_symbol", "entrez_id",
    "refseq_ids")
)
}
\arguments{
\item{dataSet}{a CellScabbard or SummarizedExperiment object}

\item{index}{the integer index of the sample of dataSet to be used}

\item{AIBSARNA}{an instance of the AIBSARNA dataset, built using the 
\code{buildAIBSARNA()} function}

\item{dataSetId}{the name of the column of gene identifiers in 
rowData(dataSet) to be used to compare dataSet to AIBSARNA.}

\item{AIBSARNAid}{the name of the column of rowData(AIBSARNA) that is
comparable to dataSetId.  One of "gene_id", "ensembl_gene_id",
"gene_symbol", "entrez_id", "refseq_ids"}
}
\value{
a named vector of gene expression values
}
\description{
Get a named vector of gene expression values from a single sample in an
outside SummarizedExperiment, for use in creating subsets of AIBSARNA with
 \code{getRelevantGenes} and comparison with that subset with
 \code{getSimScores}
}
\examples{
miniAIBSARNA <- buildAIBSARNA(mini = TRUE)
myGenes <- c(4.484885, 0.121902, 0.510035)
names(myGenes) <- c("TSPAN6", "DPM1", "C1orf112")
myGeneSet <- getRelevantGenes(myGenes, "gene_symbol", miniAIBSARNA,
    AIBSARNAid = "gene_symbol")
myGeneSampleVector <- getExternalVector(myGeneSet, index = 1, miniAIBSARNA,
    dataSetId = "gene_symbol", AIBSARNAid = "gene_symbol")

}

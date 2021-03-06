% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getExampleVector.R
\name{getExampleVector}
\alias{getExampleVector}
\title{Get an example vector for specified genes}
\usage{
getExampleVector(genes, AIBSARNA = NULL)
}
\arguments{
\item{genes}{a character vector of HGNC-compliant gene names}

\item{AIBSARNA}{an instance of the AIBSARNA dataset, built using the 
\code{buildAIBSARNA()} function}
}
\value{
a named character vector of gene-expression values
}
\description{
This function returns a named example vector of gene expression values for
the specified genes, taken from the 1st row of AIBSARNA, for use in
demonstrating \code{getSimScores}.
}
\examples{
AIBSARNA <- buildAIBSARNA(mini = TRUE)
myGenes <- c("TSPAN6", "DPM1", "C1orf112")
myExampleVector <- getExampleVector(myGenes, AIBSARNA)
}

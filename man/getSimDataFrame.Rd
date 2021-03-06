% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSimDataFrame.R
\name{getSimDataFrame}
\alias{getSimDataFrame}
\title{Get Age, Structure Acronym, and Similarity Scores Data Frame}
\usage{
getSimDataFrame(
  data = NULL,
  sim_score = NULL,
  relevantGenes = NULL,
  similarity_method = "cosine"
)
}
\arguments{
\item{data}{a CellScabbard object with non-empty relevantGenes and
similarityScores slots, or a SummarizedExperiment created using the
    \code{getRelevantGenes()} function}

\item{sim_score}{a vector or data frame of similarity scores}

\item{relevantGenes}{a SummarizedExperiment object created using the
\code{getRelevantGenes()} function}

\item{similarity_method}{currently supported similarity methods are "cosine"
and "euclidean", defaults to "cosine"}
}
\value{
a three-column data.frame or list of data frames
}
\description{
This function takes in a CellScabbard object, or both a similarity vector
or data frame as returned by \code{getSimScores} and a subset of AIBSARNA
as returned by \code{getRelevantGenes}. Constructs a \code{data.frame}
with columns age, structure_acroynym, and either cosine_similarity or
euclidean_similarity, sorted by similarity score in decreasing order.
In the case of a similarity data frame, a list of data frames is returned.
}
\examples{
AIBSARNA <- buildAIBSARNA(mini = TRUE)
# Example 1 - using CellScabbard class
# get a random sample of 3 genes
totalGenes <- nrow(AIBSARNA)
gene_idx <- sample.int(totalGenes, 3)
sample_idx <- c(1,3,5)
# Subset AIBSARNA
exprs <- assay(AIBSARNA)[gene_idx, sample_idx]
fd <- rowData(AIBSARNA)[gene_idx, ]
pd <- colData(AIBSARNA)[sample_idx, ]
# build a trimmed data set
myGenes <- CellScabbard(exprsData = exprs, phenoData = pd, featureData = fd,
                        AIBSARNA = AIBSARNA, autoTrim = TRUE)
# cosine similarity method
similarityScores(myGenes) <- getSimScores(data = myGenes,
                                          similarity_method = "cosine")
similarityDFs(myGenes) <- getSimDataFrame(data = myGenes,
                                          similarity_method = "cosine")
similarityDFs(myGenes)
# euclidean similarity method
similarityScores(myGenes) <- getSimScores(data = myGenes,
                                          similarity_method = "euclidean")
similarityDFs(myGenes) <- getSimDataFrame(data = myGenes,
                                          similarity_method = "euclidean")
similarityDFs(myGenes)

# Example 2 - manual gene selection and relevant gene extraction
myGenes <- c(4.484885, 0.121902, 0.510035)
names(myGenes) <- c("TSPAN6", "DPM1", "C1orf112")
myGeneSet <- getRelevantGenes(myGenes, AIBSARNA = AIBSARNA,
    AIBSARNAid = "gene_symbol")
myCosScore <- getSimScores(myGenes, relevantGenes = myGeneSet,
                           similarity_method = "cosine")
myEucScore <- getSimScores(myGenes, relevantGenes = myGeneSet,
    similarity_method = "euclidean")
myCosineDF <- getSimDataFrame(sim_score = myCosScore,
                              relevantGenes = myGeneSet,
                              similarity_method = "cosine")
myEuclideanDF <- getSimDataFrame(sim_score = myEucScore,
                                 relevantGenes = myGeneSet,
                                 similarity_method = "euclidean")
}

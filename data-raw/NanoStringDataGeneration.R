#NanoString data generation
library(nanoSet)
data("NanoSet")
# extract sample 1 as a named vector
NanoSample1 <- exprs(NanoSet)[, 1]
names(NanoSample1) <- as.character(fData(NanoSet)$gene_symbol)
# get relevant subset of AIBSARNA
NanoSampleSet <- getRelevantGenes(NanoSample1)
#get indices of NanoSample1 that are present in NanoSampleSet
genesToKeep <- which(names(NanoSample1) %in% fData(NanoSampleSet)$gene_symbol,
              +               arr.ind = TRUE)
#update NanoSample1 to only include comparable genes (genes in NanoSampleSet)
NanoSample1 <- NanoSample1[genesToKeep]
# remove any duplicate genes
genesToKeep <- unique(names(NanoSample1))
NanoSample1 <- NanoSample1[genesToKeep]

NanoSampleSet <- getSimScores(NanoSample1, NanoSampleSet,
                              similarity_method = "cosine")
NanoSampleSet <- getSimScores(NanoSample1, NanoSampleSet,
                              similarity_method = "euclidean")
NanoSample1cosMat <- getSimMatrix(NanoSampleSet, similarity_method = "cosine")
NanoSample1eucMat <- getSimMatrix(NanoSampleSet, similarity_method = "euclidean")

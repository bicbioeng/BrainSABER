#NanoString data generation
library(nanoSet)
data("NanoSet")
# extract sample 1 as a named vector
NanoSample1 <- exprs(NanoSet)[, 1]
names(NanoSample1) <- fData(NanoSet)$RefSeq
# get relevant subset of AIBSARNA
NanoSampleSet <- getRelevantGenes(NanoSample1, gene_names = "RefSeq")

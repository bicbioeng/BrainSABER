#' ATRT groups data generation script
#'
#' This function is the code used to generate the data sets used in the
#' ATRT heatmaps document.
#'
#' @import Biobase
#' @return void

generateATRTdata <- function(){
  #define gene lists
  MYCgenes <- c("BCL3", "COL1A2", "COL3A1", "COL5A2", "ERBB2", "EZH1", "EZH2",                     "MYC", "MYD88", "NEFH", "SHC1", "TNFRSF1A")
  SHHgenes <- c("AADAT", "ACVR2B", "APC2", "ARNT2", "CDK6", "CRNKL1", "DPYSL5",                    "EFNA2", "EZH1", "EZH2", "GLI1", "GLI2", "GNAI1", "HNRNPA3",
                "ID4", "LIFR", "MLL", "MYCN", "PBX1", "PTCH1", "RBM8A",
                "SEMA5B", "SEMA6A", "SF3A3", "SHH", "SRGAP1", "SRSF10", "TCF3",                    "WHSC1L1")
  TYRgenes <- c("ACOX2", "ACSL1", "AGPAT2", "ALDH9A1", "BAMBI", "BMP4", "CAT",                     "CCND1", "CYP7B1", "DGAT2", "DLG1", "ECHS1", "EDN1", "EPHX2",                     "ERBB2", "EZH1", "EZH2", "FZD6", "GALNT10", "GALNT12", "GALNT2",                  "GALNT5", "GDF5", "HSD17B4", "LPIN1", "LPIN2", "MITF", "NFATC1",                  "NUDT7", "OGDHL", "PSEN1", "SFRP5", "SOST", "TCF7", "TYR",
                 "VEGFA", "WNT11", "WNT7A", "WTIP", "WWC1")

  #build example vectors
  MYC8pcw <- getExampleVector(MYCgenes)
  SHH8pcw <- getExampleVector(SHHgenes)
  TYR8pcw <- getExampleVector(TYRgenes)
  ATRT8pcw <- getExampleVector(unique(c(MYCgenes, SHHgenes, TYRgenes)))

  #get relevant gene sets
  MYCset <- getRelevantGenes(MYC8pcw)
  SHHset <- getRelevantGenes(SHH8pcw)
  TYRset <- getRelevantGenes(TYR8pcw)
  ATRTset <- getRelevantGenes(ATRT8pcw)

  #get cosine similarity scores
  MYCset <- getSimScores(MYC8pcw, MYCset, similarity_method = "cosine")
  SHHset <- getSimScores(SHH8pcw, SHHset, similarity_method = "cosine")
  TYRset <- getSimScores(TYR8pcw, TYRset, similarity_method = "cosine")
  ATRTset <- getSimScores(ATRT8pcw, ATRTset, similarity_method = "cosine")

  #get euclidean distances
  MYCset <- getSimScores(MYC8pcw, MYCset, similarity_method = "euclidean")
  SHHset <- getSimScores(SHH8pcw, SHHset, similarity_method = "euclidean")
  TYRset <- getSimScores(TYR8pcw, TYRset, similarity_method = "euclidean")
  ATRTset <- getSimScores(ATRT8pcw, ATRTset, similarity_method = "euclidean")

  #build cosine similarity matrices
  MYCcosMat <- getSimMatrix(MYCset, similarity_method = "cosine")
  SHHcosMat <- getSimMatrix(SHHset, similarity_method = "cosine")
  TYRcosMat <- getSimMatrix(TYRset, similarity_method = "cosine")
  ATRTcosMat <- getSimMatrix(ATRTset, similarity_method = "cosine")

  #build euclidean distance matrices
  MYCeucMat <- getSimMatrix(MYCset, similarity_method = "euclidean")
  SHHeucMat <- getSimMatrix(SHHset, similarity_method = "euclidean")
  TYReucMat <- getSimMatrix(TYRset, similarity_method = "euclidean")
  ATRTeucMat <- getSimMatrix(ATRTset, similarity_method = "euclidean")
}

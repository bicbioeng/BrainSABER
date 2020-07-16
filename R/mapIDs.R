#' mapIDs
#'
#' @param dataSet a CellScabbard or SummarizedExperiment object containing
#' a gene expression matrix along with phenotype and sample feature data
#' @param AIBSARNA an instance of the AIBSARNA dataset
#' 
#' @return a CellScabbard object containing a subset of relevant genes from the
#' Allen institute data
#' 
#' @export

mapIDs <- function(dataSet, AIBSARNA){

  cs <- CellScabbard(exprsData = as.matrix(dataSet[,-1]),
               phenoData = data.frame(sample = colnames(dataSet[,-1]),
                                      sample1 = colnames(dataSet[,-1])),
               featureData = data.frame(gene = dataSet[,1],
                                        gene1 = dataSet[,1]),
               AIBSARNA = AIBSARNA, autoTrim = TRUE)
  allenVariantGenesPath <- file.path(system.file(package = "BrainSABER"),
                                     'variant_genes_list',
                                     'allenVariantGenes.csv')
  allenVariantGenes <- read.csv(allenVariantGenesPath, stringsAsFactors = FALSE, header = TRUE)[,2]
  allenVariantGeneID <- rowData(AIBSARNA)@listData[AIBSARNAid(cs)][[1]][allenVariantGenes]

  cs <- cs[rowData(cs)$gene %in% allenVariantGeneID,]
}


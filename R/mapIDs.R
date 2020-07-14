#' mapIDs
#'
#' placeholder
#' 
#' @export

mapIDs <- function(dataSet, AIBSARNA){

  cs <- CellScabbard(exprsData = as.matrix(dataSet[,-1]),
               phenoData = data.frame(sample = colnames(dataSet[,-1]),
                                      sample1 = colnames(dataSet[,-1])),
               featureData = data.frame(gene = dataSet[,1],
                                        gene1 = dataSet[,1]),
               AIBSARNA = AIBSARNA, autoTrim = T)

  allenVariantGenes <- read.csv(file.path(system.file('variant_genes_list', package = "shinyBrainSABER"),
                                          'allenVariantGenes.csv'), stringsAsFactors = F, header = T)[,2]
  allenVariantGeneID <- rowData(AIBSARNA)@listData[AIBSARNAid(cs)][[1]][allenVariantGenes]

  cs <- cs[rowData(cs)$gene %in% allenVariantGeneID,]
}


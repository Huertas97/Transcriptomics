setwd("~/Escritorio/Máster/Segundo_cuatrimestre/1_Transcriptómica/3_Entregas/Enrique/Datos/RESULTS/Modelo_11_estados")

# 1. Libraries' installation and importation

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager, ")
packages <- matrix(data = c("ChIPpeakAnno",
                            "TxDb.Hsapiens.UCSC.hg19.knownGene"))
for(p in packages) {
  #if(!require(p)) {BiocManager::install(p)}
  library(p, character.only=TRUE)
}

gr1 <- toGRanges("./Monocytes_E1.bed", format="BED", header=FALSE)
gr2 <- toGRanges("./Monocytes_E2.bed", format="BED", header=FALSE)
aCR1 <- assignChromosomeRegion(gr1, nucleotideLevel=FALSE,
                               precedence=c("Promoters", "immediateDownstream",
                                            "fiveUTRs", "threeUTRs",
                                            "Exons", "Introns"),
                               TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
aCR2 <- assignChromosomeRegion(gr2, nucleotideLevel=FALSE,
                               precedence=c("Promoters", "immediateDownstream",
                                            "fiveUTRs", "threeUTRs",
                                            "Exons", "Introns"),
                               TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
barplot(aCR1$percentage, las=3)
barplot(aCR2$percentage, las=3)

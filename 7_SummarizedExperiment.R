library(TCGAbiolinks)
load("G:/TFM UOC/datos/etiqueta.Rda")
load("G:/TFM UOC/datos/matriz_sna.Rda")


########################################
# Anotaciones de las sondas
#######################################
# anotaciones de las sondas 450k

library(FDb.InfiniumMethylation.hg19)
library(dplyr)

sondas <- get450k()
sondas <- sondas[rownames(matriz_sna)]
head(sondas)

##################################################
# Se construye el objeto SummarizedExperiment
############################################


data <- SummarizedExperiment::SummarizedExperiment(
  assays=S4Vectors::SimpleList(counts=matriz_sna),
  rowRanges = sondas,
  colData = etiqueta
)

save(data, file="G:/TFM UOC/datos/SummarizedExperiment.Rda")

data

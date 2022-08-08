library(curatedTCGAData)
library(TCGAutils)
library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation27kanno.ilmn12.hg19)

data("diseaseCodes", package="TCGAutils")
diseaseCodes

data("sampleTypes")
sampleTypes


###################################################################################
curatedTCGAData("*", "Methyl*", dry.run = TRUE, version = "1.1.38")
estudios <- curatedTCGAData(
  diseaseCode = c("ACC", "BLCA"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)

# el número de muestras
sampleTables(estudios)

save(estudios, file= "./estudios/estudios_1")
estudios

################################################################################
###############################################################################
estudios_2 <- curatedTCGAData(
  diseaseCode = c("BRCA"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_2)

save(estudios_2, file= "./estudios/estudios_2")
estudios_2

##############################################################################
###############################################################################
estudios_3 <- curatedTCGAData(
  diseaseCode = c("CESC", "CHOL","COAD","DLBC","ESCA"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_3)

save(estudios_3, file= "./estudios/estudios_3")
estudios_3

#############################################################################
###############################################################################
estudios_4 <- curatedTCGAData(
  diseaseCode = c("GBM", "HNSC", "KICH", "KIRC"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_4)

save(estudios_4, file= "./estudios/estudios_4")
estudios_4

#############################################################################
###############################################################################
estudios_5 <- curatedTCGAData(
  diseaseCode = c("KIRP", "LAML", "LGG", "LIHC"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_5)

save(estudios_5, file= "./estudios/estudios_5")
estudios_5

#############################################################################
###############################################################################
estudios_6 <- curatedTCGAData(
  diseaseCode = c("LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_6)

save(estudios_6, file= "./estudios/estudios_6")
estudios_6

#############################################################################
###############################################################################
estudios_7 <- curatedTCGAData(
  diseaseCode = c("PRAD", "READ", "SARC", "SKCM"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_7)

save(estudios_7, file= "./estudios/estudios_7")
estudios_7

#############################################################################
###############################################################################
estudios_8 <- curatedTCGAData(
  diseaseCode = c("TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_8)

save(estudios_8, file= "./estudios/estudios_8")
estudios_8

#############################################################################
###############################################################################
estudios_9 <- curatedTCGAData(
  diseaseCode = c("STAD"), assays = "Methylation_methyl27*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_9, vial=TRUE)

save(estudios_9, file= "./estudios/estudios_9")
estudios_9

#############################################################################
###############################################################################
estudios_10 <- curatedTCGAData(
  diseaseCode = c("STAD"), assays = "Methylation_methyl450*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_10)

save(estudios_10, file= "./estudios/estudios_10")
estudios_10

#############################################################################


load(file="./estudios/estudios_1")
load(file="./estudios/estudios_2")
load(file="./estudios/estudios_3")
load(file="./estudios/estudios_4")
load(file="./estudios/estudios_5")
load(file="./estudios/estudios_6")
load(file="./estudios/estudios_7")
load(file="./estudios/estudios_8")
load(file="./estudios/estudios_9")


###############################################################################

# obtener los barcodes
(xbarcode <- head(colnames(estudios_2)[[1]], 1000L))

# accedo a los datos finales de un estudio
assay(estudios[[1]], "counts") %>% max(na.rm=TRUE)


matriz <- assay(estudios_3[[2]], "counts") 
matriz[1:5,1:10]
dim(matriz)
sum(is.na(matriz))*100 / (dim(matriz)[1] * dim(matriz)[2])


# detalle de los experimentos en un mUltiassay
experimentos <- experiments(estudios)
experimentos
experimentos$`ACC_Methylation-20160128`
experimentos$`BLCA_Methylation-20160128`

# la matriz de fenotipos como un dataframe
fenotipos_1 <- colData(estudios)@listData %>% as.data.frame()
fenotipos_2 <- colData(estudios_2)@listData %>% as.data.frame()
fenotipos_3 <- colData(estudios_3)@listData %>% as.data.frame()
fenotipos_4 <- colData(estudios_4)@listData %>% as.data.frame()
fenotipos_5 <- colData(estudios_5)@listData %>% as.data.frame()
fenotipos_6 <- colData(estudios_6)@listData %>% as.data.frame()
fenotipos_7 <- colData(estudios_7)@listData %>% as.data.frame()
fenotipos_8 <- colData(estudios_8)@listData %>% as.data.frame()
fenotipos_9 <- colData(estudios_9)@listData %>% as.data.frame()

######################################################################
#####################################################################
# anotaciones de las sondas 450k

data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation.table <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation.table[1,]
dim(annotation.table)

ww <- assay(estudios[[2]], "counts")[1:10,] %>% as.data.frame()
colnames(ww) %>% substring(14,16) %>% table()
nombres <-row.names(ww)

annotation.table[nombres, c(1,2,3,24,25,26)]

#######################################################################
#######################################################################
# anotaciones de las sondas 27k

data("IlluminaHumanMethylation27kanno.ilmn12.hg19")
annotation27.table <- getAnnotation(IlluminaHumanMethylation27kanno.ilmn12.hg19)
annotation27.table
dim(annotation27.table)

ww <- assay(estudios_2[[1]], "counts")[1:10,] %>% as.data.frame()
colnames(ww) %>% substring(14,16) %>% table()
nombres <-row.names(ww)

annotation27.table[nombres, c(1,2,3,30,32)]

#############################################################################
################################################################################

library(dplyr)

matriz1 <- assays(estudios)
matriz1 <- matriz1[["BLCA_Methylation-20160128"]]
dim(matriz1)
matriz1[1:3, 1:3]
matriz1 <- as.matrix(matriz1)

fenotipos <- colData(estudios)
fenotipos$histological_type %>% table()
variables <- as.data.frame(fenotipos[1:5,])
write.table(variables, "./estudios/variables.txt")

etiqueta1 <- data.frame(bar_code = colnames(matriz1))
etiqueta1$sujeto <- substring(colnames(matriz1), 1, 12)
etiqueta1$estudio <- "BLCA"
etiqueta1$categoria <- substring(colnames(matriz1), 14,15) %>% as.factor()
tipos <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta1 <- merge(etiqueta1, tipos, by="sujeto", all.x=TRUE)
etiqueta1$type[etiqueta1$categoria == "11" ] <- "Control"
etiqueta1$type <- as.factor(etiqueta1$type)
etiqueta1$label <- etiqueta1$estudio
etiqueta1$label[etiqueta1$categoria == "11"] <- "Control"

##############################################################################
################################################################################

matriz2 <- assays(estudios)
matriz2 <- matriz2[["ACC_Methylation-20160128"]]
dim(matriz2)
matriz2[1:3, 1:3]
matriz2 <- as.matrix(matriz2)

etiqueta2 <- data.frame(bar_code = colnames(matriz2))
etiqueta2$sujeto <- substring(colnames(matriz2), 1, 12)
etiqueta2$estudio <- "ACC" 
etiqueta2$categoria <- substring(colnames(matriz2), 14,15) %>% as.factor()
tipos2 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta2 <- merge(etiqueta2, tipos2, by="sujeto", all.x=TRUE)
etiqueta2$type[etiqueta2$categoria == "11" ] <- "Control"
etiqueta2$type <- as.factor(etiqueta2$type)
etiqueta2$label <- etiqueta2$estudio
etiqueta2$label[etiqueta2$categoria == "11"] <- "Control"

#################################################################################
#################################################################################

matriz3 <- assays(estudios_2)
matriz3 <- matriz3[["BRCA_Methylation_methyl450-20160128"]]
dim(matriz3)
matriz3[1:3, 1:3]
matriz3 <- as.matrix(matriz3)
fenotipos <- colData(estudios_2)
fenotipos$histological_type %>% table()
etiqueta3 <- data.frame(bar_code = colnames(matriz3))
etiqueta3$sujeto <- substring(colnames(matriz3), 1, 12)
etiqueta3$estudio <- "BRCA" 
etiqueta3$categoria <- substring(colnames(matriz3), 14,15) %>% as.factor()

tipos3 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta3 <- merge(etiqueta3, tipos3, by="sujeto", all.x=TRUE)
etiqueta3$type[etiqueta3$categoria == "11" ] <- "Control"
etiqueta3$type <- as.factor(etiqueta3$type)
etiqueta3$label <- etiqueta3$estudio
etiqueta3$label[etiqueta3$categoria == "11"] <- "Control"

#################################################################################
###################################################################################
load(file="./estudios/estudios_3")
experiments(estudios_3)

matriz4 <- assays(estudios_3)
matriz4 <- matriz4[["CESC_Methylation-20160128"]]
dim(matriz4)
matriz4[1:3, 1:3]
matriz4 <- as.matrix(matriz4)
fenotipos <- colData(estudios_3)
fenotipos$histological_type %>% table()
etiqueta4 <- data.frame(bar_code = colnames(matriz4))
etiqueta4$sujeto <- substring(colnames(matriz4), 1, 12)
etiqueta4$estudio <- "CESC" 
etiqueta4$categoria <- substring(colnames(matriz4), 14,15) %>% as.factor()

tipos4 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta4 <- merge(etiqueta4, tipos4, by="sujeto", all.x=TRUE)
etiqueta4$type[etiqueta4$categoria == "11" ] <- "Control"
etiqueta4$type <- as.factor(etiqueta4$type)
etiqueta4$label <- etiqueta4$estudio
etiqueta4$label[etiqueta4$categoria == "11"] <- "Control"

#################################################################################
###################################################################################
load(file="./estudios/estudios_3")
experiments(estudios_3)

matriz5 <- assays(estudios_3)
matriz5 <- matriz5[["CHOL_Methylation-20160128"]]
dim(matriz5)
matriz5[1:3, 1:3]
matriz5 <- as.matrix(matriz5)
fenotipos <- colData(estudios_3)
etiqueta5 <- data.frame(bar_code = colnames(matriz5))
etiqueta5$sujeto <- substring(colnames(matriz5), 1, 12)
etiqueta5$estudio <- "CHOL" 
etiqueta5$categoria <- substring(colnames(matriz5), 14,15) %>% as.factor()

tipos5 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta5 <- merge(etiqueta5, tipos5, by="sujeto", all.x=TRUE)
etiqueta5$type[etiqueta5$categoria == "11" ] <- "Control"
etiqueta5$type <- as.factor(etiqueta5$type)
etiqueta5$label <- etiqueta5$estudio
etiqueta5$label[etiqueta5$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_3")
experiments(estudios_3)

matriz6 <- assays(estudios_3)
matriz6 <- matriz6[["COAD_Methylation_methyl450-20160128"]]
dim(matriz6)
matriz6[1:3, 1:3]
matriz6 <- as.matrix(matriz6)
fenotipos <- colData(estudios_3)
etiqueta6 <- data.frame(bar_code = colnames(matriz6))
etiqueta6$sujeto <- substring(colnames(matriz6), 1, 12)
etiqueta6$estudio <- "COAD" 
etiqueta6$categoria <- substring(colnames(matriz6), 14,15) %>% as.factor()

tipos6 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta6 <- merge(etiqueta6, tipos6, by="sujeto", all.x=TRUE)
etiqueta6$type[etiqueta6$categoria == "11" ] <- "Control"
etiqueta6$type <- as.factor(etiqueta6$type)
etiqueta6$label <- etiqueta6$estudio
etiqueta6$label[etiqueta6$categoria == "11"] <- "Control"


#################################################################################
###################################################################################
# FUSIONADO DE ESTUDIOS
###############################################################################

load("C:/TFM UOC/datos/matriz")
load("C:/TFM UOC/datos/etiqueta")

sum(row.names(matriz) != row.names(matriz6))

etiqueta <- rbind(etiqueta, etiqueta6)
matriz <- cbind(matriz, matriz6)

save(matriz, file="./datos/matriz")
save(etiqueta, file="./datos/etiqueta")

rm(list=ls())

##################################################################################
####################################################################################

library(Rtsne)
library(ggplot2)

load("C:/TFM UOC/datos/matriz")
load("C:/TFM UOC/datos/etiqueta")

matriz_sna <- matriz[complete.cases(matriz), ]
dim(matriz_sna)
matriz_sna <- t(matriz_sna)
matriz_sna[1:3,1:3]
sum(is.na(matriz_sna))
dim(matriz_sna)

matriz_sna_n <- normalize_input(matriz_sna)
matriz_sna_n[1:3,1:3]

tsne <- Rtsne(matriz_sna_n, dims = 2, perplexity=50, verbose=TRUE, max_iter =500)
tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2], col = etiqueta$label)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))


















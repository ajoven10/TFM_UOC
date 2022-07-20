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
assay(estudios[[2]], "counts") %>% min(na.rm=TRUE)
assay(estudios_2[[1]], "counts") %>% max(na.rm=TRUE)
assay(estudios_2[[2]], "counts") %>% max(na.rm=TRUE)
assay(estudios_3[[1]], "counts") %>% max(na.rm=TRUE)
assay(estudios_3[[2]], "counts") %>% max(na.rm=TRUE)
assay(estudios_3[[3]], "counts") %>% max(na.rm=TRUE)
assay(estudios_3[[4]], "counts") %>% max(na.rm=TRUE)
assay(estudios_3[[5]], "counts") %>% max(na.rm=TRUE)
assay(estudios_3[[6]], "counts") %>% max(na.rm=TRUE)
assay(estudios_4[[1]], "counts") %>% max(na.rm=TRUE)
assay(estudios_4[[2]], "counts") %>% max(na.rm=TRUE)
assay(estudios_4[[3]], "counts") %>% max(na.rm=TRUE)
assay(estudios_4[[4]], "counts") %>% max(na.rm=TRUE)
assay(estudios_4[[5]], "counts") %>% max(na.rm=TRUE)
assay(estudios_4[[6]], "counts") %>% max(na.rm=TRUE)

assay(estudios_5[[1]], "counts") %>% max(na.rm=TRUE)
assay(estudios_5[[2]], "counts") %>% min(na.rm=TRUE)
assay(estudios_5[[3]], "counts") %>% max(na.rm=TRUE)
assay(estudios_5[[4]], "counts") %>% max(na.rm=TRUE)
assay(estudios_5[[5]], "counts") %>% max(na.rm=TRUE)
assay(estudios_5[[6]], "counts") %>% max(na.rm=TRUE)
assay(estudios_6[[1]], "counts") %>% max(na.rm=TRUE)
assay(estudios_6[[2]], "counts") %>% max(na.rm=TRUE)
assay(estudios_6[[3]], "counts") %>% max(na.rm=TRUE)
assay(estudios_6[[4]], "counts") %>% max(na.rm=TRUE)
assay(estudios_6[[5]], "counts") %>% max(na.rm=TRUE)
assay(estudios_6[[6]], "counts") %>% max(na.rm=TRUE)
assay(estudios_6[[7]], "counts") %>% max(na.rm=TRUE)
assay(estudios_6[[8]], "counts") %>% max(na.rm=TRUE)
assay(estudios_6[[9]], "counts") %>% max(na.rm=TRUE)
assay(estudios_7[[1]], "counts") %>% max(na.rm=TRUE)
assay(estudios_7[[2]], "counts") %>% max(na.rm=TRUE)
assay(estudios_7[[3]], "counts") %>% max(na.rm=TRUE)

assay(estudios_7[[4]], "counts") %>% max(na.rm=TRUE)
assay(estudios_7[[5]], "counts") %>% max(na.rm=TRUE)
assay(estudios_8[[1]], "counts") %>% max(na.rm=TRUE)
assay(estudios_8[[2]], "counts") %>% max(na.rm=TRUE)
assay(estudios_8[[3]], "counts") %>% max(na.rm=TRUE)
assay(estudios_8[[4]], "counts") %>% max(na.rm=TRUE)
assay(estudios_8[[5]], "counts") %>% max(na.rm=TRUE)
assay(estudios_8[[6]], "counts") %>% max(na.rm=TRUE)
assay(estudios_8[[7]], "counts") %>% max(na.rm=TRUE)
assay(estudios_9[[1]], "counts") %>% max(na.rm=TRUE)

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
annotation.table
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




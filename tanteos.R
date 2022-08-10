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
#############################################################################
library(dplyr)
load(file="./estudios/estudios_3")
experiments(estudios_3)

matriz7 <- assays(estudios_3)
matriz7 <- matriz7[["DLBC_Methylation-20160128"]]
dim(matriz7)
matriz7[1:3, 1:3]
matriz7 <- as.matrix(matriz7)
fenotipos <- colData(estudios_3)
etiqueta7 <- data.frame(bar_code = colnames(matriz7))
etiqueta7$sujeto <- substring(colnames(matriz7), 1, 12)
etiqueta7$estudio <- "DLBC" 
etiqueta7$categoria <- substring(colnames(matriz7), 14,15) %>% as.factor()

tipos7 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta7 <- merge(etiqueta7, tipos7, by="sujeto", all.x=TRUE)
etiqueta7$type[etiqueta7$categoria == "11" ] <- "Control"
etiqueta7$type <- as.factor(etiqueta7$type)
etiqueta7$label <- etiqueta7$estudio
etiqueta7$label[etiqueta7$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_3")
experiments(estudios_3)

matriz8 <- assays(estudios_3)
matriz8 <- matriz8[["ESCA_Methylation-20160128"]]
dim(matriz8)
matriz8[1:3, 1:3]
matriz8 <- as.matrix(matriz8)
fenotipos <- colData(estudios_3)
etiqueta8 <- data.frame(bar_code = colnames(matriz8))
etiqueta8$sujeto <- substring(colnames(matriz8), 1, 12)
etiqueta8$estudio <- "ESCA" 
etiqueta8$categoria <- substring(colnames(matriz8), 14,15) %>% as.factor()

tipos8 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta8 <- merge(etiqueta8, tipos8, by="sujeto", all.x=TRUE)
etiqueta8$type[etiqueta8$categoria == "11" ] <- "Control"
etiqueta8$type <- as.factor(etiqueta8$type)
etiqueta8$label <- etiqueta8$estudio
etiqueta8$label[etiqueta8$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_4")
experiments(estudios_4)

matriz9 <- assays(estudios_4)
matriz9 <- matriz9[["GBM_Methylation_methyl450-20160128"]]
dim(matriz9)
matriz9 <- as.matrix(matriz9)
fenotipos <- colData(estudios_4)
etiqueta9 <- data.frame(bar_code = colnames(matriz9))
etiqueta9$sujeto <- substring(colnames(matriz9), 1, 12)
etiqueta9$estudio <- "GBM" 
etiqueta9$categoria <- substring(colnames(matriz9), 14,15) %>% as.factor()

tipos9 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta9 <- merge(etiqueta9, tipos9, by="sujeto", all.x=TRUE)
etiqueta9$type[etiqueta9$categoria == "11" ] <- "Control"
etiqueta9$type <- as.factor(etiqueta9$type)
etiqueta9$label <- etiqueta9$estudio
etiqueta9$label[etiqueta9$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_4")
experiments(estudios_4)

matriz10 <- assays(estudios_4)
matriz10 <- matriz10[["HNSC_Methylation-20160128"]]
dim(matriz10)
matriz10 <- as.matrix(matriz10)
fenotipos <- colData(estudios_4)
etiqueta10 <- data.frame(bar_code = colnames(matriz10))
etiqueta10$sujeto <- substring(colnames(matriz10), 1, 12)
etiqueta10$estudio <- "HNSC" 
etiqueta10$categoria <- substring(colnames(matriz10), 14,15) %>% as.factor()

tipos10 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta10 <- merge(etiqueta10, tipos10, by="sujeto", all.x=TRUE)
etiqueta10$type[etiqueta10$categoria == "11" ] <- "Control"
etiqueta10$type <- as.factor(etiqueta10$type)
etiqueta10$label <- etiqueta10$estudio
etiqueta10$label[etiqueta10$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_4")
experiments(estudios_4)

matriz11 <- assays(estudios_4)
matriz11 <- matriz11[["KICH_Methylation-20160128"]]
dim(matriz11)
matriz11 <- as.matrix(matriz11)
fenotipos <- colData(estudios_4)
etiqueta11 <- data.frame(bar_code = colnames(matriz11))
etiqueta11$sujeto <- substring(colnames(matriz11), 1, 12)
etiqueta11$estudio <- "KICH" 
etiqueta11$categoria <- substring(colnames(matriz11), 14,15) %>% as.factor()

tipos11 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta11 <- merge(etiqueta11, tipos11, by="sujeto", all.x=TRUE)
etiqueta11$type[etiqueta11$categoria == "11" ] <- "Control"
etiqueta11$type <- as.factor(etiqueta11$type)
etiqueta11$label <- etiqueta11$estudio
etiqueta11$label[etiqueta11$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_4")
experiments(estudios_4)

matriz12 <- assays(estudios_4)
matriz12 <- matriz12[["KIRC_Methylation_methyl450-20160128"]]
dim(matriz12)
matriz12 <- as.matrix(matriz12)
fenotipos <- colData(estudios_4)
etiqueta12 <- data.frame(bar_code = colnames(matriz12))
etiqueta12$sujeto <- substring(colnames(matriz12), 1, 12)
etiqueta12$estudio <- "KIRC" 
etiqueta12$categoria <- substring(colnames(matriz12), 14,15) %>% as.factor()

tipos12 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta12 <- merge(etiqueta12, tipos12, by="sujeto", all.x=TRUE)
etiqueta12$type[etiqueta12$categoria == "11" ] <- "Control"
etiqueta12$type <- as.factor(etiqueta12$type)
etiqueta12$label <- etiqueta12$estudio
etiqueta12$label[etiqueta12$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_5")
experiments(estudios_5)

matriz13 <- assays(estudios_5)
matriz13 <- matriz13[["KIRP_Methylation_methyl450-20160128"]]
dim(matriz13)
matriz13 <- as.matrix(matriz13)
fenotipos <- colData(estudios_5)
etiqueta13 <- data.frame(bar_code = colnames(matriz13))
etiqueta13$sujeto <- substring(colnames(matriz13), 1, 12)
etiqueta13$estudio <- "KIRP" 
etiqueta13$categoria <- substring(colnames(matriz13), 14,15) %>% as.factor()

tipos13 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta13 <- merge(etiqueta13, tipos13, by="sujeto", all.x=TRUE)
etiqueta13$type[etiqueta13$categoria == "11" ] <- "Control"
etiqueta13$type <- as.factor(etiqueta13$type)
etiqueta13$label <- etiqueta13$estudio
etiqueta13$label[etiqueta13$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_5")
experiments(estudios_5)

matriz14 <- assays(estudios_5)
matriz14 <- matriz14[["LAML_Methylation_methyl450-20160128"]]
dim(matriz14)
matriz14 <- as.matrix(matriz14)
fenotipos <- colData(estudios_5)
etiqueta14 <- data.frame(bar_code = colnames(matriz14))
etiqueta14$sujeto <- substring(colnames(matriz14), 1, 12)
etiqueta14$estudio <- "LAML" 
etiqueta14$categoria <- substring(colnames(matriz14), 14,15) %>% as.factor()

tipos14 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta14 <- merge(etiqueta14, tipos14, by="sujeto", all.x=TRUE)
etiqueta14$type[etiqueta14$categoria == "11" ] <- "Control"
etiqueta14$type <- as.factor(etiqueta14$type)
etiqueta14$label <- etiqueta14$estudio
etiqueta14$label[etiqueta14$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_5")
experiments(estudios_5)

matriz15 <- assays(estudios_5)
matriz15 <- matriz15[["LGG_Methylation-20160128"]]
dim(matriz15)
matriz15 <- as.matrix(matriz15)
fenotipos <- colData(estudios_5)
etiqueta15 <- data.frame(bar_code = colnames(matriz15))
etiqueta15$sujeto <- substring(colnames(matriz15), 1, 12)
etiqueta15$estudio <- "LGG" 
etiqueta15$categoria <- substring(colnames(matriz15), 14,15) %>% as.factor()

tipos15 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta15 <- merge(etiqueta15, tipos15, by="sujeto", all.x=TRUE)
etiqueta15$type[etiqueta15$categoria == "11" ] <- "Control"
etiqueta15$type <- as.factor(etiqueta15$type)
etiqueta15$label <- etiqueta15$estudio
etiqueta15$label[etiqueta15$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_5")
experiments(estudios_5)

matriz16 <- assays(estudios_5)
matriz16 <- matriz16[["LIHC_Methylation-20160128"]]
dim(matriz16)
matriz16 <- as.matrix(matriz16)
fenotipos <- colData(estudios_5)
etiqueta16 <- data.frame(bar_code = colnames(matriz16))
etiqueta16$sujeto <- substring(colnames(matriz16), 1, 12)
etiqueta16$estudio <- "LIHC" 
etiqueta16$categoria <- substring(colnames(matriz16), 14,15) %>% as.factor()

tipos16 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta16 <- merge(etiqueta16, tipos16, by="sujeto", all.x=TRUE)
etiqueta16$type[etiqueta16$categoria == "11" ] <- "Control"
etiqueta16$type <- as.factor(etiqueta16$type)
etiqueta16$label <- etiqueta16$estudio
etiqueta16$label[etiqueta16$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_6")
experiments(estudios_6)

matriz17 <- assays(estudios_6)
matriz17 <- matriz17[["LUAD_Methylation_methyl450-20160128"]]
dim(matriz17)
matriz17 <- as.matrix(matriz17)
fenotipos <- colData(estudios_6)
etiqueta17 <- data.frame(bar_code = colnames(matriz17))
etiqueta17$sujeto <- substring(colnames(matriz17), 1, 12)
etiqueta17$estudio <- "LUAD" 
etiqueta17$categoria <- substring(colnames(matriz17), 14,15) %>% as.factor()

tipos17 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta17 <- merge(etiqueta17, tipos17, by="sujeto", all.x=TRUE)
etiqueta17$type[etiqueta17$categoria == "11" ] <- "Control"
etiqueta17$type <- as.factor(etiqueta17$type)
etiqueta17$label <- etiqueta17$estudio
etiqueta17$label[etiqueta17$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_6")
experiments(estudios_6)

matriz18 <- assays(estudios_6)
matriz18 <- matriz18[["LUSC_Methylation_methyl450-20160128"]]
dim(matriz18)
matriz18 <- as.matrix(matriz18)
fenotipos <- colData(estudios_6)
etiqueta18 <- data.frame(bar_code = colnames(matriz18))
etiqueta18$sujeto <- substring(colnames(matriz18), 1, 12)
etiqueta18$estudio <- "LUSC" 
etiqueta18$categoria <- substring(colnames(matriz18), 14,15) %>% as.factor()

tipos18 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta18 <- merge(etiqueta18, tipos18, by="sujeto", all.x=TRUE)
etiqueta18$type[etiqueta18$categoria == "11" ] <- "Control"
etiqueta18$type <- as.factor(etiqueta18$type)
etiqueta18$label <- etiqueta18$estudio
etiqueta18$label[etiqueta18$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_6")
experiments(estudios_6)

matriz19 <- assays(estudios_6)
matriz19 <- matriz19[["MESO_Methylation-20160128"]]
dim(matriz19)
matriz19 <- as.matrix(matriz19)
fenotipos <- colData(estudios_6)
etiqueta19 <- data.frame(bar_code = colnames(matriz19))
etiqueta19$sujeto <- substring(colnames(matriz19), 1, 12)
etiqueta19$estudio <- "MESO" 
etiqueta19$categoria <- substring(colnames(matriz19), 14,15) %>% as.factor()

tipos19 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta19 <- merge(etiqueta19, tipos19, by="sujeto", all.x=TRUE)
etiqueta19$type[etiqueta19$categoria == "11" ] <- "Control"
etiqueta19$type <- as.factor(etiqueta19$type)
etiqueta19$label <- etiqueta19$estudio
etiqueta19$label[etiqueta19$categoria == "11"] <- "Control"

#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_6")
experiments(estudios_6)

matriz20 <- assays(estudios_6)
matriz20 <- matriz20[["OV_Methylation_methyl450-20160128"]]
dim(matriz20)
matriz20 <- as.matrix(matriz20)
fenotipos <- colData(estudios_6)
etiqueta20 <- data.frame(bar_code = colnames(matriz20))
etiqueta20$sujeto <- substring(colnames(matriz20), 1, 12)
etiqueta20$estudio <- "OV" 
etiqueta20$categoria <- substring(colnames(matriz20), 14,15) %>% as.factor()

tipos20 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta20 <- merge(etiqueta20, tipos20, by="sujeto", all.x=TRUE)
etiqueta20$type[etiqueta20$categoria == "11" ] <- "Control"
etiqueta20$type <- as.factor(etiqueta20$type)
etiqueta20$label <- etiqueta20$estudio
etiqueta20$label[etiqueta20$categoria == "11"] <- "Control"


#################################################################################
#############################################################################
library(dplyr)
load(file="./estudios/estudios_6")
experiments(estudios_6)

matriz21 <- assays(estudios_6)
matriz21 <- matriz21[["PAAD_Methylation-20160128"]]
dim(matriz21)
matriz21 <- as.matrix(matriz21)
fenotipos <- colData(estudios_6)
etiqueta21 <- data.frame(bar_code = colnames(matriz21))
etiqueta21$sujeto <- substring(colnames(matriz21), 1, 12)
etiqueta21$estudio <- "PAAD" 
etiqueta21$categoria <- substring(colnames(matriz21), 14,15) %>% as.factor()

tipos21 <- data.frame(sujeto = rownames(fenotipos), type = fenotipos$histological_type)
etiqueta21 <- merge(etiqueta21, tipos21, by="sujeto", all.x=TRUE)
etiqueta21$type[etiqueta21$categoria == "11" ] <- "Control"
etiqueta21$type <- as.factor(etiqueta21$type)
etiqueta21$label <- etiqueta21$estudio
etiqueta21$label[etiqueta21$categoria == "11"] <- "Control"


#################################################################################
###################################################################################
# FUSIONADO DE ESTUDIOS
###############################################################################

load("C:/TFM UOC/datos/matriz")
load("C:/TFM UOC/datos/etiqueta")

sum(row.names(matriz) != row.names(matriz21))

etiqueta <- rbind(etiqueta, etiqueta20, etiqueta21)
matriz <- cbind(matriz, matriz21)

dim(matriz)
save(matriz, file="./datos/matriz")
save(etiqueta, file="./datos/etiqueta")



##################################################################################
####################################################################################

library(Rtsne)
library(ggplot2)

load("C:/TFM UOC/datos/matriz")
load("C:/TFM UOC/datos/etiqueta")

matriz_sna <- na.omit(matriz)
dim(matriz_sna)
matriz_sna <- t(matriz_sna)
matriz_sna[1:3,1:3]
sum(is.na(matriz_sna))
dim(matriz_sna)

matriz_sna_n <- normalize_input(matriz_sna)
matriz_sna_n[1:3,1:3]
rm(matriz_sna)

tsne <- Rtsne(matriz_sna_n, partial_pca=TRUE, dims = 2, perplexity=30, verbose=TRUE, max_iter =500)
tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2], col = etiqueta$categoria)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))


      

















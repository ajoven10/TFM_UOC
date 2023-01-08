library(dplyr)
library(curatedTCGAData)
library(TCGAutils)

###############################################
# se extraen las matrices con los valores betas de cada SummarizedExperiment
# y se unen mediante cbind
#################################################

load(file="G:/TFM UOC/estudios/estudios_1")
experiments(estudios)

matriz1 <- assays(estudios)
matriz1 <- matriz1[["BLCA_Methylation-20160128"]]
matriz1 <- as.matrix(matriz1)

matriz2 <- assays(estudios)
matriz2 <- matriz2[["ACC_Methylation-20160128"]]
matriz2 <- as.matrix(matriz2)
##############################################

load(file="G:/TFM UOC/estudios/estudios_2")
experiments(estudios_2)

matriz3 <- assays(estudios_2)
matriz3 <- matriz3[["BRCA_Methylation_methyl450-20160128"]]
matriz3 <- as.matrix(matriz3)
####################################################

load(file="G:/TFM UOC/estudios/estudios_3")
experiments(estudios_3)

matriz4 <- assays(estudios_3)
matriz4 <- matriz4[["CESC_Methylation-20160128"]]
matriz4 <- as.matrix(matriz4)

matriz5 <- assays(estudios_3)
matriz5 <- matriz5[["CHOL_Methylation-20160128"]]
matriz5 <- as.matrix(matriz5)

matriz6 <- assays(estudios_3)
matriz6 <- matriz6[["COAD_Methylation_methyl450-20160128"]]
matriz6 <- as.matrix(matriz6)

matriz7 <- assays(estudios_3)
matriz7 <- matriz7[["DLBC_Methylation-20160128"]]
matriz7 <- as.matrix(matriz7)

matriz8 <- assays(estudios_3)
matriz8 <- matriz8[["ESCA_Methylation-20160128"]]
matriz8 <- as.matrix(matriz8)
########################################################

load(file="G:/TFM UOC/estudios/estudios_4")
experiments(estudios_4)

matriz9 <- assays(estudios_4)
matriz9 <- matriz9[["GBM_Methylation_methyl450-20160128"]]
matriz9 <- as.matrix(matriz9)

matriz10 <- assays(estudios_4)
matriz10 <- matriz10[["HNSC_Methylation-20160128"]]
matriz10 <- as.matrix(matriz10)

matriz11 <- assays(estudios_4)
matriz11 <- matriz11[["KICH_Methylation-20160128"]]
matriz11 <- as.matrix(matriz11)

matriz12 <- assays(estudios_4)
matriz12 <- matriz12[["KIRC_Methylation_methyl450-20160128"]]
matriz12 <- as.matrix(matriz12)
#########################################################

load(file="G:/TFM UOC/estudios/estudios_5")
experiments(estudios_5)

matriz13 <- assays(estudios_5)
matriz13 <- matriz13[["KIRP_Methylation_methyl450-20160128"]]
matriz13 <- as.matrix(matriz13)

matriz14 <- assays(estudios_5)
matriz14 <- matriz14[["LAML_Methylation_methyl450-20160128"]]
matriz14 <- as.matrix(matriz14)

matriz15 <- assays(estudios_5)
matriz15 <- matriz15[["LGG_Methylation-20160128"]]
matriz15 <- as.matrix(matriz15)

matriz16 <- assays(estudios_5)
matriz16 <- matriz16[["LIHC_Methylation-20160128"]]
matriz16 <- as.matrix(matriz16)
##############################################################

load(file="G:/TFM UOC/estudios/estudios_6")
experiments(estudios_6)

matriz17 <- assays(estudios_6)
matriz17 <- matriz17[["LUAD_Methylation_methyl450-20160128"]]
matriz17 <- as.matrix(matriz17)

matriz18 <- assays(estudios_6)
matriz18 <- matriz18[["LUSC_Methylation_methyl450-20160128"]]
matriz18 <- as.matrix(matriz18)

matriz19 <- assays(estudios_6)
matriz19 <- matriz19[["MESO_Methylation-20160128"]]
matriz19 <- as.matrix(matriz19)

matriz20 <- assays(estudios_6)
matriz20 <- matriz20[["OV_Methylation_methyl450-20160128"]]
matriz20 <- as.matrix(matriz20)

matriz21 <- assays(estudios_6)
matriz21 <- matriz21[["PAAD_Methylation-20160128"]]
matriz21 <- as.matrix(matriz21)

matriz22 <- assays(estudios_6)
matriz22 <- matriz22[["PCPG_Methylation-20160128"]]
matriz22 <- as.matrix(matriz22)
######################################################

load(file="G:/TFM UOC/estudios/estudios_7")
experiments(estudios_7)

matriz23 <- assays(estudios_7)
matriz23 <- matriz23[["PRAD_Methylation-20160128"]]
matriz23 <- as.matrix(matriz23)

matriz24 <- assays(estudios_7)
matriz24 <- matriz24[["READ_Methylation_methyl450-20160128"]]
matriz24 <- as.matrix(matriz24)

matriz25 <- assays(estudios_7)
matriz25 <- matriz25[["SARC_Methylation-20160128"]]
matriz25 <- as.matrix(matriz25)

matriz26 <- assays(estudios_7)
matriz26 <- matriz26[["SKCM_Methylation-20160128"]]
matriz26 <- as.matrix(matriz26)
##########################################################

load(file="G:/TFM UOC/estudios/estudios_8")
experiments(estudios_8)

matriz27 <- assays(estudios_8)
matriz27 <- matriz27[["TGCT_Methylation-20160128"]]
matriz27 <- as.matrix(matriz27)

matriz28 <- assays(estudios_8)
matriz28 <- matriz28[["THCA_Methylation-20160128"]]
matriz28 <- as.matrix(matriz28)

matriz29 <- assays(estudios_8)
matriz29 <- matriz29[["THYM_Methylation-20160128"]]
matriz29 <- as.matrix(matriz29)

matriz30 <- assays(estudios_8)
matriz30 <- matriz30[["UCEC_Methylation_methyl450-20160128"]]
matriz30 <- as.matrix(matriz30)

matriz31 <- assays(estudios_8)
matriz31 <- matriz31[["UCS_Methylation-20160128"]]
matriz31 <- as.matrix(matriz31)

matriz32 <- assays(estudios_8)
matriz32 <- matriz32[["UVM_Methylation-20160128"]]
matriz32 <- as.matrix(matriz32)
################################################
###################################################

##########################################################

load(file="G:/TFM UOC/estudios/estudios_10")
experiments(estudios_10)

matriz33 <- assays(estudios_10)
matriz33 <- matriz33[["STAD_Methylation_methyl450-20160128"]]
matriz33 <- as.matrix(matriz33)


# FUSIONADO DE ESTUDIOS
###############################################################################

rm(estudios, estudios_2, estudios_3, estudios_4, estudios_5, estudios_6,
   estudios_7, estudios_8, estudios_10)

matriz <- cbind(matriz1, matriz2)
rm(matriz1, matriz2)

matriz <- cbind(matriz, matriz3)
rm(matriz3)

matriz <- cbind(matriz, matriz4)
rm(matriz4)

matriz <- cbind(matriz, matriz5)
rm(matriz5)

matriz <- cbind(matriz, matriz6)
rm(matriz6)

matriz <- cbind(matriz, matriz7)
rm(matriz7)
matriz <- cbind(matriz, matriz8)
rm(matriz8)
matriz <- cbind(matriz, matriz9)
rm(matriz9)
matriz <- cbind(matriz, matriz10)
rm(matriz10)
matriz <- cbind(matriz, matriz11)
rm(matriz11)
matriz <- cbind(matriz, matriz12)
rm(matriz12)
matriz <- cbind(matriz, matriz13)
rm(matriz13)
matriz <- cbind(matriz, matriz14)
rm(matriz14)
matriz <- cbind(matriz, matriz15)
rm(matriz15)
matriz <- cbind(matriz, matriz16)
rm(matriz16)
matriz <- cbind(matriz, matriz17)
rm(matriz17)
matriz <- cbind(matriz, matriz18)
rm(matriz18)
matriz <- cbind(matriz, matriz19)
rm(matriz19)
matriz <- cbind(matriz, matriz20)
rm(matriz20)
matriz <- cbind(matriz, matriz21)
rm(matriz21)
matriz <- cbind(matriz, matriz22)
rm(matriz22)
matriz <- cbind(matriz, matriz23)
rm(matriz23)
matriz <- cbind(matriz, matriz24)
rm(matriz24)
matriz <- cbind(matriz, matriz25)
rm(matriz25)
matriz <- cbind(matriz, matriz26)
rm(matriz26)
matriz <- cbind(matriz, matriz27)
rm(matriz27)
matriz <- cbind(matriz, matriz28)
rm(matriz28)
matriz <- cbind(matriz, matriz29)
rm(matriz29)
matriz <- cbind(matriz, matriz30)
rm(matriz30)
matriz <- cbind(matriz, matriz31)
rm(matriz31)
matriz <- cbind(matriz, matriz32)
rm(matriz32)
matriz <- cbind(matriz, matriz33)
rm(matriz33)



dim(matriz)
save(matriz, file="G:/TFM UOC/datos/matriz.Rda")





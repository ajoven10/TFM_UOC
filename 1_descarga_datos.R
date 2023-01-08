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

save(estudios, file= "G:/TFM UOC/estudios/estudios_1")
estudios

################################################################################
###############################################################################
estudios_2 <- curatedTCGAData(
  diseaseCode = c("BRCA"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_2)

save(estudios_2, file= "G:/TFM UOC/estudios/estudios_2")
estudios_2

##############################################################################
###############################################################################
estudios_3 <- curatedTCGAData(
  diseaseCode = c("CESC", "CHOL","COAD","DLBC","ESCA"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_3)

save(estudios_3, file= "G:/TFM UOC/estudios/estudios_3")
estudios_3

#############################################################################
###############################################################################
estudios_4 <- curatedTCGAData(
  diseaseCode = c("GBM", "HNSC", "KICH", "KIRC"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_4)

save(estudios_4, file= "G:/TFM UOC/estudios/estudios_4")
estudios_4

#############################################################################
###############################################################################
estudios_5 <- curatedTCGAData(
  diseaseCode = c("KIRP", "LAML", "LGG", "LIHC"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_5)

save(estudios_5, file= "G:/TFM UOC/estudios/estudios_5")
estudios_5

#############################################################################
###############################################################################
estudios_6 <- curatedTCGAData(
  diseaseCode = c("LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_6)

save(estudios_6, file= "G:/TFM UOC/estudios/estudios_6")
estudios_6

#############################################################################
###############################################################################
estudios_7 <- curatedTCGAData(
  diseaseCode = c("PRAD", "READ", "SARC", "SKCM"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_7)

save(estudios_7, file= "G:/TFM UOC/estudios/estudios_7")
estudios_7

#############################################################################
###############################################################################
estudios_8 <- curatedTCGAData(
  diseaseCode = c("TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"), assays = "Methyl*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_8)

save(estudios_8, file= "G:/TFM UOC/estudios/estudios_8")
estudios_8

#############################################################################
###############################################################################
estudios_9 <- curatedTCGAData(
  diseaseCode = c("STAD"), assays = "Methylation_methyl27*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_9, vial=TRUE)

save(estudios_9, file= "G:/TFM UOC/estudios/estudios_9")
estudios_9

#############################################################################
###############################################################################
estudios_10 <- curatedTCGAData(
  diseaseCode = c("STAD"), assays = "Methylation_methyl450*", version = "1.1.38", dry.run = FALSE
)
# el número de muestras
sampleTables(estudios_10, vial=TRUE)

save(estudios_10, file= "G:/TFM UOC/estudios/estudios_10")
estudios_10

#############################################################################


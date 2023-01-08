library(readr)

betas <- read_csv("G:/TFM UOC/datos/test_TCGA_train/betas_11.txt")
betas <- read_csv("G:/TFM UOC/datos/test_TCGA_train/betas_12.txt")
betas <- read_csv("G:/TFM UOC/datos/test_TCGA_train/betas_13.txt")
betas <- read_csv("G:/TFM UOC/datos/test_colonomics/betas1_11.txt")
betas <- read_csv("G:/TFM UOC/datos/test_colonomics/betas1_1_2.txt")
betas <- read_csv("G:/TFM UOC/datos/test_GSE203061/betas1_1.txt")
betas <- read_csv("G:/TFM UOC/datos/test_GSE127985/betas1_1.txt")

library(dplyr)
library(keras)

load("C:/Users/usuario/TFM/sondas_dep_1000.Rda")
sondas_1000 <- sondas_dep_1000

tipos_sondas <- rowData(sondas_1000)$channel %>%  as.factor()
label <- sondas_1000$label %>% factor()
sondas_1000 <- row.names(sondas_1000)

modelo <- load_model_hdf5("C:/Users/usuario/TFM/modelo_1000.h5")

sondas_c <- betas[ , 1]
betas <- as.matrix(betas[, -1])
row.names(betas) <- sondas_c$cpg
rm(sondas_c)
gc()

betas_s <- data.frame(row.names = sondas_1000)
betas <- merge(betas_s, betas, by.x=0, by.y=0, all.x=TRUE)
row.names(betas) <- betas$Row.names
betas <- betas[ -1]
id <- colnames(betas)[1]
betas <- betas[sondas_1000, ]

betas <- as.data.frame(betas)
row.names(betas) <- sondas_1000
colnames(betas) <- id

m <- mean(betas[ , 1], na.rm=TRUE)
betas[is.na(betas)] <- m


# betas = betas / 1000
betas <- t(betas)

prediccion2 <- modelo %>% predict(betas, verbose=0)
clase <- which.max(prediccion2)
prob_eleccion <- max(prediccion2)


anota <- rowRanges(sondas_dep_1000) %>%  as.data.frame()
GSE42409_family_soft <- read_delim("G:/TFM UOC/datos/anotacion adicional sondas 450k/GSE42409_family.soft.txt", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)

anotaciones <- merge(anota, GSE42409_family_soft, all.x=TRUE, by.x=0, by.y="ID")
row.names(anotaciones) <- anotaciones$Row.names
anotaciones <- anotaciones[sondas_1000, ]

#########################################
# sensibilidad clases sondas
####################################

HIL_CpG <- as.factor(anotaciones$HIL_CpG_class)
HI <- levels(HIL_CpG)


salida <- data.frame(Probe_type=HI, betas_to_0=NA, new_class_0=NA,
                     betas_to_1=NA, new_class_1=NA)


for (i in 1:length(HI)) {
  rec <- betas
  
  betas[HIL_CpG == HI[i]] <- 0
  prediccion_i <- modelo %>% predict(betas, verbose=0)
  salida[i, "betas_to_0"] <- prediccion_i[clase]
  eleccion <- which.max(prediccion_i) - 1
  salida[i, "new_class_0"] <- levels(label)[eleccion]

  betas[HIL_CpG == HI[i]] <- 1
  prediccion_i <- modelo %>% predict(betas, verbose=0)
  salida[i, "betas_to_1"] <- prediccion_i[clase]
  eleccion <- which.max(prediccion_i) - 1
  salida[i, "new_class_1"] <- levels(label)[eleccion]

  betas <- rec
}

###########################################
# sensibilidad variaciones sondas individuales
#######################################


salida2 <- data.frame(Probe_type=colnames(betas), betas_to_0=NA, new_class_0=NA,
                     betas_to_1=NA, new_class_1=NA)


for (i in 1:1000) {
  rec <- betas

  betas[1,i] <- 0
  prediccion_i <- modelo %>% predict(betas, verbose=0)
  salida2[i, "betas_to_0"] <- prediccion_i[clase]
  eleccion <- which.max(prediccion_i) - 1
  salida2[i, "new_class_0"] <- levels(label)[eleccion]

  betas[1,i] <- 1
  prediccion_i <- modelo %>% predict(betas, verbose=0)
  salida2[i, "betas_to_1"] <- prediccion_i[clase]
  eleccion <- which.max(prediccion_i) - 1
  salida2[i, "new_class_1"] <- levels(label)[eleccion]

  betas <- rec
}

###############################################
# Sensibilidad cromosoma
######################################

cromosomas <- as.factor(anotaciones$seqnames)
crom <- levels(cromosomas)


salida3 <- data.frame(Probe_type=crom, betas_to_0=NA, new_class_0=NA,
                     betas_to_1=NA, new_class_1=NA)


for (i in 1:length(crom)) {
  rec <- betas
  
  betas[cromosomas == crom[i]] <- 0
  prediccion_i <- modelo %>% predict(betas, verbose=0)
  salida3[i, "betas_to_0"] <- prediccion_i[clase]
  eleccion <- which.max(prediccion_i) - 1
  salida3[i, "new_class_0"] <- levels(label)[eleccion]
  
  betas[cromosomas == crom[i]] <- 1
  prediccion_i <- modelo %>% predict(betas, verbose=0)
  salida3[i, "betas_to_1"] <- prediccion_i[clase]
  eleccion <- which.max(prediccion_i) - 1
  salida3[i, "new_class_1"] <- levels(label)[eleccion]
  
  betas <- rec
}

#####################################################






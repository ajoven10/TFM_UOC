load("G:/TFM UOC/datos/Clasificador_variabilidad/sondas_all_all_sd_20000.Rda")
sondas_20000

########################################################################
# Selección de sondas a suprimir
#####################################################################

################################################
# Carga GSE42409 anotaciones adicionales de sondas
################################################
library(GEOquery)
library(dplyr)
library(SummarizedExperiment)
elist <- getGEO("GSE42409")
GSE42409 <- elist[[1]] %>%  featureData()
GSE42409_df <- GSE42409@data %>% as.data.frame()


###################################################
# Selección de las sondas a suprimir
#######################################
sondas <- data.frame(sondas = row.names(sondas_20000))
table(substring(sondas$sondas, 1,2))

##############################################
# se quitan las que no comienzan por "cg"
##############################################
sondas <- filter(sondas, substring(sondas, 1 ,2) =="cg")

#################################################
# se añaden a sondas las anotaciones de GSE42409
####################################################
sondas <- merge(sondas, GSE42409_df, by.x="sondas", by.y="ID", all.x=TRUE)

###################################################
# se suprimen las que son problemáticas o no hay datos
# en las anotaciones adicionales de GSE42409
################################################
sondas_dep <- filter(sondas, (`Target CpG SNP`=="" | is.na(`Target CpG SNP`)) &
                       AlleleA_Hits == 1 &
                       AlleleB_Hits == 0  &
                       n_bp_repetitive == 0) 

#####################################################
# se seleccionan las 1000 sondas más variables de entre
# las que han quedado
#################################################
sondas_dep_20000 <- sondas_20000[row.names(sondas_20000) %in% sondas_dep$sondas, ]
sondas_dep_20000

sd <- rowRanges(sondas_dep_20000)$sd
sd_o <- order(sd, decreasing=TRUE)
ss <- sd[sd_o][1:1000]
sondas_dep_1000 <- sondas_dep_20000[names(ss), ]
sondas_dep_1000

#######################################################
# se seleccionan los tumores primarios y muestras de control
#################################################
codigos <- c("01", "03", "11")
sondas_dep_1000 <- sondas_dep_1000[ , sondas_dep_1000$categoria %in% codigos]
sondas_dep_1000


###################################################
# se salva el objeto final R
###################################
save(sondas_dep_1000, file = "C:/Users/usuario/TFM/sondas_dep_1000.Rda")

row.names(sondas_dep) <- sondas_dep$sondas
sondas_dep <- sondas_dep[rownames(sondas_dep) %in%  names(rowRanges(sondas_dep_1000)), ]
sondas_dep <- sondas_dep[names(rowRanges(sondas_dep_1000)),]
save(sondas_dep, file = "C:/TFM UOC/MEMORIA/anotaciones_adic_sondas_dep_1000.Rda")
load(file = "C:/TFM UOC/MEMORIA/anotaciones_adic_sondas_dep_1000.Rda")

##########################################
# tablas resumen sondas
##########################################
datos_sondas <- rowRanges(sondas_dep_1000) %>%  as.data.frame()
table(datos_sondas$seqnames)
table(substring(names(rowRanges(sondas_dep_1000)), 1,2))
table(datos_sondas$channel)
summary(datos_sondas$sd)
table(datos_sondas$platform)
table(sondas_dep$HIL_CpG_class)




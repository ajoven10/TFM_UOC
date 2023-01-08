library(dplyr)
load("G:/TFM UOC/datos/matriz.Rda")

load(file="G:/TFM UOC/estudios/estudios_1")
load(file="G:/TFM UOC/estudios/estudios_2")
load(file="G:/TFM UOC/estudios/estudios_3")
load(file="G:/TFM UOC/estudios/estudios_4")
load(file="G:/TFM UOC/estudios/estudios_5")
load(file="G:/TFM UOC/estudios/estudios_6")
load(file="G:/TFM UOC/estudios/estudios_7")
load(file="G:/TFM UOC/estudios/estudios_8")
load(file="G:/TFM UOC/estudios/estudios_9")
load(file="G:/TFM UOC/estudios/estudios_10")

etiqueta <- data.frame(bar_code = colnames(matriz), sujeto=substring(colnames(matriz), 1, 12))
etiqueta$categoria <- substring(colnames(matriz), 14,15) %>% as.factor()

# para añadir el histological_type y el disease_code

fenotipos1 <- colData(estudios)
fenotipos2 <- colData(estudios_2)
fenotipos3 <- colData(estudios_3)
fenotipos4 <- colData(estudios_4)
fenotipos5 <- colData(estudios_5)
fenotipos6 <- colData(estudios_6)
fenotipos7 <- colData(estudios_7)
fenotipos8 <- colData(estudios_8)
fenotipos9 <- colData(estudios_9)
fenotipos10 <- colData(estudios_10)
tipos1 <- data.frame(sujeto = rownames(fenotipos1), type = fenotipos1$histological_type, disease_code=fenotipos1$admin.disease_code)
tipos2 <- data.frame(sujeto = rownames(fenotipos2), type = fenotipos2$histological_type, disease_code=fenotipos2$admin.disease_code)
tipos3 <- data.frame(sujeto = rownames(fenotipos3), type = fenotipos3$histological_type, disease_code=fenotipos3$admin.disease_code)
tipos4 <- data.frame(sujeto = rownames(fenotipos4), type = fenotipos4$histological_type, disease_code=fenotipos4$admin.disease_code)
tipos5 <- data.frame(sujeto = rownames(fenotipos5), type = fenotipos5$histological_type, disease_code=fenotipos5$admin.disease_code)
tipos6 <- data.frame(sujeto = rownames(fenotipos6), type = fenotipos6$histological_type, disease_code=fenotipos6$admin.disease_code)
tipos7 <- data.frame(sujeto = rownames(fenotipos7), type = fenotipos7$histological_type, disease_code=fenotipos7$admin.disease_code)
tipos8 <- data.frame(sujeto = rownames(fenotipos8), type = fenotipos8$histological_type, disease_code=fenotipos8$admin.disease_code)
tipos9 <- data.frame(sujeto = rownames(fenotipos9), type = fenotipos9$histological_type, disease_code=fenotipos9$admin.disease_code)
tipos10 <- data.frame(sujeto = rownames(fenotipos10), type = fenotipos10$histological_type, disease_code=fenotipos10$admin.disease_code)

tipos <- rbind(tipos1, tipos2, tipos3, tipos4, tipos5, tipos6, tipos7, tipos8, tipos9, tipos10)

etiqueta <- merge(etiqueta, tipos, by="sujeto", all.x=TRUE)

# Para añadir el estudio assay

maps1 <- as.data.frame(sampleMap(estudios))
maps2 <- as.data.frame(sampleMap(estudios_2))
maps3 <- as.data.frame(sampleMap(estudios_3))
maps4 <- as.data.frame(sampleMap(estudios_4))
maps5 <- as.data.frame(sampleMap(estudios_5))
maps6 <- as.data.frame(sampleMap(estudios_6))
maps7 <- as.data.frame(sampleMap(estudios_7))
maps8 <- as.data.frame(sampleMap(estudios_8))
maps9 <- as.data.frame(sampleMap(estudios_9))
maps10 <- as.data.frame(sampleMap(estudios_10))

tipos_2 <- rbind(maps1, maps2, maps3, maps4, maps5, maps6, maps7, maps8, maps9, maps10)
tipos_2$assay <- as.character(tipos_2$assay)
ensayos <- strsplit(tipos_2$assay, split="_")
ensayos <- lapply(ensayos, "[", 1) %>% unlist()
tipos_2$assay <- ensayos

# elimino duplicados en los indicadores de muestra
tipos_2 <- tipos_2[!duplicated(tipos_2$colname), c(1,3) ]

etiqueta <- merge(etiqueta, tipos_2, by.x="bar_code", by.y= "colname", all.x=TRUE)
rownames(etiqueta) <- etiqueta$bar_code

etiqueta$label <- etiqueta$assay
etiqueta$label[etiqueta$categoria == "11"] <- "Control"

# muestras puestas en el mismo orden que la matriz de expresión
etiqueta <- etiqueta[colnames(matriz), ]

table(etiqueta$assay, useNA="always")
table(etiqueta$label, useNA="always")

sum(colnames(matriz) != rownames(etiqueta))

sujetos <- unique(etiqueta$sujeto)
table(table(etiqueta$sujeto))

save(etiqueta, file="G:/TFM UOC/datos/etiqueta.Rda")
write.table(etiqueta, "G:/TFM UOC/datos/etiqueta.txt", sep="\t", row.names=TRUE, col.names=TRUE)

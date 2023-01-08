library(Rtsne)
library(ggplot2)

load("G:/TFM UOC/datos/matriz_sna.Rda")
load("G:/TFM UOC/datos/etiqueta.Rda")

dim(matriz_sna)
sum(colnames(matriz_sna) != rownames(etiqueta))

matriz_sna_t <- t(matriz_sna)
matriz_sna_t[1:3,1:3]
sum(is.na(matriz_sna_t))

save(matriz_sna_t, file="G:/TFM UOC/datos/matriz_sna_t.Rda")

rm(matriz_sna)

# se realiza el análisis Rtnse con una muestra aleatoria de 150000 sondas,
# si se usan todas da error "complete.cases Long vectors no soported".
sed.seed=123
tsne <- Rtsne(matriz_sna_t[ , sample(dim(matriz_sna_t)[2], 150000) ], partial_pca=TRUE, 
              dims = 2, perplexity=30, verbose=TRUE, max_iter = 1000)

save(tsne, file="G:/TFM UOC/datos/tsne.Rda")

# Gráfico por patologías
tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2], col = etiqueta$label)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))

# Gráfico de las muestras de control
tsne_plot <- data.frame(x = tsne$Y[etiqueta$label=="Control", 1],
                        y = tsne$Y[etiqueta$label=="Control", 2], col = etiqueta$assay[etiqueta$label=="Control"])
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))

# Gráfico por categorías 
tsne_plot <- data.frame(x = tsne$Y[, 1], y = tsne$Y[, 2], col = etiqueta$categoria)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))

# Gráfico con solo las muestras COAD y READ
tsne_plot <- data.frame(x = tsne$Y[etiqueta$label=="READ" | etiqueta$label=="COAD", 1],
                        y = tsne$Y[etiqueta$label=="READ" | etiqueta$label=="COAD", 2], 
                        col = etiqueta$assay[etiqueta$label=="READ" | etiqueta$label=="COAD"])
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))

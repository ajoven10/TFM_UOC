library(Rtsne)
library(ggplot2)
library(gridExtra)
library(dplyr)

load("C:/Users/usuario/TFM/sondas_dep_1000.Rda")
summary(rowRanges(sondas_dep_1000)$sd)

sed.seed=123
betas <- assay(sondas_dep_1000, "counts")
label <- sondas_dep_1000$label
table(label)
betas_t <- t(betas)
dim(betas_t)

tsne <- Rtsne(betas_t, partial_pca=TRUE, dims=2, verbose=FALSE)

tsne_plot <- data.frame(x = tsne$Y[ ,1], y = tsne$Y[ ,2], col = label)

# gráfico Rtnse con todas las muestras
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col), size=1)

# gráfico Rtnse con muestras control
tsne_plot_2 <- data.frame(x = tsne$Y[sondas_dep_1000$label == "Control" , 1], 
                        y = tsne$Y[sondas_dep_1000$label == "Control" , 2], 
                        col = sondas_dep_1000$assay[sondas_dep_1000$label == "Control"])

ggplot(tsne_plot_2) + geom_point(aes(x=x, y=y, color=col))

table( sondas_dep_1000$assay[sondas_dep_1000$label == "Control"])

# histograma con número de muestras por estudio

numero_muestras <- data.frame(ensayo=sondas_dep_1000$assay, label=sondas_dep_1000$label)
numero_muestras$label[numero_muestras$label != "Control"] <- "Tumor"

p <- ggplot(data=numero_muestras, aes(x=factor(ensayo), fill=label)) + geom_bar(stat="count") +
  labs(x="Estudios", y="Número de muestras") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))     
p


# box-plot valores betas por etiqueta.
box <- data.frame(betas_t, tumor = label)

ggplot(box, aes(x=tumor, y=rowMeans(betas_t[ , 1:1000]))) + geom_boxplot() +
  labs(x="Tejidos", y="Media de valores betas") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))   


# grafico violin valores betas por etiqueta.
box <- data.frame(betas_t, tumor = label)

ggplot(box, aes(x=tumor, y=rowMeans(betas_t[ , 1:1000]), fill=tumor)) + 
  geom_violin(trim=FALSE, draw_quantiles=c(0.25, 0.5, 0.75)) +
  labs(x="Tejidos", y="Media de valores betas") +
 # coord_flip() +
  geom_boxplot(width=0.15) +
#  scale_fill_brewer() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), legend.position="none")   


# histograma con número de sondas por cromosoma
datos_sondas <- rowRanges(sondas_dep_1000) %>%  as.data.frame()
cromos_df <- data.frame(crom=datos_sondas$seqnames)

p <- ggplot(data=cromos_df, aes(crom)) + geom_bar() +
  labs(x="Cromosomas", y="Número de sondas") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))     
p

# descriptiva valores sd
load("C:/TFM UOC/MEMORIA/sd_all_all.Rda")
summary(sd)

qq <- quantile(sd, probs = seq(0, 1, by= 0.01))[100] 
sd_df_095 <- sd[sd > qq]
sd_df_095 <- data.frame(sd = sd_df_095)

sd_1000 <- rowRanges(sondas_dep_1000)$sd
sd_df_1000 <- data.frame(sd=sd_1000)

p <- ggplot(data= sd_df_095, aes(sd)) + 
  geom_histogram(data=sd_df_095, binwidth=0.0005) +
  geom_histogram(data=sd_df_1000, binwidth = 0.0005, fill="red") + 
  labs(x="desviaciones típicas", y="Número de sondas")
  
p

# normalidad de las sondas elección de 4 al azar
l <- label[label == "BRCA"]

graficos <- list()
for (i in 1:4){
  
ss <-  sample(1:dim(betas_t)[2],1)
df <- data.frame(betas = betas_t[label =="BRCA" , ss], label = l )

g <- ggplot(data = df , aes(x=betas)) +
  geom_histogram(aes(y=..density.., fill = ..count..)) +
  scale_fill_gradient(low="#DCDCDC", high = "#7C7C7C") +
  stat_function(fun=dnorm, colour="#0C3D7D9F", args=list(mean=mean(df$betas), sd = sd(df$betas))) +
  ggtitle(colnames(betas_t)[ss]) + facet_grid( . ~ df$label)

graficos[[i]] = g
}
grid.arrange(graficos[[1]], graficos[[2]], graficos[[3]], graficos[[4]], ncol=2)


################################################################
# box plot estudios pareados tumor vs control
##############################################

 pareados <- colData(sondas_dep_1000)[ ,  c("assay", "label")] %>% as.data.frame()
 pareados$label[pareados$label != "Control"] <- "Tumor"


es <- c("BRCA", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "THCA")


  bb <- betas_t[pareados$assay %in% es, ]
  ll <- pareados[pareados$assay %in% es, ]
  df2 <- data.frame(bb, ll)

   g <- ggplot(df2, aes(x=label, y=rowMeans(df2[ , 1:1000]), fill = assay)) + 
    geom_violin(trim=FALSE, draw_quantiles=c(0.25, 0.5, 0.75)) +
    labs( y="Media de valores betas") +
    geom_boxplot(width=0.15) +
    scale_fill_brewer() +
    facet_wrap( ~ assay ) 
   
g

######################################################################
pareados <- colData(sondas_dep_1000)[ ,  c("assay", "label")] %>% as.data.frame()
pareados$label[pareados$label != "Control"] <- "Tumor"


es <- c("BRCA", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "THCA")


bb <- betas_t[pareados$assay %in% es, ]
ll <- pareados[pareados$assay %in% es, ]
df2 <- data.frame(bb, ll)

g <- ggplot(df2, aes(x=label, y=rowMeans(df2[ , 1:1000]), fill = assay)) + 
 
  labs( y="Media de valores betas") +
  geom_boxplot() +
  scale_fill_brewer() +
  facet_wrap( ~ assay ) 

g



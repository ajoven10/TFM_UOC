library(sva)
library(dplyr)
library(TCGAutils)
library(SummarizedExperiment)

load("G:/TFM UOC/datos/data_1524sondas.Rda")
edata <- assay(data_1524sondas, "counts")
pdata <- colData(data_1524sondas)

nombres <- assay(data_1524sondas, "counts") %>% colnames()
plate_id <- TCGAbiospec(nombres) 
plate_id <- plate_id$plate

pdata$plate_id <- plate_id

mod0 <- model.matrix( ~ 1, data = pdata)
mod <- model.matrix( ~ as.factor(label), data = pdata)

combat_edata <- ComBat(dat = edata, batch = plate_id,
                       mod = mod0, par.prior=TRUE)



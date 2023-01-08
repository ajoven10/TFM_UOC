# GNU Free Documentation License (GNU FDL)
# Copyright ©  2022 Alberto Joven Álvarez.
# Permission is granted to copy, distribute and/or modify this document under 
# the terms of the GNU Free Documentation License, Version 1.3 or any later
# version published by the Free Software Foundation; with no Invariant Sections, 
# no Front-Cover Texts, and no Back-Cover Texts. 
# A copy of the license is included in the section entitled "GNU Free Documentation License".


library(shiny)
library(dplyr)
library(impute)
library(readr)
# library(SummarizedExperiment)
library(keras)
library(TCGAutils)
data("diseaseCodes")

library(BiocManager)
options(repos = BiocManager::repositories())
options(shiny.maxRequestSize=2000*1024^2)

# modelo sondas depuradas tipos 1-2-3-3-4 y muestra 01-03-11
# load("./sondas_dep_1000.Rda")


#label <- sondas_dep_1000$label %>% factor()
load("./label.Rda")
#sondas_1000 <- row.names(sondas_dep_1000)
load("./sondas_1000.Rda")

# gc()
modelo <- load_model_hdf5("./modelo_1000.h5", compile=FALSE)
r <- data.frame(a="")

shinyServer(function(input, output, session) {
  
  data <- reactive({
    req(input$betas)
    read_csv(input$betas$datapath)
    }
  )
  
  etiquetas <- reactive({
    req(input$etiqueta)
    read_csv(input$etiqueta$datapath)
    
  })
  
  
  divide <- reactive({
    divide <- input$divide
  })
  
  confussion <- reactive({
    confussion <- input$confussion
  })
  
  output$table1 <- renderTable({
    if (confussion())  {
      betas <- data() 
      sondas_c <- betas[ , 1]
      betas <- as.matrix(betas[, -1])
      row.names(betas) <- sondas_c$cpg
      rm(sondas_c)
      # gc()
      betas_s <- data.frame(row.names = sondas_1000)
      betas <- merge(betas_s, betas, by.x=0, by.y=0, all.x=TRUE)
      row.names(betas) <- betas$Row.names
      betas <- betas[ -1]
      id <- colnames(betas)[1]
      betas <- betas[sondas_1000, ]
      
      if (is.null(dim(betas))) {
        betas <- as.data.frame(betas)
        row.names(betas) <- sondas_1000
        colnames(betas) <- id
        m <- mean(betas[ , 1], na.rm=TRUE)
        betas[is.na(betas)] <- m
      }
      else {
        betas <- impute.knn(as.matrix(betas))$data 
      }
      
      if (divide()) {betas = betas / 1000}
      
      betas <- t(betas)
      rm(betas_s)
      # gc()
      
      fenotipos <- etiquetas()
      fenotipos <- fenotipos[ , 1:2]
      fenotipos_s <- intersect(fenotipos$id_clx, rownames(betas))
      fenotipos_r <- fenotipos[fenotipos$id_clx %in% fenotipos_s,  ]
      
      if (dim(betas)[1] != 1) {betas <- betas[fenotipos_s, ]}
      
      prediccion <- modelo %>% predict(betas, verbose=0) %>%
        k_argmax() %>%
        as.array() %>% as.integer()
      
      l <- as.list(1:34)
      names(l) <- levels(label)
      f <- names(l)[prediccion]
      
      t <- table( prediccion = f , type = paste("truth class:", fenotipos_r$type))
      tdf <- as.data.frame.matrix(t)
      prediccion <- paste("predicted class:", row.names(tdf))
      tdf <- cbind(prediccion, tdf)
    }
    else {tdf=r}
    
    tdf
  }
  )
  
  #############################################################################
  output$table2 <- renderTable({
    if (!confussion()) {
      betas <- data()
      sondas_c <- betas[ , 1]
      betas <- as.matrix(betas[, -1])
      row.names(betas) <- sondas_c$cpg
      rm(sondas_c)
      # gc()
      betas_s <- data.frame(row.names = sondas_1000)
      betas <- merge(betas_s, betas, by.x=0, by.y=0, all.x=TRUE)
      row.names(betas) <- betas$Row.names
      betas <- betas[ -1]
      id <- colnames(betas)[1]
      betas <- betas[sondas_1000, ]
      
      if (is.null(dim(betas))) {
        betas <- as.data.frame(betas)
        row.names(betas) <- sondas_1000
        colnames(betas) <- id
        m <- mean(betas[ , 1], na.rm=TRUE)
        betas[is.na(betas)] <- m
      }
      else {
        betas <- impute.knn(as.matrix(betas))$data 
      }
      
      if (divide()) {betas = betas / 1000}
      
      betas <- t(betas)
      rm(betas_s)
      # gc()
      
      #####################################################
      prediccion <- modelo %>% predict(betas, verbose=0) %>%
        k_argmax() %>%
        as.array() %>% as.integer()
      
      prediccion2 <- modelo %>% predict(betas, verbose=0)
      
      l <- as.list(1:34)
      names(l) <- levels(label)
      f <- names(l)[prediccion]
      
      p <- apply(prediccion2, 1, FUN=max)
      
      df <- data.frame(id_clx = row.names(betas),
                       predict = f,
                       probability = p)
    }
    else{df = r}
    
    df
  }
  )
  
  output$table3 <- renderTable({
    betas <- data()
    if (dim(betas)[2] == 2){
        sondas_c <- betas[ , 1]
        betas <- as.matrix(betas[, -1])
        row.names(betas) <- sondas_c$cpg
        rm(sondas_c)
       # gc()
        betas_s <- data.frame(row.names = sondas_1000)
        betas <- merge(betas_s, betas, by.x=0, by.y=0, all.x=TRUE)
        row.names(betas) <- betas$Row.names
        betas <- betas[ -1]
        id <- colnames(betas)[1]
        betas <- betas[sondas_1000, ]
        
        if (is.null(dim(betas))) {
          betas <- as.data.frame(betas)
          row.names(betas) <- sondas_1000
          colnames(betas) <- id
          m <- mean(betas[ , 1], na.rm=TRUE)
          betas[is.na(betas)] <- m
        }
        else {
          betas <- impute.knn(as.matrix(betas))$data 
        }
        
        if (divide()) {betas = betas / 1000}
        
        betas <- t(betas)
        rm(betas_s)
        # gc()
        anotaciones <- read_csv("./HIL_class.txt")
        HIL_CpG <- as.factor(anotaciones$HIL_CpG_class)
        HI <- levels(HIL_CpG)
        
        prediccion2 <- modelo %>% predict(betas, verbose=0)
        clase <- which.max(prediccion2)
        prob_eleccion <- max(prediccion2)
        
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
        salida
    } 
    else {salida = r}
    salida
  }
  ) 
  output$table4 <- renderTable({
    betas <- data()
    if (dim(betas)[2] == 2){
      sondas_c <- betas[ , 1]
      betas <- as.matrix(betas[, -1])
      row.names(betas) <- sondas_c$cpg
      rm(sondas_c)
      # gc()
      betas_s <- data.frame(row.names = sondas_1000)
      betas <- merge(betas_s, betas, by.x=0, by.y=0, all.x=TRUE)
      row.names(betas) <- betas$Row.names
      betas <- betas[ -1]
      id <- colnames(betas)[1]
      betas <- betas[sondas_1000, ]
      
      if (is.null(dim(betas))) {
        betas <- as.data.frame(betas)
        row.names(betas) <- sondas_1000
        colnames(betas) <- id
        m <- mean(betas[ , 1], na.rm=TRUE)
        betas[is.na(betas)] <- m
      }
      else {
        betas <- impute.knn(as.matrix(betas))$data 
      }
      
      if (divide()) {betas = betas / 1000}
      
      betas <- t(betas)
      rm(betas_s)
      # gc()
      
      anotaciones <- read_csv("./HIL_class.txt")
      cromosomas <- as.factor(anotaciones$seqnames)
      crom <- levels(cromosomas)
      
      prediccion2 <- modelo %>% predict(betas, verbose=0)
      clase <- which.max(prediccion2)
      prob_eleccion <- max(prediccion2)
      
      salida3 <- data.frame(Cromosome=crom, betas_to_0=NA, new_class_0=NA,
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
      
    } 
    else {salida3 = r}
    salida3
    }
  )
  
  output$table5 <- renderPrint({
   anotaciones <- read_csv("./HIL_class.txt")
   table(anotaciones$HIL_CpG_class)
    
  })
  
  output$table6 <- renderPrint({
    anotaciones <- read_csv("./HIL_class.txt")
    table(anotaciones$seqnames)
  })
  
  output$table7 <- renderDataTable({
    anotaciones <- read_csv("./HIL_class.txt")
    anotaciones
  })

  output$table8 <- renderTable({
    codigos <- diseaseCodes[ diseaseCodes$Available=="Yes", c(1,4)]
    codigos

  })
} 
)

    
 


 
    

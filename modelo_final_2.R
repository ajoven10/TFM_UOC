library(keras)
library(caret)
library(TCGAutils)
library(sva)
library(SummarizedExperiment)
library(dplyr)

########################################################
# Carga de los datos
#######################################################
load("C:/Users/usuario/TFM/sondas_dep_1000.Rda")
sondas_dep_1000

train_data <- assay(sondas_dep_1000, "counts") %>%  t()
label <- sondas_dep_1000$label %>% factor()
label_c <- to_categorical(as.integer(label))


#####################################################
# Construcción del modelo
##############################################
build_model <- function() {
  model <- keras_model_sequential() %>% 
    layer_dense(units=64, activation="relu", input_shape=dim(train_data)[[2]]) %>% 
    layer_dense(units=32, activation="relu") %>% 
    layer_dense(units=64, activation = "relu") %>% 
    layer_dense(units=35, activation = "softmax")
  
  
  model %>%  compile(
    optimizer = "rmsprop",
    loss = "categorical_crossentropy",
    metrics = c("accuracy")
  )
}

build_model() %>% summary()

####################################################
# Desglose train y test
########################################################
set.seed(123)
etiqueta <- sondas_dep_1000$label

in_train <- createDataPartition(etiqueta, p=0.75, list=FALSE) %>% as.vector()

train <- sondas_dep_1000[ , in_train]
test <- sondas_dep_1000[ , -in_train]

train_data <- assay(train, "counts")  %>% t()
label <- train$label %>%  factor()
label_c <- to_categorical(as.integer(label))

################################################################################
# Validación cruzada
########################################################

k=4
folds <- createFolds(label, k)

num_epoch = 40
all_scores = c()

for (i in 1:k) {
  cat("procesing fold #", i, "\n")
  partial_data <- train_data[ -folds[[i]] , ]
  partial_train_label <- label_c[-folds[[i]] , ]
  
  val_data <- train_data[folds[[i]] , ]
  val_label <- label_c[folds[[i]],]
  
  model <- build_model()
  
  history <- model %>% fit(partial_data, 
                           partial_train_label,
                           epoch = num_epoch
  )
  
  results <- model %>%  evaluate(val_data, val_label)
  all_scores <- c(all_scores, results[2])
}

all_scores
mean(all_scores)

############################################
# Validación del modelo usando los datos test
###########################################
test_data <- assay(test, "counts") %>%  t()
fenotipos_test <- test$label %>% factor()
test_labels <- to_categorical(as.integer(fenotipos_test))

modelo <- build_model()
num_epoch = 40
history <- modelo %>%  fit(train_data,
                           label_c,
                           epoch=num_epoch,
                           validation_data=list(test_data, test_labels))
plot(history)
metrics <- modelo %>% evaluate(test_data, test_labels)
metrics

##############################################
# Evaluación final del modelo
############################################
prediccion <- modelo %>% predict(test_data) %>% 
  k_argmax() %>% 
  as.array() %>% as.integer()

l <- as.list(1:34)
names(l) <- levels(fenotipos_test)
f <- names(l)[prediccion]
lev <- levels(fenotipos_test)
prediccion <- factor(f, levels=lev)

c3 <- confusionMatrix(data=prediccion, reference=fenotipos_test)
c3$overall


###########################################
# Entrenamiento del modelo con todos los datos sin distinción train test
###########################################
train_data <- assay(sondas_dep_1000, "counts") %>%  t()
label <- sondas_dep_1000$label %>% factor()
label_c <- to_categorical(as.integer(label))

modelo <- build_model()

history <- modelo %>% fit(train_data , 
               label_c,
               epoch = num_epoch )
plot(history)


##############################################
# Se guarda el modelo final
####################################
# save_model_hdf5(modelo, "C:/Users/usuario/TFM/modelo_1000.h5")
modelo <- load_model_hdf5("C:/Users/usuario/TFM/modelo_1000.h5")


###############################################
###############################################
##############################################
# Ajuste valores betas con Combat: covariate: plate_id
##########################################

load("C:/Users/usuario/TFM/sondas_dep_1000.Rda")
sondas_dep_1000

edata <- assay(sondas_dep_1000, "counts")

nombres <- colnames(edata)
plate_id <- TCGAbiospec(nombres)
plate_id <- plate_id$plate

sondas_dep_1000$plate_id <- plate_id
pdata <- colData(sondas_dep_1000)


mod0 <- model.matrix( ~ 1, data=pdata)
mod <- model.matrix( ~ as.factor(label), data=pdata)

combat_edata <- ComBat(dat=edata, batch=plate_id,
                       mod=mod0, mean.only=TRUE, par.prior=TRUE)

assay(sondas_dep_1000, "counts") <- combat_edata 
save(sondas_dep_1000, file = "C:/Users/usuario/TFM/sondas_dep_1000_combat.Rda")

#############################################################
# Carga de los datos
#######################################################

train_data <- assay(sondas_dep_1000, "counts") %>%  t()
label <- sondas_dep_1000$label %>% factor()
label_c <- to_categorical(as.integer(label))

####################################################
# Desglose train y test
########################################################
set.seed(123)
etiqueta <- sondas_dep_1000$label

in_train <- createDataPartition(etiqueta, p=0.75, list=FALSE) %>% as.vector()

train <- sondas_dep_1000[ , in_train]
test <- sondas_dep_1000[ , -in_train]

train_data <- assay(train, "counts")  %>% t()
label <- train$label %>%  factor()
label_c <- to_categorical(as.integer(label))


################################################################################
# Validación cruzada
########################################################

k=4
folds <- createFolds(label, k)

num_epoch = 40
all_scores = c()

for (i in 1:k) {
  cat("procesing fold #", i, "\n")
  partial_data <- train_data[ -folds[[i]] , ]
  partial_train_label <- label_c[-folds[[i]] , ]
  
  val_data <- train_data[folds[[i]] , ]
  val_label <- label_c[folds[[i]],]
  
  model <- build_model()
  
  history <- model %>% fit(partial_data, 
                           partial_train_label,
                           epoch = num_epoch
  )
  
  results <- model %>%  evaluate(val_data, val_label)
  all_scores <- c(all_scores, results[2])
}

all_scores
mean(all_scores)

############################################
# Validación del modelo usando los datos test
###########################################
test_data <- assay(test, "counts") %>%  t()
fenotipos_test <- test$label %>% factor()
test_labels <- to_categorical(as.integer(fenotipos_test))

modelo <- build_model()
num_epoch = 40
history <- modelo %>%  fit(train_data,
                           label_c,
                           epoch=num_epoch,
                           validation_data=list(test_data, test_labels))
plot(history)
metrics <- modelo %>% evaluate(test_data, test_labels)
metrics

##############################################
# Evaluación final del modelo 
############################################
prediccion <- modelo %>% predict(test_data) %>% 
  k_argmax() %>% 
  as.array() %>% as.integer()

l <- as.list(1:34)
names(l) <- levels(fenotipos_test)
f <- names(l)[prediccion]
lev <- levels(fenotipos_test)
prediccion <- factor(f, levels=lev)

c3 <- confusionMatrix(data=prediccion, reference=fenotipos_test)
c3$overall

###########################################
# Entrenamiento del modelo con todos los datos sin distinción train test
###########################################
train_data <- assay(sondas_dep_1000, "counts") %>%  t()
label <- sondas_dep_1000$label %>% factor()
label_c <- to_categorical(as.integer(label))

modelo <- build_model()

history <- modelo %>% fit(train_data, 
                          label_c,
                          epoch = num_epoch )

plot(history)


##############################################
# Se guarda el modelo final
####################################
# save_model_hdf5(modelo, "C:/Users/usuario/TFM/modelo_1000_combat.h5")





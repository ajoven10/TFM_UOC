load("G:/TFM UOC/datos/matriz.Rda")
dim(matriz)

library(dplyr)

nas <- is.na(matriz) %>% rowSums()

# Sondas que tiene NAs en todas las observaciones: 89512
sum(nas == 9707)
# sondas con más del 10% de observaciones NAs: 90258:
sum(nas > 9707 * .1)

save(nas, file="G:/TFM UOC/datos/matrices intermedias nas/nas.Rda")

# porcentaje de nas estimados
nn <- sum(nas[nas<9707*.1]) 
(100 * nn) / (395319*9707)

# se eliminan las sondas con más de 10% de las observaciones NAs, quedan 395319:
matriz2 <- matriz[nas < 9707 * 0.1, ]
dim(matriz2)
rm(matriz)

save(matriz2, file="G:/TFM UOC/datos/matriz2.Rda")

# se imputan valores a los NAs restantes con la librería impute que usa el promedio knn:
# sondas en filas y muestras en las columnas
# se parte la matriz en 9 para poder ejecutar la función y luego se vuelve a unir
library(impute)

matriz3 <- matriz2[ , 1:1000]
matriz4 <- matriz2[ , 1001:2000]
matriz5 <- matriz2[ , 2001:3000]
matriz6 <- matriz2[ , 3001:4000]
matriz7 <- matriz2[ , 4001:5000]
matriz8 <- matriz2[ , 5001:6000]
matriz9 <- matriz2[ , 6001:7000]
matriz10 <- matriz2[ , 7001:8000]
matriz11 <- matriz2[ , 8001:9000]
matriz12 <- matriz2[ , 9001:9707]

rm(matriz2)

matriz33 <- impute.knn(matriz3)
save(matriz33, file="G:/TFM UOC/datos/matrices intermedias nas/matriz33")
rm(matriz3)

matriz44 <- impute.knn(matriz4)
save(matriz44, file="G:/TFM UOC/datos/matrices intermedias nas/matriz44")
rm(matriz4)

matriz55 <- impute.knn(matriz5)
save(matriz55, file="G:/TFM UOC/datos/matrices intermedias nas/matriz55")
rm(matriz5)

matriz66 <- impute.knn(matriz6)
save(matriz66, file="G:/TFM UOC/datos/matrices intermedias nas/matriz66")
rm(matriz6)

matriz77 <- impute.knn(matriz7)
save(matriz77, file="G:/TFM UOC/datos/matrices intermedias nas/matriz77")
rm(matriz7)

matriz88 <- impute.knn(matriz8)
save(matriz88, file="G:/TFM UOC/datos/matrices intermedias nas/matriz88")
rm(matriz8)

matriz99 <- impute.knn(matriz9)
save(matriz99, file="G:/TFM UOC/datos/matrices intermedias nas/matriz99")
rm(matriz9)

matriz1010 <- impute.knn(matriz10)
save(matriz1010, file="G:/TFM UOC/datos/matrices intermedias nas/matriz1010")
rm(matriz10)

matriz1111 <- impute.knn(matriz11)
save(matriz1111, file="G:/TFM UOC/datos/matrices intermedias nas/matriz1111")
rm(matriz11)

matriz1212 <- impute.knn(matriz12)
save(matriz1212, file="G:/TFM UOC/datos/matrices intermedias nas/matriz1212")
rm(matriz12)

matriz_sna <- cbind(matriz33$data, matriz44$data, matriz55$data, matriz66$data, matriz77$data, matriz88$data,
                    matriz99$data, matriz1010$data, matriz1111$data, matriz1212$data)

save(matriz_sna, file="G:/TFM UOC/datos/matriz_sna.Rda")


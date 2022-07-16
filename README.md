PREGUNTA INICIAL DE INTER�S BIOL�GICO:
�Es posible entrenar una red neuronal para que a partir de los datos de posiciones metiladas de una muestra sea capaz de identificar a qu� categor�a de tejido pertenece (clasificarla en una de las 20 categor�as de tejidos del TCGA)?.
Variaci�n de la pregunta:
* Lo mismo pero considerando un algoritmo random forest.
* Lo mismo pero considerando un algoritmo �rbol de decisi�n. Este algoritmo tiene la ventaja de ser f�cilmente interpretable, por lo que en el supuesto de que tuviera una capacidad predictiva aceptable podr�amos buscar la coherencia de sus criterios de clasificaci�n con los resultados de estudios previos de metilaci�n diferencial entre grupos. 
* Lo mismo pero en lugar de utilizar datos de metilaci�n usando los datos (picos) de histonas modificadas.

Se partir�a de los datos depositados en TCGA (aproximadamente �? 11000 muestras existentes en TCGA) cuyas tablas de fenotipos contemplan unas �? 400 variables, entre las que se encuentran los tipos de tejido.

El planteamiento trata de ser algo diferente del que he observado en otros trabajos y que resumido consistir�a en partir de dos grupos de muestras (dos grupos fenot�picamente distintos) y obtener las �features� (expresi�n de genes, locis metilados, picos de histonas) que estad�sticamente poseen valores diferenciales.
 
En cuanto a la b�squeda de bibliograf�a:

Partiendo de B�squedas en GOOGLE scholar:
* epigenetic�methylation
* promoter�methylation
* histone�methylation
* dna methylation

Ordenar�a por fecha de publicaci�n, puesto que las librer�as R para manejar los repositorios de datos de TCGA son recientes no creo que sea preciso remontarse muchos a�os atr�s para comprobar si existe publicado alg�n trabajo con este enfoque. 

Este es el apartado en el que, de ser posible, querr�a alguna indicaci�n para poder avanzar estos dos meses de verano: disponer de alg�n criterio de selecci�n de art�culos dado que su n�mero es abrumador.



PROBLEM�TICA ESPERADA

* Obtener una Conclusi�n del trabajo de dos palabras: �NO FUNCIONA�. Si el objetivo es la comparaci�n de dos grupos de muestras fenot�picamente distintas siempre vamos a encontrar alg�n grupo de �features� con valores significativamente diferentes y por lo tanto llegaremos a algo que nos permite redactar una �Conclusi�n�. En esta propuesta no esta garantizado que sea capaz de dise�ar y entrenar un modelo con resultados predictivos aceptables.
* Que el volumen de muestras disponible no sea ni de lejos el necesario para entrenar una red neuronal y por lo tanto la matriz de confusi�n final sea decepcionante  (estar�a en la situaci�n descrita en el punto anterior). 
* Que las capacidades de c�lculo no sean suficientes y no consiga entrenar el algoritmo.
* Las habituales de todo el mundo de limitaci�n de tiempo y conocimientos a las que hay que a�adir las m�as propias.


Los documentos y scripts del trabajo Fin de M�ster de m�ster UOC se almacenan en la ubicaci�n:
C:\TFM UOC
https://github.com/ajoven10/TFM_UOC.git



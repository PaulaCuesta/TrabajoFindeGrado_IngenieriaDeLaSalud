library (GEOquery)
library (limma)
library (dplyr)
library ("hgu133a2.db")
library("hgu133a.db")



#Obtenemos el conjunto de datos "GSE14520" sobre datos de expresión de cáncer hepático en GEO
datos_cancer <- getGEO ("GSE14520")
#Mostramos por pantalla los datos que hemos obtenido
datos_cancer



#Creamos un objeto "array_GPL571" que recoge los datos que se han analizado en la plataforma GPL571
array_GPL571 <-datos_cancer[2]
#Mostramos por pantalla este objeto que hemos creado
array_GPL571


#Creamos otro objeto "array_GPL3921" que recoge los datos analizados en la plataforma GPL3921
array_GPL3921 <-datos_cancer[1]
#Mostramos por pantalla el nuevo objeto "array_GPL3921"
array_GPL3921

#Se aplica una transformación logarítmica en base 2 a la matriz con los datos de expresión para normalizar los valores y poder compararlos más fácilmente
exp_array_GPL571 <- log2 (array_GPL571$`GSE14520-GPL571_series_matrix.txt.gz`@assayData$exprs)
#Se muestran los datos de expresión normalizados por pantalla
exp_array_GPL571

dim(exp_array_GPL571)

#Se realiza la misma transformación logarítmica a la matriz con los datos de la plataforma GPL3921 para normalizarlos y facilitar la comparación entre los valores de expresión
exp_array_GPL3921 <- log2 (array_GPL3921$`GSE14520-GPL3921_series_matrix.txt.gz`@assayData$exprs)
#Se muestran los datos normalizados nuevamente por pantalla
exp_array_GPL3921


#Creamos un objeto que almacene los nombres de la columnas de los datos, es decir, los nombres de las muestras que están siendo analizadas
muestras_array_GPL571 <- colnames (exp_array_GPL571)
#Mostramos por pantalla los nombres de las columnas, es decir, las muestras analizadas
muestras_array_GPL571
#Almacenamos la descripción de cada una de estas muestras en un nuevo objeto "desc_array_GPL571"
desc_array_GPL571 <- array_GPL571$`GSE14520-GPL571_series_matrix.txt.gz`@phenoData@data$title
#Mostramos por pantalla el nuevo objeto con las descripciones que hemos creado
desc_array_GPL571

#Creamos un nuevo DataFrame en el que convinamos el identificador de la muestra (GSM...) con su descripción (Liver Tumor Tissue o Liver Non-Tumor Tissue)
df_muestras_array_GPL571 <- data.frame (Muestra = muestras_array_GPL571, Descripción = desc_array_GPL571 )
df_muestras_array_GPL571


#Creamos otro objeto en el que se almacenen las muestras que se analizan en la plataforma GPL3921
muestras_array_GPL3921 <- colnames (exp_array_GPL3921)
#Mostramos este objeto por pantalla
muestras_array_GPL3921
#Creamos un objeto que contenga la descripción de cada muestra, es decir, si se trata de un tejido tumoral o no
desc_array_GPL3921<- array_GPL3921$`GSE14520-GPL3921_series_matrix.txt.gz`@phenoData@data$title
#Mostramos este nuevo objeto que hemos creado por pantalla
desc_array_GPL3921
#Creamos un nuevo DataFrame que combine el ID de las muestras de la plataforma GPL3921 con su información de patogenicidad
df_muestras_array_GPL3921 <- data.frame (Muestra = muestras_array_GPL3921, Descripción = desc_array_GPL3921 )
#Mostramos este nuevo objeto que hemos creado por pantalla
df_muestras_array_GPL3921


# Convertimos la descripción de las muestras del nuevo DF a cadena de texto, para poder operar sobre ellas
df_muestras_array_GPL571$Descripción <- as.character(df_muestras_array_GPL571$Descripción)
#Realizamos una conversión de la descripción a cadena de caracteres, ya que estaba en forma de factor
df_muestras_array_GPL3921$Descripción <- as.character(df_muestras_array_GPL3921$Descripción)

#Buscamos la cadena de texto "Non-Tumor" en la descripción de las muestras, y cambiamos la descripción de las muestras que cumplen la condición por "Sanos"
df_muestras_array_GPL571$Descripción[grepl("Non-Tumor", df_muestras_array_GPL571$Descripción)]<- "Sanos"

#Hacemos lo mismo, pero poniendo la descripción "Donantes" a aquellas muestras en cuya descripción aparece la palabra "donors"
df_muestras_array_GPL571$Descripción[grepl("donors", df_muestras_array_GPL571$Descripción)]<- "Donantes"
#Por último, describimos las muestras como "Patogénicos" cuando en su descripción aparecen las palabras "Liver Tumor"
df_muestras_array_GPL571$Descripción[grepl("Liver Tumor", df_muestras_array_GPL571$Descripción)]<- "Patogénicos"


#Creamos un objeto que contenga las muestras cuya descripción equivalga a los pacientes sanos
sanos_GPL571 <- df_muestras_array_GPL571[df_muestras_array_GPL571$Descripción=="Sanos",]
#Obtenemos el número de muestras de pacientes sanos que encontramos en la plataforma GPL571
dim(sanos_GPL571)
#Creamos un objeto en el que almacenamos las muestras cuya descripción sea "Patogénicos"
enfermos_GPL571 <- df_muestras_array_GPL571[df_muestras_array_GPL571$Descripción=="Patogénicos",]
enfermos_GPL571[1:10,]
#Calculamos el número de muestras patogénicas que tenemos en el array GPL571
dim (enfermos_GPL571)
#Creamos otro objeto con los donantes, es decir, aquellos en cuya descripción aparece la palabra "Donantes"
donantes_GPL571 <- df_muestras_array_GPL571[df_muestras_array_GPL571$Descripción=="Donantes",]
#Calculamos el número de muestras que tenemos de donantes en la plataforma
dim(donantes_GPL571)

#Combinamos los tres DataFrames individuales en uno solo, uno detrás de otro, en el orden especificado
muestras_array_GPL571 <- rbind.data.frame (sanos_GPL571, enfermos_GPL571, donantes_GPL571)
#Mostramos el objeto que almacena los 3 DataFrames por pantalla
muestras_array_GPL571


#Cambiamos la descripción de las muestras en los que aparecen en su descripción las palabras "Non-Tumor" por la descripción "Sanos"
df_muestras_array_GPL3921$Descripción[grepl("Non-Tumor", df_muestras_array_GPL3921$Descripción)] <- "Sanos"
#Hacemos lo mismo para aquellas muestras en cuya descripción aparecen las palabras "Liver Tumor", y cambiamos esta por la palabra "Patogénicos"
df_muestras_array_GPL3921$Descripción[grepl("Liver Tumor", df_muestras_array_GPL3921$Descripción)] <- "Patogénicos"

#Creamos un objeto que contenga las muestras cuya descripción es "Sanos"
sanos_GPL3921 <- df_muestras_array_GPL3921[df_muestras_array_GPL3921$Descripción=="Sanos",]
#Calculamos el número de muestras procedentes de tejidos sanos que tenemos en la plataforma GPL3921
dim(sanos_GPL3921)
#Creamos un nuevo objeto que contenga las muestras que presentan la descripción "Patogénicos"
enfermos_GPL3921 <- df_muestras_array_GPL3921[df_muestras_array_GPL3921$Descripción=="Patogénicos",]
#Obtenemos el número de muestras de tejido patogénico que tenemos en GPL3921
dim(enfermos_GPL3921)
#Combinamos en un nuevo DataFrame los objetos DataFrame que teníamos, uno a continuación de otro
muestras_array_GPL3921 <- rbind.data.frame(sanos_GPL3921, enfermos_GPL3921)
#Mostramos este nuevo DataFrame que hemos creado por pantalla
muestras_array_GPL3921

#Creamos un factor que nos devuelva "TRUE" si la descripción de las muestras es "Patogénicos" y FALSE en caso contrario, es decir, si la descripción es "Sanos" o "Donantes"
comp_GPL571 <- as.factor(muestras_array_GPL571$Descripción == "Patogénicos")
#Creamos una matriz de diseño en función del factor, para poder ajustar los modelos lineales que tenemos
matriz_GPL571 <- model.matrix(~comp_GPL571)
#Creamos un objeto en el que se almacena el ajuste del modelo lineal que hacemos para cada gen para poder comparar los valores de expresión entre los genes
ajuste_GPL571 <- lmFit (exp_array_GPL571, matriz_GPL571)
#Aplicamos un ajuste bayesiano a los resúmenes estadísticos ajustados de la expresión de los genes
res_bayesiano_GPL571 <- eBayes (ajuste_GPL571)
#Generamos una tabla ordenada con los valores de expresión diferencial de todos los genes que aparecen en las muestras que hemos analizado
resultados_GPL571 <- topTable (res_bayesiano_GPL571, coef=2, number = nrow(res_bayesiano_GPL571))
#Mostramos la tabla ordenada con los resultados de la expresión diferencial de los genes
resultados_GPL571
#Seleccionamos a partir de la tabla, una subtabla que contenga unos datos significativos, es decir, aquellas muestras cuyo p-valor sea menor de 0,05. Es decir, aquellas muestras de las que realmente tengamos evidencia científica de que se encuentran diferencialmente expresadas
pvalor_GPL571 <- subset (resultados_GPL571, resultados_GPL571$adj.P.Val<0.05)
#Mostramos esta nueva tabla que hemos creado por pantalla
pvalor_GPL571


#Creamos un objeto en el que almacenamos valores de TRUE cuando la descripción de la muestra sea "Patogénicos" y FALSE en el caso contrario, cuando dicha descripción sea "Sanos"
comp_GPL3921 <- as.factor(muestras_array_GPL3921$Descripción == "Patogénicos")
#Creamos una matriz de diseño en base al factor que hemos calculado previamente
matriz_GPL3921 <- model.matrix(~comp_GPL3921)
#Realizamos el ajuste lineal de los datos de la expresión al modelo lineal
ajuste_GPL3921 <- lmFit (exp_array_GPL3921, matriz_GPL3921)
#Aplicamos el ajuste bayesiano a los datos de expresión que han sido ajustados linealmente
res_bayesiano_GPL3921 <- eBayes (ajuste_GPL3921)
#Creamos una tabla con las muestras ordenadas en función del resultado de su expresión diferencial
resultados_GPL3921 <- topTable (res_bayesiano_GPL3921, coef=2, number = nrow(res_bayesiano_GPL3921))
#Mostramos esta tabla que hemos creado por pantalla
resultados_GPL3921
#De la tabla con los resultados de expresión diferencial, seleccionamos aquellas entradas que sean significativas científicamente, creando una nueva tabla con ellas. Estas entradas son aquellas cuyo p-valor es inferior a 0.05
pvalor_GPL3921 <- subset (resultados_GPL3921, resultados_GPL3921$adj.P.Val<0.05)
pvalor_GPL3921

dim(pvalor_GPL3921)


#Creamos un objeto en el que almacenamos los identificadores de cada una de las filas de la tabla que hemos creado con los datos relativos a la expresión diferencial
cols_GPL571 <- rownames (pvalor_GPL571)
#Convertimos nuestra tabla a DataFrame añadiendo una columna más que equivalga a los identificadores de cada una de las filas, es decir, que equivalga a los identificadores de cada una de las sondas utilizadas. Esta columna la colocamos además como la primera en el DataFrame
pvalor_GPL571 <- cbind.data.frame(id_sonda = cols_GPL571, pvalor_GPL571)
#Asociamos el nombre del gen correspondiente a cada una de las sondas utilizadas, en base a la biblioteca que se ha utilizado, que en este caso ha sido "HG-U133A_2"
genes_GPL571 <- select (hgu133a2.db, as.vector(pvalor_GPL571$id_sonda), c("SYMBOL"))
#Mostramos el objeto que acabamos de crear por pantalla
genes_GPL571


#Creamos un nuevo objeto con los identificadores de las sondas que hemos utilizado en el estudio de la expresión diferencial de las muestras
cols_GPL3921 <- rownames (pvalor_GPL3921)
#Creamos un DataFrame en el que añadimos al principio la columna con los identificadores de las sondas, en combinación de los datos de expresión de cada una de las muestras que tenemos
pvalor_GPL3921 <- cbind.data.frame(id_sonda = cols_GPL3921, pvalor_GPL3921)
#Asociamos el identificador de las sondas con el nombre del gen al que hacen referencia, en base a la biblioteca que se ha utilizado para renombrar las sondas, que en este caso ha sido "HT_HG-U133A"
genes_GPL3921 <- select (hgu133a.db, as.vector(pvalor_GPL3921$id_sonda), c("SYMBOL"))
#Mostramos por pantalla el nuevo objeto que ha sido creado
genes_GPL3921



# Función que combina dos DataFrames, df1 y df2, por una de sus columnas, la cual presenta el mismo nombre en los 2 DataFrames y hace referencia al nombre de la sonda, combinando así el nombre de cada gen con el resultado de los datos de expresión
combinar_df <- function(df1,df2){
  # Guardamos el número de filas (longitud) del primer DataFrame
  filas_df1 <- nrow(df1)
  #Guardamos el número de filas (longitud) del segundo DataFrame
  filas_df2 <- nrow(df2)
  #Inicializamos un vector vacío al que denominamos "logFC"
  logFC <- c()
  #Inicializamos un vector vacío al que denominamos "Symbol"
  Symbol <- c()
  #Recorremos todas las filas tanto del primer como del segundo vector, de forma que si el nombre de la sonda coincide en ambos se asigna al vector "logFC" el valor de los datos de expresión y al vector "Symbol" el gen al que corresponde la muestra
  for(i in 1:filas_df1){
    for(j in 1:filas_df2){
      if(df1[i,1]==df2[j,1]){
        logFC <- append(logFC,df1[i,2])
        Symbol <- append(Symbol,df2[j,2])
      }
    }
  }
  #Combinamos en un DataFrame estas nuevas variables que hemos creado
  result <- cbind.data.frame(SYMBOL=Symbol,logFC=logFC)
  #Imprimimos por pantalla el resultado
  print(result)
}

#Utilizamos la función "combinar_df" para aplicar esta combinación a los valores de expresión de la plataforma GPL571
datos_genes_GPL571 <- combinar_df (pvalor_GPL571, genes_GPL571)
#Eliminamos las apariciones en las que se desconoce el nombre del gen
datos_genes_GPL571 <- na.omit (datos_genes_GPL571)
#Como en esta plataforma solo se había analizado una muestra con un valor científico significativo, pero se desconoce el gen de la sonda, lo vamos a eliminar, y ya no vamos a continuar con el análisis de esta plataforma

# Como el DataFrame que tenemos de la plataforma es muy grande lo vamos a dividir en trozos por el método de Slicing o rebanado
pvalor_GPL3921_1 <- pvalor_GPL3921 [1:1500, ]
pvalor_GPL3921_2 <- pvalor_GPL3921 [1501:3000, ]
pvalor_GPL3921_3 <- pvalor_GPL3921 [3001:4500, ]
pvalor_GPL3921_4 <- pvalor_GPL3921 [4501:6000, ]
pvalor_GPL3921_5 <- pvalor_GPL3921 [6001:7500, ]
pvalor_GPL3921_6 <- pvalor_GPL3921 [7501:9000, ]
pvalor_GPL3921_7 <- pvalor_GPL3921 [9001:10329, ]


# Vamos a realizar la combinación de los DataFrames, realizando la combinación de cada subDataFrame por separado
comb1 <- combinar_df (pvalor_GPL3921_1, genes_GPL3921)
comb2 <- combinar_df (pvalor_GPL3921_2, genes_GPL3921)
comb3 <- combinar_df (pvalor_GPL3921_3, genes_GPL3921)
comb4 <- combinar_df (pvalor_GPL3921_4, genes_GPL3921)
comb5 <- combinar_df (pvalor_GPL3921_5, genes_GPL3921)
comb6 <- combinar_df (pvalor_GPL3921_6, genes_GPL3921)
comb7 <- combinar_df (pvalor_GPL3921_7, genes_GPL3921)

#Combinamos en un nuevo DataFrame los resultados obtenidos de la combinación de DataFrames
datos_genes_GPL3921 <- rbind.data.frame (comb1, comb2, comb3, comb4, comb5, comb6, comb7)
#Eliminamos las muestras para las que se desconoce la sonda, y, por lo tanto, nuestro gen de interés
datos_genes_GPL3921 <- na.omit (datos_genes_GPL3921)
#Ordenamos alfabéticamente nuestros genes expresados diferencialmente
datos_genes_GPL3921 <- datos_genes_GPL3921[order(datos_genes_GPL3921$SYMBOL), ]

#Eliminamos los nombres de las filas del DataFrame, es decir, los nombres de las sondas que utiliza el array
rownames (datos_genes_GPL3921) <- c()
datos_genes_GPL3921[1:10, ]

#Añadimos una nueva columna al DataFrame que contiene los genes junto a sus valores de expresión para indicar si estos se encuentran sobre o subexpresados. 
datos_genes_GPL3921 <- cbind.data.frame(datos_genes_GPL3921,status=c("Expresión"))

#Si el valor de logFC es inferior a 0 significa que los datos se encuentran subexpresados
datos_genes_GPL3921$status[datos_genes_GPL3921$logFC<0] <- "Subexpresados"
#Por el contrario, si logFC es positivo, dichos genes se encuentran subexpresados
datos_genes_GPL3921$status[datos_genes_GPL3921$logFC>0] <- "Sobrexpresados"

#Mostramos el objeto que hemos creado con la nueva columna que indica si se encuentran sobre o subexpresados.
datos_genes_GPL3921
#Calculamos el número de genes diferencialmente expresados
dim(datos_genes_GPL3921)

#Creamos un objeto que almacene los genes sobreexpresados con un logFC mayor a 0.15
genes_sobreexp_GPL3921 <- subset(datos_genes_GPL3921, datos_genes_GPL3921$logFC >0.15)
#Mostramos por pantalla estos genes sobreexpresados
genes_sobreexp_GPL3921

#Recuperamos el número de genes que se encuentran sobreexpresados en base a nuestro criterio del valor de logFC superior a 0.15
dim(genes_sobreexp_GPL3921)

#Creamos un objeto con los genes subexpresados con un valor de logFC inferior a 0.15
genes_subexp_GPL3921 <- subset(datos_genes_GPL3921, datos_genes_GPL3921$logFC < -0.15)
#Mostramos el objeto que hemos creado por pantalla
genes_subexp_GPL3921

#Calculamos el número de genes que se encuentran subexpresados en nuestra plataforma
dim(genes_subexp_GPL3921)

#Obtenemos un boxplot que represente la distribución de las muestras entre "Sanas" y "Enfermas", así como su valor de logFC para poder estudiar la distribución de la expresión de los genes en estas muestras.
boxplot(datos_genes_GPL3921$logFC,col=c("blue"),ylab="logFC",main="Comparation of logFC 
values",cex.lab=1.5,cex.main=1.5)
abline(h=0,col=1,lty=2)
axis(1,at=1,labels=c("Sanos - Patogénicos"),cex.axis=1.5)

#Creamos una gráfica que compara las muestras de forma individual en base a su nivel de expresión
par(mfrow=c(1,1),mar=c(2, 4, 2, 2))
plot(datos_genes_GPL3921$logFC,col=3,pch=18,ylab="logFC",main="Valores de logFC en pacientes sanos y enfermos")
abline(h=0,col=2)






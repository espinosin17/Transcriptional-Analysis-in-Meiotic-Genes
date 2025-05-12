# Instalar y cargar el paquete jsonlite
#install.packages("jsonlite")

# Cargar el paquete jsonlite
library(jsonlite) #Para utilizar el comando fromJSON

# Leer el archivo como texto
path <- "RNAseq/Data/Klepikova eFP (RNA-Seq data)-AT1G04880-AtFM5.txt" #Dirección del archivo
file_content <- readLines(path, encoding = "UTF-8") #Leer el archivo txt línea por línea

# Encontrar la línea que contiene la data JSON "\\[" y unirla en un solo texto sin espacios
json_text <- paste(file_content[grep("\\[", file_content):length(file_content)], collapse = "")

# Convertir el JSON a un data frame
data <- fromJSON(json_text)

# Mostrar los primeros registros para verificar
head(data)
nrow(data)

# Crear dos vectores con los nombres de los estadios y sus respectivos códigos de medición
fase <- c("F1", "F2", "F3", "F4", "F5", "F6_8", "F9_11", "F12_14", "F15_18", "F19_more")
name <-  c("SRR3581693","SRR3581694", #Estos códigos pertenecen a la primera réplica
           "SRR3581695","SRR3581696",
           "SRR3581697","SRR3581698",
           "SRR3581699","SRR3581700",
           "SRR3581701","SRR3581702")

#Generar un dataframe con las fases y códigos de la primera réplica
RNAdata <- as.data.frame(cbind(fase, name))
RNAdata

#Generar un dataframe que recupere los valores de RNAdata según los códigos de réplica
ejemplo <- merge(RNAdata, data, by = "name", all.x = TRUE)
ejemplo

#Actualizar RNAdata con las columnas fase y valor
#Eliminar la columna con el código de la primera réplica
RNAdata <- ejemplo[,c(2,3)]
RNAdata

#Asignar el nombre Rep1 para los valores de la primera réplica (columna "value")
colnames(RNAdata)[colnames(RNAdata) == "value"] <- "Rep1"
RNAdata

#Actualizo el objeto "name" con los códigos de los valores de la segunda réplica
name <- c("SRR3581859",
          "SRR3581860",
          "SRR3581861",
          "SRR3581862",
          "SRR3581863",
          "SRR3581864",
          "SRR3581865",
          "SRR3581866",
          "SRR3581867",
          "SRR3581868")
#Añadir la columna "name" al dataframe RNAdata
RNAdata <- cbind(RNAdata, name)
RNAdata

#Seleccionar los valores según los códigos de la segunda réplica
ejemplo <- merge(RNAdata, data, by = "name", all.x = TRUE)
ejemplo

#Elegir las columnas "fase", "Rep1" y "value"
#value almacena las mediciones de la segunda réplica
RNAdata <- ejemplo[,c(2,3,4)]
RNAdata

#Cambiar "value" por "Rep2"
colnames(RNAdata)[colnames(RNAdata) == "value"] <- "Rep2"
RNAdata

#Realizar un promedio de la primera y segunda medición
value <- (as.numeric(RNAdata$Rep1) + as.numeric(RNAdata$Rep2))/2
value

RNAdata <- cbind(RNAdata, value)
RNAdata

RNAdata <- RNAdata[, c(1,4)]
RNAdata

gene_name <- "AtUBC22" ###################################################### nuevo nombre
gene <- c(rep(gene_name, 10))
gene

# install.packages("tibble")
library(tibble)
RNAdata <- add_column(RNAdata, gene, .after = "fase") 
RNAdata

file_name <- "../Data/UBC22_RNAdata" ####################################### nuevo nombre de archivo
write.table(x = RNAdata, file = file_name, sep = ',',col.names = TRUE, quote = FALSE, row.names = FALSE)





# Fijar directorio de trabajo
setwd(dir = "./RNAseq")
getwd()

# Importar datos de expresión de RNAm (Klepikova et al., 2016)
DMC1RNAdata <- read.table(file = "./Data/FlowerBudsExpression/DMC1_RNAdata", header = TRUE, sep = ",")
SPO11_1RNAdata <- read.table(file = "./Data/FlowerBudsExpression/SPO11-1_RNAdata", header = TRUE, sep = ",")
DUETRNAdata <- read.table(file = "./Data/FlowerBudsExpression/DUET_RNAdata", header = TRUE, sep = ",")
FSTRNAdata <- read.table(file = "./Data/FlowerBudsExpression/FST_RNAdata", header = TRUE, sep = ",")
UBC22RNAdata <- read.table(file = "./Data/FlowerBudsExpression/UBC22_RNAdata", header = TRUE, sep = ",")
FM5RNAdata <- read.table(file = "./Data/FlowerBudsExpression/FM5_RNAdata", header = TRUE, sep = ",")
FM13RNAdata <- read.table(file = "./Data/FlowerBudsExpression/FM13_RNAdata", header = TRUE, sep = ",")
FM12RNAdata <- read.table(file = "./Data/FlowerBudsExpression/FM12_RNAdata", header = TRUE, sep = ",")
ARP6RNAdata <- read.table(file = "./Data/FlowerBudsExpression/ARP6_RNAdata", header = TRUE, sep = ",")
ASK1RNAdata <- read.table(file = "./Data/FlowerBudsExpression/ASK1_RNAdata", header = TRUE, sep = ",")

# Definir función para obtener los valores relativos para cada gen
RelativeMaths <- function(datos){
  RelativeValue <- datos$value/max(datos$value) #Divide el valor de expresión entre el máximo de los datos en los 10 estadios
  datos <- cbind(datos, RelativeValue) 
  return(datos)
}

# Obtener los valores relativos para los 10 genes meióticos
DMC1RNAdata <- RelativeMaths(DMC1RNAdata)
SPO11_1RNAdata <- RelativeMaths(SPO11_1RNAdata)
DUETRNAdata <- RelativeMaths(DUETRNAdata)
FSTRNAdata <- RelativeMaths(FSTRNAdata)
UBC22RNAdata <- RelativeMaths(UBC22RNAdata)
FM5RNAdata <- RelativeMaths(FM5RNAdata)
FM13RNAdata <- RelativeMaths(FM13RNAdata)
FM12RNAdata <- RelativeMaths(FM12RNAdata)
ARP6RNAdata <- RelativeMaths(ARP6RNAdata)
ASK1RNAdata <- RelativeMaths(ASK1RNAdata)

# Combinar los datos de los 10 genes en una sola tabla
data_combinada <- rbind(DMC1RNAdata, SPO11_1RNAdata, DUETRNAdata, FSTRNAdata, UBC22RNAdata, FM5RNAdata, FM13RNAdata, FM12RNAdata, ARP6RNAdata, ASK1RNAdata)
data_combinada # Se observaon 100 filas (10 genes x 10 estadios)

# Importar paqueterías
install.packages("ggplot2")

###################################### VISUALIZACIÓN DE DATOS
# Cargar paquetes
library("ggplot2")

# Definir objetos con el nombre de los promotores y los estadios del desarrollo floral
orden_fase <- c("F19_more", "F15_18", "F12_14", "F9_11", "F6_8", "F5", "F4", "F3", "F2", "F1")  # 10 estadios en orden
orden_promoters <- c("AtARP6", "AtASK1", "AtFM13", "AtFST", "AtUBC22", "AtDMC1", "AtSPO11-1", "AtFM5", "AtDUET", "AtFM12") # 10 genes

# Convertir orden_fase y orden_promotores en factores de data combinada
data_combinada$fase <- factor(exp_long$fase, levels = orden_fase) # Factor aplicado a la columna fase
data_combinada$gene <- factor(exp_long$gene, levels = orden_promoters) # Factor aplicado a la columna gene

# Graficar mapa de calor de expresión relativa de ARNm
exp_heatmap <- ggplot(data = data_combinada, mapping = aes(x = fase,
                                                           y = gene,
                                                           fill = RelativeValue)) + # Crear el canvas
  geom_tile() + # Mostrar los valores de expresión relativa en "tiles"
  ylab(label = "Genes meióticos") + # Nombre del eje y
  xlab(label = "Estadio del desarrollo floral") + # Nombre del eje x
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + # Rota los nombres del eje x en 45°
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red", # Rellena el mosaico considerando estos colores y sus valores
                       midpoint = 0.5,  # Ajustar el punto medio a 0.5
                       limits = c(min(exp_long$RelativeValue), max(exp_long$RelativeValue))) +# Ajusta los límites del gráfico según mis datos
  scale_y_discrete(labels = c("ARP6", "ASK1", "FM13", "FST", "UBC22", "DMC1", "SPO11-1", "FM5", "DUET", "FM12")) + # Nombra los genes en el eje y
  scale_x_discrete(labels = c("F19+","F15-18", "F12-14", "F9-11", "F6-8", "F5", "F4","F3", "F2", "F1")) + # Nombra los estadios en el eje x
  labs(fill = "Expresión relativa") # Añade el nombre de la leyenda
exp_heatmap # Mostrar el mapa de calor de expresión relativa de ARNm

####################################### ANÁLISIS ESTADÍSTICO
# Definir función para agrupar los valores
meiotic_group <- function(datos){
  relativos_meiosis <- datos$RelativeValue[c(7:10)] # Incluye a los estadios F19+, F15-18, F12-14 y F6-8
  return(relativos_meiosis)
}

postmeiotic_group <- function(datos){
  relativos_postmeioticos <- datos$RelativeValue[c(1:6)] # Incluye desde F1 hasta F6-8
  return(relativos_postmeioticos)
}

# ¿Estarán los genes meióticos más expresados durante la meiosis?
# Test de Wilcoxon para dos muestras independientes
# Ho: la diferencia de medianas es igual 0 
# Ha: la diferencia de medianas es mayor a 0 (mediana meiotica – mediana postmeiotica)
wilcox.test(meiotic_group(DMC1RNAdata), postmeiotic_group(DMC1RNAdata), alternative = "greater")
wilcox.test(meiotic_group(SPO11_1RNAdata), postmeiotic_group(SPO11_1RNAdata), alternative = "greater")
wilcox.test(meiotic_group(DUETRNAdata), postmeiotic_group(DUETRNAdata), alternative = "greater")
wilcox.test(meiotic_group(FM12RNAdata), postmeiotic_group(FM12RNAdata), alternative = "greater")
wilcox.test(meiotic_group(FM5RNAdata), postmeiotic_group(FM5RNAdata), alternative = "greater")
wilcox.test(meiotic_group(UBC22RNAdata), postmeiotic_group(UBC22RNAdata), alternative = "greater")
wilcox.test(meiotic_group(FSTRNAdata), postmeiotic_group(FSTRNAdata), alternative = "greater")
wilcox.test(meiotic_group(FM13RNAdata), postmeiotic_group(FM13RNAdata), alternative = "greater")
wilcox.test(meiotic_group(ASK1RNAdata), postmeiotic_group(ASK1RNAdata), alternative = "greater")
wilcox.test(meiotic_group(ARP6RNAdata), postmeiotic_group(ARP6RNAdata), alternative = "greater")

# ¿Estarán los otros 9 genes igualmente expresados que DMC1 durante los estadios meióticos del desarrollo floral? 
# Test de Wilcoxon para una muestra
# Utilicé la mediana de expresión relativa de DMC1 como una constante comparable entre genes
# Ho: la diferencia de medianas es igual 0
# Ha: la diferencia de medianas es mayor a 0 (mediana meiótica – mediana postmeiótica)
mediana_meiotica_DMC1 <- median(meiotic_group(DMC1RNAdata))
wilcox.test(meiotic_group(FM12RNAdata), mu = mediana_meiotica_DMC1, alternative = "two.sided")
wilcox.test(meiotic_group(DUETRNAdata), mu = mediana_meiotica_DMC1, alternative = "two.sided")
wilcox.test(meiotic_group(FM5RNAdata), mu = mediana_meiotica_DMC1, alternative = "two.sided")
wilcox.test(meiotic_group(SPO11_1RNAdata), mu = mediana_meiotica_DMC1, alternative = "two.sided")
wilcox.test(meiotic_group(UBC22RNAdata), mu = mediana_meiotica_DMC1, alternative = "two.sided")
wilcox.test(meiotic_group(FSTRNAdata), mu = mediana_meiotica_DMC1, alternative = "two.sided")
wilcox.test(meiotic_group(FM13RNAdata), mu = mediana_meiotica_DMC1, alternative = "two.sided")
wilcox.test(meiotic_group(ASK1RNAdata), mu = mediana_meiotica_DMC1, alternative = "two.sided")
wilcox.test(meiotic_group(ARP6RNAdata), mu = mediana_meiotica_DMC1, alternative = "two.sided")

#  A nivel del ARN, para un gen meiótico estricto, podríamos establecer un 35% de expresión permisible en los estadios post-meióticos.  
# ¿Se tendrán genes meióticos estrictos en el grupo de genes meióticos?
# Wilcoxon de una muestra
# Ho: la mediana de expresio n post-meiótica relativa del gen es igual a 0.35
# Ha: la mediana de expresio n post-meiótica relativa del gen es menor que 0.35
wilcox.test(postmeiotic_group(DMC1RNAdata), mu = 0.35, alternative = "less")
wilcox.test(postmeiotic_group(FM12RNAdata), mu = 0.35, alternative = "less")
wilcox.test(postmeiotic_group(DUETRNAdata), mu = 0.35, alternative = "less")
wilcox.test(postmeiotic_group(FM5RNAdata), mu = 0.35, alternative = "less")
wilcox.test(postmeiotic_group(SPO11_1RNAdata), mu = 0.35, alternative = "less")
wilcox.test(postmeiotic_group(UBC22RNAdata), mu = 0.35, alternative = "less")
wilcox.test(postmeiotic_group(FSTRNAdata), mu = 0.35, alternative = "less")
wilcox.test(postmeiotic_group(FM13RNAdata), mu = 0.35, alternative = "less")
wilcox.test(postmeiotic_group(ASK1RNAdata), mu = 0.35, alternative = "less")
wilcox.test(postmeiotic_group(ARP6RNAdata), mu = 0.35, alternative = "less")















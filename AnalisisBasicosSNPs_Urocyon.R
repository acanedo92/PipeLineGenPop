####################################################################################
# Análisis de Diversidad y Estructura Genética a partir de SNP's de un archivo VCF.#
####################################################################################
# Instalación de paquetes:
install.packages(c("adegenet", "hierfstat", "ape", "vcfR"))

# Cargar paqueterías:

library(vcfR)
library(adegenet)
library(hierfstat)
library(ape)


setwd ("/media/gabriela/ADATA_HD710/PipeLineGenPop/")

# Cargar un archivo tipo vcf y convertirlo en objeto genind para analizar con adegente. En este caso nuestro archivo esta en formato Genpop y no es necesrio transformarlo.
# a <- read.vcfR("File.vcf")
# b <- vcfR2genind(a)
# c <- vcfR2genlight(a)

# Cargar archivo tipo genpop
a <- read.genepop("GP_withgrays.gen")
class(a)
##################################
# Explorar nuestra base de datos #
##################################

a # Ver que tipo de objetos creamos
a@ploidy # Ver la haploidia de los datos
a@pop # Ver a las poblaciones en el objeto genind (@pop debe ser un factor)
a@tab # Es la matriz de todos los datos o genotipos
a@tab[c(1:5),c(1:5)] # Ver los datos en las primeras 5 filas y 5 columnas
a@tab[,1] # Ver el nombre de los individuos

###########################################
# Analisis basicos de diversidad genetica #
###########################################

divers <- basic.stats(a) # Analisis global, todas las pobalcaiones
divers$overall # Promedio de todos los datos. Todas las medidas de diversidad

# Para obtener los promedios por poblacionn
colnames(divers$Ho) # Para ver el nombre de las poblaciones
levels(a@pop) 

nombresPobs <- colnames(divers$Ho) # Hacer un objeto con el nombre de las poblaciones
HetObs <-  colMeans(divers$Ho)   # Heterocigosis observada promedio de la poblacion 1
Richness <- allelic.richness(a) # Estima la riqueza alélica, los recuentos alélicos enrarecidos, por locus y población

# Hacer un grafico de barras con los valores de heterocigosis observada.
barplot (HetObs, xlab = "Populations", ylab="HO", col="red", names.arg = nombresPobs)
abline (h=0.1253) # Si queremos agregar una linea horizontal que se??ale el valor promedio, tomando el pormedio de la tabla divers$overall

##########################
# Podemos hacer lo mismo para obtener los valores promedio de las demas medidas de diversidad.
##########################

# dev.new()

par(mfrow=c(1,3)) 
myCol <- c("darkblue","purple","green","orange","red","blue", "black")
barplot(colMeans(divers$Hs, na.rm = T), col = myCol , names.arg = nombresPobs, las=2, axis.lty = 6, axisnames = T, ylab = "H.S.", xlab = "Populations")
barplot(colMeans(divers$Ho, na.rm = T), col = myCol, names.arg = nombresPobs, las=2, axis.lty = 6, axisnames = T, ylab = "H.O.",  xlab = "Populations")
barplot(colSums(Richness$Ar),  las=2, col = myCol, ylab = "N. Alleles")
                  
###################################
# Analisis de estructura genetica #
###################################

# Podemos estimar distancias para hacer heatmap de distancia gen??tica entre individuos
x.dist <- dist(a)
heatmap(as.matrix(x.dist))

# Estimar distancia genetica de Nei entre poblaciones
b <- genind2genpop(a) # Requerimos primero trasnformar nuestros datos a formato genpop
y.dist<- dist.genpop(b, method = 1)
heatmap(as.matrix(y.dist))

# Hacer un Neighbor joining
plot(nj(y.dist), type="fan")

# Hacer un UPGMA
c <- hclust(y.dist, method = "average")
plot(c, main= "UPGMA dendrogram", xlab = "Urocyon litoralis")

# Hacer un DAPC con los colores para identificar grupos geneticos. Se puede hacer con las poblaciones ya definidas
dpca1 <- dapc(a)
scatter.dapc(dpca1, cex=4, scree.da = FALSE)
compoplot.dapc(dpca1)

# Podemos buscar el numero de grupos geneticos presentes en los datos
grp <- find.clusters(a, max.n.clust=15)

myCol <- c("darkblue","purple","green","orange","red","blue", "black")
scatter(dpca1, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, col=myCol, solid=.4,cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:7))

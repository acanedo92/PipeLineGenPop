library (vcfR)
library (adegenet)
library (hierfstat)
library (ape)

setwd ("Documents/ECOSUR/Cursos/GeneticaConservacion2018/PracticaSNPs/")

# Cargar un archivo tipo vcf y convertirlo en objeto genind para analizar con adegente. En este caso nuestro archivo esta en formato Genpop y no es necesrio transformarlo.
# a <- read.vcfR("FiltrosPlink/Filtro_LD.vcf")
# b <- vcfR2genind(a)
# c <- vcfR2genlight(a)

# Cargar archivo tipo genpop
a <- read.genepop("GP_withgrays.gen")

##################################
# Explorar nuestra base de datos #
##################################

a # Ver que tipo de objetos creamos
a@ploidy # Ver la haploidia de los datos
a@pop # Ver a las poblaciones en el objeto genind (@pop debe ser un factor)
a@tab # Es la matriz de todos los datos o genotipos
a@tab [c(1:5),c(1:5)] # Ver los datos en las primeras 5 filas y 5 columnas
a@tab [,1] # Ver el nombre de los individuos

###########################################
# Analisis basicos de diversidad genetica #
###########################################

divers <- basic.stats(a) # Analisis global, todas las pobalcaiones
divers$overall # Promedio de todos los datos. Todas las medidas de diversidad

# Para obtener los promedios por poblaci??n
colnames(divers$Ho) # Para ver el nombre de las poblaciones
nombresPobs <- colnames(divers$Ho) # Hacer un objeto con el nombre de las poblaciones
HO <- c(mean(divers$Ho[ ,1]), # Heterocigosis observada promedio de la poblaci??n 1
        mean(divers$Ho[ ,2]),
        mean(divers$Ho[ ,3]),
        mean(divers$Ho[ ,4]),
        mean(divers$Ho[ ,5]),
        mean(divers$Ho[ ,6]),
        mean(divers$Ho[ ,7]))

HO
divers$overall$Ho

# Hacer un grafico de barras con los valores de heterocigosis observada.
barplot (HO, xlab = "Populations", ylab="HO", col="red", names.arg = nombresPobs)
abline (h=0.1253) # Si queremos agregar una linea horizontal que se??ale el valor promedio, tomando el pormedio de la tabla divers$overall

##########################
# Podemos hacer lo mismo para obtener los valores promedio de las demas medidas de diversidad.
##########################

###################################
# Analisis de estructura genetica #
###################################

# Podemos estimar distancias para hacer heatmap de distancia gen??tica entre individuos
x.dist <- dist (a)
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

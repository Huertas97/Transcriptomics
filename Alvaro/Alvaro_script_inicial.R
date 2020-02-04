
setwd("D:/Programs/MEGA/Universidad/UAM/Segundo_semestre/Transcriptomica/Trabajo") 

# From Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("affy")
BiocManager::install("limma")
BiocManager::install("genefilter")

# Remove previous loaded objects to avoid errors
rm(list=ls())

# Load required libraries
library("affy")
library("limma")
library("genefilter")

# Target file 
targets <- readTargets("target.txt", row.names="FileName")

# Check that the .CEL files available in the directory are the one expected
CELfiles <- list.celfiles()
CELfiles
length(CELfiles) == 12

# Read .CEL files 
# data <- ReadAffy(filenames=targets$FileName)
rawdata <- ReadAffy(filenames = CELfiles)
rawdata

# Exploratory graphs that provide information about the quality of the arrays

library("RColorBrewer")

# Select 4 different color in the palette Set1. Return them in a vector
usr.col <- brewer.pal(4, "Set1")

# Repeat each color three times, such there are three arrays per sample type
mycols <- rep(usr.col, each = 3)

#########################
# Intensity distributions
#########################

# Show intensity distributrions with hist function 
hist(rawdata, lty = rep(1, length(CELfiles)), col = mycols,
     main = "Intensity distributions")

# Add legend in the top right corner 
legend("topright", rownames(pData(rawdata)),
       lty = rep(1, length(CELfiles)), 
       col = mycols,
       cex = 0.6)

# Tema 19 página 637. Vemos que no hay problemas y que no tenemos que descartar
# ningún array, aunque sí es importante normalizar. No hay ninguno con un valor
# muy alto en log intensity. También vemos que todos tienen una distribución 
# más o menos igual. 

#############################
# boxplot antes de normalizar
#############################

boxplot(rawdata,
        main="Boxplot Before Normalization",
        col = mycols,
        names = sampleNames(rawdata),
        las = 2,  # vertical labels
        cex.axis = 0.6, # size labels
        yaxt = "n") # do not show y labels
axis(2, cex.axis = 1) # plot y labels with higher size


# aquí ya quita el background. 
##################################################
# Ver si puedo eliminar yo el background yo solo
##################################################
eset <- expresso(rawdata,
                 bg.correct = TRUE, 
                 bgcorrect.method="rma", # corregir fondo
                 normalize = TRUE, 
                 normalize.method="quantiles", 
                 pmcorrect.method="pmonly", # los que machean perfectamente
                 summary.method="medianpolish", # resumeme por la mediana par aun mismo gen
                 verbose = TRUE,
) 

#########################
# Intensity distributions
#########################

# Show intensity distributrions
dens <- apply(data.frame(exprseset), 2, density)

plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")))
mapply(lines, dens, col=mycols)

# Add legend in the top right corner
legend("topright", legend=names(dens), lty = rep(1, length(CELfiles)), col = mycols, cex = 0.6)

## Hay que cambiar la separacion de los ejes

# boxplot de todos los genes después de normalizar
exprseset <- as.data.frame(exprs(eset))	# cojo las intensidades normalizadas	
boxplot(data.frame(exprseset),
        main="Boxplot After Normalization (log scale)",
        col = mycols,
        names = sampleNames(rawdata),
        las = 2,  # vertical labels
        cex.axis = 0.6, # size labels
        yaxt = "n") # do not show y labels
axis(2, cex.axis = 1) # plot y labels with higher size


## Comprobamos el filtro del rango intercuartílico
dim(exprs(eset))
dim(exprs(esetIQR))

# boxplot de solo los genes normalizados con un ranog intercuartílico superior a 0.5 para que tengan variabilidad
esetIQR <- varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)

exprseset <- as.data.frame(exprs(esetIQR))	# cojo las intensidades normalizadas	
boxplot(data.frame(exprseset),
        main="Boxplot After Normalization (log scale)",
        col = "red")


##################################################
## Añadir nuevos parámetros de control de calidad
# Vulcano Plot
##################################################


# design<-cbind(C1=c(rep(1, 3), rep(0, 9)), C2=c(rep(0,3), rep(1, 3), rep(0, 6)), S1=c(rep(0, 6), rep(1, 3), rep(0, 3)), S2 = c(rep(0, 9), rep(1, 3)) )
design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2,3,3,3,4,4,4)))
rownames(design)<-targets$FileName
colnames(design) <- c("C1", "C2", "S1", "S2")

## ajustamos cada nube de puntos correspondiente a cada grupo a una recta
fit<-lmFit(esetIQR,design)

## establecemos las comparaciones que queremos realizar
## Tenemos que dividir entre 2 para que tengan todos el mismo peso. Las dos últimas no tienen ningún sentido
contrast.matrix <- makeContrasts(C1-C2, S1-C1, S2-C2, ((S1+S2)/2)-((C1+C2)/2), ((S1+S2)/2)-C1, ((S1+S2)/2)-C2, levels=design)

## comparamos las regresiones lineales según el patrón de contraste pedido
fit2 <- contrasts.fit(fit, contrast.matrix)

fit2<-eBayes(fit2)

#Table with DEGs results number te cogerá todos los genes (filas) que tenia esetIQR

toptableIQR<-topTable(fit2, number=dim(exprs(esetIQR))[1], adjust.method="BH", coef = 3,  sort.by="p") # para poder ordenar por "p" necesito poner el coeficiente, es decir, que sepa que comparación debe coger
# toptableIQR<-toptable(fit2, number=dim(exprs(esetIQR))[1], adjust.method="BH", sort.by="p", coef = 1) DESAPARECERÁ

head(toptableIQR)



###
# script for plotting output from tuxedo suite ab initio RNA seq
# uses cuffdiff output

# clear memory

rm(list=ls())

#turn off graphics
graphics.off()

#source("https://bioconductor.org/biocLite.R")
#biocLite("cummeRbund")

#load packages
library(cummeRbund)
library(RSQLite)
library(ggplot2)
library(reshape2)
library(plyr)
library(fastcluster)
library(rtracklayer)
library(Gviz)
library(BiocGenerics)
library(Hmisc)


#loads in the cuffdiff data
myDir <- paste(getwd(), "/DE/", sep="")
cuff <- readCufflinks(dir=myDir)


#geneAnno <- annotation(genes(cuff))
#geneFPKM <- fpkm(genes(cuff))
#generepFPKM <- repFpkm(genes(cuff))
#isoformsFPMK <- fpkm(isoforms(cuff))
#fpkmMatrix(genes(cuff))
#repFpkmMatrix(genes(cuff))
#countMatrix(genes(cuff))
#featureNames(genes(cuff))
repInfo <- replicates(cuff)
geneDiff <- diffData(genes(cuff))
sampleNames <- samples(genes(cuff))

myReps <- TRUE
if (dim(replicates(cuff))[1] == dim(samples(cuff))[1] ){
  myReps <- FALSE
}

####
#script for DE of pipeline samples
#Samples <- genes(cuff)
#myDendro <- csDendro(Samples, replicates=T)
#rm(Samples)
dispPlot <- dispersionPlot(genes(cuff))
cvPlot <- fpkmSCVPlot(genes(cuff))
densPlot <- csDensity(genes(cuff))
densPlotrep <- csDensity(genes(cuff), replicates=T)
boxPlt <- csBoxplot(genes(cuff))
boxPltrep <- csBoxplot(genes(cuff), replicates=T)
scatPlot <- csScatterMatrix(genes(cuff), replicates=T) #csScatter(genes, s1, s2, smooth=T)

maList <- list()
for (i in 1:(length(sampleNames)-1)){
  for (j in (i+1):length(sampleNames)){
    maList[[length(maList)+1]] <- MAplot(genes(cuff), sampleNames[i], sampleNames[j])
  }
}

macountList <- list()
for (i in 1:(length(sampleNames)-1)){
  for (j in (i+1):length(sampleNames)){
    macountList[[length(macountList)+1]] <- MAplot(genes(cuff), sampleNames[i], sampleNames[j], useCount=T)
  }
}

volcList <- list(0)
for (i in 1:(length(sampleNames)-1)){
  for(j in (i+1):length(sampleNames)){
    volcList[[length(volcList)+1]] <- csVolcano(genes(cuff), sampleNames[i], sampleNames[j], alpha=0.05, showSignificant=T)
  }
}

###create gene subset of signficant genes from cuff data

sigGeneDiff <- geneDiff[ which(geneDiff$significant =="yes"), ]
sortedSGD <- sigGeneDiff[order(-abs(sigGeneDiff$log2_fold_change)),]
sortedSGD200 <- sortedSGD[1:200, ]
mySGD200geneIDs <- sortedSGD200$gene_id
myGenes200 <- getGenes(cuff, mySGD200geneIDs)

myHeatmap <- csHeatmap(myGenes200, cluster='row')
myHeatmaprep <-  csHeatmap(myGenes200, cluster='row', replicates=T)

#prints all plots to pdf
pdf(file="Differential_Expression_plots.pdf")

#print(myDendro)
csDendro(genes(cuff), replicates=T)
#print(repInfo)
print(dispPlot)

if(myReps == TRUE){
  print(cvPlot)
  par(mfrow=c(1,2))
  print(densPlot)
  print(densPlotrep)
  print(boxPlt)
  print(boxPltrep)
  par(mfrow=c(1,1))
}else{
  print(densPlot)
  print(boxPlt)
}

print(scatPlot)

for (i in 1:length(maList)){
  print(maList[[i]])
}

for (i in 1:length(macountList)){
  print(macountList[[i]])
}

for (i in 1:length(volcList)){
  print(volcList[[i]])
}

print(myHeatmap)
if(myReps == TRUE){
  print(myHeatmaprep)
}

dev.off()


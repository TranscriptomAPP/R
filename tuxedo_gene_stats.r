###
#

#clear memory
rm(list=ls())

#turn off graphics
graphics.off()

#options
options(scipen=999)

sampleFiles <- vector()
sampleFiles <- list.files(recursive=T, pattern="isoforms.fpkm_tracking")
#sampleFiles <- c("testing/cufflinks_out_2cells/isoforms.fpkm_tracking", "testing/cufflinks_out_6h/isoforms.fpkm_tracking")
sampleNames <- gsub("^.*out_", "", sampleFiles)
sampleNames <- gsub("\\/.*$", "", sampleNames)
plotCol <- rainbow(length(sampleNames))

sampleList <- list()
sampleSizes <- vector()
uniqSampleList <- list()
uniqSampleSizes <- vector()
for(i in 1:length(sampleFiles)){
  genes <- read.table(file=sampleFiles[i], sep="\t", stringsAsFactors = F, header=T)
  myTable <- genes[, c("gene_id", "length", "FPKM")]
  #create additional table with the longest isoform for each gene
  uniq_names <- unique(myTable[myTable$gene_id != "CUFF","gene_id"])
  
  myUniqGenesTable <- data.frame(matrix(nrow=0, ncol=3))
  colnames(myUniqGenesTable) <- colnames(myTable)
  for(j in 1:length(uniq_names)){
    tmpTable <- myTable[myTable$gene_id == uniq_names[j],]
    tmpmax <- max(tmpTable$length)
    tmpRow <- tmpTable[tmpTable$length==tmpmax,]
    myUniqGenesTable[j,] <- tmpRow
  }
  sampleList[[i]] <- myTable
  uniqSampleList[[i]] <- myUniqGenesTable
  sampleSizes <- c(sampleSizes, nrow(myTable))
  uniqSampleSizes <- c(uniqSampleSizes, nrow(myUniqGenesTable))
  names(sampleList)[i] <- sampleNames[i]
  names(uniqSampleList)[i] <- sampleNames[i]
}
maxSize <- max(sampleSizes)

lengthTable <- matrix(nrow=maxSize, ncol=length(sampleList))
colnames(lengthTable) <- sampleNames
fpkmTable <- matrix(nrow=maxSize, ncol=length(sampleList))
colnames(fpkmTable) <- sampleNames

pdf(file="RTranscripts.pdf")
for(i in 1:length(sampleList)){
  cols <- cm.colors(4)
  anno <- sampleList[[i]][which(sampleList[[i]]$gene_id != "CUFF"),"length"]
  length(anno) <- nrow(sampleList[[i]])
  noanno <- sampleList[[i]][which(sampleList[[i]]$gene_id == "CUFF"),"length"]
  length(noanno) <- nrow(sampleList[[i]])
  tmpTable <- cbind(sampleList[[i]]$length, anno, noanno)
  colnames(tmpTable) <- c("All", "Annotated", "Unannotated")
  boxplot(tmpTable, main=paste(sampleNames[i], "Transcript Length Distribution", sep=" "), col=cols, log="y", ylab="Length/bp - Log Scale")
  lengthAll <- sampleList[[i]]$length
  length(lengthAll) <- maxSize
  lengthTable[,i] <- lengthAll
  
  anno <- sampleList[[i]][which(sampleList[[i]]$gene_id != "CUFF"),"FPKM"]
  length(anno) <- nrow(sampleList[[i]])
  noanno <- sampleList[[i]][which(sampleList[[i]]$gene_id == "CUFF"),"FPKM"]
  length(noanno) <- nrow(sampleList[[i]])
  tmpTable <- cbind(sampleList[[i]]$FPKM, anno, noanno)
  colnames(tmpTable) <- c("All", "Annotated", "Unannotated")
  boxplot(tmpTable, main=paste(sampleNames[i], "Transcript FPKM Distribution", sep=" "), col=cols, ylab="FPKM")
  fpkmAll <- sampleList[[i]]$FPKM
  length(fpkmAll) <- maxSize
  fpkmTable[,i] <- fpkmAll
}


boxplot(lengthTable, col=plotCol, main="Transcript Length Distribution - All Samples", log="y", ylab="Length/bp - Log Scale")
boxplot(fpkmTable, cols=plotCol, main="Transcript FPKM Distribution - All Samples", ylab="FPKM")

# Plot number of isoforms per annotated gene
for(i in 1:length(sampleList)){
  genesVec <- sampleList[[i]][which(sampleList[[i]]$gene_id != "CUFF"),"gene_id"]
  genesTable <- table(genesVec)
  genesFreqVec <- as.vector(genesTable)
  genesFreqVec6max <- replace(genesFreqVec, genesFreqVec>5, 6)
  genesFreqPlot <- table(genesFreqVec6max)
  names(genesFreqPlot) <- replace(names(genesFreqPlot), names(genesFreqPlot)=="6", "6+")
  barplot(genesFreqPlot, main=paste("Annotated Genes - Isoform Count",sampleNames[i], sep="\n"), xlab="Number of Isoforms", ylab="Frequency", col=plotCol[i])
}


ltyvec <- rep(c(1,2), length(sampleNames))[1:length(sampleNames)]
plot(density(fpkmTable[,1]), col=plotCol[1], log="x", xlab="FPKM - Log Scale", main="Transcript FPKM - All Samples")
for(i in 2:length(sampleNames)){
  lines(density(fpkmTable[,i]), col=plotCol[i], log="x", lty=ltyvec[i])
}
legend("topright", sampleNames, col=plotCol, lty=ltyvec)

plot(density(lengthTable[,1]), col=plotCol[1], xlab="Length/bp", main="Transcript Length - All Samples")
for(i in 2:length(sampleNames)){
  lines(density(lengthTable[,i]), col=plotCol[i], lty=ltyvec[i])
}
legend("topright", sampleNames, col=plotCol, lty=ltyvec)

#plot fpkm distribution of longest transcript for each annotated gene
uniqMaxSize <- max(uniqSampleSizes)
uniqFpkmTable <- matrix(nrow=uniqMaxSize, ncol=length(uniqSampleList))
colnames(uniqFpkmTable) <- sampleNames
for(i in 1:length(uniqSampleList)){
  fpkmAll <- uniqSampleList[[i]]$FPKM
  length(fpkmAll) <- uniqMaxSize
  uniqFpkmTable[,i] <- fpkmAll
}  
plot(density(uniqFpkmTable[,1]), col=plotCol[1], log="x", xlab="FPKM - Log Scale", main="Longest Transcript FPKM - All Samples")
for(i in 2:length(sampleNames)){
  lines(density(uniqFpkmTable[,i]), col=plotCol[i], lty=ltyvec[i], log="x")
}
legend("topright", sampleNames, col=plotCol, lty=ltyvec)

#plot length distribution of longest transcript for each annotated gene
uniqLengthTable <- matrix(nrow=uniqMaxSize, ncol=length(uniqSampleList))
colnames(uniqLengthTable) <- sampleNames
for(i in 1:length(uniqSampleList)){
  lengthAll <- uniqSampleList[[i]]$length
  length(lengthAll) <- uniqMaxSize
  uniqLengthTable[,i] <- lengthAll
}  
plot(density(uniqLengthTable[,1]), col=plotCol[1], xlab="Length/bp", main="Longest Transcript Length - All Samples")
for(i in 2:length(sampleNames)){
  lines(density(uniqLengthTable[,i]), col=plotCol[i], lty=ltyvec[i])
}
legend("topright", sampleNames, col=plotCol, lty=ltyvec)



dev.off()

#x11()
par(mfrow=c(1,2))
boxplot(lengthTable, col=plotCol, main="Transcript Length Distribution\nAll Samples", log="y", ylab="Length/bp - Log Scale")
boxplot(fpkmTable, cols=plotCol, main="Transcript FPKM Distribution\nAll Samples", ylab="FPKM")



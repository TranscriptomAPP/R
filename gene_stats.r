###
#

#clear memory
rm(list=ls())

#turn off graphics
graphics.off()

#options
options(scipen=999)

sampleFiles <- vector()
#sampleFiles <- list.files(recursive=T, pattern="isoforms.fpkm_tracking")
sampleFiles <- c("testing/cufflinks_out_2cells/isoforms.fpkm_tracking", "testing/cufflinks_out_6h/isoforms.fpkm_tracking")
sampleNames <- gsub("^.*out_", "", sampleFiles)
sampleNames <- gsub("\\/.*$", "", sampleNames)
sampleList <- list()
plotCol <- rainbow(length(sampleNames))

sampleSizes <- vector()
for(i in 1:length(sampleFiles)){
  genes <- read.table(file=sampleFiles[i], sep="\t", stringsAsFactors = F, header=T)
  myTable <- genes[, c("gene_id", "length", "FPKM")]
  sampleList[[i]] <- myTable
  sampleSizes <- c(sampleSizes, nrow(myTable))
  names(sampleList)[i] <- sampleNames[i]
}
maxSize <- max(sampleSizes)

lengthTable <- matrix(nrow=maxSize, ncol=length(sampleList))
colnames(lengthTable) <- sampleNames
fpkmTable <- matrix(nrow=maxSize, ncol=length(sampleList))
colnames(fpkmTable) <- sampleNames

pdf(file="Transcript_stats.pdf")
for(i in 1:length(sampleList)){
  cols <- cm.colors(4)
  anno <- sampleList[[i]][which(sampleList[[i]]$gene_id != "CUFF"),"length"]
  length(anno) <- nrow(sampleList[[i]])
  noanno <- sampleList[[i]][which(sampleList[[i]]$gene_id == "CUFF"),"length"]
  length(noanno) <- nrow(sampleList[[i]])
  tmpTable <- cbind(sampleList[[i]]$length, anno, noanno)
  colnames(tmpTable) <- c("All", "Annotated", "Unannotated")
  boxplot(tmpTable, main=paste(sampleNames[i], "Transcript Length Distribution", sep=" "), col=cols, log="y", ylab="Log Scale")
  lengthAll <- sampleList[[i]]$length
  length(lengthAll) <- maxSize
  lengthTable[,i] <- lengthAll
  
  anno <- sampleList[[i]][which(sampleList[[i]]$gene_id != "CUFF"),"FPKM"]
  length(anno) <- nrow(sampleList[[i]])
  noanno <- sampleList[[i]][which(sampleList[[i]]$gene_id == "CUFF"),"FPKM"]
  length(noanno) <- nrow(sampleList[[i]])
  tmpTable <- cbind(sampleList[[i]]$FPKM, anno, noanno)
  colnames(tmpTable) <- c("All", "Annotated", "Unannotated")
  boxplot(tmpTable, main=paste(sampleNames[i], "Transcript FPKM Distribution", sep=" "), col=cols)
  fpkmAll <- sampleList[[i]]$FPKM
  length(fpkmAll) <- maxSize
  fpkmTable[,i] <- fpkmAll
}


boxplot(lengthTable, col=plotCol, main="Transcript Length Distribution - All Samples", log="y", ylab="Log Scale")
boxplot(fpkmTable, cols=plotCol, main="Transcript FPKM Distribution - All Samples")

ltyvec <- rep(c(1,2), length(sampleNames))
plot(density(fpkmTable[,1]), col=plotCol[1], log="x", xlab="FPKM - Log Scale", main="Transcript FPKM - All Samples")
for(i in 2:length(sampleNames)){
  lines(density(fpkmTable[,i]), col=plotCol[i], log="x", lty=ltyvec[i])
}
legend("topright", sampleNames, col=plotCol, lty=ltyvec)

plot(density(lengthTable[,1]), col=plotCol[1], xlab="Length", main="Transcript Length - All Samples")
for(i in 2:length(sampleNames)){
  lines(density(lengthTable[,i]), col=plotCol[i], lty=ltyvec[i])
}
legend("topright", sampleNames, col=plotCol, lty=ltyvec)


dev.off()

x11()
par(mfrow=c(1,2))
boxplot(lengthTable, col=plotCol, main="Transcript Length Distribution - All Samples", log="y", ylab="Log Scale")
boxplot(fpkmTable, cols=plotCol, main="Transcript FPKM Distribution - All Samples")



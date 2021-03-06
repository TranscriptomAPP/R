###
#script checking read mapping rate single end data.

#clear memory
rm(list=ls())

#close graphics
graphics.off()

#options
options(scipen=999)

#load in align_summary.txt data
sampleFiles <- vector()
sampleFiles <- list.files(recursive=T, pattern="align_summary.txt")
sampleNames <- gsub("^.*out_", "", sampleFiles)
sampleNames <- gsub("\\/.*$", "", sampleNames)

samplesMatrix <- matrix(ncol=3, nrow=length(sampleFiles))
colnames(samplesMatrix) <- c("Input Reads", "Mapped Reads", "Multiple Mapped")
rownames(samplesMatrix) <- sampleNames
plotCol <- rainbow(nrow(samplesMatrix))

sampMatrixPerc <- samplesMatrix


readMapPerc <- vector()
for(i in 1:length(sampleFiles)){
  readMap <- read.table(file=sampleFiles[i], sep="\t", stringsAsFactors = F)
  readMap2 <- gsub("^.*:", "", readMap[-c(1, 5), 1])
  readMap2 <- as.numeric(gsub("\\(.*$", "", readMap2))
  readMapPerc <- gsub("%.*$", "%", readMap[-c(1, 5), 1])
  readMapPerc <- gsub("^.*\\(", "", readMapPerc)
  readMapPerc[1] <- "100%"
  
  samplesMatrix[i, ] <- readMap2
  sampMatrixPerc[i, ] <- readMapPerc
}

pdf(file="Rmapping.pdf")

for(i in 1:length(sampleNames)){
  par(mar=c(9, 4, 4, 2)+0.1)
  coords <- barplot(samplesMatrix[i,], main=paste(sampleNames[i], "Read mapping", sep=" - "), col=plotCol[i], las=2, cex.names=0.9, ylim=c(0, max(samplesMatrix)*1.2), las=3, ylab="Count")
  text(coords, y=samplesMatrix[i,]+max(samplesMatrix)/25, labels=sampMatrixPerc[i,], cex=0.7)
}

par(mar=c(5, 4, 4, 7)+0.1, xpd=T)
mappingCols <- cm.colors(4)
coords <-barplot(t(samplesMatrix), beside=T, main="Total Read Pair Mapping", col=mappingCols, ylim=c(0, max(samplesMatrix)*1.2), xlab="Sample ID", ylab="Count")
legend("topright", inset=c(-0.25,0), colnames(samplesMatrix), fill=mappingCols)
text(t(coords), y=samplesMatrix+max(samplesMatrix)/25, labels=sampMatrixPerc, cex=0.7)

dev.off()

x11()
par(mar=c(5, 4, 4, 7)+0.1, xpd=T)
mappingCols <- cm.colors(4)
coords <-barplot(t(samplesMatrix), beside=T, main="Total Read Pair Mapping", col=mappingCols, ylim=c(0, max(samplesMatrix)*1.2), xlab="Sample ID", ylab="Count")
legend("topright", inset=c(-0.25,0), colnames(samplesMatrix), fill=mappingCols)
text(t(coords), y=samplesMatrix+max(samplesMatrix)/25, labels=sampMatrixPerc, cex=0.7)


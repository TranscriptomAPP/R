#clear memory

rm(list=ls())

#close graphics
graphics.off()

#options
options(scipen=999)

#load in read_stats.txt data
sampleFiles <- vector()
sampleFiles <- list.files(recursive=T, pattern="read_stats.txt")
#sampleFiles <- c("testing/trinity_out_2cells6h/rsem_out_2cells/read_stats.txt", "testing/trinity_out_2cells6h/rsem_out_6h/read_stats.txt")
sampleNames <- gsub("^.*out_", "", sampleFiles)
sampleNames <- gsub("\\/.*$", "", sampleNames)

trimFile <- list.files(recursive=T, pattern="trimmomatic_forR.txt")
#trimFile <- "trimmomatic_forR.txt"
trimTable <- read.table(trimFile, header=T, row.names=1)

samplesMatrix <- matrix(ncol=4, nrow=length(sampleFiles))
colnames(samplesMatrix) <- c("Total", "Mapped", "Concordant", "Disconcordant")
rownames(samplesMatrix) <- sampleNames
plotCol <- rainbow(length(sampleNames))
samplesPercMatrix <-samplesMatrix

pdf(file="Denovo_read_mapping_summary.pdf")

for(i in 1:length(sampleFiles)){
  readMap <- read.table(file=sampleFiles[i], sep="\t", stringsAsFactors = F, fill=T, blank.lines.skip=T, row.names=1)
  readMap[5, 1] <- sum(readMap[-5, 1])
  totalMap <- readMap[5, 1]
  readMap <- readMap[-5,]
  readMap <- rbind(readMap, c(trimTable[sampleNames[i], "BothSurv"]-totalMap, NA))
  readMap[, 2] <- prop.table(as.matrix(readMap[,1]), margin = 2)*100
  colnames(readMap) <- c("Count", "Percentage")
  readMap <- rbind(colSums(readMap[-5,]), readMap)
  readMap <- rbind(colSums(readMap[-1,]), readMap)
  rownames(readMap) <- c("Total", "Mapped", "Concordant", "Left only", "Right only", "Disconcordant", "Unmapped")
  
  samplesMatrix[i, ] <- readMap[c(1, 2, 3, 6), "Count"]
  samplesPercMatrix[i,] <- readMap[c(1, 2, 3, 6), "Percentage"]

  par(mar=c(9, 4, 4, 2)+0.1)
  coords <- barplot(readMap$Count, main=paste(sampleNames[i], "Read Pair Mapping", sep=" - "), col=plotCol[i], cex.names=0.9, ylim=c(0, max(readMap$Count)*1.2), las=2, las=3, names.arg=rownames(readMap))
  text(coords, y=readMap$Count+max(readMap$Count)/25, labels=paste(round(readMap$Percentage, digits=2), "%", sep=""), cex=0.7)
}

par(mar=c(9, 4, 4, 7)+0.1, xpd=T)
mappingCols <- cm.colors(4)
coords <-barplot(t(samplesMatrix), beside=T, main="Read Pair Mapping", col=mappingCols, ylim=c(0, max(samplesMatrix)*1.2), xlab="Sample ID", ylab="Count", las=2, las=3)
legend("topright", inset=c(-0.25,0), colnames(samplesMatrix), fill=mappingCols, cex=0.7)
text(t(coords), y=samplesMatrix+max(samplesMatrix)/25, labels=paste(round(samplesPercMatrix, digits=2), "%", sep=""), cex=0.7)

dev.off()



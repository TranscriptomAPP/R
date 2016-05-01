####
# script to generate plots for trinity transcript stats

#clear memory
rm(list=ls())

#cloes graphics
graphics.off()

#options
options(scipen=999)

sampleFiles <- vector()
sampleFiles <- list.files(recursive=T, pattern="RSEM.genes.results")
#sampleFiles <- c("testing/trinity_out_2cells6h/rsem_out_2cells/RSEM.genes.results", "testing/trinity_out_2cells6h/rsem_out_6h/RSEM.genes.results")
sampleNames <- gsub("^.*out_", "", sampleFiles)
sampleNames <- gsub("\\/.*$", "", sampleNames)
plotCol <- rainbow(length(sampleNames))

sampleList <- list()
sampleSizes <- vector()

for(i in 1:length(sampleFiles)){
  genes <- read.table(file=sampleFiles[i], sep="\t", stringsAsFactors = F, header=T)
  myTable <- genes[which(genes$FPKM > 0), ]
  sampleList[[i]] <- myTable
  names(sampleList)[i] <- sampleNames[i]
  sampleSizes <- c(sampleSizes, nrow(myTable))
  if(i==1){
    transcriptCount <- data.frame(matrix(nrow=length(genes$length), ncol=length(sampleNames)+2))
    colnames(transcriptCount) <- c("length", "Unexpressed", sampleNames)
    transcriptCount$length <- genes$length
    transcriptCount$Unexpressed <- rep(0, length(transcriptCount$length))
  }
  transcriptCount[, sampleNames[i]] <- genes$FPKM
  transcriptCount[, sampleNames[i]] <- replace(transcriptCount[, sampleNames[i]], transcriptCount[, sampleNames[i]]>0, 1)
}

pdf(file="Trinity_gene_stats.pdf")

#comment
transcriptCount[rowSums(transcriptCount[, sampleNames])==0, "Unexpressed"] <-1
totalCounts <- colSums(transcriptCount)
totalCounts[1] <- nrow(transcriptCount)
names(totalCounts)[1] <- "All"
barplot(totalCounts, main="Assembled Trinity Genes", ylab="Total Count", col=c(terrain.colors(2), plotCol))

maxSize <- max(sampleSizes)
lengthPlot <- data.frame(matrix(nrow=maxSize, ncol=length(sampleNames)))
colnames(lengthPlot) <- sampleNames
fpkmPlot <- lengthPlot
tpmPlot <- lengthPlot

for(i in 1:length(sampleFiles)){
  len <- sampleList[[i]]$length
  fpkm <- sampleList[[i]]$FPKM
  tpm <- sampleList[[i]]$TPM
  length(len) <- maxSize
  length(fpkm) <- maxSize
  length(tpm) <- maxSize
  lengthPlot[,i] <- len
  fpkmPlot[,i] <- fpkm
  tpmPlot[,i] <- tpm
}

boxplot(lengthPlot, col=plotCol, main="Transcript Length Distribution - All Samples", log="y", ylab="Log Scale")
boxplot(fpkmPlot, col=plotCol, main="Transcript FPKM Distribution - All Samples", log="y", ylab="Log Scale")
boxplot(tpmPlot, col=plotCol, main="Transcript TPM Distribution - All Samples", log="y", ylab="Log Scale")

maxDen <- vector()
myDens <- list()
myDens[[1]] <- density(transcriptCount$length)
maxDen <- c(maxDen, max(myDens[[1]]$y))
for(i in 2:ncol(transcriptCount)){
  myDens[[i]] <- density(transcriptCount[transcriptCount[,i]==1, "length"])
  maxDen <- c(maxDen, max(myDens[[i]]$y))
}
maxY <- max(maxDen)
minX <- min(transcriptCount$length)
maxX <- max(transcriptCount$length)
denCol <- c(terrain.colors(2), plotCol)
plot(myDens[[1]], xlim=c(minX, maxX), ylim=c(0, maxY), col=denCol[1], log="x", main="Trinity Genes Length Distribution")
for(i in 2:length(myDens)){
  lines(myDens[[i]], col=denCol[i])
}
legend("topright", colnames(transcriptCount), lty=1, col=denCol)

#calculate isoforms
for(i in 1: length(sampleFiles)){
  isoformList <- strsplit(sampleList[[i]]$transcript_id.s., ",", fixed=T)
  isoformCount <- unlist(lapply(isoformList, function(x) length(x)))
  isoformCountMax6 <- replace(isoformCount, isoformCount>5, 6)
  isoformFreqPlot <- table(isoformCountMax6)
  names(isoformFreqPlot) <- replace(names(isoformFreqPlot), names(isoformFreqPlot)=="6", "6+")
  barplot(isoformFreqPlot, main=paste("Assembled Trinity Genes - Isoform Count",sampleNames[i], sep="\n"), xlab="Number of Isoforms", ylab="Frequency", col=plotCol[i])
}
dev.off()

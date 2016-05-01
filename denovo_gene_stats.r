####
# script to generate plots for denovo transcript stats

#clear memory
rm(list=ls())
args <- commandArgs(TRUE)
print(commandArgs(TRUE)[1])

#cloes graphics
graphics.off()

#options
options(scipen=999)

args <- commandArgs(trailingOnly=T)

#if(length(args)==0){
  userFPKM <- 0
#}else{
#  userFPKM <- as.numeric(args[1])
#}

	sampleFiles <- vector()
	sampleFiles <- list.files(recursive=T, pattern=".genes.results")
	#sampleFiles <- c("testing/trinity_out_2cells6h/rsem_out_2cells/RSEM.genes.results", "testing/trinity_out_2cells6h/rsem_out_6h/RSEM.genes.results")
	sampleNames <- gsub("^.*out_", "", sampleFiles)
	sampleNames <- gsub("\\/.*$", "", sampleNames)
	plotCol <- rainbow(length(sampleNames))

sampleList <- list()
sampleSizes <- vector()

for(i in 1:length(sampleFiles)){
  genes <- read.table(file=sampleFiles[i], sep="\t", stringsAsFactors = F, header=T)
  myTable <- genes[which(genes$FPKM > userFPKM), ]
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

pdf(file="Denovo_gene_stats.pdf")

#comment
par(mar=c(9, 4, 4, 2)+0.1)
transcriptCount[rowSums(transcriptCount[, sampleNames])==0, "Unexpressed"] <-1
totalCounts <- colSums(transcriptCount)
totalCounts[1] <- nrow(transcriptCount)
names(totalCounts)[1] <- "All"
barplot(totalCounts, main=paste("Denovo Assembled Genes\n",commandArgs(TRUE)[1]), ylab="Total Count", col=c(terrain.colors(2), plotCol),las=2, las=3)

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

boxplot(lengthPlot, col=plotCol, main=paste("Transcript Length Distribution - All Samples\nFPKM filter = ", userFPKM, sep=""), log="y", ylab="Length/bp - Log Scale",las=3, las=2)
boxplot(fpkmPlot, col=plotCol, main=paste("Transcript FPKM Distribution - All Samples\nFPKM filter = ", userFPKM, sep=""), log="y", ylab="FPKM - Log Scale", las=3, las=2)
boxplot(tpmPlot, col=plotCol, main=paste("Transcript TPM Distribution - All Samples\nFPKM filter = ", userFPKM, sep=""), log="y", ylab="TPM - Log Scale", las=3, las=2)

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
plot(myDens[[1]], xlim=c(minX, maxX), ylim=c(0, maxY), col=denCol[1], log="x", main="Denovo Genes Length Distribution", xlab="Length/bp - Log Scale")
for(i in 2:length(myDens)){
  lines(myDens[[i]], col=denCol[i])
}
colnames(transcriptCount)[1] <- "All"
legend("topright", colnames(transcriptCount), lty=1, col=denCol)

#calculate isoforms
for(i in 1: length(sampleFiles)){
  isoformList <- strsplit(sampleList[[i]]$transcript_id.s., ",", fixed=T)
  isoformCount <- unlist(lapply(isoformList, function(x) length(x)))
  isoformCountMax6 <- replace(isoformCount, isoformCount>5, 6)
  isoformFreqPlot <- table(isoformCountMax6)
  names(isoformFreqPlot) <- replace(names(isoformFreqPlot), names(isoformFreqPlot)=="6", "6+")
  barplot(isoformFreqPlot, main=paste( sampleNames[i], " Isoform Count, FPKM ",userFPKM, sep=""), xlab="Number of Isoforms", ylab="Frequency", col=plotCol[i])
}
dev.off()



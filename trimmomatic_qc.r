#R script to plot trimmomatic QC metrics

#clears memory
rm(list=ls())

#turns off graphics
graphics.off()

options(scipen=999)

#reads in file and converts to matrix
qc <- read.table(file="trimmomatic_forR.txt", header=F, sep="\t",stringsAsFactors=F,row.names=1, check.names = FALSE)
rn <- rownames(qc)
qc <-as.matrix(qc)
qc <- c(qc[1], qc[2],qc[3],qc[4],qc[5])
pr <- c((qc[1]/qc[1])*100, (qc[2]/qc[1])*100,(qc[3]/qc[1])*100, (qc[4]/qc[1])*100, (qc[5]/qc[1])*100)
pr <- round (pr, digits=2)

#opens pdf for printing
pdf(file=paste(rn, ".pdf", sep = ""))
par(mar=c(9, 4, 4, 2)+0.1)
  coords <- barplot(qc, main=paste("Trimming results", rn, sep =" - "),names.arg=c("Raw Reads Pairs", "Trimmed Pairs", "Trimmed Fwd. Only", "Trimmed Rev. Only", "Dropped"), ylab="Read Pair Count", ylim=c(0,max(qc)*1.2), las=2, las=3)
  text(coords, y=qc+max(qc)/35 , labels=paste(pr, "%", sep=""), cex=0.7)



#closes pdf
dev.off()



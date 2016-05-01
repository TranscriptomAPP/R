args <- commandArgs(TRUE)
# folder where the script is executed
# plot the nb of transcripts vs TPM 
setwd(args[1])
png('Nb_Transcripts_per_TPM.png')
data = read.table("Nb_Transcripts_per_TPM.txt", header=T)
plot(data, xlim=c(-100,0), ylim=c(0,100000), t='b',main="Expressed transcripts above a minimum TPM threshold ")
dev.off()

#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

filename <- args[1]
data = read.table(filename, header=T, stringsAsFactors=F)
idx = order(data[,"TPM"], decreasing=T)
table=data[idx[1:10], c("gene_id", "expected_count", "TPM")]
write.csv(table, file = "top10_gene.csv")



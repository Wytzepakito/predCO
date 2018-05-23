#!/usr/bin/Rscript

library(DNAshapeR)
library(seqinr)

args <- commandArgs(trailingOnly=T)  

file = args[1] #fasta file

pred = getShape(file)

pdf(paste(file, 'DNAshapes_plot.pdf', sep="."))

par(mfrow=c(2,2))

for (i in 1:4) {plotShape(pred[[i]], colLine= 'red', main=names(pred)[i])} 

dev.off()


fas = read.fasta(file, seqtype="DNA", as.string=T)

seqlen= c()

for (i in 1:4) { 

	rownames(pred[[i]]) = names(fas)

	pred[[i]] = data.frame(pred[[i]])

	seqlen[i] = ncol(pred[[i]]) # sequence length 
}


for (i in 1:4) {

	pred[[i]] = cbind(pred[[i]], 	mean=apply(pred[[i]], 1, function(x) mean(x, na.rm=T)), 
					max= apply(pred[[i]], 1, function(x) max(x, na.rm=T)), 
					min= apply(pred[[i]], 1, function(x) min(x, na.rm=T)))

}

for (i in 1:4) { 

	n = ncol(pred[[i]])

	for (j in 1:ncol(pred[[i]])) { colnames(pred[[i]])[j] = paste(names(pred)[i], colnames(pred[[i]])[j], sep=".") }

	write.table(pred[[i]][,(n-2):n], file=paste(file, names(pred)[i], "txt" ,sep="."), sep="\t", col.names=T, row.names=F,append=FALSE, quote=FALSE) 
}

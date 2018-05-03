#!/usr/bin/Rscript

### works in R3.3

library(DNAshapeR)
library(seqinr)

args <- commandArgs(trailingOnly=T)  

## args[1] = fasta file 

file = args[1]

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


#write.table(pred, file=paste(file,"dnashaper.out",sep="."), append=FALSE, quote=FALSE) 



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

#### BINNING

#bi= 10  # bin interval length


#for (i in 1:4) { 

#	bn= as.integer(seqlen[i] / bi)  # number of bins

#	n = ncol(pred[[i]])

#	for (k in 0:(bn-1)) { 
#		pred[[i]][,n+k+1] = apply(pred[[i]][,(k*bi+1):((k+1)*bi)], 1, function(x) mean(x, na.rm=T))
#	}
#	
#	for (j in 1:ncol(pred[[i]])) { colnames(pred[[i]])[j] = paste(names(pred)[i], colnames(pred[[i]])[j], sep=".") }

#	write.table(pred[[i]][,(n-2):(n+bn)], file=paste(file,"binned", bi ,names(pred)[i], "txt" ,sep="."), sep="\t", col.names=T, row.names=T,append=FALSE, quote=FALSE) 

#} ## it is from 1:10 .. 11:20 






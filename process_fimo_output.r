#!/usr/bin/Rscript

# process fimo.txt from fimo output 

## assumptions: 1) the motifs are string and not a pwm which gives same score for every match. 2) the motifs can overlap then weighted score is reported.

options(warn=1)

args <- commandArgs(trailingOnly=T) 

data = read.table(args[1], header = F)

pos = read.table(args[2], header= F)

colnames(pos) = c("chrom", "start", "end")

pos[,"id"] = paste0(pos$chrom, ":", pos$start, "-", pos$end )

colnames(data) = c("pattern_name", "sequence_name", "start", "stop", "strand", "score", "p-value", "q-value", "matched_sequence")

motifs = as.character(unique(data[,"pattern_name"]))

output = data.frame()  # make an N by 3 dataframe, where N is the number of cases in the dataset fed to fimo. 

output = data.frame(matrix(vector(), nrow(pos), 0), stringsAsFactors=F)

rownames(output) = pos$id  # maybe convert to list

names = list()

for (i in 1:length(motifs)){

	names = c(names, paste0(motifs[i], "_presence"), paste0(motifs[i], "_count"), paste0(motifs[i], "_score_weighted"))
}

#colnames(output) = names

for (col in names){

	output[,col] = 0

}

for (mot in motifs){

	instances = as.character(unique(data[which(data$pattern_name == mot),"sequence_name"]))

	for (ins in instances){

		insDf = data[which(data$sequence_name == ins & data$pattern_name == mot),]

		insDf = insDf[order(insDf$start),]  ## sort by start 

		output[ins, paste0(mot, "_presence")] = 1 

		count= 1 # number of non-overlapping motifs

		insDf_idx =c(1)  # non-ov motif start location per instance. 

		if (nrow(insDf) > 1){

			for (i in 1:(nrow(insDf)-1)){
		  
				# non-overlap 
				if (insDf[i, "stop"] < insDf[(i+1), "start"]){
		  
					count = count +1
			  
					insDf_idx = c(insDf_idx, i+1)

				}   
			}
		}
		
		ov_size = 0
		
		for (m in 1:length(insDf_idx)){ # for each non-ov motif
		  
			first_loc = insDf[insDf_idx[m], "start"]
		        
			if (m == length(insDf_idx)){
		  
				last_loc = insDf[nrow(insDf), "stop"]
		  
			} else {
		  
				last_loc = insDf[(insDf_idx[m+1]-1), "stop"]
		  
			}  
		  
			ov_size = ov_size + (last_loc - first_loc +1) # total size where motifs covered. 

			motif_score = insDf[1,"score"]

			motif_size = length(mot)

			unit_score = motif_score / motif_size 
		  
			weighted_score = unit_score * ov_size 
		  
			# assume all the scores are same due to stricht motif search (i.e. not the pwm but search the string)
	   
	  	}
	  
		output[ins, paste0(mot, "_count")] = count
	  
		output[ins, paste0(mot, "_score_weighted")] = weighted_score

	}

}

write.table(output, args[3], sep="\t", quote=F, row.names=T, col.names=T)



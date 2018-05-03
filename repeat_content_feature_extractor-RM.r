
#!/usr/bin/Rscript

#setwd("/home/demir004/scratch/prediction/data/rice/")

args <- commandArgs(trailingOnly=T)  

## for rice: 
#sed '1,3d' ref_genome/masked/IRGSP-1.0_genome.fasta.out | awk -v OFS="\t" '{print $5, $6-1, $7, $11}' > ref_genome/IRGSP-1.0_repeat_aggressive.bed
#sed -i 's/?//g' IRGSP-1.0_repeat_aggressive.bed
#cut -f4 ref_genome/IRGSP-1.0_repeat_aggressive.bed | sort | uniq > ref_genome/aggresive_repeats_list.txt

##--- input preparation 
##for file in positive.all.4000.max20NN.maxoverlap1kb.bed negative.all.4000.max20NN.maxoverlap1kb.bed
##do 
##    bedtools intersect -wo -a $file -b ref_genome/*repeats_aggressive.bed > $file.repeats.overlap.out



#args = c("positive.all.4000.max20NN.maxoverlap1kb.bed.repeats.overlap.out",
#        "positive.all.4000.max20NN.maxoverlap1kb.bed",
#        "positive.all.4000.max20NN.maxoverlap1kb.bed.repeat.content.txt",
#        "ref_genome/aggresive_repeats_list.txt")




file = args[1]

ori_file = args[2]

outfile = args[3]

namesfile = args[4]


df = read.table(file, header=F, sep="\t", stringsAsFactors=F)

colnames(df) = c("region_chrom", "region_start", "region_end", 
                                            "chrom", "start", "end", "type", "overlap")
                                             

ori = read.table(ori_file, header=F, sep="\t", stringsAsFactors=F)

colnames(ori) = c("chrom", "start", "end")


names = read.table(namesfile, header=F, sep="\t", stringsAsFactors=F)

names = names[,1]



for (i in 1:nrow(df)){
    
    df[i,'id'] = paste(df[i,'region_chrom'], ":", df[i, 'region_start'], "-", df[i, 'region_end'], sep="")
    
}

for (i in 1:nrow(ori)){
    
    ori[i,'id'] = paste(ori[i,'chrom'], ":", ori[i, 'start'], "-", ori[i, 'end'], sep="")
    
}



#idlist = unique(df[,'id']) 

idlist = ori$id  ### take the id list from original file 



dfout = data.frame(matrix(NA, ncol = 3, nrow = length(idlist)))

dfout[,'id'] = idlist

count = 0 

for (i in 1:nrow(dfout)){
#for (i in 1:2){
    
    
    id = dfout[i,'id']
    
    dfI = df[which(id == df[,'id']),]
    
    #//  this section counts the overlapping hits whose repeat class is different 
    #    for example if LTR and LTR/Copia overlaps for the same hit.. 
    
    if (nrow(dfI) > 1){ 
    
        dfI = dfI[order(dfI$start),]

        for (k in 1:(nrow(dfI)-1)){

            if (dfI[k, 'end'] >= dfI[k+1,'start'] ){

                if (dfI[k,'type'] != dfI[k+1, 'type']){
                    
                    count = count + 1

#                    print("overlap found")
#                    print(dfI[k,])
#                    print(dfI[k+1,])

                }
            }
        } 
    }
    
    #// 

    
#    types = unique(dfI[,'type'])  
    types = names
    
    for (type in types){
        
        dfIt = dfI[which(dfI[,'type'] == type),]
        
        if (nrow(dfIt) == 0) {
            
            dfout[i, type] = 0
            
            next
        }
        
        if (nrow(dfIt) == 1){
            
            ratio = dfIt[1,'overlap'] / (dfIt[1,'region_end'] - dfIt[1,'region_start']) 

            dfout[i, type] = ratio

            
        } else {  ## MERGING 
            
           # sort by start 
            # check if it is overlaps to the next one
            # if not overlaps, save the first ones  'overlap' value to somewhere
            # if it is overlaps, save the distance from first row start to second row start 
        
           
            
             s = dfIt[order(dfIt$start),] # sorted 
            
             
            
            for (k in 1:nrow(s)){   # merging 
                
                n = nrow(s)
            
                j=1 
                
                I = c()
                
            while (j < n){
                
               if( s[j, 'end'] >= s[j+1, 'start']){ # overlapping
                   
                   s[j, 'start'] = min(s[j, 'start' ], s[j+1, 'start'])
                   
                   s[j, "end"] = max(s[j, 'end'], s[j+1, 'end'])
                
                   I = c(I, j)
                   
                   j = j+ 2
               } else {
                   
                   I =c(I, j)
                   
                   j = j+1 
               } 
                
                if (j == n)  I = c(I, j)
                
                
            }
                
            if (length(I) != 0) s = s[I,]   # update s 
                
            }
            
            
   
                
            # find the new overlap of each region in s  
            
            for (j in 1:nrow(s)){
                
                st = s[j, 'start'] -1 
                end = s[j, 'end'] -1
                rst = s[j, 'region_start']
                rend = s[j, 'region_end']
                
                
                l = sort(c(st, end, rst, rend))
                
                s[j, 'overlap'] = l[3] - l[2]
                
            }

            
            dfout[i, type] = sum(s[,'overlap']) / (dfIt[1,'region_end'] - dfIt[1,'region_start']) 
            
        }
          
        
    }
    
}

dfout = dfout[,4:ncol(dfout)]






## for NA count

#for (j in 1:ncol(df)){ print(colnames(dfout)[j]); print(length(which(is.na(dfout[,j]) ==T)))}


## for 0 count 

for (j in 1:ncol(dfout)) { print(colnames(dfout)[j]); 
                          print(length(which(dfout[,j] == 0)) / nrow(dfout) * 100) }






# dfout[is.na(dfout)] = 0 ### replace NA's with 0

# out = dfout[, c("id", "LTR", "LTR/Copia", "LTR/Gypsy", "Simple_repeat", "Low_complexity" )]

out = dfout # the column names are same order as in names list. 



write.table(out, file=outfile, 
            col.names=T, row.names=F, sep="\t", quote=F, )










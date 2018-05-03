#!/usr/bin/bash 

### DATA PREPARATION ####

do_genome=false

prepare_positive_bed=false

negative_set="pericentromere"  #sample_gene // genome // pericentromere

## Required scripts: 

# find_NNs_v3_print_bed.py	## takes *genome.fasta and prepares NNs.bed file 

# write_genome_lengths.py	## takes *genome.fasta and prepares sizes of genome file 


SIZE=4000

src="/path/to/scripts"

REF=$species"/ref_genome"

RscriptPath="/path/to/bin/R-3.3.2/bin/Rscript"


if true ; then 

    newFolder=$2

    if [ ! -d "$newFolder" ]; then

        mkdir $newFolder

    fi

    cd $newFolder

fi


###############################
## PREPARE GENOME files #######
###############################

## only once 

if $do_genome

    then 

    python $src/find_NNs_v3_print_bed.py $REF/*genome.fasta > $REF/NNs.bed

    python $src/write_genome_lengths.py $REF/*genome.fasta > $REF/lengths.genome

    ## change lengths.genome to lengths.chr.genome to use only chromosomes in bedtools 

fi 

######################
## POSITIVE SET ######
######################

if $prepare_positive_bed

then 

#for maize add Chr to the begining of chr numbers
  #sed '1,3d' maize_ncomms7648-s3.csv | cut -f5,2-3 | sed 's/chr//g' |awk '{printf "Chr%02i\t", $3} {print $1"\t"$2}' > maize.ncomms7648.CO.csv


  ###for bed style but directly from excel file. Assumed to be 1-based ## for arabidopsis

  awk '$3-$2 <= '"$SIZE"'-2000' *CO.bed | awk -F'\t' -v OFS='\t' '{mid=int(($3-$2)/2)+$2; $2=mid-('"$SIZE"'/2)-1; $3=mid+('"$SIZE"'/2); print $0}' | sort -k1,1 -k2,2n | uniq > CO.$SIZE.bed

  bedtools intersect -v -a CO.$SIZE.bed -b $REF/*NNs.bed > positive.bed


  bedtools merge -d -1000 -i positive.bed | awk '$3-$2 > 4001' > positive.bed.overlapping.bed

  bedtools merge -d -1000 -i positive.bed | awk '$3-$2 == 4001' > positive.bed.non-overlapping.mrthn1kb.bed

  awk -v OFS='\t' 'function roll(n) {return 1+ int(rand()*n)} {if (roll(2)==1) print $1, $2, $2+4001; else print $1, $3-4001, $3}' positive.bed.overlapping.bed > positive.bed.overlapping_selected.bed

  cat positive.bed.overlapping_selected.bed positive.bed.non-overlapping.mrthn1kb.bed | sort -k1,1 -k2,2n > positive.all.4000.maxoverlap1kb.bed

else 

  cp $REF/../positive.all.4000.maxoverlap1kb.bed ./

fi 

## get the sequences

bedtools getfasta -fi $REF/*genome.fasta -bed positive.all.4000.maxoverlap1kb.bed -fo positive.all.4000.maxoverlap1kb.fasta

#for genomes that are masked and replaced with lower-case letters
sed 's/g/G/g' positive.all.4000.maxoverlap1kb.fasta | sed 's/a/A/g' | sed 's/t/T/g' | sed 's/c/C/g' > positive.all.4000.maxoverlap1kb.tmp

mv positive.all.4000.maxoverlap1kb.tmp positive.all.4000.maxoverlap1kb.fasta




######################
## NEGATIVE SET ######
######################


## Negative set from random regions:

if [ $negative_set != "none" ]

then 

  ### !!!! EXLUDE FROM ALL CO AND THEIR NEAR +-!kb 

  cat positive.all.4000.maxoverlap1kb.bed $REF/*NNs.bed | sort -k1,1 -k2,2n > exclude_region.bed

  bedtools merge -i exclude_region.bed | sort -k1,1 -k2,2n > exclude_region.merged.bed


  if [ $negative_set == "genome" ]  # -eq

  then 

    bedtools shuffle -excl exclude_region.merged.bed -chrom -noOverlapping -g $REF/*.genome -maxTries 10000000 -seed 1234567 -i positive.all.4000.maxoverlap1kb.bed | sort -k1,1 -k2,2n > negative.all.4000.maxoverlap1kb.bed

    # Try  differrent seed to have all instances in the genome so that negative set and positive set has the same size

    while [ "$(wc -l negative.all.4000.maxoverlap1kb.bed | cut -d" " -f1)" -ne  "$(wc -l positive.all.4000.maxoverlap1kb.bed | cut -d" " -f1)" ]

    do 

      bedtools shuffle -excl exclude_region.merged.bed -chrom -noOverlapping -g $REF/*.genome -maxTries 10000000 -i positive.all.4000.maxoverlap1kb.bed | sort -k1,1 -k2,2n > negative.all.4000.maxoverlap1kb.bed

    done 
    
    
  
  elif [ $negative_set == "sample_gene" ] # -eq

  then 

    python3 $src/density_genes_general_with_opt_v3.py positive.all.4000.maxoverlap1kb.bed $REF/lengths.genome exclude_region.merged.bed $REF/gene.mid.txt ./  # Bandwith is optimized inside the script

    sort -k1,1 -k2,2n sampled.bed > negative.all.4000.maxoverlap1kb.bed
    
    
  elif [ $negative_set == "pericentromere" ]
  
  then 
  
    grep pericentromere $REF/genome.bed | bedtools merge -i - > $REF/pericentromere.bed
    
    bedtools subtract -a $REF/pericentromere.bed -b $REF/*NNs.bed > iclude_regions.bed

    bedtools shuffle -incl iclude_regions.bed -chrom -noOverlapping -g $REF/*.genome -maxTries 10000000 -i positive.all.4000.maxoverlap1kb.bed | sort -k1,1 -k2,2n > negative.all.4000.maxoverlap1kb.bed
  

  fi 

  ### If negative set is already there

  bedtools getfasta -fi $REF/*genome.fasta -bed negative.all.4000.maxoverlap1kb.bed -fo negative.all.4000.maxoverlap1kb.fasta


	  if [ "$?" != "0" ]; then

		  echo "Error in bedtools getfasta" 1>&2
		  exit 1
	  fi

  sed 's/g/G/g' negative.all.4000.maxoverlap1kb.fasta | sed 's/a/A/g' | sed 's/t/T/g' | sed 's/c/C/g' > negative.all.4000.maxoverlap1kb.tmp

  mv negative.all.4000.maxoverlap1kb.tmp negative.all.4000.maxoverlap1kb.fasta

fi 

if [ $REF != "arabidopsis/ref_genome" ]

then 

    for file in *tive.all.4000.maxoverlap1kb.fasta

    do 
            sed -i 's/Ch/ch/g' $file # for tomato and rice ; maize doesn't have chr in front of chromosomes

    done
    
fi

######################
## FEATURES  #########
######################


####SEQUENCE BASED
##################


## Dinucleotide freqs   ## NNs are not taken into account

# gives 17 column file: sequence name + 16 dinuc features. 

for file in negative.all.4000.maxoverlap1kb.fasta positive.all.4000.maxoverlap1kb.fasta

do 
	python $src/calc_dinuc_freq.py $file > $file.dinuc.txt 

done 


## DNA structure information


for file in *tive.all.4000.maxoverlap1kb.fasta 

do 

	$RscriptPath $src/DNAshapeR.r $file # gives 4 files each 

	paste $file.dinuc.txt $file.HelT.txt $file.MGW.txt $file.ProT.txt $file.Roll.txt > $file.fromSeq.txt

done 


### ANNOTATION BASED
####################

## distance to the nearest transcription start site


## for rice 

#awk -v OFS="\t" '{if ($7=="+") print $1,$4-1,$4,$3,$6,$7; if ($7=="-") print $1,$5,$5+1,$3,$6,$7}' $REF/*five_prime_UTR.gff3 | sort -k1,1 -k2,2n > $REF/five_prime_UTR.bed

for file in *tive.all.4000.maxoverlap1kb.bed

do 

	awk -v OFS='\t' '{diff=int(($3-$2)/2); mid=$2+diff; print $1, mid-1, mid, $1":"$2"-"$3}' $file | sort -k1,1 -k2,2n > $file.mids.bed
	
        # for tomato : 
        #cp tomato/ref_genome/ITAG2.4_gene_models.mRNA.start.bed tomato/ref_genome/five_prime_UTR.bed


	bedtools closest -nonamecheck -D b -t first -a $file.mids.bed -b $REF/five_prime_UTR.bed | awk '{print $4"\t"$NF}' | uniq > $file.temp

	echo -e "#name\tdistanceTSS" | cat - $file.temp > $file.mids.bed.closest.tss.txt

	rm $file.temp

done 

## gene coverage & repeat coverage


## for rice: prepare the repeatmasker out file 
## go to maize/ref_genome 
#	sed '1,3d' masked/*genome.fasta.out | awk -v OFS="\t" '{print $5, $6-1, $7, $11}' > repeats_aggressive.bed
#	grep -v scaffold repeats_aggressive.bed > repeats_aggressive_chr.bed
#	sed -i 's/?//g' repeats_aggressive_chr.bed
#	cut -f4 repeats_aggressive_chr.bed | sort | uniq > aggressive_repeats_list.txt
#	nano aggressive_repeats_list.txt # to remove artefact and Unknown if present


#cut -f3 $REF/Zea_mays.AGPv3.21.gff3 | sort | uniq | grep -v '##' > genelist

#for maize:
#cp maize/ref_genome/Zea_mays.AGPv3.21.gff3 maize/ref_genome/AGPv3.21_gene_models.gff3


for file in *tive.all.4000.maxoverlap1kb.bed

do 

	bedtools intersect -nonamecheck -wo -a $file -b $REF/*gene_models.gff3 > $file.genes.overlap.out  ## -wao option will print all a entries: A features w/o overlap are also reported with a NULL B feature and overlap = 0.

	sed -i 's/#/_/g' $file.genes.overlap.out  ### these changes are necessary to load the data to R
	sed -i 's/%/_/g' $file.genes.overlap.out

	Rscript $src/gene_content_feature_extractor.r $file.genes.overlap.out $file $file.gene.content.txt ~/pred/genelist

	# also for tomato:
	# mv ITAG2.4_repeats_aggressive_type.bed ITAG2.4_repeats_aggressive.bed
	
	bedtools intersect -nonamecheck -wo -a $file -b $REF/*repeats_aggressive.bed > $file.repeats.overlap.out
	
	## repeat_content_feature_extractor.r accepts 4 arguments: 
	## - 1 out file (bed + gff from bedtools), for tomato ## for rice: out file is bed + bed from bedtools
	## - 2 bed file, 
	## - 3 outfile name, 
	## - 4 list of all posible repeat classes

	Rscript $src/repeat_content_feature_extractor-RM.r $file.repeats.overlap.out $file $file.repeat.content.txt $REF/aggressive_repeats_list.txt

	rm $file.repeats.overlap.out $file.genes.overlap.out 

done

#for file in *.content.txt ; do awk '{if (NR==1) print "#"$0; else print $0}' $file > $file.header; mv $file.header $file; done



for file in *.all.4000.maxoverlap1kb.bed 
do 

paste $file.mids.bed.closest.tss.txt $file.gene.content.txt $file.repeat.content.txt | cut -f1,2,4-8,10- > $file.fromAnn.txt

done



## prepare motifs

for file in positive negative

do

fimo --oc ${file}_fimo/ --verbosity 1 --thresh 1.0E-4 ~/pred/motifs.meme $file.all.4000.maxoverlap1kb.fasta

Rscript $src/process_fimo_output.r ${file}_fimo/fimo.txt $file.all.4000.maxoverlap1kb.bed ${file}_fimo/$file.all.4000.maxoverlap1kb.bed.fimo.txt

awk '{if (NR==1) print "name\t"$0; else print $0}' ${file}_fimo/$file.all.4000.maxoverlap1kb.bed.fimo.txt > ${file}_fimo/$file.all.4000.maxoverlap1kb.bed.fimo.txt.header 
mv ${file}_fimo/$file.all.4000.maxoverlap1kb.bed.fimo.txt.header $file.all.4000.maxoverlap1kb.bed.fimo.txt

done

## combine all features 

for file in positive negative  
do 

cut -f1 $file.all.4000.maxoverlap1kb.bed.fimo.txt > $file.index

cut -f2- $file.all.4000.maxoverlap1kb.fasta.fromSeq.txt > $file.Seq.txt 

cut -f2- $file.all.4000.maxoverlap1kb.bed.fromAnn.txt > $file.Ann.txt

cut -f2- $file.all.4000.maxoverlap1kb.bed.fimo.txt > $file.fimo.txt

paste $file.index $file.Seq.txt $file.Ann.txt $file.fimo.txt > $file.features.txt

rm $file.Seq.txt $file.Ann.txt $file.fimo.txt 

done


python3 $src/prediction_from_features.py ./



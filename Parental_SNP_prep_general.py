#!/usr/bin/pyhton3

## takes 2 vcf files sorted!!!! 
## combines the vcf files and prints the SNPs between 2 files as 

##EXAMPLE INPUT:  should be no header!!!!

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  calmd

#SL2.50ch00      3235    .       A       G       228     .       
#DP=38;VDB=0.856678;SGB=-0.693143;MQSB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,19,19;MQ=60
#    GT:PL:DP        1/1:255,114,0:38

#SL2.50ch00      56776   .       AGGGGGGGGGGGG   AGGGGGGGGGGG    22.2442 .       
#INDEL;IDV=13;IMF=0.65;DP=20;VDB=0.431879;SGB=-0.662043;MQSB=0.974597;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=60 
#GT:PL:DP        1/1:49,27,0:9


##EXAMPLE OUTPUT: 
#tab-separated bed file: col4: homozygous count; col5: heterozygous count

#SL2.50ch00	578	579	1	0
#SL2.50ch00	1795	1796	1	0
#SL2.50ch00	2013	2014	1	0
#SL2.50ch00	2771	2772	1	0
#SL2.50ch00	2845	2846	1	0


import sys 
import os
import subprocess


#file1 = path + "ERR418039.bam.filtered.bam.dedup.bam.bcf.vcf.SNPs.10"  # sys.argv[1]   ## preferably smaller file
 
#file2 = path + "CGN15528.bam.filtered.bam.dedup.bam.bcf.vcf.SNPs.10"   # sys.argv[2]

file1 = sys.argv[1]

file2 = sys.argv[2]


def print_line(haplotype, chrom, pos):
    
    "This prints bed formatted entry with its haplotype"
    
    homo = str(1) + "\t" + str(0)  ## you can define this as global on top, but now it is not necessary :p

    hetero = str(0) + "\t" + str(1)
    
    if haplotype == "1/1": 
            
        print(chrom + "\t" + str(pos-1) + "\t" + str(pos) + '\t' + homo)
            
    else: 
            
        print(chrom + "\t" + str(pos-1) + "\t" + str(pos) + '\t' + hetero)
        
    return


def process_line(line):

    "this function returns the haplotype, .. from a gff file line"

    l1 = line.strip().split('\t')

    hap = l1[9].split(":")[0] # could be 1/1, 0/1, 1/2

    chrom = l1[0]

    pos = int(l1[1])

    ALT = l1[4]

    return hap, chrom, pos, ALT 


f1 = open(file1, 'r')

f2 = open(file2, 'r')


line1 = f1.readline()

line2 = f2.readline()

homo = str(1) + "\t" + str(0)

hetero = str(0) + "\t" + str(1)


while line1 and line2:

    hap1, chrom1, pos1, ALT1 = process_line(line1)
    
    hap2, chrom2, pos2, ALT2 = process_line(line2)
    
    
    
    if chrom1 == chrom2 and pos1 == pos2:
        
       
        ## process the alleles, write the bed formatted file 
        
#        if hap1 == hap2 and hap2 == "1/1" and ALT1 == ALT2:

            ## this is not a snp 
            
        if hap1 == hap2 and hap2 == "1/1" and ALT1 != ALT2:
            
            ## this is homozygous snp 
            
            print(chrom1 + "\t" + str(pos1-1) + "\t" + str(pos1) + '\t' + homo) 
            
        elif hap1 != hap2:
            
            ## this is heterozygous snp 
            
            print(chrom1 + "\t" + str(pos1-1) + "\t" + str(pos1) + '\t' + hetero)
            
        
        line1 = f1.readline()

        line2 = f2.readline()
        
    elif chrom1 == chrom2 and pos1 < pos2:
        
        # write the small one = line1
        
        print_line(hap1, chrom1, pos1)
        
        line1 = f1.readline()
            
    elif chrom1 == chrom2 and pos1 > pos2:
        
        # write the small one = line2
        
        print_line(hap2, chrom2, pos2)

        line2 = f2.readline()
        
    elif int(chrom1) < int(chrom2): 
         
        # write line 1 
        
        print_line(hap1, chrom1, pos1)
        
        line1 = f1.readline()
        
    else: 
        
        # write line 2 
        
        print_line(hap2, chrom2, pos2)
        
        line2 = f2.readline()
        
        
while line1:
    
    hap1, chrom1, pos1, ALT1 = process_line(line1)
    
    #write line1 
    
    print_line(hap1, chrom1, pos1)

    line1 = f1.readline()
    
while line2:
    
    hap2, chrom2, pos2, ALT2 = process_line(line2)

    #write line2

    print_line(hap2, chrom2, pos2)

    line2 = f2.readline()

f1.close()
f2.close()


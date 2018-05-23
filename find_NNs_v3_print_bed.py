#!/usr/bin/python

## this script finds and reports the location of N's in fasta sequence and prints 0-based bed file format 

import sys, os
from Bio import SeqIO
import re

f = open(sys.argv[1], 'rU')  #fasta file

for record in SeqIO.parse(f, "fasta"):

	##record.seq = the sequence ; record.id = name of sequence 
	
	nnindices = [(m.start(0), m.end(0)) for m in re.finditer('N+', str(record.seq))] # nonoverlapping search gives # [(1,5),(9,13)] # nnindices[k][0] is start # nnindices[k][1] is end

	nnlength = list()

	for i, x in enumerate(nnindices):

		nnlength.append(nnindices[i][1] - nnindices[i][0])
		
		chrx = record.id.split()[0]
		
		print chrx + '\t' + str(nnindices[i][0]) + '\t' + str(nnindices[i][1])

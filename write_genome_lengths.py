#!/usr/bin/python

import sys, os
from Bio import SeqIO
import re

## Input = fasta file 

## output = 2-column sequence lengths in fasta file i.e .genome file 

f = open(sys.argv[1], 'rU')

for record in SeqIO.parse(f, "fasta"):

	##record.seq = the sequence ; record.id = name of sequence 
	
	l = len(record.seq)

	chrx = record.id.split()[0]
		
	print chrx + '\t' + str(l)

		 


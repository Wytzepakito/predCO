#!/usr/bin/python

## dinucleotide frequencies from seqeunces, input is multi sequence fasta format; 

import gzip
import os
import sys
import re


def dinuc_list(sequence):

	l1 = re.findall('.{1,2}', sequence)

	l2 = re.findall('.{1,2}', sequence[1:])

	sequence_list = l1 + l2

	return sequence_list

def dinuc_counts(dnt_list, sequence_list): 

	freq = [0] * len(dnt_list)

	for i, pattern in enumerate(dnt_list):

		count = sequence_list.count(pattern)  ## list.count(x)

		freq[i] = float(count) # / float(len(seq)-1)

	return freq 



f = open(sys.argv[1], 'r')


## all possible sets:

nt = ['A', 'G', 'T', 'C']

dnt = [n + nt[i] for n in nt for i in range(4) ]


## write the header

dntw = [d + '\t' for d in dnt]

dntw = ''.join(dntw)

newline = '#name' + '\t' + dntw.strip()

print newline




for line in f: 
	
	line = line.strip();

	if line.startswith('>'):

		seqname = line[1:];
		continue

	else:

		seq = line	


	seqlist = dinuc_list(seq)

	dcounts = dinuc_counts(dnt, seqlist)

	freq = [float(d) / float(len(seq)-1) for d in dcounts ]
		
	# print freq

	freq = ['%.4f' % i + '\t' for i in freq]

	freq = ''.join(freq)

	newline = seqname + '\t' + freq.strip() 

	print newline
	

#!/usr/bin/python3

import numpy as np
import pandas as pd

import sys

from subprocess import Popen, PIPE, STDOUT
import shlex

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from sklearn.neighbors import KernelDensity


## check if they overlap 

def remove_overlaps(samples, overlap):
    
    "This function takes a list and removes the elements if they overlap with defined 'overlap'"
    
    idX = list()
    
    samples = [int(i) for i in samples]
    
    samples = np.sort(samples)

    d = np.diff(samples)
   
    if min(d) < overlap:
        
        for i, num in enumerate(d):
            if num < overlap:
                idX.append(i)   # i is the index of the first of the overlapping pair in sorted list
       
        samples = np.delete(samples, idX, 0)
        
    return samples
    

def check_range(x, maxlen):
    
    xnew = [i for i in x if i > 2001 and i < (maxlen - 2001)]
    
    return xnew


def remove_exclude(samples):
    
    "This function takes list of positions, calls bedtools to remove the regions +-2000"
    "when they overlap with defined exclude_regions"
    
    # convert the samples to bed formatted file 

    sample_lines = ""

    for i in range(len(samples)):
        st = samples[i] -2001 
        end = samples[i] +2000
        bedline= chrname + '\t' + str(st) + '\t' + str(end) + '\n'
        sample_lines = sample_lines + bedline

    sample_lines = sample_lines.strip()

    ## call bedtools to remove exlude regions 

    command_line = "bedtools intersect -nonamecheck -v -a stdin -b " + f_exclude  
    cmd = shlex.split(command_line)

    p = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=STDOUT, universal_newlines=True)
    p_stdout, p_stderr = p.communicate(input=sample_lines) 

    ## convert bedfile back to samples list format 

    lines = p_stdout.strip().split('\n')

    newsamples = [int(line.split('\t')[1]) + 2001 for line in lines]

    return newsamples


def optimize_bd(dfGenome, dfPos, dfGene, outpath):
    
    "Bandwidth optimization by fitting the density to positive set"

    dfPos['mid'] = ((dfPos['end'] - dfPos['start'])/2) + dfPos['start']

    chrs = list(dfGenome.chrom.unique())

    bdlist = list(np.linspace(1000, 1000000, 1000))

    sc = np.array([0.0] * (len(bdlist)+1))

    for chrname in chrs: 

        chrlen = int(dfGenome[dfGenome.chrom == chrname].length)

        N = dfPos[dfPos.chrom == chrname].shape[0]

        dfchr = dfGene[dfGene.chrom == chrname]

        dfPosChr = dfPos[dfPos.chrom == chrname]

        Xp = np.array(list(dfPosChr['mid']))[:, np.newaxis]

        X = np.array(list(dfchr['mid']))[:, np.newaxis]

        ## estimate the density at each 1000 bp

        X_plot = np.linspace(0, chrlen, int(chrlen/1000))[:, np.newaxis]  

        b= np.array([[0,0]])
    
        print("optimization for", chrname)

        for bd in bdlist:

            kde = KernelDensity(kernel='gaussian', bandwidth=bd).fit(X)

            a = np.c_[bd, kde.score(Xp)]

            b = np.r_[b, a]

        sc[:] =  sc[:] + b[:,1]


    end = np.c_[bdlist, list(sc[1:,])]

    idxrow = np.argwhere(end == max(end[:,1]))[0,0]
    newbd = int(end[idxrow, 0])

    print("the bandwith is", newbd)


    #plt.plot(bdlist, list(sc[1:,]))
    #plt.title("genome")
    #plt.xlabel("bandwidth (bp)")
    #plt.ylabel("log score of positive set")
    #plt.savefig(path + 'gene_density_optimization.png')
    #plt.close()

    dfout = pd.DataFrame({'A' : bdlist, 'B' : sc[1:,]})

    dfout.to_csv(path_or_buf= outpath + "bandwidth_trials.txt", sep='\t', header=False, index=False)
    
    return newbd


if __name__=="__main__":


    f_pos = sys.argv[1] # path + "positive.all.4000.maxoverlap1kb.bed"

    f_genome = sys.argv[2] # path + "ref_genome/lengths.genome"  ## changes for each species

    f_exclude = sys.argv[3] # path + "exclude_region.merged.bed"

    f_gene = sys.argv[4] # path + "ref_genome/gene.mid.txt"
    
    outpath = sys.argv[5]


    dfGene = pd.read_csv(f_gene, sep="\t", header=None, names=["chrom", "mid", "length"],
                     dtype = {'chrom': object, 'mid': np.int64, 'length' : np.int64})

    #dfNN = pd.read_csv(f_NN, sep="\t", header=None, names=["chrom", "start", "end"])

    dfGenome = pd.read_csv(f_genome, sep="\t", header=None, names=["chrom", "length"],
                    dtype = {'chrom': object, 'length' : np.int64})

    dfPos = pd.read_csv(f_pos, sep="\t", header=None, names=["chrom", "start", "end"],
                    dtype = {'chrom': object, 'start': np.int64, 'end' : np.int64})



    chrs = list(dfGenome.chrom.unique())

    df = pd.DataFrame()

    bd = optimize_bd(dfGenome, dfPos, dfGene, outpath)



    for chrname in chrs: 

        chrlen = int(dfGenome[dfGenome.chrom == chrname].length)

        N = dfPos[dfPos.chrom == chrname].shape[0]

        dfchr = dfGene[dfGene.chrom == chrname]

        X = np.array(list(dfchr['mid']))[:, np.newaxis]
        
        ## estimate the density at each 1000 bp

        X_plot = np.linspace(0, chrlen, int(chrlen/1000))[:, np.newaxis]  
        
        print("sampling for", chrname)
        
        kde = KernelDensity(kernel='gaussian', bandwidth=bd).fit(X)

        samples = kde.sample(N)
        
        ov = 4000

        samples = check_range(samples, chrlen)
        samples = remove_overlaps(samples, ov)
        samples = remove_exclude(samples)

        while len(samples) < N:
            subs = N - len(samples)
            newsamples = kde.sample(subs)
            samples = np.append(samples, newsamples)
            samples = check_range(samples, chrlen)
            samples = remove_overlaps(samples, ov)
            samples = remove_exclude(samples)

        ## plot the gene density 
        
        log_dens = kde.score_samples(X_plot)

        plt.plot(X_plot, np.exp(log_dens), color="lightblue", lw=1)

        plt.plot(samples, np.zeros(len(samples)), '+k', color="red")

        plt.title(chrname)
        
        plt.savefig(outpath + chrname + '_gene_density_bw_' + str(bd) + '.png')

        plt.close()
        
        ## append the sampled chromosomes 

        dfout = pd.DataFrame({'A' : chrname, 
                  'B' : np.array([i - 2001 for i in samples]),
                  'C' : np.array([i + 2000 for i in samples])})
        
        if chrname == chrs[0]:
        
            df = dfout
        
        else: 
        
            df = df.append(dfout, ignore_index=True)

    df.to_csv(path_or_buf= outpath + "sampled.bed", sep='\t', header=False, index=False)



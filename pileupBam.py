#!/bin/env python

import pysam
import sys
import numpy as np
import re
import argparse
from multiprocessing import Pool

def getOption():
    """
    Get all option fomr command lines
    """
    parser = argparse.ArgumentParser(description='Pile up bam file and get mapped bases (excluded clipped bases and low quality base).')
    parser.add_argument('-i','--bamfile',required=True,
            help='position sorted bam file')
    parser.add_argument('-q','--qualThresh',type=int, default = 33,
            help='Base calling quality threshold (default: 33)')
    parser.add_argument('-r','--refFasta', required=True,
            help='reference fasta')
    parser.add_argument('-d','--depth', type=int, 
            default = 8000, help='Maximum depth (default: 8000)')
    parser.add_argument('-p','--threads', type=int, 
            default = 1, help='Threads to used (default: 1)')
    args = parser.parse_args()
    bamfile = args.bamfile
    qualThresh = args.qualThresh
    refFasta = args.refFasta
    depth =  args.depth
    threads =  args.threads
    return bamfile, qualThresh, refFasta, depth, threads
    

def cigarToSeq(cigar):
    """
    Input a cigar string (eg. 1S16M5I10M4S)
    output a line: SMMMMMMMMMMMMMMMMIIIIIMMMMMMMMMMSSSS
    """
    cigarNum = np.array(re.findall('[0-9]+',cigar),dtype='int64')
    cigarStr = np.array(re.findall('[A-Z]',cigar),dtype='string')
    usable = np.in1d(cigarStr,np.array(['S','M','I'],dtype='string'))
    cigarStr = cigarStr[usable]
    cigarNum = cigarNum[usable]
    cigarSeq = ''
    for s, n in zip(cigarStr, cigarNum):
        cigarSeq += int(n)*str(s)
    return cigarSeq

def extractQual(args):
    """
    for each alignment that mapped to the position,
    extract the base and the cigar annotation
    if the base is a Mapped position ('M') and have quality higher than the 
    given threshold, return the base
    """
    sequence, quality, cigar, matchPos, indel, qualThresh = args
    if indel == 0 and matchPos != 0 and matchPos != (len(sequence)-1):
        cigarSeq = cigarToSeq(cigar)
        assert len(cigarSeq) == len(sequence), '\n%s\n%s' %(cigarSeq,sequence)
        qual = quality[matchPos]
        base = sequence[matchPos]
        if ord(qual) - 33 > qualThresh and cigarSeq[matchPos] == 'M':
            return base
    
def printLine(referenceBase, bases, count , position, coverage):
    """
    for each column position (pileup column), output the 
    1. position on the chromosome
    2. reference base on the ref-fasta
    3. count of A, C, T, G
    """
    regularBase = np.array(['A','T','C','G'],dtype='string')
    baseDict = {}
    coverage = np.sum(count[np.in1d(bases,regularBase)]) #only high qual base were counted in coverage
    for base in regularBase:
        if base not in bases:
            baseDict[base] = 0
        else: 
            baseDict[base] = count[bases==base][0]
    line = str(position) +'\t' + referenceBase  + '\t' + str(coverage)
    for base in regularBase:
        line += '\t' + str(baseDict[base])
    print line
    return 0

def analyzePosition(pileupColumn, refBase, threads, position, qualThresh):
    """
    for each pileup position, extracted all alignments using pysam
    processing each alignment is computationally heavy.
    used multiprocessing in here
    """
    cov = pileupColumn.n
    pool = Pool(processes=threads, maxtasksperchild = 1)
    bases = pool.map(extractQual,[(aln.alignment.seq, aln.alignment.qual, 
                                    aln.alignment.cigarstring, aln.qpos,
                                    aln.indel, qualThresh) for aln in pileupColumn.pileups \
                                            if (not aln.alignment.is_secondary)])
    pool.close()
    pool.join()
    bases = filter(None,bases)
    bases = np.array(bases,dtype='string')
    bases, count = np.unique(bases,return_counts=True)
    printLine(refBase, bases, count, position, cov) 
    return 0

def printHeader():
    print 'refpos\trefBase\tcoverage\tA\tT\tC\tG'
    return 0

def main():
    """
    Using pysam to pileup bam file
    1. extract aligned reads at each position
    2. extract the base on each mapped aligned reads on the position (computationally heavy)
    3. write out base count lines.
    """
    programnam = sys.argv[0]
    bamfile, qualThresh, ref, depth, threads = getOption()
    refFasta = pysam.Fastafile(ref) 
    with pysam.Samfile(bamfile,'rb') as bam:
        if not bam._hasIndex() :
            sys.stderr.write('[%s] Indexing %s \n' %(programnam,bamfile))
            pysam.index(bamfile)
    with pysam.Samfile(bamfile,'rb') as bam:
        sys.stderr.write('[%s] Pileup BAM file: %s \n' %(programnam,bamfile))
        refName = bam.references[0]
        refLength = int(bam.lengths[0])
        printHeader()
        for pileupColumn in bam.pileup(refName,
                                    fastafile=refFasta,
                                    max_depth=depth):
            position = pileupColumn.pos
            refBase = refFasta.fetch(refName,position, position+1)
            analyzePosition(pileupColumn, refBase, threads, position, qualThresh)
            if position % 10 == 0:
                sys.stderr.write('[%s] Piled up %s: %i \n' %(programnam,refName,position))
    sys.stderr.write('[%s] Finished extracting %s.\n' %(programnam,bamfile))
    return 0
        
if __name__ == '__main__':
    main()

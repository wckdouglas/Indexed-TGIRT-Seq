#!/usr/bin/env python

import pysam
import re
import numpy as np
import sys
import argparse
from numba import jit

def mismatchMatrix():
    """
        read
    ref    A C T G
    A      0 0 0 0
    C      0 0 0 0
    T      0 0 0 0
    G      0 0 0 0
    """
    return np.zeros(16).reshape(4,4), {'A':0,'C':1,'T':2,'G':3}

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

def MDToSeq(mdtag):
    """
    input mdTag: 31^CA31
    output: -------------------------------DDD-------------------------------

    input mdTag: 8G11^TG14T0T28T3
    output: --------G-----------DD--------------TT----------------------------T---
    """
    mdNum = np.array(re.findall('[0-9]+',mdtag),dtype='int64')
    mdStr = np.array(re.split('[0-9]+',mdtag),dtype='string')[1:]
    mdSeq = ''
    for s, n in zip(mdStr,mdNum):
        string = n * '-' + (len(s)-1) * 'D' if '^' in s else n * '-' + s
        mdSeq += string
    return mdSeq

def processAln(aln, misMat, tTable, qualThresh, lenCut):
    startpos = aln.qstart
    sequence = aln.seq
    MDtag = [tag[1] for tag in aln.tags if tag[0] == 'MD'][0]
    cigarSeq = cigarToSeq(aln.cigarstring) #(MDtags have no soft clipped)
    qual = np.array([ord(q)-33 for q in aln.qual],dtype='int64')
    newSeq, newCigar, newQual = '', '', []
    for s, c, q in zip(sequence, cigarSeq, qual):
        if c == 'M':
            newSeq += s
            newCigar += c
            newQual.append(q)
    mdSeq = re.sub('D','',MDToSeq(MDtag)) #(cigarSeq do not contain 'D' anymore)
    assert len(newSeq) == len(mdSeq), newSeq + '\n' + mdSeq + '\n' + newCigar
    assert len(newSeq) == len(newCigar), newSeq + '\n' + mdSeq + '\n' + newCigar + '\n' 
    assert len(newSeq) == len(newQual), newSeq + '\n' + mdSeq + '\n' + newCigar + '\n' + ''.join(chr(q+33) for q in newQual)
    length = len(newSeq)
    pos = np.arange(len(newSeq))
    for b, m, c, q, p in zip(newSeq, mdSeq, newCigar, newQual, pos):
        if q > qualThresh and c == 'M' and p > lenCut and p < length - lenCut:
            assert m != b, '\nread: '+ m +\
                           '\n ref: ' + b +\
                           '\n'+ newSeq +\
                           '\n' + mdSeq +\
                           '\n' + MDtag +\
                           '\n' + sequence +\
                           '\n' + cigarSeq
            if m == '-':
                misMat[tTable[b],tTable[b]] += 1
            else:
                misMat[tTable[b],tTable[m]] += 1
    return misMat

def main(inBam, qualThresh, lenCut, programname):
    misMat, tTable = mismatchMatrix()
    with pysam.Samfile(inBam,'rb') as bam:
        count = 0
        for aln in bam:
            if not any(np.in1d(['I','D'],list(aln.cigarstring))):
                processAln(aln, misMat, tTable, qualThresh, lenCut)
            count += 1
            if count % 100000 == 0:
                sys.stderr.write('[%s] Processed %i alignments\n' %(programname,count))
    print 'ref\tA\tC\tT\tG'
    rownum = 0
    bases = 'ACTG'
    for row in misMat:
        print bases[rownum] + '\t' + '\t'.join(map(str,row))
        rownum += 1
    sys.stderr.write('[%s] Finished parsing %s\n' %(programname, inBam))
    return 0

if __name__=='__main__':
    programname = sys.argv[0]
    parser = argparse.ArgumentParser(description='Extract mismatches from alignments in a bam file')
    parser.add_argument('-i','--inBam',required=True, help='input a bam file ( -  if stdout)')
    parser.add_argument('-q','--qual',default=33, type = int, help='BAQ threshold (default: 33)')
    parser.add_argument('-s','--skip',default=3, type = int,  help='Skipping n bases from head and tail (default: 3)')
    args = parser.parse_args()
    inBam =  args.inBam
    qualThresh = args.qual
    lenCut = args.skip
    main(inBam, qualThresh, lenCut, programname)

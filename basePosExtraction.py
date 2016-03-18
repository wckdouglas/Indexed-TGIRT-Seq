#!/usr/bin/env python

import pysam
import re
import numpy as np
import sys
import argparse
from pyfaidx import Fasta

def mismatchMatrix(seqLength):
    """
        read
    pos    A C T G
    0      0 0 0 0
    1      0 0 0 0
    2      0 0 0 0
    3      0 0 0 0
    .      0 0 0 0
    .      0 0 0 0
    .      0 0 0 0
    seqlength 0 0 0 0
    """
    assert isinstance(seqLength,int), seqLength 
    return np.zeros(seqLength * 4).reshape(seqLength,4), {'A':0,'C':1,'T':2,'G':3}

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
    clipped = sum(n for n, s in zip(cigarNum, cigarStr) if s == 'S')
    headClipped = 0 if cigarStr[0] != 'S' else cigarNum[0]
    tailClipped = 0 if cigarStr[-1] != 'S' else cigarNum[-1]
    return cigarSeq, headClipped, tailClipped

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

def processAln(aln, misMat, tTable, qualThresh, lenCut, refSeq, headCut1, tailCut1, headCut2, tailCut2):
    passed = 1
    startpos = aln.pos
    sequence = aln.seq
    MDtag = [tag[1] for tag in aln.tags if tag[0] == 'MD'][0]
    cigarSeq, headClipped, tailClipped = cigarToSeq(aln.cigarstring) #(MDtags have no soft clipped)
    qual = np.array([ord(q)-33 for q in aln.qual],dtype='int64')
    newSeq, newCigar, newQual = '', '', []
    for s, c, q in zip(sequence, cigarSeq, qual):
        if c == 'M':
            newSeq += s
            newCigar += c
            newQual.append(q)
    mdSeq = re.sub('D','',MDToSeq(MDtag)) #(cigarSeq no longer contain 'D' anymore)
    seqlength = len(newSeq)
    assert seqlength == len(mdSeq) and seqlength == len(newCigar) and seqlength == len(newQual), \
            newSeq + '\n' + mdSeq + '\n' + newCigar + '\n' + ''.join(chr(q+33) for q in newQual)
    length = len(newSeq)
    pos = np.arange(len(newSeq))
    for b, m, c, q, p in zip(newSeq, mdSeq, newCigar, newQual, pos):
        if q > qualThresh and c == 'M' and p > lenCut and p < (length - lenCut) and ((tailClipped < tailCut2 and aln.is_read2 and headClipped == headCut2 ) or (headClipped == headCut1 and tailClipped < tailCut1 and aln.is_read1)):
            position = startpos + p
            assert m != b, '\nread: '+ m +\
                           '\n ref: ' + b +\
                           '\n'+ newSeq +\
                           '\n' + mdSeq +\
                           '\n' + MDtag +\
                           '\n' + sequence +\
                           '\n' + cigarSeq
#            refBase = str(refSeq[position])
#            assert (m == '-' and b == refBase) or (m != '-' and m == refBase), 'Wrong interpretation of refbase and md'
            misMat[position,tTable[b]] += 1
            passed = 0
    return passed

def parseRegion(refID, bam, misMat, tTable, qualThresh, lenCut, refSeq, indel, headCut1, tailCut1, headCut2, tailCut2):
    count = 0
#    outBam = ''
#    out = pysam.Samfile(outBam, 'wb',template=bam)
    for aln in bam.fetch(refID) :
        if (not indel and not any(np.in1d(['I','D'],list(aln.cigarstring)))) or (indel):
            passed = processAln(aln, misMat, tTable, qualThresh, lenCut, refSeq, headCut1, tailCut1, headCut2, tailCut2)
 #           if passed == 0:
 #               out.write(aln)
        count += 1
        if count % 100000 == 0:
            sys.stderr.write('[%s] Processed %i alignments in %s\n' %(programname,count, refID))
    rownum = 0
    out.close()
    for row in misMat:
        print '\t'.join([refID, str(rownum), str(refSeq[rownum]), str(sum(row)) ,'\t'.join(map(str,row))])
        rownum += 1
    return 0


def main(inBam, refFa, qualThresh, lenCut, indel, programname, headCut1, tailCut1, headCut2, tailCut2):
    sys.stderr.write('[%s] Starting %s\n' %(programname, inBam))
    print 'refID\trefPos\trefBase\tcoverage\tA\tC\tT\tG'
    indexing =  pysam.index(inBam)
    with pysam.Samfile(inBam,'rb') as bam, Fasta(refFa) as fa:
        refnames = bam.references
        refsizes = bam.lengths
        for refID, refsize in zip(refnames, refsizes):
            misMat, tTable = mismatchMatrix(refsize)
            refSeq = fa[refID]
            parseRegion(refID, bam, misMat, tTable, qualThresh, lenCut, refSeq, indel, headCut1, tailCut1, headCut2, tailCut2)
            sys.stderr.write('[%s] Parsed %s in %s\n' %(programname, refID, inBam))
    sys.stderr.write('[%s] Finished parsing %s\n' %(programname, inBam))
    return 0

if __name__=='__main__':
    programname = sys.argv[0]
    parser = argparse.ArgumentParser(description='Extract mismatches from alignments in a bam file')
    parser.add_argument('-i','--inBam',required=True, help='input a bam file ( -  if stdout)')
    parser.add_argument('-f','--refFasta', required=True,  help='Reference fasta')
    parser.add_argument('-q','--qual',default=33, type = int, help='BAQ threshold (default: 33)')
    parser.add_argument('-s','--skip',default=3, type = int,  help='Skipping n bases from head and tail (default: 3)')
    parser.add_argument('-d','--indel', action='store_true', help='Allow indel alignments')
    parser.add_argument('-a','--h1', default = 0, type = int, help="Need extact n softclipped base on read1 5' (default: 0)")
    parser.add_argument('-b','--h2', default = 0, type = int, help="Need extact n softclipped base on read2 5' (default: 0)")
    parser.add_argument('-y','--t1', default = 0, type = int, help="Allow n softclipped base on read1 3' (defualt: 0)")
    parser.add_argument('-z','--t2', default = 0, type = int, help="Allow n softclipped base on read2 3' (default: 0)")
    args = parser.parse_args()
    args = parser.parse_args()
    inBam =  args.inBam
    qualThresh = args.qual
    lenCut = args.skip
    refFa = args.refFasta
    indel = args.indel
    tailCut1 = args.t1
    tailCut2 = args.t2
    headCut1 = args.h1
    headCut2 = args.h2
    main(inBam, refFa, qualThresh, lenCut, indel, programname, headCut1, tailCut1, headCut2, tailCut2)

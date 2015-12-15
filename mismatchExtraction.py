#!/bin/env python

from __future__ import division
import subprocess
import argparse
from sys import stderr
import sys
import time
import os
import pileupBamToBase
import readClusterPairs

index= '/corral-repl/utexas/2013lambowitz/Ref/syntheticOligos/celSpikeinLong/spikein.fa'
programname = sys.argv[0]

def getOpt():
    parser =  argparse.ArgumentParser(description='pipeline for mismatch analysis.')
    parser.add_argument('-1', '--fq1', required=True, help='read 1 fastq')
    parser.add_argument('-2', '--fq2', required=True, help='read 2 fastq')
    parser.add_argument('-o', '--outputPath', required=True, help='Output directoty')
    parser.add_argument('-t', '--threads', default=1,  help='Threads',type=int)
    args = parser.parse_args()
    fq1 = args.fq1
    fq2 = args.fq2
    outputPath = args.outputPath
    threads = args.threads
    return fq1, fq2, outputPath, threads

def runCommand(command, samplename):
    start = time.time()
    stderr.write('[%s] %s: %s\n' %(programname,samplename, command))
    subprocess.call([command],shell=True)
    timelapse = time.time()-start
    stderr.write('[%s] %s: Time lapse %.3f min\n' %(programname,samplename,timelapse/60))
    return 0


def makeErrorFree(fq1, fq2, outputPath, samplename, threads):
    outputprefix = '%s/rawData/%s-ErrorFree' %(outputPath, samplename)
    inFastq1 = fq1
    inFastq2 = fq2
    idxBase = 13
    minReadCount = 4
    retainN = False
    barcodeCutOff = 30
    voteCutOff = 0.9
    printScore = False
    readClusterPairsmain(outputprefix, inFastq1, inFastq2, idxBase, threads, minReadCount, retainN, barcodeCutOff, voteCutOff, printScore)
    return 0

def mapping(fq1, fq2, outputPath, samplename, threads):
    bampath = outputPath + '/bamFile'
    bamFile = '%s/%s.bam ' %(bampath, samplename)
    command = 'bwa mem -t %i ' %(threads)   +\
           '%s %s %s ' %(index ,fq1, fq2)+\
           '| samtools view -b@ %i -F 4 -F 256 -F 2048 -' %(threads) +\
           '| samtools sort -@ %i -O bam -T %s/%s - ' %(threads, bampath, samplename) + \
           '> %s' %(bamFile)
    runCommand(command, samplename)
    return bamFile

def pileup(bamFile, outputPath, samplename, threads):
    depth = 30000000
    bampath = outputPath + '/bamFile'
    mismatchFile = '%s/mismatchData/%s.tsv' %(outputPath,samplename)
    print samplename
    qualThresh = 33
    ref = index
    skipBases = 3
    outFile = open(mismatchFile, 'w')
    pileupBamToBase.main(bamfile, qualThresh, ref, depth, threads, skipBases, outFile)
    return mismatchFile

def makedirs(dir):
    if not os.path.isdir(dir):
        os.mkdir(dir)
    else:
        stderr.write('[%s] %s existed! \n' %(programname, dir))
    return 0 

def main():
    fq1, fq2, outputPath, threads = getOpt()
    folders = ['bamFile','mismatchData','rawData']
    [makedirs(outputPath + '/' + folder) for folder in folders]
    samplename = os.path.basename(fq1).split('_')[0]
    makeErrorFree(fq1,fq2,outputPath, samplename, threads)
    newfq1 = '%s/rawData/%s-ErrorFree_R1_001.fastq.gz' %(outputPath,samplename)
    newfq2 = '%s/rawData/%s-ErrorFree_R2_001.fastq.gz' %(outputPath,samplename)
    newsamplename = newfq1.split('/')[-1].split('_')[0]
    bamFiles = [mapping(f1, f2, outputPath, name, threads) for f1, f2, name in zip([fq1,newfq1], [fq2,newfq2], [samplename,newsamplename])]
    mismatchFiles = [pileup(bamFile, outputPath, name, threads) for bamFile, name in zip(bamFiles, [samplename, newsamplename])]
    #bamFile = mapping(newfq1, newfq2, outputPath, newsamplename) 
    #mismatchFiles = pileup(bamFile, outputPath, newsamplename) 
    return 0

if __name__ == '__main__':
    main()

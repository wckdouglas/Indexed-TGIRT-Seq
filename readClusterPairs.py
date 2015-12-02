#!/bin/env python

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from sys import stderr
import sys
from shove import Shove
from multiprocessing import Pool, Manager, Process
import argparse
import numpy as np
import glob
import gzip
import time
import math

#    ==================      Sequence class sotring left right record =============
class seqRecord:
    def __init__(self):
        self.seqListRight = np.array([],dtype='string')
        self.qualListRight = np.array([],dtype='string')
        self.seqListLeft = np.array([],dtype='string')
        self.qualListLeft = np.array([],dtype='string')

    def addRecord(self, seqRight, qualRight, seqLeft, qualLeft):
        self.seqListRight = np.append(self.seqListRight,seqRight)
        self.qualListRight = np.append(self.qualListRight,qualRight)
        self.seqListLeft = np.append(self.seqListLeft,seqLeft)
        self.qualListLeft = np.append(self.qualListLeft,qualLeft)

    def readCounts(self):
        assert len(self.seqListLeft) == len(self.seqListRight), 'Not equal pairs'
        return len(self.seqListRight)

    def readLengthRight(self):
        return np.array([len(seq) for seq in self.seqListRight],dtype='int64')

    def readLengthLeft(self):
        return np.array([len(seq) for seq in self.seqListLeft],dtype='int64')

#======================  starting functions =============================
def getOptions():
    '''reading input 
    '''
    descriptions = 'Clustering fastq reads to fasta reads with the first $IDXBASE bases as cDNA-synthesis barcode. ' +\
                'Concensus bases are sleected by the highest vote, ' +\
                'quality scores are generated by the average score of the concensus bases at the given position' +\
                'ties are settled by base calling quality if vote threshold is not set'
    parser = argparse.ArgumentParser(description=descriptions)
    parser.add_argument('-o', '--outputprefix', required=True,
        help='Paired end Fastq files with R1_001.fastq.gz as suffix for read1, and R2_001.fastq.gz as suffix for read2')
    parser.add_argument('-1', '--fastq1', required=True,
        help='Paired end Fastq file 1 with four line/record')
    parser.add_argument('-2', '--fastq2',required=True,
        help='Paired end Fastq file 2 with four line/record')
    parser.add_argument('-m', '--cutoff', type=int,default=4,
        help="minimum read count for each read cluster (default: 4)")
    parser.add_argument("-t", "--threads", type=int, default=1,
        help="number of threads to use (default: 1)")
    parser.add_argument("-x", "--idxBase", type=int, default=13,
        help="how many base in 5' end as index? (default: 13)")
    parser.add_argument("-n", "--retainN", action='store_true',
        help="Use N-containing sequence for concensus base vote and output sequences containing N (defulat: False)")
    parser.add_argument('-q', '--barcodeCutOff', type=int, default=30,
        help="Average base calling quality for barcode sequence (default=30)")
    parser.add_argument('-l', '--loglikThreshold', type=float, default=100,
            help="log likelihood threshold, if the p(max likelihood base)/p(2nd max likelihood base) do not exceed threshold, the position would have 'N' as base. (default: 100)")
    parser.add_argument('-s', '--seqError', type=float, default=0.01,
            help="Sequencing error probability (default: 0.01)")
    args = parser.parse_args()
    outputprefix = args.outputprefix
    inFastq1 = args.fastq1
    inFastq2 = args.fastq2
    idxBase = args.idxBase
    threads = args.threads
    minReadCount = args.cutoff
    retainN = args.retainN
    barcodeCutOff = args.barcodeCutOff
    seqErr = args.seqError
    loglikThreshold = args.loglikThreshold
    return outputprefix, inFastq1, inFastq2, idxBase, threads, \
            minReadCount, retainN, barcodeCutOff, seqErr, loglikThreshold


def getLikelihood(countDict, seqErr, base):
    """
    Using sequencing error rate as conditional probability: e.g. P(A|base=C) = 0.01
    calculate the likelihood of the given base is the concensus base
    """
    probabilities = np.array([np.power((1-seqErr),countDict[basekey]) \
                                if basekey == base \
                                else np.power(seqErr,countDict[basekey]) \
                                for basekey in countDict.keys()],dtype='float64')
    likelihood = np.prod(probabilities)
    return likelihood

def likelihoodSelection(uniqueBases, baseCount, seqErr, loglikThreshold):
    """
    for each base in {A,T,C,G}
    1. get the count of each base at the position
    2. get the likelihood function assuming the certain base is the concensus base
    3. perform log likelihood test and compare to threshold
    4. output concensus base or 'N'
    """
    regBase = np.array(['A','T','C','G'],dtype='string')
    countDict = {} 
    for base in regBase:
        countDict[base] = baseCount[uniqueBases==base][0] \
                            if base in uniqueBases \
                            else 0 
    likelihood = np.array([getLikelihood(countDict, seqErr, base) \
                            for base in regBase],dtype='float64')
    sortedLikihood = np.sort(likelihood)[::-1]
    maxLikelihood = sortedLikihood[0]
    otherLikelihood = np.sum(sortedLikihood[1:])
    loglikRatio = np.log(np.divide(maxLikelihood,otherLikelihood))\
                    if otherLikelihood != 0 \
                    else loglikThreshold + 1
    concensusBase = regBase[likelihood == maxLikelihood][0] \
                    if loglikRatio > loglikThreshold \
                    else 'N'
    return concensusBase

def calculateConcensusBase(arg):
    """Given a list of sequences, 
        a list of quality line and 
        a position, 
    return the maximum likelihood base at the given position,
        along with the mean quality of these concensus bases.
    """
    seqList, qualList, pos, seqErr, loglikThreshold = arg
    columnBases = np.array([],dtype='string')
    qualities = np.array([],dtype='int64')
    for seq, qual in zip(seqList, qualList):
        columnBases = np.append(columnBases,seq[pos])
        qualities = np.append(qualities,ord(qual[pos]))
    uniqueBases, baseCount = np.unique(columnBases, return_counts=True)
    concensusBase = likelihoodSelection(uniqueBases, baseCount, seqErr, loglikThreshold)
    # offset -33
    quality = np.mean(qualities[columnBases==concensusBase]) \
            if concensusBase in columnBases \
            else 33
    return concensusBase, quality

def concensusSeq(seqList, qualList, positions, seqErr, loglikThreshold):
    """given a list of sequences, a list of quality and sequence length. 
        assertion: all seq in seqlist should have same length (see function: selectSeqLength)
    return a consensus sequence and the mean quality line (see function: calculateConcensusBase)
    """
    concensusPosition = map(calculateConcensusBase,[(seqList, qualList, pos, seqErr, loglikThreshold) \
                            for pos in positions])
    bases, quals = zip(*concensusPosition)
    sequence = ''.join(list(bases))
    quality = ''.join([chr(int(q)) for q in list(quals)])
    return sequence, quality


def concensusPairs(reads, seqErr, loglikThreshold):
    """ given a pair of reads as defined as the class: seqRecord
    return concensus sequence and mean quality of the pairs, 
        as well as the number of reads that supports the concnesus pairs
    see function: concensusSeq, calculateConcensusBase
    """
    # get concensus left reads first
    sequenceLeft, qualityLeft = concensusSeq(reads.seqListLeft, reads.qualListLeft,  
                                            range(np.unique(reads.readLengthLeft())[0]),  
                                            seqErr, loglikThreshold)
    assert len(sequenceLeft) == len(qualityLeft), 'Wrong concensus sequence and quality!'
    # get concensus right reads first
    sequenceRight, qualityRight = concensusSeq(reads.seqListRight, reads.qualListRight,  
                                                range(np.unique(reads.readLengthRight())[0]), 
                                                seqErr, loglikThreshold)
    assert len(sequenceRight) == len(qualityRight), 'Wrong concensus sequence and quality!'
    return sequenceLeft, qualityLeft, len(reads.seqListLeft), \
            sequenceRight, qualityRight, len(reads.seqListRight)

def selectSeqLength(readLengthArray):
    """
    Given a list of sequence length of a read cluster from either side of the pair,
    select the sequence length with highest vote
    """
    seqlength, count = np.unique(readLengthArray, return_counts=True)
    return seqlength[count==max(count)][0]

def getCompatibleSeqs(seqlist, qualList, readLengthList):
    """ filter the read clusters to retrain all the sequence 
    in the list of sequence with the highest vote of read length
    see function: selectSeqLength
    """
    selectedLength =  selectSeqLength(readLengthList)
    seqs = seqlist[readLengthList == selectedLength]
    quals = qualList[readLengthList == selectedLength]
    return seqs, quals

def filterRead(readCluster):
    """
    Given a read cluster (seqRecord object),
    filter out all sequences that have a distinctive length (for getting concensus sequence from same length reads)
    and create a seqRecord object with the filtered sequences only
    """
    leftSeqs, leftQuals = getCompatibleSeqs(readCluster.seqListLeft, 
                                            readCluster.qualListLeft, 
                                            readCluster.readLengthLeft())
    rightSeqs, rightQuals = getCompatibleSeqs(readCluster.seqListRight, 
                                            readCluster.qualListRight, 
                                            readCluster.readLengthRight())
    reads = seqRecord()
    [reads.addRecord(rs, rq, ls, lq) for rs, rq, ls, lq in zip(rightSeqs, rightQuals, leftSeqs, leftQuals)]
    return reads

def errorFreeReads(args):
    """
    main function for getting concensus sequences from read clusters.
    return  a pair of concensus reads with a 4-line fastq format
    see functions: 1. filterRead, 
                  2. concensusPairs,
                  3. calculateConcensusBase
    """
    readCluster, index, counter, lock, minReadCount, \
            retainN, seqErr, loglikThreshold = args
    # skip if not enough sequences to perform voting
    reads = filterRead(readCluster)
    if reads.readCounts() > minReadCount:
        sequenceLeft, qualityLeft, supportedLeftReads, \
        sequenceRight, qualityRight, supportedRightReads = concensusPairs(reads, 
                                                                        seqErr, 
                                                                        loglikThreshold)
        if (retainN == False and 'N' not in sequenceRight and 'N' not in sequenceLeft) \
                or (retainN == True and set(sequenceLeft)!={'N'}):
            lock.acquire()
            counter.value += 1
            clusterCount = counter.value
            leftFile = '@cluster_%i %s %i reads\n%s\n+\n%s\n' \
                %(counter.value, index, supportedLeftReads, sequenceLeft, qualityLeft)
            rightFile = '@cluster_%i %s %i reads\n%s\n+\n%s\n' \
                %(counter.value, index, supportedRightReads, sequenceRight, qualityRight)
            if clusterCount % 100000 == 0:
                stderr.write('[readCluster] Processed %i read clusters.\n' %clusterCount)
            lock.release()
            return leftFile,rightFile

def readClustering(args):
    """
    generate read cluster with a dictionary object and seqRecord class.
    index of the dictionary is the barcode extracted from first /idxBases/ of read 1 
    """
    read1, read2, barcodeDict, idxBase, barcodeCutOff, retainedN = args
    idLeft, seqLeft, qualLeft = read1
    idRight, seqRight, qualRight = read2
    barcode = seqLeft[:idxBase]
    barcodeQualmean = int(np.mean([ord(q) for q in qualLeft[:idxBase]]) - 33)
    if ((retainedN==True and 'N' not in seqLeft and 'N' not in seqRight) or retainedN==False) \
            and ('N' not in barcode and barcodeQualmean > barcodeCutOff):
        seqLeft = seqLeft[idxBase:]
        barcodeDict.setdefault(barcode,seqRecord()) 
        barcodeDict[barcode].addRecord(seqRight, qualRight, seqLeft, qualLeft)
    return 0

def writeFile(stringList,outputFile):
    """
    write fastq lines to gzip files
    """
    with gzip.open(outputFile,'wb') as f:
        [f.write(string) for string in stringList]
    return 0

def main():
    """
    main function:
        controlling work flow
        1. generate read clusters by reading from fq1 and fq2
        2. obtain concensus sequence from read clusters
        3. writing concensus sequence to files
    """
    outputprefix, inFastq1, inFastq2, idxBase, threads, minReadCount, \
            retainN, barcodeCutOff, seqErr, loglikThreshold = getOptions()
    start = time.time()

    #print out parameters
    stderr.write( '[readCluster] Using parameters: \n')
    stderr.write( '[readCluster]     threads:                           %i\n' %threads)
    stderr.write( '[readCluster]     indexed bases:                     %i\n' %idxBase)
    stderr.write( '[readCluster]     minimum coverage:                  %i\n' %minReadCount)
    stderr.write( '[readCluster]     outputPrefix:                      %s\n' %outputprefix)
    stderr.write( '[readCluster]     retaining N-containing sequence:   %r\n' %retainN)
    stderr.write( '[readCluster]     Sequencing error rate:             %.3f\n' %seqErr)
    stderr.write( '[readCluster]     log likelihood threhold:           %.3f\n' %loglikThreshold)
    
    #barcodeDict = Shove('file://%s' %outputprefix) # reduced memory by: shove package
    barcodeDict = {}
    # generate read clusters
    with gzip.open(inFastq1,'rb') as fq1, gzip.open(inFastq2,'rb') as fq2:
        map(readClustering,[(read1,read2,barcodeDict, idxBase, barcodeCutOff, retainN)  \
            for read1,read2 in zip(FastqGeneralIterator(fq1),FastqGeneralIterator(fq2))])
    stderr.write('[readCluster] Extracted: %i barcodes sequence\n' %len(barcodeDict.keys()))

    # From index library, generate error free reads
    # using multicore to process read clusters
    counter = Manager().Value('i',0)
    lock = Manager().Lock()
    pool = Pool(processes=threads, maxtasksperchild = 1)
    results = pool.map(errorFreeReads, [(barcodeDict[index],index, counter, lock, \
                        minReadCount, retainN, seqErr, loglikThreshold) \
                        for index in barcodeDict.keys()])
    pool.close()
    pool.join()
    # since some cluster that do not have sufficient reads
    # will return None, results need to be filtered
    results = filter(None,results)
    if (len(results) == 0):
        sys.exit('No concensus clusters!! \n')
    left, right = zip(*results)
    stderr.write('[readCluster] Extracted error free reads\n')
    # output file name
    read1File = outputprefix + '_R1_001.fastq.gz'
    read2File = outputprefix + '_R2_001.fastq.gz'
    # use two cores for parallel writing file
    processes = [Process(target = writeFile, args = (seqRecords,file)) \
                for seqRecords, file in zip([list(left),list(right)], [read1File,read2File])]
    [p.start() for p in processes]
    [p.join() for p in processes]
    # all done!
    stderr.write('[readCluster] Finished writing error free reads\n')
    stderr.write('[readCluster]     read1:            %s\n' %(read1File))
    stderr.write('[readCluster]     read2:            %s\n' %(read2File))
    stderr.write('[readCluster]     output clusters:  %i\n' %(len(left)))
    stderr.write('[readCluster]     time lapsed:      %2.3f min\n' %(np.true_divide(time.time()-start,60)))
    return 0
        
if __name__ == '__main__':
    main()

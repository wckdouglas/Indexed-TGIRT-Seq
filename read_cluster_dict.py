#!/bin/env python

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from sys import stderr
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import argparse
import glob
import gzip
import time
import os
from itertools import izip
from multiprocessing import Pool, Manager
from scipy.spatial.distance import hamming
from collections import defaultdict
sns.set_style('white')
programname = os.path.basename(sys.argv[0]).split('.')[0]
minQ = 33
maxQ = 73

#    ==================      Sequence class sotring left right record =============
class seqRecord:
    def __init__(self):
        self.seqListRight = []
        self.qualListRight = []
        self.seqListLeft = []
        self.qualListLeft = []

    def addRecord(self, seqRight, qualRight, seqLeft, qualLeft):
        self.seqListRight.append(seqRight)
        self.qualListRight.append(qualRight)
        self.seqListLeft.append(seqLeft)
        self.qualListLeft.append(qualLeft)

    def readCounts(self):
        assert len(self.seqListLeft) == len(self.seqListRight), 'Not equal pairs'
        return len(self.seqListRight)

    def readLengthRight(self):
        return np.array([len(seq) for seq in self.seqListRight],dtype=np.int64)

    def readLengthLeft(self):
        return np.array([len(seq) for seq in self.seqListLeft],dtype=np.int64)

#======================  starting functions =============================
def getOptions():
    '''reading input
    '''
    descriptions = 'Clustering fastq reads to fasta reads with the first $IDXBASE bases as cDNA-synthesis barcode. ' +\
                'Concensus bases are called only when the fraction of reads that contain the concensus base exceed some threshold. '+ \
                'Quality scores are generated by the average score for the bases that matched concensus base. '
    parser = argparse.ArgumentParser(description=descriptions)
    parser.add_argument('-o', '--outputprefix', required=True,
        help='Paired end Fastq files with R1_001.fastq.gz as suffix for read1, and R2_001.fastq.gz as suffix for read2')
    parser.add_argument('-1', '--fastq1', required=True,
        help='Paired end Fastq file 1 with four line/record')
    parser.add_argument('-2', '--fastq2',required=True,
        help='Paired end Fastq file 2 with four line/record')
    parser.add_argument('-m', '--cutoff', type=int,default=4,
        help="minimum read count for each read cluster (default: 4)")
    parser.add_argument("-x", "--idxBase", type=int, default=13,
        help="how many base in 5' end as index? (default: 13)")
    parser.add_argument('-q', '--barcodeCutOff', type=int, default=30,
        help="Average base calling quality for barcode sequence (default=30)")
    parser.add_argument("-t", "--threads", type=int, default = 1,
        help="Threads to use (default: 1)")
    parser.add_argument("-c", "--constant_region", default='',
            help="Constant sequence after tags (default: '')")
    args = parser.parse_args()
    outputprefix = args.outputprefix
    inFastq1 = args.fastq1
    inFastq2 = args.fastq2
    idxBase = args.idxBase
    minReadCount = args.cutoff
    barcodeCutOff = args.barcodeCutOff
    threads = args.threads
    constant = args.constant_region
    return outputprefix, inFastq1, inFastq2, idxBase, minReadCount, barcodeCutOff, threads, constant

def qual2Prob(q):
    '''
    Given a q list,
    return a list of prob
    '''
    return np.power(10, np.true_divide(-q,10))

def calculatePosterior(guessBase, columnBases, qualities):
    qualHit = qualities[columnBases==guessBase]
    qualMissed = qualities[columnBases!=guessBase]
    if len(qualMissed) > 0:
        hit = np.prod(1- qual2Prob(qualHit)) if len(qualHit) > 0 else 0
        missed = np.prod(np.true_divide(qual2Prob(qualMissed),3))
        posterior = missed * hit
    else:
        posterior = 1
    return posterior

def calculateConcensusBase(arg):
    """Given a list of sequences,
        a list of quality line and
        a position,
    return the maximum likelihood base at the given position,
        along with the mean quality of these concensus bases.
    """
    seqList, qualList, pos = arg
    no_of_reads = len(seqList)
    acceptable_bases = np.array(['A','C','T','G'], dtype='string')
    columnBases = np.empty(no_of_reads,dtype='string')
    qualities = np.empty(no_of_reads,dtype=np.int64)
    for seq, qual, i  in zip(seqList, qualList, range(no_of_reads)):
        columnBases[i] = seq[pos]
        qualities[i] = ord(qual[pos]) -33
    posteriors = [calculatePosterior(guessBase, columnBases, qualities) for guessBase in acceptable_bases]
    posteriors = np.true_divide(posteriors, np.sum(posteriors))
    maxLikHood = np.argmax(posteriors)
    concensusBase = acceptable_bases[maxLikHood]
    posterior = posteriors[maxLikHood]
    quality = -10 * np.log10(1 - posterior) if posterior < 1 else maxQ
    return concensusBase, quality

def concensusSeq(seqList, qualList, positions):
    """given a list of sequences, a list of quality and sequence length.
        assertion: all seq in seqlist should have same length (see function: selectSeqLength)
    return a consensus sequence and the mean quality line (see function: calculateConcensusBase)
    """
    if len(seqList) > 1:
        concensusPosition = map(calculateConcensusBase,[(seqList, qualList, pos) for pos in positions])
        bases, quals = zip(*concensusPosition)
        quality = np.array(quals,dtype=np.int64)
        quality[quality<minQ] = minQ
        quality[quality > maxQ] = maxQ
        sequence = ''.join(list(bases))
        quality = ''.join(map(chr,quality))
    else:
        sequence = seqList[0]
        quality = qualList[0]
    return sequence, quality

def concensusPairs(reads):
    """ given a pair of reads as defined as the class: seqRecord
    return concensus sequence and mean quality of the pairs,
        as well as the number of reads that supports the concnesus pairs
    see function: concensusSeq, calculateConcensusBase
    """
    # get concensus left reads first
    sequenceLeft, qualityLeft = concensusSeq(reads.seqListLeft, reads.qualListLeft, range(np.unique(reads.readLengthLeft())[0]))
    assert len(sequenceLeft) == len(qualityLeft), 'Wrong concensus sequence and quality!'
    # get concensus right reads first
    sequenceRight, qualityRight = concensusSeq(reads.seqListRight, reads.qualListRight, range(np.unique(reads.readLengthRight())[0]))
    assert len(sequenceRight) == len(qualityRight), 'Wrong concensus sequence and quality!'
    return sequenceLeft, qualityLeft, len(reads.seqListLeft), sequenceRight, qualityRight, len(reads.seqListRight)

def selectSeqLength(readLengthArray):
    """
    Given a list of sequence length of a read cluster from either side of the pair,
    select the sequence length with highest vote
    """
    seqlength, count = np.unique(readLengthArray, return_counts=True)
    return seqlength[count==max(count)][0]

def errorFreeReads(args):
    """
    main function for getting concensus sequences from read clusters.
    return  a pair of concensus reads with a 4-line fastq format
    see functions: 1. filterRead,
                  2. concensusPairs,
                  3. calculateConcensusBase
    """
    #if readCluster.readCounts() > minReadCount:
    #    reads = filterRead(readCluster)
    # skip if not enough sequences to perform voting
    readCluster, index, counter, minReadCount, lock = args
    if readCluster is not None and readCluster.readCounts() > minReadCount:
        sequenceLeft, qualityLeft, supportedLeftReads, sequenceRight, qualityRight, supportedRightReads = concensusPairs(readCluster)
        lock.acquire()
        count = counter.value
        count += 1
        counter.value = count
        lock.release()
        leftRecord = '@cluster_%i %s %i readCluster\n%s\n+\n%s\n' \
            %(count, index, supportedLeftReads, sequenceLeft, qualityLeft)
        rightRecord = '@cluster_%i %s %i readCluster\n%s\n+\n%s\n' \
            %(count, index, supportedRightReads, sequenceRight, qualityRight)
        if count % 100000 == 0:
            stderr.write('[%s] Processed %i read clusters.\n' %(programname, count))
        return(leftRecord,rightRecord)

def hammingDistance(expected_constant, constant_region):
    dist = hamming(list(expected_constant),list(constant_region))
    return dist

def readClustering(read1, read2, barcodeDict, idxBase, barcodeCutOff, constant, constant_length, hamming_threshold, usable_seq):
    """
    generate read cluster with a dictionary object and seqRecord class.
    index of the dictionary is the barcode extracted from first /idxBases/ of read 1
    """
    idLeft, seqLeft, qualLeft = read1
    idRight, seqRight, qualRight = read2
    assert idLeft.split(' ')[0] == idRight.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(idRight, idLeft)
    barcode = seqLeft[:idxBase]
    constant_region = seqLeft[idxBase:usable_seq]
    barcodeQualmean = int(np.mean(map(ord,qualLeft[:idxBase])) - 33)
    if ('N' not in barcode \
            and barcodeQualmean > barcodeCutOff \
            and not any(pattern in barcode for pattern in ['AAAAA','CCCCC','TTTTT','GGGGG']) \
            and hammingDistance(constant, constant_region) <= hamming_threshold):
        seqLeft = seqLeft[usable_seq:]
        qualLeft = qualLeft[usable_seq:]
        barcodeDict[barcode].addRecord(seqRight, qualRight, seqLeft, qualLeft)
    return 0

def writeFile(outputprefix, leftReads, rightReads):
    """
    write fastq lines to gzip files
    """
    # output file name
    read1File = outputprefix + '_R1_001.fastq.gz'
    read2File = outputprefix + '_R2_001.fastq.gz'
    with gzip.open(read1File,'wb') as read1, gzip.open(read2File,'wb') as read2:
        for left, right in zip(leftReads,rightReads):
            assert left.split(' ')[0] == right.split(' ')[0], 'Wrong order pairs!!'
            read1.write(left)
            read2.write(right)
    return read1File, read2File

def plotBCdistribution(barcodeDict, outputprefix):
    #plotting inspection of barcode distribution
    barcodeCount = map(lambda x: barcodeDict[x].readCounts(), barcodeDict.keys())
    barcodeCount = np.array(barcodeCount, dtype=np.int64)
    hist, bins = np.histogram(barcodeCount[barcodeCount<50],bins=50)
    centers = (bins[:-1] + bins[1:]) / 2
    width = 0.7 * (bins[1] - bins[0])
    hist = np.true_divide(hist,np.sum(hist))
    figurename = '%s.png' %(outputprefix)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(centers,hist,align='center',width=width)
    ax.set_xlabel("Number of occurence")
    ax.set_ylabel("Count of tags")
#    ax.set_yscale('log',nonposy='clip')
    ax.set_title(outputprefix.split('/')[-1])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.savefig(figurename)
    stderr.write('Plotted %s.\n' %figurename)
    return 0

def clustering(outputprefix, inFastq1, inFastq2, idxBase, minReadCount, barcodeCutOff, threads, constant):
    manager = Manager()
    barcodeDict = defaultdict(seqRecord)
    read_num = 0
    constant_length = len(constant)
    hamming_threshold = float(1)/constant_length
    usable_seq = idxBase + constant_length
    with gzip.open(inFastq1,'rb') as fq1, gzip.open(inFastq2,'rb') as fq2:
        for read1,read2 in izip(FastqGeneralIterator(fq1),FastqGeneralIterator(fq2)):
            readClustering(read1,read2,barcodeDict, idxBase, barcodeCutOff,
                    constant, constant_length, hamming_threshold, usable_seq)
            read_num += 1
            if read_num % 1000000 == 0:
                stderr.write('[%s] Parsed: %i sequence\n' %(programname,read_num))
    stderr.write('[%s] Extracted: %i barcodes sequence\n' %(programname,len(barcodeDict.keys())))
    p = plotBCdistribution(barcodeDict, outputprefix)

    # From index library, generate error free reads
    # using multicore to process read clusters
    counter = manager.Value('i',0)
    pool = Pool(threads)
    lock = manager.Lock()
    dict_iter = barcodeDict.iteritems()
    iterator = iter([(seq_record, index, counter, minReadCount, lock) for index, seq_record in dict_iter])
    processes = pool.imap_unordered(errorFreeReads, iterator)
    results = [p for p in processes]
    pool.close()
    pool.join()
    results = filter(None,  results)
    # since some cluster that do not have sufficient reads
    # will return None, results need to be filtered
    if (len(results) == 0):
        sys.exit('[%s] No concensus clusters!! \n' %(programname))
    left, right = zip(*results)
    stderr.write('[%s] Extracted error free reads\n' %(programname))
    # use two cores for parallel writing file
    read1File, read2File = writeFile(outputprefix, list(left), list(right))

    # all done!
    stderr.write('[%s] Finished writing error free reads\n' %programname)
    stderr.write('[%s]     read1:            %s\n' %(programname, read1File))
    stderr.write('[%s]     read2:            %s\n' %(programname, read2File))
    stderr.write('[%s]     output clusters:  %i\n' %(programname, len(left)))
    return 0

def main(outputprefix, inFastq1, inFastq2, idxBase, minReadCount,
        barcodeCutOff, threads, constant):
    """
    main function:
        controlling work flow
        1. generate read clusters by reading from fq1 and fq2
        2. obtain concensus sequence from read clusters
        3. writing concensus sequence to files
    """
    start = time.time()

    #print out parameters
    stderr.write( '[%s] Using parameters: \n' %(programname))
    stderr.write( '[%s]     indexed bases:                     %i\n' %(programname,idxBase))
    stderr.write( '[%s]     threads:                           %i\n' %(programname, threads))
    stderr.write( '[%s]     minimum coverage:                  %i\n' %(programname,minReadCount))
    stderr.write( '[%s]     outputPrefix:                      %s\n' %(programname,outputprefix))
    stderr.write( '[%s]     using constant regions:   %s\n' %(programname,constant))

    # divide reads into subclusters
    clustering(outputprefix, inFastq1, inFastq2, idxBase, minReadCount, barcodeCutOff, threads, constant)
    stderr.write('[%s]     time lapsed:      %2.3f min\n' %(programname, np.true_divide(time.time()-start,60)))
    return 0

if __name__ == '__main__':
    outputprefix, inFastq1, inFastq2, idxBase, minReadCount, barcodeCutOff, threads, constant = getOptions()
    main(outputprefix, inFastq1, inFastq2, idxBase, minReadCount, barcodeCutOff, threads, constant)

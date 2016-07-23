
#/usr/bin/env python

from scipy.spatial.distance import hamming
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab
import matplotlib.pyplot as plt
import seaborn as sns
from sys import stderr
import h5py
import gzip
from multiprocessing import Pool, Manager
from itertools import imap
import shelve
sns.set_style('white')
minQ = 33
maxQ = 73
maxProb = 0.999999

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
    qualities = np.empty(no_of_reads,dtype=np.int8)
    for seq, qual, i  in zip(seqList, qualList, range(no_of_reads)):
        columnBases[i] = seq[pos]
        qualities[i] = ord(qual[pos]) -33
    posteriors = [calculatePosterior(guessBase, columnBases, qualities) for guessBase in acceptable_bases]
    likelihoods = np.true_divide(posteriors, np.sum(posteriors))
    argMaxLikHood = np.argmax(likelihoods)
    concensusBase = acceptable_bases[argMaxLikHood]
    posterior = posteriors[argMaxLikHood]
    return concensusBase, posterior

def qualString(posteriors):
    posteriors = np.array(posteriors, dtype=np.float32)
    posteriors[posteriors > maxProb] = maxProb
    quality =  -10 * np.log10(1 - posteriors)
    quality = np.array(quality,dtype=np.int8) + 33
    quality[quality<minQ] = minQ
    quality[quality > maxQ] = maxQ
    quality = ''.join(map(chr,quality))
    return quality

def concensusSeq(seqList, qualList, positions):
    """given a list of sequences, a list of quality and sequence length.
        assertion: all seq in seqlist should have same length (see function: selectSeqLength)
    return a consensus sequence and the mean quality line (see function: calculateConcensusBase)
    """
    if len(seqList) > 1:
        concensusPosition = map(calculateConcensusBase,[(seqList, qualList, pos) for pos in positions])
        bases, posteriors = zip(*concensusPosition)
        sequence = ''.join(list(bases))
        quality = qualString(posteriors)
    else:
        sequence = seqList[0]
        quality = qualList[0]
    return sequence, quality

def hammingDistance(expected_constant, constant_region):
    dist = hamming(list(expected_constant),list(constant_region))
    return dist

def plotBCdistribution(barcodeCount, outputprefix):
    #plotting inspection of barcode distribution
    barcodeCount = np.array(barcodeCount, dtype=np.int64)
    hist, bins = np.histogram(barcodeCount[barcodeCount<50],bins=50)
    centers = (bins[:-1] + bins[1:]) / 2
    width = 0.7 * (bins[1] - bins[0])
    figurename = '%s.pdf' %(outputprefix)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(centers,hist,align='center',width=width)
    ax.set_xlabel("Number of occurence")
    ax.set_ylabel("Count of tags")
    ax.set_yscale('log',nonposy='clip')
    ax.set_title(outputprefix.split('/')[-1])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.savefig(figurename)
    stderr.write('Plotted %s.\n' %figurename)
    return 0

def dictToh5File(barcodeDict, h5_file):
    """
    converting sequence dict to a h5 file for minimizing memory use
    """
    with h5py.File(h5_file,'w') as h5:
        group = h5.create_group('barcodes')
        [group.create_dataset(index, data = index_family) for index, index_family in barcodeDict.iteritems()]
    print 'Finished writting %s'  %h5_file
    return 0

def concensusPairs(table):
    """ given a pair of reads as defined as the class: seqRecord
    return concensus sequence and mean quality of the pairs,
        as well as the number of reads that supports the concnesus pairs
    see function: concensusSeq, calculateConcensusBase
    """
    # get concensus left reads first
    sequenceLeft, qualityLeft = concensusSeq(table[:,0], table[:,2], range(len(table[:,0][0])))
    assert len(sequenceLeft) == len(qualityLeft), 'Wrong concensus sequence and quality!'
    # get concensus right reads first
    sequenceRight, qualityRight = concensusSeq(table[:,1], table[:,3], range(len(table[:,1][0])))
    assert len(sequenceRight) == len(qualityRight), 'Wrong concensus sequence and quality!'
    return sequenceLeft, qualityLeft, sequenceRight, qualityRight

def errorFreeReads(args):
    """
    main function for getting concensus sequences from read clusters.
    return  a pair of concensus reads with a 4-line fastq format
    see functions: 1. filterRead,
                  2. concensusPairs,
                  3. calculateConcensusBase
    """
    # skip if not enough sequences to perform voting
    index, h5_file, minReadCount = args
    leftRecord, rightRecord = 0, 0
    with h5py.File(h5_file,'r') as h5:
        table = h5['barcodes'][index]
        member_count = table.shape[1]
        if member_count >= minReadCount:
            sequenceLeft, qualityLeft, sequenceRight, qualityRight = concensusPairs(table)
            leftRecord = '%s_%i_readCluster\n%s\n+\n%s\n' %(index, member_count, sequenceLeft, qualityLeft)
            rightRecord = '%s_%i_readCluster\n%s\n+\n%s\n' %(index, member_count, sequenceRight, qualityRight)
    return leftRecord, rightRecord

@profile
def writingAndClusteringReads(outputprefix, minReadCount, h5_file, threads):
    # From index library, generate error free reads
    # using multicore to process read clusters
    counter = 0
    output_cluster_count = 0
    read1File = outputprefix + '_R1_001.fastq.gz'
    read2File = outputprefix + '_R2_001.fastq.gz'
    with h5py.File(h5_file) as h5:
        indexes = h5['barcodes'].keys()
    with gzip.open(read1File,'wb') as read1, gzip.open(read2File,'wb') as read2:
        args = ((str(index).strip(), h5_file, minReadCount) for index in indexes)
        pool = Pool(threads)
        processes = pool.imap_unordered(errorFreeReads, args, chunksize = 10000)
        #processes = imap(errorFreeReads, args)
        for p in processes:
            counter += 1
            if counter % 1000000 == 0:
                stderr.write('Processed %i read clusters.\n' %(counter))
            if p != (0,0):
                leftRecord, rightRecord = p
                read1.write('@cluster%i_%s' %(output_cluster_count, leftRecord))
                read2.write('@cluster%i_%s' %(output_cluster_count, rightRecord))
                output_cluster_count += 1
    pool.close()
    pool.join()
    return output_cluster_count, read1File, read2File

def errorFreeReadsDict(args):
    """
    main function for getting concensus sequences from read clusters.
    return  a pair of concensus reads with a 4-line fastq format
    see functions: 1. filterRead,
                  2. concensusPairs,
                  3. calculateConcensusBase
    """
    # skip if not enough sequences to perform voting
    table, index, minReadCount = args
    leftRecord, rightRecord = 0, 0
    member_count = table.shape[1]
    if member_count >= minReadCount:
        sequenceLeft, qualityLeft, sequenceRight, qualityRight = concensusPairs(table)
        leftRecord = '%s_%i_readCluster\n%s\n+\n%s\n' %(index, member_count, sequenceLeft, qualityLeft)
        rightRecord = '%s_%i_readCluster\n%s\n+\n%s\n' %(index, member_count, sequenceRight, qualityRight)
    return leftRecord, rightRecord

@profile
def writingAndClusteringReadsDict(outputprefix, minReadCount, barcode_dict, threads):
    # From index library, generate error free reads
    # using multicore to process read clusters
    counter = 0
    output_cluster_count = 0
    read1File = outputprefix + '_R1_001.fastq.gz'
    read2File = outputprefix + '_R2_001.fastq.gz'
    args =((np.array(table), index, minReadCount) for index,table in barcode_dict.iteritems())
    with gzip.open(read1File,'wb') as read1, gzip.open(read2File,'wb') as read2:
        pool = Pool(threads)
        processes = pool.imap_unordered(errorFreeReadsDict, args)
        #processes = imap(errorFreeReads, args)
        for p in processes:
            counter += 1
            if counter % 1000000 == 0:
                stderr.write('Processed %i read clusters.\n' %(counter))
            if p != (0,0):
                leftRecord, rightRecord = p
                read1.write('@cluster%i_%s' %(output_cluster_count, leftRecord))
                read2.write('@cluster%i_%s' %(output_cluster_count, rightRecord))
                output_cluster_count += 1
    pool.close()
    return output_cluster_count, read1File, read2File

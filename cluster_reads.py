
#/usr/bin/env python

from scipy.spatial.distance import hamming
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab
import matplotlib.pyplot as plt
import seaborn as sns
from sys import stderr
import cjson
import gzip
from multiprocessing import Pool, Manager
from itertools import imap,izip
import shelve
from functools import partial
sns.set_style('white')
min_q = 33
max_q = 73
max_prob = 0.999999
acceptable_bases = np.array(['A','C','T','G'], dtype='string')

def qualToString(posteriors):
    posteriors = np.array(posteriors, dtype=np.float64)
    posteriors[posteriors > max_prob] = max_prob
    quality =  -10 * np.log10(1 - posteriors)
    quality = np.array(quality,dtype=np.int8) + 33
    quality[quality<min_q] = min_q
    quality[quality > max_q] = max_q
    quality = ''.join(map(chr,quality))
    return quality

def qualToStringToInt(q):
    return ord(q)-33

def qual2Prob(base_qual):
    '''
    Given a q list,
    return a list of prob
    '''
    return np.power(10, np.true_divide(-base_qual,10))

def calculatePosterior(column_bases, column_qualities, guess_base):
    qual_missed = column_qualities[column_bases!=guess_base]
    qual_hit = column_qualities[column_bases==guess_base]
    hit = np.prod(1- qual2Prob(qual_hit))
    missed = np.prod(np.true_divide(qual2Prob(qual_missed),3))
    posterior = missed * hit
    return posterior

def calculateConcensusBase(arg):
    """Given a list of sequences,
        a list of quality line and
        a position,
    return the maximum likelihood base at the given position,
        along with the mean quality of these concensus bases.
    """
    seq_list, qual_list, pos = arg
    no_of_reads = len(seq_list)
    column_bases = np.empty(no_of_reads,dtype='string')
    column_qualities = np.zeros(no_of_reads,dtype=np.int8)
    column_bases = seq_list[:,pos]
    column_qualities = np.array(map(qualToStringToInt,qual_list[:,pos]))
    bases = np.unique(column_bases)
    if len(bases) == 1:
        posterior_correct_probability = 1
        concensus_base = bases[0]
    else:
        posteriors = [calculatePosterior(column_bases, column_qualities, guess_base) for guess_base in bases]
        likelihoods = np.true_divide(posteriors, np.sum(posteriors))
        arg_max_likelihood = np.argmax(likelihoods)
        concensus_base = bases[arg_max_likelihood]
        posterior_correct_probability = likelihoods[arg_max_likelihood]
    return concensus_base, posterior_correct_probability

def concensusSeq(seq_list, qual_list):
    """given a list of sequences, a list of quality and sequence length.
        assertion: all seq in seqlist should have same length (see function: selectSeqLength)
    return a consensus sequence and the mean quality line (see function: calculateConcensusBase)
    """
    if len(seq_list) > 1:
        seq_len = len(seq_list[0])
        seq_list = np.array(map(list,seq_list))
        qual_list = np.array(map(list, qual_list))
        iter_list = ((seq_list, qual_list, pos) for pos in xrange(seq_len))
        concensus_position = map(calculateConcensusBase, iter_list)
        bases, posterior_error_probs = zip(*concensus_position)
        sequence = ''.join(list(bases))
        quality = qualToString(posterior_error_probs)
    else:
        sequence = seq_list[0]
        quality = qual_list[0]
    return sequence, quality

def hammingDistance(expected_constant, constant_region):
    dist = hamming(list(expected_constant),list(constant_region))
    return dist

def plotBCdistribution(barcode_family_count, outputprefix):
    #plotting inspection of barcode distribution
    barcode_family_count = np.array(barcode_family_count, dtype=np.int64)
    hist, bins = np.histogram(barcode_family_count[barcode_family_count<50],bins=50)
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

def concensusPairs(table):
    """ given a pair of reads as defined as the class: seqRecord
    return concensus sequence and mean quality of the pairs,
        as well as the number of reads that supports the concnesus pairs
    see function: concensusSeq, calculateConcensusBase
    """
    #extract table
    seq_left_list, qual_left_list = table[:,0], table[:,2]
    seq_right_list, qual_right_list = table[:,1], table[:,3]

    # get concensus left reads first
    sequence_left, quality_left = concensusSeq(seq_left_list, qual_left_list)
    # get concensus right reads
    sequence_right, quality_right = concensusSeq(seq_right_list, qual_right_list)
    return sequence_left, quality_left, sequence_right, quality_right

def dictToJson(barcode_dict, json_file):
    with open(json_file,'w') as f:
        [f.write(cjson.encode(items) + '\n') for items in barcode_dict.iteritems()]
    stderr.write('written %s' %(json_file) + '\n')
    return 0

def errorFreeReads(min_family_member_count, record):
    """
    main function for getting concensus sequences from read clusters.
    return  a pair of concensus reads with a 4-line fastq format
    see functions: 1. filterRead,
                  2. concensusPairs,
                  3. calculateConcensusBase
    """
    # skip if not enough sequences to perform voting
    record = cjson.decode(record)
    index = record[0]
    table = np.array(record[1])
    member_count = table.shape[0]
    if member_count >= min_family_member_count:
        sequence_left, quality_left, sequence_right, quality_right = concensusPairs(table)
        left_record = '%s_%i_readCluster\n%s\n+\n%s\n' %(index, member_count, sequence_left, quality_left)
        right_record = '%s_%i_readCluster\n%s\n+\n%s\n' %(index, member_count, sequence_right, quality_right)
        return left_record, right_record
    else:
        return 'No'

def writingAndClusteringReads(outputprefix, min_family_member_count, json_file, barcode_count,threads):
    # From index library, generate error free reads
    # using multicore to process read clusters
    counter = 0
    output_cluster_count = 0
    read1File = outputprefix + '_R1_001.fastq.gz'
    read2File = outputprefix + '_R2_001.fastq.gz'
    with gzip.open(read1File,'wb') as read1, gzip.open(read2File,'wb') as read2,open(json_file,'r') as infile:
        func = partial(errorFreeReads, min_family_member_count)
        pool = Pool(threads,maxtasksperchild=100)
        for result in pool.imap_unordered(func, infile, chunksize = 1000):
            if result != 'No':
                read1.write('@cluster%i_%s' %(output_cluster_count, result[0]))
                read2.write('@cluster%i_%s' %(output_cluster_count, result[1]))
                output_cluster_count += 1
            counter += 1
            if counter % 1000000 == 0:
                stderr.write('Processed %i read clusters.\n' %(counter))
        pool.close()
        pool.join()
    return output_cluster_count, read1File, read2File

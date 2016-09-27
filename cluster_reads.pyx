

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
from itertools import imap, izip
from functools import partial
from numpy cimport ndarray
import re
sns.set_style('white')

cdef:
    int min_q = 33
    int max_q = 73
    float max_prob = 0.999999
    ndarray acceptable_bases = np.array(['A','C','T','G'], dtype='string')

np_ord = np.vectorize(ord,otypes=[np.int16])

cpdef str qualToString(ndarray posteriors):
    cdef:
        ndarray quality
        str quality_str

    posteriors = np.clip(posteriors, 0, max_prob)
    quality =  -10 * np.log10(1 - posteriors)
    quality = np.array(quality,dtype=np.int16) + 33
    quality = np.clip(quality, min_q, max_q)
    quality_str = ''.join(map(chr,quality))
    return quality_str


cpdef ndarray qualToInt(ndarray qs):
    cdef:
        ndarray out_qs
    out_qs = np_ord(qs) - 33
    return out_qs

cpdef ndarray qual2Prob(ndarray base_qual):
    '''
    Given a q list,
    return a list of prob
    '''
    return np.power(10, np.true_divide(-base_qual,10))


cpdef float calculatePosterior(ndarray column_bases, ndarray column_qualities, guess_base):
    cdef:
        ndarray qual_missed, qual_hit
        float prod, missed, posterior

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
    cdef:
        ndarray column_bases
        ndarray in_column_qualities
        ndarray column_qualities
        ndarray bases, likelihoods, posteriors
        float posterior_correct_probability
        int arg_max_likelihood

    column_bases, in_column_qualities = arg
    column_qualities = qualToInt(in_column_qualities)
    bases = np.unique(column_bases)
    if len(bases) == 1:
        posterior_correct_probability = 1
        concensus_base = bases[0]
    else:
        posteriors = np.array([calculatePosterior(column_bases, column_qualities, guess_base) for guess_base in bases])
        likelihoods = np.true_divide(posteriors, np.sum(posteriors))
        arg_max_likelihood = np.argmax(likelihoods)
        concensus_base = bases[arg_max_likelihood]
        posterior_correct_probability = likelihoods[arg_max_likelihood]
    return concensus_base, posterior_correct_probability

def voteConcensusBase(arg):
    """Given a list of sequences,
        a list of quality line and
        a position,
    return the maximum likelihood base at the given position,
        along with the mean quality of these concensus bases.
    """
    cdef:
        ndarray column_bases, column_qualities, in_column_qualities
        ndarray bases, counts, posteriors
        float prob, likelihood
        int depth

    column_bases, in_column_qualities = arg
    depth = len(column_bases)
    column_qualities = qualToInt(in_column_qualities)
    bases, counts = np.unique(column_bases, return_counts = True)
    if np.true_divide(max(counts), depth) > 0.66:
        base = bases[np.argmax(counts)]
        if len(bases) == 1:
            prob = 1 - np.prod(qual2Prob(column_qualities))
        else:
            posteriors = np.array([calculatePosterior(column_bases, column_qualities, guess_base) for guess_base in bases])
            likelihoods = np.true_divide(posteriors, np.sum(posteriors))
            prob = likelihoods[bases == base]
        #sum_qual = np.sum(column_qualities[column_bases == base])
        #prob = 1 - qual2Prob(np.array([sum_qual]))[0]
    else:
        base = 'N'
        prob = 1
    return base, prob

def concensusSeq(ndarray in_seq_list, ndarray in_qual_list):
    """given a list of sequences, a list of quality and sequence length.
        assertion: all seq in seqlist should have same length (see function: selectSeqLength)
    return a consensus sequence and the mean quality line (see function: calculateConcensusBase)
    """
    cdef:
        int seq_len, pos
        ndarray seq_list, qual_list
        str sequence, quality

def concensusSeq(ndarray in_seq_list, ndarray in_qual_list):
    """given a list of sequences, a list of quality and sequence length.
        assertion: all seq in seqlist should have same length (see function: selectSeqLength)
    return a consensus sequence and the mean quality line (see function: calculateConcensusBase)
    """
    cdef:
        int seq_len, pos
        ndarray seq_list, qual_list
        str sequence, quality

    if len(in_seq_list) > 1:
        seq_len = len(in_seq_list[0])
        seq_list = np.array(map(list, in_seq_list))
        qual_list = np.array(map(list, in_qual_list))
        iter_list = ((seq_list[:,pos], qual_list[:,pos]) for pos in xrange(seq_len))
        #concensus_position = map(calculateConcensusBase, iter_list)
        concensus_position = map(voteConcensusBase, iter_list)
        bases, posterior_error_probs = zip(*concensus_position)
        sequence = ''.join(list(bases))
        quality = qualToString(np.array(posterior_error_probs, dtype=np.float64))
    else:
        sequence = str(in_seq_list[0])
        quality = str(in_qual_list[0])
    return sequence, quality

cpdef float hammingDistance(str expected_constant, str constant_region):
    cdef float dist = hamming(list(expected_constant),list(constant_region))
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

def concensusPairs(ndarray table):
    """ given a pair of reads as defined as the class: seqRecord
    return concensus sequence and mean quality of the pairs,
        as well as the number of reads that supports the concnesus pairs
    see function: concensusSeq, calculateConcensusBase
    """
    #extract table
    cdef:
        ndarray seq_left_list, qual_left_list
        ndarray seq_right_list, qual_right_list
        str sequence_left, quality_left
        str sequence_right, quality_right

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

def errorFreeReads(int min_family_member_count, str json_record):
    """
    main function for getting concensus sequences from read clusters.
    return  a pair of concensus reads with a 4-line fastq format
    see functions: 1. filterRead,
                  2. concensusPairs,
                  3. calculateConcensusBase
    """
    # skip if not enough sequences to perform voting
    cdef:
        str index
        int member_count
        ndarray table
        str sequence_left, quality_left
        str sequence_right, quality_right
        str left_record, right_record

    record = cjson.decode(json_record)
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

def writeSeqToFiles(read1, read2, output_cluster_count, result):
    if result!='No':
        read1.write('@cluster%i_%s' %(output_cluster_count, result[0]))
        read2.write('@cluster%i_%s' %(output_cluster_count, result[1]))
        return 1
    else:
        return 0

def writingAndClusteringReads(outputprefix, min_family_member_count, json_file, threads):
    # From index library, generate error free reads
    # using multicore to process read clusters
    cdef:
        int counter = 0
        int output_cluster_count = 0

    read1File = outputprefix + '_R1_001.fastq.gz'
    read2File = outputprefix + '_R2_001.fastq.gz'
    with gzip.open(read1File,'wb') as read1, gzip.open(read2File,'wb') as read2, open(json_file,'r') as infile:
        error_func = partial(errorFreeReads, min_family_member_count)
        write_func = partial(writeSeqToFiles,read1, read2)
        pool = Pool(threads,maxtasksperchild=1000)
        processes = pool.imap_unordered(error_func, infile, chunksize = 1000)
        #processes = imap(error_func, infile)
        for result in processes:
            output_cluster_count += write_func(output_cluster_count, result)
            counter += 1
            if counter % 1000000 == 0:
                stderr.write('Processed %i read clusters.\n' %(counter))
        pool.close()
        pool.join()
    return output_cluster_count, read1File, read2File



############### clustering #####################
def readClusteringR2(barcode_dict, idx_base, barcode_cut_off, constant,
                   constant_length, hamming_threshold, usable_seq, failed_file,
                   low_complexity_composition, read1, read2):
    """
    generate read cluster with a dictionary object and seqRecord class.
    index of the dictionary is the barcode extracted from first /idx_bases/ of read 1
    """
    id_left, seq_left, qual_left = read1
    id_right, seq_right, qual_right = read2
    assert id_left.split(' ')[0] == id_right.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(id_right, id_left)
    barcode = seq_right[:idx_base]
    constant_region = seq_right[idx_base:usable_seq]
    barcodeQualmean = int(np.mean(map(ord,qual_right[:idx_base])) - 33)

    no_N_barcode = 'N' not in barcode
    is_low_complexity_barcode = bool(re.search(low_complexity_composition, barcode))
    hiQ_barcode = barcodeQualmean > barcode_cut_off
    accurate_constant = hammingDistance(constant, constant_region) <= hamming_threshold

    if no_N_barcode and hiQ_barcode and accurate_constant: #and not is_low_complexity_barcode):
        seq_right = seq_right[usable_seq:]
        qual_right = qual_right[usable_seq:]
        barcode_dict[barcode].append([seq_left,seq_right,qual_left, qual_right])
        return 0
    else:
        failed_file.write('\t'.join([id_left, seq_left, qual_left, seq_right, qual_right]) + '\n')
        return 1

def readClusteringR1(barcode_dict, idx_base, barcode_cut_off, constant,
                     constant_length, hamming_threshold, usable_seq, failed_file,
                     low_complexity_composition, read1, read2):
    """
    generate read cluster with a dictionary object and seqRecord class.
    index of the dictionary is the barcode extracted from first /idx_bases/ of read 1
    """
    id_left, seq_left, qual_left = read1
    id_right, seq_right, qual_right = read2
    assert id_left.split(' ')[0] == id_right.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(id_right, id_left)
    barcode = seq_left[:idx_base]
    constant_region = seq_left[idx_base:usable_seq]
    barcodeQualmean = int(np.mean(map(ord,qual_left[:idx_base])) - 33)

    no_N_barcode = 'N' not in barcode
    is_low_complexity_barcode = bool(re.search(low_complexity_composition, barcode))
    hiQ_barcode = barcodeQualmean > barcode_cut_off
    accurate_constant = hammingDistance(constant, constant_region) <= hamming_threshold


    if no_N_barcode and hiQ_barcode and accurate_constant: #and not low_complexity_barcode:
        seq_left = seq_left[usable_seq:]
        qual_left = qual_left[usable_seq:]
        barcode_dict[barcode].append([seq_left, seq_right, qual_left, qual_right])
        return 0
    return 1

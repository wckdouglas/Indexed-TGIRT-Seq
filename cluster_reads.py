#/usr/bin/env python

from scipy.spatial.distance import hamming
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab
import matplotlib.pyplot as plt
import seaborn as sns
from sys import stderr
sns.set_style('white')
minQ = 33
maxQ = 73

#    ==================      Sequence class sotring left right record =============
class seqRecord:
    def __init__(self):
        self.seqListRight = []
        self.qualListRight = []
        self.seqListLeft = []
        self.qualListLeft = []
        self.member_count = 0

    def addRecord(self, seqRight, qualRight, seqLeft, qualLeft):
        self.seqListRight.append(seqRight)
        self.qualListRight.append(qualRight)
        self.seqListLeft.append(seqLeft)
        self.qualListLeft.append(qualLeft)
        self.member_count += 1

    def readLengthRight(self):
        return np.array([len(seq) for seq in self.seqListRight],dtype=np.int64)

    def readLengthLeft(self):
        return np.array([len(seq) for seq in self.seqListLeft],dtype=np.int64)

#======================  starting concensus functions =============================

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
#    quality = -10 * np.log10(1 - posterior) if posterior < 1 else maxQ
    return concensusBase, posterior #, quality

def qual_string(posteriors):
    posteriors = np.array(posteriors, dtype=np.float64)
    posteriors[posteriors == 1] = 0.999999
    quality =  -10 * np.log10(1 - posteriors)
    quality = np.array(quality,dtype=np.int64)  + 33
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
	quality = qual_string(posteriors)
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

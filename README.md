# Indexed-TGIRT-Seq

This folder contain scripts for tag-based error correction of TGIRT-seq developed by [Douglas Wu](wckdouglas@gmail.com).

Scripts:
	readClusterPairs.py		using maximum likelihood to infer the concensus base from a cluster of reads. A [note](http://rawgit.com/wckdouglas/Indexed-TGIRT-Seq/master/notes/tagBased-Error.html) on developing the tool is store in the repo as well (depends on biopython, numpy)
	pileupBam.py			pileup bam file and  count the bases at each position column (depends on pysam) 

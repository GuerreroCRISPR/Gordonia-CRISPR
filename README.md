# Gordonia-CRISPR
Scripts used to the reconstruction of Gordonia CRISPR

This pipeline allow to reconstruct a specific CRISPR array diversity from several samples using illumina pair-end reads. Reads length must be >=100 bp or longer to get a better result. As a rule the read length >= 2 spacers + 1 repeat.

May be possible to use filtered files by read length and quality 

Requirements:

bbduck (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide)

Mothur (https://mothur.org)

Python modules:

Biopython (pip install biopython)

Use some of the available software to identify possible CRISPR arrays and extract the repeat sequence. Consider the number of mismatches.

Example:




Using this information run bbduk over all samples in order to extract the reads with the repeat sequence:

bbduk.sh
in=R1_sample_1.fastq
in2=R2_sample_1.fastq
outm=R1_matched_sample_1.fq
outm2=R2_matched_sample_1.fq
k=			# k-mer length (max 31)
mm=f
literal= 		# repeat sequence   
hdist= 		# allowed mismatches
rcomp=T


Make a reads.names.txt file of the reads samples fastq files as shown:

sample1_R1.fastq


Make a source file to be used by the script findRepeatCRISPR.py:

The source file is a tab delimited file with a sample identifier and the reads file names:

sample1	R1.sample1.fastq	R2.sample1.fastq
sample2	R1.sample2.fastq	R2.sample2.fastq
sample3	R1.sample3.fastq	R2.sample3.fastq
…
samplen	R1.samplen.fastq	R2.samplen.fastq

# Gordonia-CRISPR
Scripts used to the reconstruction of Gordonia CRISPR. All these scripts run under Linux operating system.

This pipeline allow to reconstruct a specific CRISPR array diversity from several samples using illumina pair-end reads. Reads length must be >=100 bp or longer to get a better result. As a rule the read length >= 2 spacers + 1 repeat.

May be possible to use filtered files by read length and quality 

Requirements:

bbduck (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide)

Mothur (https://mothur.org)

Blast 

Python modules:

Biopython (pip install biopython)

Use some of the available software to identify possible CRISPR arrays and extract the repeat sequence. Consider the number of mismatches.


-Using this information run bbduk over all samples in order to extract the reads with the repeat sequence:

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



-Make a source file to be used by the script findRepeatCRISPR.py:

The source file is a tab delimited file with a sample identifier and the reads file names:

sample1	R1.sample1.fastq	R2.sample1.fastq

sample2	R1.sample2.fastq	R2.sample2.fastq

sample3	R1.sample3.fastq	R2.sample3.fastq

…

samplen	R1.samplen.fastq	R2.samplen.fastq


-Run the following Mothur scripts in order to process the sequences. Choose the desired parameters:

screen.seqs(fasta=, maxlength=130)

unique.seqs(fasta=)

pre.cluster(fasta=, name=, diffs=3)

cluster.fragments(fasta=, name=, diffs=3)

-Mothur can't resolve the merge of identical sequences of diferent length. To overcome this problem run the mergeLengthCRISPR.py script.

-Use the generated files to make the input nodes and edges files using the networkCRISPR.py script. These two files can be uploaded to gephi to visualize the raw network.

-In order to filter the network eliminating low quality connections and nodes, the qualityCRISPR.py script uses a series of files to generate a file which can be used into gephi:

1) The last .fasta file created (mergeLengthCRISPR.py step).

2) The merged reads (fastq) files for all samples (e.g.: all.R1.fastq,all.R2.fastq).

3) The last .name file created (mergeLengthCRISPR.py step).

4) The original .seq file created at the begining (findRepeatCRISPR.py step).




-Spacers fasta file creation

The easiest way to extract the spacers as fasta file is from the network.nodes.csv file:


sed 's/ /\t/g' $bin.network.nodes.csv | cut -f 1,3 | sed 's/^/>/' | sed 's/\t/\n/' | tail -n+3 > $bin.spacersL.fa

sed 's/ /\t/g' $bin.network.nodes.csv | cut -f 1,5 | sed 's/^/>/' | sed 's/\t/\n/' | tail -n+2 > $bin.spacersR.fa






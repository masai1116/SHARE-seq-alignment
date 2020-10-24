# SHARE-seq-alignment
Pipeline for demultiplex and align both ATAC and RNA data generated in SHARE-seq

Dependencies: GNU parallel, Bcl2fastq, fastp, zcat, STAR, bowtie2, python2, umi_tools, samtools, picard, R, featureCounts, read_distribution.py from RSeQC, bedtools

After installing, update the following paths in main script "Split_seq_example.sh":
1) myPATH # change it to where the SHARE-seq scripts are installed. e.g. myPATH='/mnt/users/Script/share-seq-github-v1/'
2) pythohPATH # change it to where python2 is installed e.g. pythohPATH='/usr/bin/python' 
3) picardPATH # change it to where picard is installed e.g. picardPATH='/mnt/bin/picard/picard.jar'

# How to run the script?
A small set of fastq files for testing are in the test_fastq_nova/ folder

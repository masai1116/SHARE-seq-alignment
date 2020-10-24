# SHARE-seq-alignment
Pipeline for demultiplex and align both ATAC and RNA data generated in SHARE-seq

Dependencies: GNU parallel, Bcl2fastq, fastp, zcat, STAR, bowtie2, python2, umi_tools, samtools, picard, R, featureCounts, read_distribution.py from RSeQC, bedtools

After installing, update the general configuration section in main script "Split_seq_example.sh":
1) myPATH # change it to where the SHARE-seq scripts are installed. e.g. myPATH='/mnt/users/Script/share-seq-github-v1/'
2) pythohPATH # change it to where python2 is installed e.g. pythohPATH='/usr/bin/python' 
3) picardPATH # change it to where picard is installed e.g. picardPATH='/mnt/bin/picard/picard.jar'

# How to run the script?
A small set of fastq files for testing are in the test_fastq_nova/ folder
Before running, three sections in the main script "Split_seq_example.sh" need to be updated for each run, inlcuding 
A) paths B) sample configuration C) fastq configuration

## sample configuration
1) Project=(sp.rna sp.atac.first) # use differnt name for each sample, the sample 
2) Type=(RNA ATAC)  # ATAC or RNA
3) Genomes=(hg19 both) # both mm10 hg19
RawReadsPerBarcode and ReadsPerBarcode options are designed to remove barcode with few reads and speed up processing. 
4) RawReadsPerBarcode=(10 10) # reads cutoff for the unfiltered bam file. Recommend to use 100 for full run; 10 for QC run
5) ReadsPerBarcode=(1 1) # reads cutoff for the filtered bam file. Recommend to use 100 for full run, 1 for QC run
6) keepMultiMapping=(F F)  # T or F; default is F. Keep or discard multi-mapping reads

## fastq configuration
1) Indexed=F # T or F; defaul is F. Indicate if the index reads are already attached to biological reads. Use F, when started with BCL file.
2) Start=Fastq_Merge # Bcl or Fastq_Merge (when fastq were generated per run) or Fastq_SplitLane (when fastq were generated per sequencing lane)
3) Runtype=full # QC or full;  QC only analyze 12M reads to get a quick sense of data
4) Sequencer=Novaseq # Novaseq or Nextseq;  miseq or nova-seq has the same sequencing direction, use "Novaseq" for either

## paths
rawdir=./test_fastq_nova/ # there the raw data is e.g. ~/xxx/201021_SL-NVQ_0277_AHTHLLDRXX/ or ./test_fastq_nova/
dir=~/test/ # where output data will be stored
ymal=./config_example.ymal # where the ymal configuration file is

## RNA-seq options
The pipeline also offers flexible RNA-seq specific options for advanced users. 
1) removeSingelReadUMI=F # T or F; default is F. If T, UMIs with single read will be removed.
2) keepIntron=T # T or F; default is T. If F, intronic RNA reads will be discarded.
3) matchPolyT=F # T or F; default is F. If T, will try to find match of TTTTTT (allowing 1 mis-match) in 11-16 bp position of biological read2. If match is not identified, read will be disgarded. Only works if Read2 is longer than 16 bp.
4) SkipPolyGumi=F # T or F; default is F, pipeline will remove polyG UMIs. If T, pipeline will keep polyG UMIs.
5) genename=gene_name # gene_name (official gene symbol) or gene_id (ensemble gene name), gene_name is default
6) refgene=gencode # gencode or genes; gencode is default; genes is UCSC refseq genes; gencode also annotates nc-RNA

# sample barcode table
SHARE-seq allows mutiplexing samples in one run.. This alignment pipeline use ymal file to store two levels of sample barcode information, including Round1 hybridization barcode (R1.xx), and PCR barcode (P1.xx). See "config_example.ymal" as an example. 

# Pipeline output
After successfully running the pipeline, data for each project will be kept in the folder with the same name as the project. 
Several QC plots will be generated. 
For ATAC, the most important files are the fragment file (*.rmdup.cutoff.bed.gz) and a summarize report (*.counts.csv.gz).
For RNA, the most important files are UMIxCell matrix (*.UMIcounts.csv.gz) and a summarize report (*.counts.csv.gz).

This pipeline currently keeps many intermedia files. If preferred, they can be manually removed afterwards.

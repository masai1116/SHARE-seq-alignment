<html>
  <head>
    <meta name="google-site-verification" content="FC8xEq44yQepn6GHyTtgwE8snz8H4lnnfgcJE16CIHY" />
  </head>
</html>

# SHARE-seq-alignment v1.0
Pipeline for demultiplexing and aligning both ATAC and RNA data generated in SHARE-seq\
The pipeline is updated to a newer version. Please refer to V2 for future usage.\
[https://github.com/masai1116/SHARE-seq-alignmentV2/](https://github.com/masai1116/SHARE-seq-alignmentV2/)
**Author: Sai Ma. sai@broadinstitute.org**

# Important note
The ATAC and RNA barcodes would in the format of R1.xx,R2.xx,R3.xx,P1.xx, where xx is in the range of 01-96.\
The P1.xx indicates the primers that were used to amplify each ATAC and RNA sub-libraries. It is expected to be different between RNA and ATAC assay. The rest part of barcode (R1.xx,R2.xx,R3.xx) would be the same for ATAC and RNA library.\
Depending on the downstream analysis pipeline, sometimes the comma in the barcode would be converted to period (R1.xx,R2.xx,R3.xx,P1.xx --> R1.xx.R2.xx.R3.xx.P1.xx).\
A barcode translation table indicating the P1.xx used in the assay would be beneficial, when submitting data to GEO.

# Installation
This pipeline requires following packages to be properly installed and added to system path: GNU parallel, Bcl2fastq, fastp, zcat, STAR, bowtie2, python2, umi_tools, samtools, picard (2.14.1, newer version may result in error), R, featureCounts, read_distribution.py from RSeQC, bedtools

The SHARE-seq-alignment scripts can be directly downloaded from the github website.\
[https://github.com/masai1116/SHARE-seq-alignment/](https://github.com/masai1116/SHARE-seq-alignment/)

After downloading all scripts, update the general configuration section in main script "Split_seq_example.sh":
1) myPATH # where the SHARE-seq scripts are installed. e.g. myPATH='/mnt/users/Script/share-seq-github-v1/'
2) pythohPATH # where python2 is installed e.g. pythohPATH='/usr/bin/python' 
3) picardPATH # where picard is installed e.g. picardPATH='/mnt/bin/picard/picard.jar'

The pipeline also requres gtf files and aligner index files to be download and unziped into the right location.\
GTF files can be downloaded [here](https://drive.google.com/file/d/1HuGLf0vSHO58Ek5HibTRiwXWBn9fBMTz/view?usp=sharing).\
Bowtie2 index files (Hg19 and mm10) can be downloaded [here](https://drive.google.com/file/d/1bXIxznwirsZ6DZhqK1gw6ZKlj-UjFRhn/view?usp=sharing).\
Assuming SHARE-seq aligment scripts are installed to "/home/SHARE-seq-alignment/", the gtf files should be placed in the "/home/SHARE-seq-alignment/gtf" folder.\
The bowtie2 index files should be placed in the "/home/SHARE-seq-alignment/refGenome/bowtie2" folder.\
Three sets of index files (hg19, mm10 and hg19-mm10 combined genome) for star aligner should be prepared according to star aligner [manual](https://github.com/alexdobin/STAR), or downloded from here: [hg19](https://drive.google.com/file/d/1IXI4DP-mjh2qc-EQe1WnWJQOCVEI4KVX/view?usp=sharing), [mm10](https://drive.google.com/file/d/1n0UwzOeUbX7TIBOrcBbXjgH3i0UH-Ka5/view?usp=sharing), [combined genome](https://drive.google.com/file/d/15Z2YMUDiavYG0s9zLFAbbwqA0VhVNu-f/view?usp=sharing).\
The unziped index files should be placed in the "/home/SHARE-seq-alignment/refGenome/star/hg19", "/home/SHARE-seq-alignment/refGenome/star/mm10", and "/home/SHARE-seq-alignment/refGenome/star/both", respectively.\
The index file for hg19-mm10 combined genome can be downloaded from [10x Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest).

# How to run the script?
A small set of fastq files for testing are in the test_fastq_nova/ folder
Before running, three sections in the main script "Split_seq_example.sh" need to be updated for each run, inlcuding 
A) paths B) sample configuration C) fastq configuration. After update all specific information in Split_seq_example.sh and config_example.ymal, run the script by ```./Split_seq_example.sh```
## A) paths
1) rawdir=./test_fastq_nova/ # where the raw data is
2) dir=~/test/ # where output data will be stored
3) ymal=./config_example.ymal # where the ymal configuration file is

## B) sample configuration
1) Project=(sp.rna sp.atac.first) # use differnt name for each sample 
2) Type=(RNA ATAC)  # ATAC or RNA
3) Genomes=(hg19 both) # both mm10 hg19 \
RawReadsPerBarcode and ReadsPerBarcode options are designed to remove barcodes with too few reads and speed up processing. 
4) RawReadsPerBarcode=(10 10) # reads cutoff for the unfiltered bam file. Recommend to use 100 for full run; 10 for QC run
5) ReadsPerBarcode=(1 1) # reads cutoff for the filtered bam file. Recommend to use 100 for full run, 1 for QC run
6) keepMultiMapping=(F F)  # T or F; default is F. Keep or discard multi-mapping reads

## C) fastq configuration
1) Indexed=F # T or F; defaul is F. Indicate if the index reads are already attached to biological reads. Use F, when started with BCL file.
2) Start=Fastq_Merge # Bcl or Fastq_Merge (when fastq were generated per run) or Fastq_SplitLane (when fastq were generated per sequencing lane)
3) Runtype=full # QC or full;  QC only analyze 12M reads to get a quick sense of data
4) Sequencer=Novaseq # Novaseq or Nextseq;  miseq or nova-seq has the same sequencing direction, use "Novaseq" for either

## RNA-seq options
The pipeline also offers flexible RNA-seq specific options for advanced users. 
1) removeSingelReadUMI=F # T or F; default is F. If T, UMIs with single read will be removed.
2) keepIntron=T # T or F; default is T. If F, intronic RNA reads will be discarded.
3) matchPolyT=F # T or F; default is F. If T, will try to find TTTTTT (allowing 1 mis-match) in 11-16 bp position of biological read2. If TTTTTT is not identified, read will be disgarded. Only works if Read2 is longer than 16 bp.
4) SkipPolyGumi=F # T or F; default is F, pipeline will remove polyG UMIs. If T, pipeline will keep polyG UMIs.
5) genename=gene_name # gene_name (official gene symbol) or gene_id (ensemble gene name), gene_name is default
6) refgene=gencode # gencode or genes; gencode is default; genes is UCSC refseq genes; gencode also indcludes nc-RNA

# Sample barcode table
SHARE-seq allows mutiplexing samples in one run. We use ymal file to store two levels of sample barcode information, including Round1 hybridization barcode (R1.xx), and PCR barcode (P1.xx, refers to the Ad1.xx primers used in the PCR step). See ```config_example.ymal``` as an example. This file needs to be updated for each sample and each sequencing run. When all the 96 barcodes in the Round1 plate are used during the experiment, see ```config_example2.yaml``` as an example.
```
---
Project1:
    Name: sp.atac.first
    Primer:
        - P1.01
        - P1.02
    Round1:
        - R1.05
        - R1.13
        - R1.21
        - R1.29
...
```        
R1.xx can be R1.01, R1.02, ..., R1.96.\
P1.xx can be P1.01, P1.02, ..., P1.96.\
The detialed information about these barcode can be found in [SHARE-seq manuscript](https://www.sciencedirect.com/science/article/pii/S0092867420312538).

# Pipeline output
After successfully running the pipeline, data for each project will be kept in the folder with the same name as the project. 
Several QC plots will be generated.\
For ATAC, the most important files are the fragment file after removing duplicates and mito reads (.rmdup.cutoff.bed.gz) and a summarize report (.counts.csv.gz).\
For RNA, the most important files are UMIxCell matrix (.UMIcounts.csv.gz) and a summarize report (.counts.csv.gz).\
This pipeline currently keeps many intermedia files. If preferred, they can be manually removed afterwards.

# Read data example
Public SHARE-seq datasets are available at [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140203).\
A set of fastq files for species mixing experiment can be downloaded [here](https://drive.google.com/drive/folders/19HdjJuWrpRJz8OeB6YNUMSojPyTMV7OP?usp=sharing).

# Cite us
For more details, please refer to [Ma et al. Chromatin Potential Identified by Shared Single-Cell Profiling of RNA and Chromatin, Cell 2020](https://www.sciencedirect.com/science/article/pii/S0092867420312538)

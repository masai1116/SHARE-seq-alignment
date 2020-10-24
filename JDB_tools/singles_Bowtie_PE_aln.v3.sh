#!/bin/bash

# read from command line which files to align
ref=$1 # Reference 
list=$2  # Sample list
dir=$3  # Directory

# help
if [ -z "$dir" ]
then
    echo "This will align a directory of single cells"
    echo "First input is the reference"
    echo "Second input is the sample list"
    echo "Third input is the directory"
    exit
fi

# make fastq dir
mkdir fastqs
mkdir tmp
#set -o errexit
set -E

# set directories
toolPATH='/mnt/Apps/JDB_tools/'
tssFilesPATH='/mnt/Apps/JDB_tools/01_additionalData/TSSfiles/'
picardPATH='/mnt/bin/picard/picard.jar'

# loop through
while read sample
do
    # get data names
    file=`echo $sample | cut -d" " -f1`
    
    # pass if header
    if [ $file = "Name" ]
    then
	continue
    fi

    # get files
    echo "Processing: " $file
    p1=`ls $dir/$file"_"*R1*.fastq.gz | sed 's/.trim.fastq/.fastq/g'`
    p2=`echo $p1 | sed 's/R1/R2/g'`
	
    # trim
    trimp1=`basename $p1 | sed 's/.fastq/.trim.fastq/g'`
    trimp2=`basename $p2 | sed 's/.fastq/.trim.fastq/g'`
    if [ -f fastqs/$trimp1 ]   # if already in dir
    then
	echo "Found: " $trimp1 "and" $trimp2
    elif [ -f $dir/$trimp1 ]   # if trimmed then move to fastqs dir
    then
	echo "Found: " $trimp1 "and" $trimp2 " in:" $dir
	mv $dir/$trimp1 fastqs/; mv $dir/$trimp2 fastqs/
	#cp $dir/$trimp1 fastqs/; cp $dir/$trimp2 fastqs/
    else   # trim then move to fastqs dir
	echo "Trimming: " $p1 "and" $p2
	$toolPATH/pyadapter_trim.py -a $p1 -b $p2
	mv $trimp1 fastqs/; mv $trimp2 fastqs/
	#cp $trimp1 fastqs/; cp $trimp2 fastqs/ 
    fi
    
    # make out dir
    if [ -d $file ]
    then
	echo "Directory exists: " $file
    else
	mkdir $file
    fi
    
    # align
    out1=$file/`basename $p1 | sed 's/.fastq.gz/.bam/g' | sed 's/_R1//g'`
    out2=$file/`basename $p1 | sed 's/.fastq.gz/.align.log/g' | sed 's/_R1//g'`
    out3=`echo $out1 | sed 's/.bam/.st/'`
    if [ -f $out1 ]
    then
	echo "Found: " $out1
    elif [ -f $out3.bam ]
    then
	echo "Found: " $out3.bam
    else
	echo "Aligning: " $trimp1 "and" $trimp2
	(bowtie2 -X2000 -p4 --rg-id $file -x $ref -1 <(gunzip -c fastqs/$trimp1) -2 <(gunzip -c fastqs/$trimp2) | samtools view -bS - -o $out1) 2>$out2
    fi
    
    # sort
    out3=`echo $out1 | sed 's/.bam/.st/'`
    if [ -f $out3.bam ]
    then
	echo "Found: " $out3.bam
    else
	echo "Sorting: " $out1
	java -jar -Djava.io.tmpdir=`pwd`/tmp $picardPATH SortSam SO=coordinate I=$out1 O=$out3.bam VALIDATION_STRINGENCY=SILENT
	samtools index $out3.bam
	rm $out1  #### remove unsorted bam file
    fi
    
    # remove dups and mito
    out5=$file/`basename $p1 | sed 's/.fastq.gz/.dups.log/g' | sed 's/_R1//g'`
    out6=$out3".rmdup.flt.bam"
    chrs=`samtools view -H $out3.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | grep -v Y | awk '{if(length($0)<6)print}'`
    if [ -f $out6 ]
    then
	echo "Found: " $out6
    else
	echo "Removing duplicates and unwanted chrs: " $out3.bam
	samtools view -b -q 30 -f 0x2 $out3.bam -o temp.bam `echo $chrs`
	java -jar -Djava.io.tmpdir=`pwd`/tmp $picardPATH MarkDuplicates INPUT=temp.bam OUTPUT=$out6 METRICS_FILE=$out5 REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT
	samtools index $out6
    	rm temp.bam
    fi
    
    # get final quality stats
    out7=$file/`basename $p1 | sed 's/.fastq.gz/.stats.log/g' | sed 's/_R1//g'`
    echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Raw" > $out7
    samtools idxstats $out3.bam >> $out7
    echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Filtered" >> $out7
    samtools idxstats $out6 >> $out7
    
    # get iSize histogram
    out8=`echo $out6 | sed 's/.bam/.hist_data.log/'`
    out9=`echo $out6 | sed 's/.bam/.hist_data.pdf/'`
    if [ -f $out8 ]
    then
	echo "Found: " $out8
    else
    	# get insert-sizes
	echo '' > $out8
	java -jar $picardPATH CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT I=$out6 O=$out8 H=$out9 W=1000	
    fi
    
    # make TSS pileup fig
    genome=`basename $ref`
    if [ -f $file/$file.RefSeqTSS ]
    then
	echo "Found: " $file.RefSeqTSS
    else
	echo "Creating TSS pileup"; set +e
	$toolPATH/pyMakeVplot.py -a $out6 -b $tssFilesPATH/$genome.TSS.bed -e 2000 -p ends -v -u -o $file/$file.RefSeqTSS
	set -e
    fi
done < <(grep '' $list)

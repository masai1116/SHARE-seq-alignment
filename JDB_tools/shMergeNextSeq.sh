#!/bin/bash

# read from command line which files to align
p1=$1  # first input is the sample sheet
p2=$2  # second input is the fastq dir

# help
if [ -z "$p1" ]
then
    echo "This will merge fastqs"
    echo "First input is the sample sheet"
    echo "Second input is the fastq dir"
    exit
fi

# process
process='0'

# iterate through samples
while read line
do
    # define output dir
    sample=`echo $line | cut -f1 -d"," | sed 's/_/-/g' | sed 's/\./-/g'`
    #sample=`echo $line | cut -f1 -d","`
    if [[ $process == '1' ]]
    then
        echo "Processing... "$sample
        #nohup cat `ls -v $p2/$sample"_"*R1*.fastq.gz | grep trim` > $sample"_R1.trim.fastq.gz" &
        #nohup cat `ls -v $p2/$sample"_"*R2*.fastq.gz | grep trim`> $sample"_R2.trim.fastq.gz" &
        nohup cat `ls -v $p2/$sample"_"*R1*.fastq.gz | grep -v trim` > $sample"_R1.fastq.gz" &
        cat `ls -v $p2/$sample"_"*R2*.fastq.gz | grep -v trim`> $sample"_R2.fastq.gz"
    fi  

    # keep track
    if [[ $sample == 'Sample-ID' ]]
    then
	process='1'
    fi
    
done < <(grep '' $p1)
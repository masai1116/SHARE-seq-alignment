#!/bin/bash

# The following code will trim and concatenate a dir of fastqs
# example run: shMXp1processor_processing.sh L3_AGGCAG_L003_R1 L3_AGGCAG_L003_R2 L3_AGGCAG_L003

# read from command line which files to trim
read1=$1
read2=$2
name=$3

# initialize
rm temp*
rm ad_trim_log.sh

# initialize
set -e
dir=`pwd`
r1=`ls *$1*`
r2=`ls *$2*`

if [ -z $read1 ]
then
    r1=`ls *R1*`
    r2=`ls *R2*`
    name=`ls *R1* | head -1 | sed 's/_R1//g' | sed 's/.fastq//g' | sed 's/.gz//g'`
fi

# generate script to be run to get all read 1 files
for i in $r1
do
   echo "/mnt/Apps/JDB_tools/pyadapter_trim.py -a " $dir"/"$i >> tempr1
done

# get all read2 files
for j in $r2
do
    echo " -b " $dir"/"$j " &" >> tempr2
done

# paste command file
paste tempr1 tempr2 > ad_trim_log.sh

# run command file and save nohup as log
val=`cat ad_trim_log.sh`
s_paren="("
end_paren=") 2>&1 | cat -u >> log_pyadaptertrim.out"
echo $s_paren $val $end_paren > ad_trim_log.sh
sh ad_trim_log.sh

# after nohup is completed
#r1_files=`ls *R1* | grep trim`
#cat $r1_files > $name"_R1.trim.fastq.gz"
#r2_files=`ls *R2* | grep trim`
#cat $r2_files > $name"_R2.trim.fastq.gz"

set +e
# remove trim files
#rm $r1_files
#rm $r2_files
#mkdir 00_source
#mv *.fastq 00_source/
#mv *.fastq.gz 00_source/
#mv 00_source/$name"_R1.trim.fastq.gz" .
#mv 00_source/$name"_R2.trim.fastq.gz" .

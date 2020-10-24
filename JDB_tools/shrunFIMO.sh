#!/bin/bash
set -o errexit

# read from command line which files to align
fastaFile=$1
outDir=$2

# help
if [ -z "$fastaFile" ]
then
    echo "This will fimo in parrallel"
    echo "shRunFimo.sh <.fasta> <outDir>"
    exit
fi

# if no outDir then make it up
if [ -z "$p2" ]
then
    p2=`echo $fastaFile | sed 's/.fasta//g'`
fi

# make out directory
mkdir $outDir
bg='../01_additionalData/motifFiles/tier1_markov1.norc.txt'
motifFile='../01_additionalData/motifFiles/JASPAR_CORE_2014_vertebrates_withChen.meme'
memePATH='/home/wjg/Applications/MEME/meme_4.10.0/bin/'

# loop through motifs
for i in `grep MOTIF $motifFile | sed 's/:/-/g' | sed 's/--/-/g' | sed 's/ /:/g'`;do 
    echo "Processing: " $i
    motif=`echo $i | cut -f2 -d":"`
    outName=`echo $i | cut -f3 -d":"`
    if [[ $outName == '' ]];then
	outName='Motif'
    fi
    
    # relax threshold is .0001, stringet = .00001
    nohup $memePATH/fimo --bgfile $bg --text --thresh .00005 --motif $motif $motifFile $fastaFile 2>status.log |  awk '{OFS="\t"}{split($2,a,":");split(a[2],b,"-");if(NR>1)print a[1],$3+b[1],$4+b[1],$5,$6,$7,$1}' | gzip > $outDir/$outName"-"$motif.bed.gz &
    while [ ` ps -ef | grep -c $memePATH ` -gt "20" ]; do
	sleep 10s
    done
done

# print
echo "Motif finding completed..."


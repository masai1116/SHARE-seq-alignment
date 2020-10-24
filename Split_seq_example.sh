#!/bin/bash
# Author: Sai Ma <sai@broadinstitute.org>
# Last modified date: 2020/10/20
# Designed for processing SHARE-atac, rna-seq data

## paths
rawdir=/mnt/users/sai/Script/share-seq-github-v1/test_fastq_nova/ # e.g. ~/xxx/201021_SL-NVQ_0277_AHTHLLDRXX/ or ./test_fastq_nova/
dir=~/test/
ymal=/mnt/users/sai/Script/share-seq-github-v1/config_example.ymal

## sample configuration
Project=(sp.rna sp.atac.first) # use differnt name for each sample
Type=(RNA ATAC)  # ATAC or RNA 
Genomes=(hg19 both) # both mm10 hg19
RawReadsPerBarcode=(10 10) # reads cutoff for the unfiltered bam, use 100 to remove barcode with few reads for “full" run and speed up processing ; use 10 for QC runs
ReadsPerBarcode=(1 1) # reads cutoff for the filtered bam 100 for full run of ATAC and RNA, 1 for QC run 
keepMultiMapping=(F F)  # T or F; default is F

## fastq configuration
Indexed=F # T or F; defaul is F. Indicate if the index reads are already attached to biological reads. Use F, when started with BCL file.
Start=Fastq_Merge # Bcl or Fastq_Merge (when fastq were generated per run) or Fastq_SplitLane (when fastq were generated per sequencing lane)
Runtype=full # QC or full;  QC only analyze 12M reads to get a quick sense of data
Sequencer=Novaseq # Novaseq or Nextseq;  miseq or nova-seq should use "Novaseq"

# RNA-seq options (for advanced users)
removeSingelReadUMI=F # T or F; default is F. If T, UMIs with single read will be removed.
keepIntron=T # T or F; default is T. If F, intronic RNA reads will be discarded.
matchPolyT=F # T or F; default is F. If T, will try to find match of TTTTTT (allowing 1 mis-match) in 11-16 bp position of biological read2. If match is not identified, read will be disgarded. Only works if Read2 is longer than 16 bp.
SkipPolyGumi=F # T or F; default is F, pipeline will remove polyG UMIs. If T, pipeline will keep polyG UMIs.
genename=gene_name # gene_name (official gene symbol) or gene_id (ensemble gene name), gene_name is default
refgene=gencode # gencode or genes; gencode is default; genes is UCSC refseq genes; gencode also annotates nc-RNA

## general configuration
myPATH='/mnt/users/sai/Script/share-seq-github-v1/'
pythohPATH='/usr/bin/python'
picardPATH='/mnt/bin/picard/picard.jar'

### don't modify the following
toolPATH=$myPATH/JDB_tools/
tssFilesPATH=$myPATH/TSSfiles/ # cleaned TSSs on unknown chr for mm10
bowtieGenome=$myPATH/refGenome/bowtie2/
starGenome=$myPATH/refGenome/star/
genomeBed=$myPATH/genomeBed/
cores=18 # multi-threads

export SHELL=$(type -p bash)
source ~/.bashrc

echo "the number of projects is" ${#Project[@]}
echo "Running $Runtype pipeline"

if [ ! -d $dir ]; then mkdir $dir; fi
if [ ! -d $dir/fastqs ]; then mkdir $dir/fastqs ; fi
if [ ! -d $dir/temp ]; then mkdir $dir/temp ; fi

cp $myPATH/Split*.sh $dir/
cp $ymal $dir/
cd $dir 
if [ -f $dir/Run.log ]; then rm $dir/Run.log; fi

export PATH="/mnt/users/sai/miniconda2/bin:$PATH"

# if start with bcl file
if [ "$Start" = Bcl ]; then
    echo "Bcl2fastq"
    if [ -f $rawdir/fastqs/Undetermined_S0_R1_001.fastq.gz ] || [ -f $rawdir/fastqs/Undetermined_S1_R1_001.fastq.gz ]; then
	echo "Found Undetermined_S0_L001_I1_001.fastq.gz, skip Bcl2Fastq"
    else
	echo "Converting bcl to fastq"
	mkdir $rawdir/fastqs/
	bcl2fastq -p $cores -R $rawdir --no-lane-splitting --mask-short-adapter-reads 0 -o $rawdir/fastqs/ --create-fastq-for-index-reads  2>>$dir/Run.log
    fi
    cd $rawdir/fastqs/
    
    mv Undetermined_S0_R1_001.fastq.gz Undetermined_S1_R1_001.fastq.gz
    mv Undetermined_S0_I1_001.fastq.gz Undetermined_S1_R2_001.fastq.gz
    mv Undetermined_S0_I2_001.fastq.gz Undetermined_S1_R3_001.fastq.gz
    mv Undetermined_S0_R2_001.fastq.gz Undetermined_S1_R4_001.fastq.gz
    
    if [ -f $dir/fastqs/filesi1.xls ]; then rm $dir/fastqs/filler*; fi
    
    rawdir=$rawdir/fastqs/
    Start=Fastq_Merge
fi

if [ "$Start" = Fastq_SplitLane ] || [ "$Start" = Fastq_Merge ]; then
    echo "Skip bcltofastq"
    if [ -f $dir/fastqs/filesi1.xls ]; then rm $dir/fastqs/filler*; fi
    cd $rawdir
    if [ "$Start" = Fastq_SplitLane ]; then
	temp=$(ls 1_*.1.1.fastq.gz)
	Run=$(echo $temp | sed -e 's/^1_//' | sed -e 's/.1.1.fastq.gz//')
	echo "Run number is:" $Run
	echo "Spliting to 3 million reads files"
	mkdir $dir/smallfastqs/
	dosplit(){
            fastp -i $3/$1_$4.$1.barcode_1.fastq.gz -o $2/smallfastqs/$4.$1.I1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 
	    fastp -i $3/$1_$4.$1.barcode_2.fastq.gz -o $2/smallfastqs/$4.$1.I2.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 
	    fastp -i $3/$1_$4.$1.1.fastq.gz -o $2/smallfastqs/$4.$1.R1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q
	    fastp -i $3/$1_$4.$1.2.fastq.gz -o $2/smallfastqs/$4.$1.R2.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 
	}
	export -f dosplit
	
	dosplitQC(){
            zcat $3/$1_$4.$1.1.fastq.gz | head -n 12100000 | gzip > $2/temp1.$1.fastq.gz &
	    zcat $3/$1_$4.$1.2.fastq.gz | head -n 12100000 | gzip > $2/temp2.$1.fastq.gz &
	    zcat $3/$1_$4.$1.barcode_1.fastq.gz | head -n 12100000 | gzip > $2/temp3.$1.fastq.gz &
	    zcat $3/$1_$4.$1.barcode_2.fastq.gz | head -n 12100000 | gzip > $2/temp4.$1.fastq.gz &
	    wait
	    
	    fastp -i $2/temp1.$1.fastq.gz -o $2/smallfastqs/$4.$1.R1.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 
	    fastp -i $2/temp2.$1.fastq.gz -o $2/smallfastqs/$4.$1.R2.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q
	    fastp -i $2/temp3.$1.fastq.gz -o $2/smallfastqs/$4.$1.I1.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 
	    fastp -i $2/temp4.$1.fastq.gz -o $2/smallfastqs/$4.$1.I2.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 
	}
	export -f dosplitQC
	if [ -f $dir/smallfastqs/0001.$Run.1.I1.fastq.gz ]; then
            echo "Found 0001.Undetermined.1.I1.fastq, skip split fastqs"
	else
            if [ "$Runtype" = QC ]; then
		if [ "$Sequencer" = Novaseq ]; then
                    parallel --delay 1 --jobs 2 dosplitQC {} $dir $rawdir $Run ::: 1 &
		    parallel --delay 1 --jobs 2 dosplitQC {} $dir $rawdir $Run ::: 2 &
		    wait
		else
		    parallel --delay 1 --jobs 4 dosplitQC {} $dir $rawdir $Run ::: 1 &
		    parallel --delay 1 --jobs 4 dosplitQC {} $dir $rawdir $Run ::: 2 &
		    parallel --delay 1 --jobs 4 dosplitQC {} $dir $rawdir $Run ::: 3 &
		    parallel --delay 1 --jobs 4 dosplitQC {} $dir $rawdir $Run ::: 4 &
		    
		    wait
		fi

            elif [ "$Runtype" = full ]; then
		if [ "$Sequencer" = Novaseq ]; then
		    parallel --delay 1 --jobs 2 dosplit {} $dir $rawdir $Run ::: 1 &
		    parallel --delay 1 --jobs 2 dosplit {} $dir $rawdir $Run ::: 2 &
		    wait
		else
                    parallel --delay 1 --jobs 1 dosplit {} $dir $rawdir $Run ::: 1 &
		    parallel --delay 1 --jobs 1 dosplit {} $dir $rawdir $Run ::: 2 &
		    parallel --delay 1 --jobs 1 dosplit {} $dir $rawdir $Run ::: 3 &
		    parallel --delay 1 --jobs 1 dosplit {} $dir $rawdir $Run ::: 4 &
		    wait
		fi
            else
		echo "Unknown sequencer type, exiting" && exit;
            fi
	fi
    elif [ "$Start" = Fastq_Merge ]; then
        temp=$(ls *S1_R1_001.fastq.gz)
	Run=$(echo $temp | sed -e 's/\_\S1\_\R1\_\001.fastq.gz//')
        echo "Run number is:" $Run
        echo "Spliting to 3 million reads files"
	mkdir $dir/smallfastqs/
        if [ -f $dir/smallfastqs/0001.$Run.1.I1.fastq.gz ]; then
            echo "Found 0001.$Run.1.I1.fastq, skip split fastqs"
        else
            if [ "$Runtype" = QC ]; then
		if [ ! -f $dir/temp4.fastq.gz ]; then
		    echo "Runing QC pipeline"
		    zcat $rawdir/"$Run"_S1_R1_001.fastq.gz | head -n 12100000 | sed 's/1:N:0:1/1:N:0:/g' | sed 's/1:N:0:0/1:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp1.fastq.gz &
		    zcat $rawdir/"$Run"_S1_R4_001.fastq.gz | head -n 12100000 | sed 's/4:N:0:1/2:N:0:/g' | sed 's/2:N:0:0/2:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp2.fastq.gz &
		    zcat $rawdir/"$Run"_S1_R2_001.fastq.gz | head -n 12100000 | sed 's/2:N:0:1/1:N:0:/g' | sed 's/1:N:0:0/1:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp3.fastq.gz &
		    zcat $rawdir/"$Run"_S1_R3_001.fastq.gz | head -n 12100000 | sed 's/3:N:0:1/2:N:0:/g' | sed 's/2:N:0:0/2:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp4.fastq.gz &
		    wait
		    
		    fastp -i $dir/temp1.fastq.gz -o $dir/smallfastqs/$Run.1.R1.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
		    fastp -i $dir/temp2.fastq.gz -o $dir/smallfastqs/$Run.1.R2.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
		    fastp -i $dir/temp3.fastq.gz -o $dir/smallfastqs/$Run.1.I1.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
		    fastp -i $dir/temp4.fastq.gz -o $dir/smallfastqs/$Run.1.I2.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
		    wait
		    rm $dir/temp*.fastq.gz
		    touch $dir/temp4.fastq.gz
		fi
            elif [ "$Runtype" = full ]; then
		echo "Runing full pipeline"
		if [ ! -f $dir/temp4.fastq.gz ]; then
                    echo "Modify fastqs"
		    zcat $rawdir/"$Run"_S1_R1_001.fastq.gz | sed 's/1:N:0:1/1:N:0:/g' | sed 's/1:N:0:2/1:N:0:/g' | sed 's/1:N:0:0/1:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp1.fastq.gz &
		    zcat $rawdir/"$Run"_S1_R4_001.fastq.gz | sed 's/4:N:0:1/2:N:0:/g' | sed 's/4:N:0:2/2:N:0:/g' | sed 's/2:N:0:0/2:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp2.fastq.gz &
		    zcat $rawdir/"$Run"_S1_R2_001.fastq.gz | sed 's/2:N:0:1/1:N:0:/g' | sed 's/2:N:0:2/1:N:0:/g' | sed 's/1:N:0:0/1:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp3.fastq.gz &
		    zcat $rawdir/"$Run"_S1_R3_001.fastq.gz | sed 's/3:N:0:1/2:N:0:/g' | sed 's/3:N:0:2/2:N:0:/g' | sed 's/2:N:0:0/2:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp4.fastq.gz &
		    wait
		fi
		echo "Split fastqs to small files"
		fastp -i $dir/temp1.fastq.gz -o $dir/smallfastqs/$Run.1.R1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
		fastp -i $dir/temp2.fastq.gz -o $dir/smallfastqs/$Run.1.R2.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
		fastp -i $dir/temp3.fastq.gz -o $dir/smallfastqs/$Run.1.I1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
		fastp -i $dir/temp4.fastq.gz -o $dir/smallfastqs/$Run.1.I2.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
		wait
		rm $dir/temp*.fastq.gz
		touch $dir/temp4.fastq.gz
	    else
                echo "Unknown sequencer type, exiting" && exit;
            fi
        fi
    fi
    
    if [ -f $dir/fastp.json ]; then rm $dir/fastp.json $dir/fastp.html; fi
    ls $dir/smallfastqs | grep R1 > $dir/filesr1.xls
    ls $dir/smallfastqs | grep R2 > $dir/filesr2.xls
    ls $dir/smallfastqs | grep I1 > $dir/filesi1.xls
    ls $dir/smallfastqs | grep I2 > $dir/filesi2.xls
    
    cd $dir/
    if [ "$Indexed" = F ]; then
	mkdir $dir/Indexed/
	if [ -f $dir/Indexed/Sub.0001.$Run.1.R1.fastq.gz ]; then
            echo "Found Sub.0001.$Run.1.R1.fastq.R1.fastq, skip adding index"
	else
            echo "Adding index to fastqs"
            paste filesr1.xls filesr2.xls filesi1.xls filesi2.xls | awk -v OFS='\t' '{print $1, $2, $3, $4}'> Filelist2.xls
	    parallel --jobs 16 --colsep '\t' ' if [ -f '$dir'/Indexed/Sub.{1}.R1.fastq.gz ]; then echo "found '$dir'/Indexed/Sub.{1}.R1.fastq.gz"; \
                  else      '$pythohPATH' '$myPATH'/updateIndex_Next_core_V2.py -R1 '$dir'/smallfastqs/{1} -R2 '$dir'/smallfastqs/{2} \
                  --out '$dir'/Indexed/Sub.{1} -Index1 '$dir'/smallfastqs/{3} -Index2 '$dir'/smallfastqs/{4}; fi 2>>'$dir'/Run.log' :::: Filelist2.xls
	fi
	cd $dir/Indexed
	ls *R1.fastq.gz.R1.fastq.gz | while read -r i; do temp=$(basename -s .R1.fastq.gz.R1.fastq.gz $i) && mv $i $temp.R1.fastq.gz; done
	ls *R1.fastq.gz.R2.fastq.gz | while read -r i; do temp=$(basename -s .R1.fastq.gz.R2.fastq.gz $i) && mv $i $temp.R2.fastq.gz; done
	rawdir=$dir/Indexed/
    else
	mkdir $dir/Indexed/
	cd $dir/Indexed
	if [ -f $dir/Indexed/Sub.0001.$Run.1.R1.fastq.gz ]; then
            echo "Found Sub.0001.$Run.1.R1.fastq.R1.fastq, skip adding index"
        else
            echo "Adding index to fastqs"
	    mv $dir/smallfastqs/*fastq.gz  $dir/Indexed/
	    touch $dir/smallfastqs/0001.$Run.1.I1.fastq.gz
	    ls *fastq.gz | while read -r i; do  mv $i Sub.$i; done
            mv $dir/smallfastqs/*fastq.gz  $dir/Indexed/
            rawdir=$dir/Indexed/
	fi
    fi
fi 

# trim
cd $dir/fastqs/
if [ "$Start" = Fastq_core ]; then
    rawdir=$dir/Indexed/
fi
ls $rawdir | grep R1.fastq.gz > $rawdir/Read1.xls
ls $rawdir | grep R2.fastq.gz > $rawdir/Read2.xls

count=`ls -1 $dir/fastqs/Sub*.trim.fastq.gz 2>/dev/null | wc -l`
if [ $count != 0 ]; then
    echo "Found trimmed files, skip trimming"
else
    echo "Trim fastqs"
    paste $rawdir/Read1.xls $rawdir/Read2.xls | awk -v OFS='\t' '{print $1, $2}'> $rawdir/Filelist.xls
    cd $rawdir
    parallel --will-cite --jobs 16 --colsep '\t' ''$pythohPATH' '$myPATH'/pyadapter_trim.py -a {1} -b {2} >>'$dir'/Run.log' :::: $rawdir/Filelist.xls
    mv $rawdir/*trim.fastq.gz $dir/fastqs/
fi

# split projects and align
index=0
for Name in ${Project[@]}; do
    echo "project $index : $Name" 
    if [ -d $dir/$Name ]; then
	echo "Found $Name dir, skip this project"
    else    
	# split fastas to projects
	if [ -f $dir/fastqs/$Name.R1.trim.fastq.gz ] || [ -f $dir/fastqs/$Name.R1.trim.fastq ]; then
	    echo "Found $dir/fastqs/$Name.R1.trim.fastq.gz, skip splitting project"
	else
	    echo "Spliting project"
	    cd $dir/fastqs/
	    ls Sub*R1.trim.fastq.gz  > Trimmed1.xls
            ls Sub*R2.trim.fastq.gz  > Trimmed2.xls
	    cat Trimmed1.xls | sed 's/.R1.trim.fastq.gz//g' | sed 's/'Sub'/'"$Name".Sub'/g'  > $Name.out.xls
	    paste Trimmed1.xls Trimmed2.xls $Name.out.xls | awk -v OFS='\t' '{print $1, $2, $3}'> $Name.Filelist.xls
	    rm $Name.out.xls
	    if [ "$Sequencer" = Novaseq ]; then
                parallel --jobs 12 --colsep '\t' 'if [ -f {3}.R1.fq.gz ]; then echo "Found {3}.R1.fq.gz"; \
                         else '$pythohPATH' '$myPATH'/splitProjects_V4_novaseq_V2.py -R1 {1} -R2 {2} -yaml '$ymal' -Project '$Name' -out {3}; fi' :::: $Name.Filelist.xls
            else
		parallel --jobs 12 --colsep '\t' 'if [ -f {3}.R1.fq.gz ]; then echo "Found {3}.R1.fq.gz"; \
                         else '$pythohPATH' '$myPATH'/splitProjects_V4_nextseq_V2.py -R1 {1} -R2 {2} -yaml '$ymal' -Project '$Name' -out {3}; fi' :::: $Name.Filelist.xls
	    fi
	    cat $Name.Sub*R1.fq.gz > $Name.R1.trim.fastq.gz 
            cat $Name.Sub*R2.fq.gz > $Name.R2.trim.fastq.gz 
	    
	    if [ ${Type[$index]} == "ATAC" ]; then
		zcat $dir/fastqs/$Name.R1.trim.fastq.gz | sed "s/ 1:N:0/_1:N:0/g" | gzip > temp1.fastq.gz && \
		    mv temp1.fastq.gz $dir/fastqs/$Name.R1.trim.fastq.gz &
		zcat $dir/fastqs/$Name.R2.trim.fastq.gz | sed "s/ 2:N:0/_2:N:0/g" | sed "s/ 4:N:0/_2:N:0/g" | gzip > temp2.fastq.gz && \
		    mv temp2.fastq.gz $dir/fastqs/$Name.R2.trim.fastq.gz &
		wait

            elif [ ${Type[$index]} == "RNA" ]; then	
		# add UMI to read tag for RNA-seq
		echo "Add UMI to read name"
		cd $dir/fastqs/
		ls $Name.Sub*R1.fq.gz > Trimmed3.xls
		ls $Name.Sub*R2.fq.gz > Trimmed4.xls
		# Match 6 poly T to remove low quality reads
		paste Trimmed3.xls Trimmed4.xls | awk -v OFS='\t' '{print $1, $2}'> $Name.Filelist.xls
		if [ $matchPolyT == "T" ]; then
		    parallel --jobs 8 --delay 1 --colsep '\t' 'fastp -i {1} -I {2} -o {1}.temp -O {2}.temp --thread 1 -A -G -Q -l 28 2>>'$dir'/Run.log && \
		    umi_tools extract --bc-pattern=NNNNNNNNNN --stdin {2}.temp --stdout {2}.extract --read2-in {1}.temp --read2-out {1}.extract >>'$dir'/Run.log && \
                    python '$myPATH'/matchPolyT.py -R1 {1}.extract -R2 {2}.extract --out {1} && \
                    rm {1}.extract {2}.extract {1} {2} {1}.temp {2}.temp' :::: $Name.Filelist.xls
		else
		    parallel --jobs 8 --delay 1 --colsep '\t' 'umi_tools extract --bc-pattern=NNNNNNNNNN --stdin {2} --stdout {2}.extract --read2-in {1} --read2-out {1}.extract >>'$dir'/Run.log && \
		    mv {1}.extract {1}.R1polyT && rm {2}.extract {1} {2}' :::: $Name.Filelist.xls
		fi
		mv $Name.R1.trim.fastq.gz $Name.R1.trimRaw.fastq.gz
		mv $Name.R2.trim.fastq.gz $Name.R2.trimRaw.fastq.gz
		echo "Merge trimmed fastqs"
		cat $Name.Sub*R1polyT | sed "s/ 1:N:0/_1:N:0/g" > $Name.R1.trim.fastq
		rm $Name*polyT
	    fi
	fi

	# align
	if [ ${Genomes[$index]} == "both" ]; then
	    if [ ${Type[$index]} == "ATAC" ] ; then
		Genome1=(hg19 mm10) # Genome1 for alignment, Genome2 for other steps
	    else
		Genome1=(both) # for RNA and Cellhash
	    fi
            Genome2=(hg19 mm10)
	else
	    Genome1=${Genomes[$index]}
	    Genome2=${Genomes[$index]}
	fi
	
	for Species in ${Genome1[@]}; do
	    if [ -d $dir/$Name ]; then
		echo "Found $Name folder, skip alignment"
	    else
		echo "Align $Name to $Species ${Type[$index]} library"
		cd $dir/fastqs/
		if [ -f $dir/fastqs/$Name.$Species.st.bam ]; then
		    echo "Found $Name.$Species.st.bam, skip alignment"
		else
		    if [ ${Type[$index]} == "ATAC" ]; then
			(bowtie2 -X2000 -p $cores --rg-id $Name \
				 -x $bowtieGenome/$Species/$Species \
				 -1 $dir/fastqs/$Name.R1.trim.fastq.gz \
				 -2 $dir/fastqs/$Name.R2.trim.fastq.gz | \
			     samtools view -bS -@ $cores - -o $Name.$Species.bam) 2>$Name.$Species.align.log
		    elif [ ${Type[$index]} == "RNA" ]; then
			STAR --chimOutType WithinBAM \
			    --runThreadN $cores \
			    --genomeDir $starGenome/$Species/ \
			    --readFilesIn $dir/fastqs/$Name.R1.trim.fastq  \
			    --outFileNamePrefix $dir/fastqs/$Name.$Species. \
			    --outFilterMultimapNmax 20 \
			    --outFilterMismatchNoverLmax 0.06 \
			    --limitOutSJcollapsed 2000000 \
			    --outSAMtype BAM Unsorted \
			    --limitIObufferSize 400000000 \
			    --outReadsUnmapped Fastx
			# modify chrMT to chrM
			if [ $refgene == "gencode" ]; then
			    samtools view -h $dir/fastqs/$Name.$Species.Aligned.out.bam | sed 's/chrMT/chrM/g' | samtools view -bS > $dir/fastqs/$Name.$Species.bam
			else
			    mv $dir/fastqs/$Name.$Species.Aligned.out.bam $dir/fastqs/$Name.$Species.bam
			fi
			rm -r *_STARtmp *Log.progress.out *SJ.out.tab *Unmapped.out.mate1 *_STARpass1 *_STARgenome $dir/fastqs/$Name.$Species.Aligned.out.bam
			mv $dir/fastqs/$Name.$Species.Log.final.out $dir/fastqs/$Name.$Species.align.log  
		    else
			echo "Unknown parameter"
			exit
		    fi
                        echo "Sort $Name.$Species.bam"
		        cd $dir/fastqs/
                        java -Xmx24g -Djava.io.tmpdir=$dir/temp/ -jar $picardPATH SortSam SO=coordinate I=$Name.$Species.bam O=$Name.$Species.st.bam VALIDATION_STRINGENCY=SILENT TMP_DIR=$dir/temp/ 2>>$dir/Run.log
                        samtools index -@ $cores $Name.$Species.st.bam
                        rm $Name.$Species.bam
		fi
		
		# Update RGID's and Headers
		if [ -f $dir/fastqs/$Name.$Species.rigid.reheader.st.bam.bai ]; then
		    echo "Found $Name.$Species.rigid.reheader.st.bam, skip update RGID"
		elif [ ${Type[$index]} == "ATAC" ] || [ ${Type[$index]} == "RNA" ]; then
		    echo "Update RGID for $Name.$Species.st.bam"
		    if [ -f $Name.$Species.rigid.st.bam ]; then
			echo "Found $Name.$Species.rigid.st.bam, skip update RG tag"
		    else
			if [ "$Sequencer" = Novaseq ]; then
			    $pythohPATH $myPATH/updateRGID_Singles_novaseq_V4.py --bam $Name.$Species.st.bam --out $Name.$Species.rigid.st.bam --err $Name.$Species.discard.st.bam --libtype ${Type[$index]}
			else
			    $pythohPATH $myPATH/updateRGID_Singles_nextseq_V4.py --bam $Name.$Species.st.bam --out $Name.$Species.rigid.st.bam --err $Name.$Species.discard.st.bam --libtype ${Type[$index]}
			fi
		    fi
		    samtools view -H $Name.$Species.st.bam > $Name.$Species.st.header.sam
		    samtools view -@ $cores $Name.$Species.rigid.st.bam | cut -f1 | \
			sed 's/_/\t/g' | cut -f2 | sort --parallel=$cores -S 10G | \
			uniq | awk -v OFS='\t' '{print "@RG", "ID:"$1, "SM:Barcode"NR}' > header.temp.sam
		    sed -e '/\@RG/r./header.temp.sam' $Name.$Species.st.header.sam > $Name.$Species.rigid.st.header.sam
		    samtools reheader -P $Name.$Species.rigid.st.header.sam $Name.$Species.rigid.st.bam > $Name.$Species.rigid.reheader.st.bam
		    samtools index -@ $cores $Name.$Species.rigid.reheader.st.bam
		    rm *rigid.st*
		    rm $Name.$Species.st.header.sam header.temp.sam
		fi

		
		# remove dups and mito
		if [ -f $dir/fastqs/$Name.$Species.rigid.reheader.unique.st.bam ] || [ -f $dir/fastqs/$Name.$Species.wdup.all.bam ]; then
		    echo "Found $Name.$Species.rigid.reheader.nomito.st.bam, skip mito removal"
		else
		    if   [ ${Type[$index]} == "ATAC" ]; then
			echo "Remove low quality reads and unwanted chrs:" $Name.$Species.rigid.reheader.st.bam
			chrs=`samtools view -H $Name.$Species.st.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | grep -v Y | awk '{if(length($0)<6)print}'`
			samtools view -@ $cores -b -q 30 -f 0x2 $Name.$Species.rigid.reheader.st.bam -o $Name.$Species.wdup.all.bam `echo $chrs`
			samtools index -@ $cores $Name.$Species.wdup.all.bam
		    elif [ ${Type[$index]} == "RNA" ]; then
			echo "Remove low quality reads and unwanted chrs:" $Name.$Species.rigid.reheader.st.bam
			if [ ${keepMultiMapping[$index]} == "T" ]; then
			    samtools view -@ $cores -b $Name.$Species.rigid.reheader.st.bam -o $Name.$Species.rigid.reheader.unique.st.bam            
			else    
			    samtools view -@ $cores -b -q 30 $Name.$Species.rigid.reheader.st.bam -o $Name.$Species.rigid.reheader.unique.st.bam
			fi
			samtools index -@ $cores $Name.$Species.rigid.reheader.unique.st.bam
		    fi
		fi
		
	       
		
		if [ ${Type[$index]} == "ATAC" ]; then
		    # count unfiltered reads
		    if [ ! -f $Name.$Species.unfiltered.counts.csv ]; then
			echo "Count unfiltered reads" $Name.$Species.wdup.all.bam
			if [ ! -f $Name.$Species.wdup.RG.bed ]; then
			    samtools view -@ $cores $Name.$Species.wdup.all.bam | cut -f1 |  sed 's/_/\t/g' | cut -f2 > $Name.$Species.wdup.RG.bed
			fi
			if [ ! -f $Name.$Species.wdup.RG.freq.bed ]; then
			    cat $Name.$Species.wdup.RG.bed | sort --parallel=$cores -S 10G | uniq -c > $Name.$Species.wdup.RG.freq.bed
			fi
			Rscript $myPATH/sum_reads.R $dir/fastqs/ $Name.$Species.wdup.RG.freq.bed --save
			mv $Name.$Species.wdup.RG.freq.bed.csv $Name.$Species.unfiltered.counts.csv
			rm $Name.$Species.wdup.RG.bed $Name.$Species.wdup.RG.freq.bed
		    fi
		    
		    # remove barcode combination that has less then N unfiltered reads
		    if [ -f $dir/fastqs/$Name.$Species.wdup.bam.bai ]; then
			echo "Skip removing barcode combinations"
		    else
			echo "Remove barcode combinations that have less then ${RawReadsPerBarcode[$index]} reads"
			sed -e 's/,/\t/g' $Name.$Species.unfiltered.counts.csv | awk -v OFS=',' 'NR>=2 {if($5 >= '${RawReadsPerBarcode[$index]}') print $1,$2,$3,$4}'  > $Name.$Species.barcodes.txt
			# too many barcode combinations will result in picard failure
			samtools view -@ $cores -R $Name.$Species.barcodes.txt $Name.$Species.wdup.all.bam -o $Name.$Species.wdup.bam
			samtools index -@ $cores $Name.$Species.wdup.bam	
		    fi
		    # remove duplicates
		    if [ -f $dir/fastqs/$Name.$Species.rmdup.bam ]; then
			echo "Skip split $Name.$Species.wdup.bam and mark duplicate"
		    else
			if [ ${Type[$index]} == "ATAC" ]; then
			    # split bam and remove unused barcode in the bam header, which will spead up dedupliaction a lot!
			    if [ -f SubBarcodes.$Name.$Species.0000.bam ]; then
				echo "Found SubBarcodes.$Name.$Species.0001.bam, skip splitting bam file"
			    else
				echo "Split on bam"
				samtools view -@ $cores -H $Name.$Species.wdup.bam | grep -e @HD -e @SQ -e @PG > $Name.$Species.wdup.header.sam
				split -l 1000 -d -a 4 $Name.$Species.barcodes.txt SubBarcodes.$Name.$Species.
				ls SubBarcodes.$Name.$Species.* | while read -r i; do
				    samtools view -@ $cores -R $i $Name.$Species.wdup.bam -o $i.bam &&
					cat $i | awk -v OFS='\t' '{print "@RG", "ID:"$1, "SM:Barcode"NR}' > header.temp.sam &&
					cat $Name.$Species.wdup.header.sam header.temp.sam > header.temp2.sam &&
					mv header.temp2.sam header.temp.sam &&
					samtools reheader header.temp.sam $i.bam > $i.bam2 
				    mv $i.bam2 $i.bam; done
			    fi
			    echo "Mark duplicates: $Name.$Species.wdup.bam"
			    parallel --will-cite --jobs 2 --delay 0.5 'if [ -f {.}.rmdup ]; then echo \"Found {.}.rmdup\"; else java -Xmx4g -Djava.io.tmpdir='$dir'/temp/ -jar '$picardPATH' MarkDuplicates INPUT={} OUTPUT={.}.rmdup METRICS_FILE={.}.dups.log REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT BARCODE_TAG=RG TMP_DIR='$dir'/temp/ 2>>'$dir'/Run.log && echo \"Marked duplicates for {}\"; fi' ::: SubBarcodes.$Name.$Species.*.bam
			
			    samtools merge -@ 6 -c $Name.$Species.rmdup.bam SubBarcodes.$Name.$Species.*.rmdup
			    cp SubBarcodes.$Name.$Species.0000.dups.log $Name.$Species.dups.log
			    samtools index -@ $cores $Name.$Species.rmdup.bam
			    rm SubBarcodes.$Name.$Species.* $Name.$Species.wdup.header.sam header.temp.sam
			fi
		    fi
		fi
		
		#  split into mouse and human
		if [ ${Type[$index]} == "RNA" ] && [ ! -f $Name.mm10.rigid.reheader.unique.st.bam ]; then
		    cd $dir/fastqs/
		    if [ ${Genomes[$index]} == "both" ]; then
			echo "Split into hg and mm"
			chrs1=`samtools view -H $Name.$Species.rigid.reheader.unique.st.bam | grep hg19 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
			chrs2=`samtools view -H $Name.$Species.rigid.reheader.unique.st.bam | grep mm10 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
			samtools view -@ $cores -b $Name.$Species.rigid.reheader.unique.st.bam -o temp1.bam `echo ${chrs1[@]}`
			samtools view -@ $cores -b $Name.$Species.rigid.reheader.unique.st.bam -o temp2.bam `echo ${chrs2[@]}`
			samtools view -@ $cores -h temp1.bam | sed 's/hg19_/chr/g' | sed 's/chrMT/chrM/g'| samtools view -@ $cores -b -o $Name.hg19.rigid.reheader.unique.st.bam
			samtools view -@ $cores -h temp2.bam | sed 's/mm10_/chr/g' | sed 's/chrMT/chrM/g'| samtools view -@ $cores -b -o $Name.mm10.rigid.reheader.unique.st.bam
			samtools index -@ $cores $Name.hg19.rigid.reheader.unique.st.bam &
			samtools index -@ $cores $Name.mm10.rigid.reheader.unique.st.bam &
			wait
			rm temp1.bam temp2.bam
		    else
			chrs=`samtools view -H $Name.$Species.rigid.reheader.unique.st.bam | grep chr | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
			samtools view -@ $cores -b $Name.$Species.rigid.reheader.unique.st.bam -o temp1.bam `echo ${chrs[@]}` 
			mv temp1.bam $Name.$Species.rigid.reheader.unique.st.bam
		    fi		
		fi

		

		# assign feature to reads
		cd  $dir/fastqs/
		if [ ${Type[$index]} == "RNA" ]; then
		    for Species2 in ${Genome2[@]}; do
			if [ -f $dir/fastqs/$Name.$Species2.exon.wdup.bam.bai ]; then
			    echo "Skip exon feasure count"
			else
		    	    # excliude multimapping, uniquely mapped reads only, -Q 30, for real sample, might consider include multi-mapping
			    echo "Feature counting on exons"
			    # count exon
                            if [ ${keepMultiMapping[$index]} == "T" ]; then
				featureCounts -T $cores -Q 0 -M -a $myPATH/gtf/$Species2.$refgene.gtf -t exon -g $genename \
                                    -o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.rigid.reheader.unique.st.bam >>$dir/Run.log
                            else
				featureCounts -T $cores -Q 30 -a $myPATH/gtf/$Species2.$refgene.gtf -t exon -g $genename \
                                    -o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.rigid.reheader.unique.st.bam >>$dir/Run.log
                            fi
			    # Extract reads that assigned to genes
			    mv $Name.$Species2.rigid.reheader.unique.st.bam.featureCounts.bam $Name.$Species2.exon.featureCounts.bam
			    samtools view -H $Name.$Species2.exon.featureCounts.bam > $Name.$Species2.header.sam
			    samtools view -@ $cores $Name.$Species2.exon.featureCounts.bam | grep XT > $Name.$Species2.temp.sam
			    cat $Name.$Species2.header.sam $Name.$Species2.temp.sam | samtools view -b -@ $cores -o $Name.$Species2.exon.wdup.bam
			    samtools sort -@ $cores -o $Name.$Species2.exon.wdup.st.bam $Name.$Species2.exon.wdup.bam 
			    mv $Name.$Species2.exon.wdup.st.bam $Name.$Species2.exon.wdup.bam 
			    samtools index -@ $cores $Name.$Species2.exon.wdup.bam
			    rm $Name.$Species2.header.sam $Name.$Species2.temp.sam
			fi
			if [ -f $dir/fastqs/$Name.$Species2.gene.wdup.bam.bai ]; then
			    echo "Skip intron and exon feasure count"
			else
			    # count both intron and exon
			    echo "Count feature on both intron and exon"
			    if [ $keepIntron == "T" ]; then
				if [ ${keepMultiMapping[$index]} == "T" ]; then
				    featureCounts -T $cores -Q 0 -M -a $myPATH/gtf/$Species2.$refgene.gtf -t gene -g $genename \
					-o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.exon.featureCounts.bam >>$dir/Run.log
				else
				    featureCounts -T $cores -Q 30 -a $myPATH/gtf/$Species2.$refgene.gtf -t gene -g $genename \
					-o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.exon.featureCounts.bam >>$dir/Run.log
				fi
				mv $Name.$Species2.exon.featureCounts.bam.featureCounts.bam $Name.$Species2.gene.featureCounts.bam
			    else
				cp $Name.$Species2.exon.featureCounts.bam $Name.$Species2.gene.featureCounts.bam
			    fi

			    # Extract reads that assigned to genes
			    samtools view -H $Name.$Species2.gene.featureCounts.bam > $Name.$Species2.header.sam
			    samtools view -@ $cores $Name.$Species2.gene.featureCounts.bam | grep XT > $Name.$Species2.temp.sam
			    cat $Name.$Species2.header.sam $Name.$Species2.temp.sam | \
				samtools view -b -@ $cores -o $Name.$Species2.gene.wdup.bam
			    samtools sort -@ $cores -o $Name.$Species2.wdup.st.bam $Name.$Species2.gene.wdup.bam
			    mv $Name.$Species2.wdup.st.bam $Name.$Species2.gene.wdup.bam
			    samtools index -@ $cores $Name.$Species2.gene.wdup.bam
			    rm $Name.$Species2.header.sam $Name.$Species2.temp.sam
			fi
		    done
		fi

                cd  $dir/fastqs/
		# count UMIs
		if [ ${Type[$index]} == "RNA" ]; then
		    for Species2 in ${Genome2[@]}; do
			if [ -f $dir/fastqs/$Name.$Species2.UMIcounts.csv ]; then
			    echo "Found $Name.$Species.UMIcounts.csv, skip counting UMIs"
			else
			    echo "Count UMIs and generate $Name.$Species.UMIcounts.csv"
			
			if [ ${Type[$index]} == "RNA" ] && [ ! -f $dir/fastqs/$Name.$Species2.assigned.exon.st.bam.bai ]; then
			    samtools sort -@ $cores $dir/fastqs/$Name.$Species2.gene.featureCounts.bam -o $dir/fastqs/$Name.$Species2.assigned.st.bam
			    samtools sort -@ $cores $dir/fastqs/$Name.$Species2.exon.featureCounts.bam -o $dir/fastqs/$Name.$Species2.assigned.exon.st.bam
			    rm $dir/fastqs/$Name.$Species2.gene.featureCounts.bam $dir/fastqs/$Name.$Species2.exon.featureCounts.bam
			    
			    samtools index -@ $cores $dir/fastqs/$Name.$Species2.assigned.st.bam &
			    samtools index -@ $cores $dir/fastqs/$Name.$Species2.assigned.exon.st.bam &
			    wait
			fi
			# Counting molecules
			echo "Count molecules in " $Name.$Species2.gene.wdup.bam
			echo "Group reads to unique UMIs"
			if [ -f $dir/fastqs/$Name.$Species2.groups.tsv ]; then
			    echo "Skip grouping UMIs"
			else
			    umi_tools group --extract-umi-method=read_id \
				--per-gene --gene-tag=XT --per-cell \
				-I $dir/fastqs/$Name.$Species2.gene.wdup.bam  \
				--output-bam -S $dir/fastqs/$Name.$Species2.grouped.bam \
				--group-out=$dir/fastqs/$Name.$Species2.groups.tsv --skip-tags-regex=Unassigned >>$dir/Run.log
			fi
			# filter UMI that has only 1 read
			if [ -f $dir/fastqs/$Name.$Species2.wdup.all.bam.bai ]; then
			    echo "Skip removing single-read UMIs and GGGG"
			else
			    if [ $SkipPolyGumi == "F" ]; then
				
				echo "Filter UMI that has only 1 read and GGGGGGGGGG read"
				# note: if there is severe tilt of species mixing, use $8>1
				if [ $removeSingelReadUMI == "T" ]; then
				    less $dir/fastqs/$Name.$Species2.groups.tsv | awk '$8 > 1' | \
					awk '$5 != "GGGGGGGGGG" {print $1}' > $dir/fastqs/$Name.$Species2.UMIsave.txt                
				else
				    less $dir/fastqs/$Name.$Species2.groups.tsv | awk '$8 > 0' | \
					awk '$5 != "GGGGGGGGGG" {print $1}' > $dir/fastqs/$Name.$Species2.UMIsave.txt
				fi
				samtools view -@ $cores $dir/fastqs/$Name.$Species2.grouped.bam > $dir/fastqs/$Name.$Species2.temp.sam
				split -l 1000000 -d -a 4 $dir/fastqs/$Name.$Species2.UMIsave.txt $dir/fastqs/SubUMI.$Name.$Species2.
				ls $dir/fastqs/SubUMI.$Name.$Species2.* | \
				    while read -r i; do fgrep -w -f $i $dir/fastqs/$Name.$Species2.temp.sam >> $dir/fastqs/$Name.$Species2.temp2.sam; done
				
				cat <(samtools view -H $dir/fastqs/$Name.$Species2.grouped.bam) $dir/fastqs/$Name.$Species2.temp2.sam | \
				    samtools view -b -@ $cores | samtools sort -@ $cores -o $dir/fastqs/$Name.$Species2.grouped.fil.bam
				rm $dir/fastqs/$Name.$Species2.temp.sam  $dir/fastqs/$Name.$Species2.temp2.sam  $dir/fastqs/SubUMI.$Name.$Species2.*
				mv $dir/fastqs/$Name.$Species2.grouped.fil.bam $dir/fastqs/$Name.$Species2.wdup.all.bam
				samtools index -@ $cores $dir/fastqs/$Name.$Species2.wdup.all.bam
			    elif [ $SkipPolyGumi == "T" ]; then
				cp $dir/fastqs/$Name.$Species2.grouped.bam $dir/fastqs/$Name.$Species2.wdup.all.bam
				samtools index -@ $cores $dir/fastqs/$Name.$Species2.wdup.all.bam
			    fi
			fi
			# count unfiltered reads
			if [ ! -f $Name.$Species2.unfiltered.counts.csv ]; then
			    echo "Count unfiltered reads" $Name.$Species2.wdup.all.bam
			    if [ ! -f $Name.$Species2.wdup.RG.bed ]; then
				samtools view -@ $cores $Name.$Species2.wdup.all.bam | cut -f1 | sed 's/_/\t/g' | cut -f2 > $Name.$Species2.wdup.RG.bed
			    fi
			    if [ ! -f $Name.$Species2.wdup.RG.freq.bed ]; then
				cat $Name.$Species2.wdup.RG.bed | sort --parallel=$cores -S 10G | uniq -c > $Name.$Species2.wdup.RG.freq.bed
			    fi
			    Rscript $myPATH/sum_reads.R $dir/fastqs/ $Name.$Species2.wdup.RG.freq.bed --save
			    mv $Name.$Species2.wdup.RG.freq.bed.csv $Name.$Species2.unfiltered.counts.csv
			    rm $Name.$Species2.wdup.RG.bed $Name.$Species2.wdup.RG.freq.bed
			fi
			
			# remove barcode combination that has less then N reads
			if [ -f $dir/fastqs/$Name.$Species2.wdup.bam.bai ]; then
			    echo "Skip removing barcode combinations that have less then ${RawReadsPerBarcode[$index]} reads"
			else
			    echo "Remove barcode combinations that have less then ${RawReadsPerBarcode[$index]} reads"
			    sed -e 's/,/\t/g' $Name.$Species2.unfiltered.counts.csv | \
				awk -v OFS=',' 'NR>=2 {if($5 >= '${RawReadsPerBarcode[$index]}') print $1,$2,$3,$4}'  > $Name.$Species2.barcodes.txt			    
			    samtools view -@ $cores -R $Name.$Species2.barcodes.txt $Name.$Species2.wdup.all.bam -o $Name.$Species2.wdup.bam
			    samtools index -@ $cores $Name.$Species2.wdup.bam
			    # rm $Name.$Species2.barcodes.txt
			fi
			# deduplication of umi
			if [ -f $dir/fastqs/$Name.$Species2.rmdup.bam.bai ]; then
			    echo "Skip removing duplicates"
			else
			    cd $dir/fastqs/
			    echo "Remove duplicates"
			    if [ -f SubBarcodes.$Name.$Species2.0000.bam ]; then
				echo "Found SubBarcodes.$Name.$Species2.0000.bam, skip splitting bam file"
                            else
				echo "Split on bam"
				samtools view -@ $cores -H $Name.$Species2.wdup.bam | grep -e @HD -e @SQ -e @PG > $Name.$Species2.wdup.header.sam
				split -l 1000 -d -a 4 $Name.$Species2.barcodes.txt SubBarcodes.$Name.$Species2.
				ls SubBarcodes.$Name.$Species2.* | while read -r i; do samtools view -@ $cores -R $i $Name.$Species2.wdup.bam -o $i.bam &&
				cat $i | awk -v OFS='\t' '{print "@RG", "ID:"$1, "SM:Barcode"NR}' > header.temp.sam &&
                                cat $Name.$Species2.wdup.header.sam header.temp.sam > header.temp2.sam &&
                                mv header.temp2.sam header.temp.sam &&
                                samtools reheader header.temp.sam $i.bam > $i.bam2
				mv $i.bam2 $i.bam
				samtools index -@ $cores $i.bam; done
                            fi
                            echo "Mark duplicates: $Name.$Species2.wdup.bam"
			    parallel --will-cite --jobs 6 --delay 1 'if [ -f {.}.rmdup ]; then echo \"Found {.}.rmdup\"; else umi_tools dedup --extract-umi-method=read_id --per-gene --gene-tag=XT --per-cell -I {} --output-stats='$dir'/fastqs/'$Name'.'$Species2'.dedup -S {.}.rmdup --skip-tags-regex=Unassigned >>'$dir'/Run.log && echo \"Marked duplicates for {}\"; fi' ::: SubBarcodes.$Name.$Species2.*.bam
			    
                            samtools merge -@ 6 -c $Name.$Species2.rmdup.bam SubBarcodes.$Name.$Species2.*.rmdup
                            samtools index -@ $cores $dir/fastqs/$Name.$Species2.rmdup.bam
                            rm SubBarcodes.$Name.$Species2.* $Name.$Species2.wdup.header.sam header.temp.sam	    
			fi
			fi
		    done
		fi

		## calculate read distribution
		if [ ${Type[$index]} == "RNA" ]; then
                    # split bam to hg and mm
                    cd $dir/fastqs/
                    for Species2 in ${Genome2[@]}; do
                         if [ -f $dir/fastqs/$Name.$Species2.read_distribution.txt ]; then
                            echo "Skip calculate read disbution"
                        else
                            echo "Calculate read distribution"
                            if [ $Runtype = QC ]; then
                                read_distribution.py -i $dir/fastqs/$Name.$Species2.assigned.st.bam -r $myPATH/$Species2.UCSC_RefSeq.bed > $Name.$Species2.read_distribution.txt 2>>$dir/Run.log
                            else
                                # only use 1% of reads
                                samtools view -s 0.01 -o $dir/fastqs/temp.bam $dir/fastqs/$Name.$Species2.assigned.st.bam
                                read_distribution.py -i $dir/fastqs/temp.bam -r $myPATH/$Species2.UCSC_RefSeq.bed > $Name.$Species2.read_distribution.txt 2>>$dir/Run.log
                                rm $dir/fastqs/temp.bam
                            fi
                            # plot reads disbution
                            tail -n +5 $Name.$Species2.read_distribution.txt | head -n -1 > temp1.txt
                            head -n 3  $Name.$Species2.read_distribution.txt | grep Assigned | sed 's/Total Assigned Tags//g' | sed 's/ //g' > temp2.txt
                            Rscript $myPATH/Read_distribution.R $dir/fastqs/ $Name.$Species2 --save
                            rm temp1.txt temp2.txt
                        fi
                    done
		fi
		
		if [ ${Type[$index]} == "ATAC" ]; then
		    # get final quality stats
		    echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Raw" > $Name.$Species.stats.log
		    samtools idxstats $Name.$Species.st.bam >> $Name.$Species.stats.log
		    echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Filtered" >> $Name.$Species.stats.log
		    samtools idxstats $Name.$Species.rmdup.bam >> $Name.$Species.stats.log
		    
		    # get insert-sizes
		    if [ -f $Name.$Species.rmdup.hist_data.pdf ]; then
			echo "Skip checking insert size"
		    else
			echo '' > $Name.$Species.rmdup.hist_data.log
			java -jar $picardPATH CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT I=$Name.$Species.rmdup.bam O=$Name.$Species.rmdup.hist_data.log H=$Name.$Species.rmdup.hist_data.pdf W=1000  2>>$dir/Run.log
		    fi
		fi
		
		# make TSS pileup fig
		if [ ${Type[$index]} == "ATAC" ] && [ ! -f $Name.$Species.RefSeqTSS ]; then
		    echo "Create TSS pileup"; set +e
		    $toolPATH/pyMakeVplot.py -a $Name.$Species.rmdup.bam -b $tssFilesPATH/$Species.TSS.bed -e 2000 -p ends -v -u -o $Name.$Species.RefSeqTSS
		fi

						
       		if [ -d "$dir/fastqs/tmp" ]; then rm -r $dir/fastqs/tmp; fi
	    fi
	done
	
	# count reads for ATAC and RNA
	for Species2 in ${Genome2[@]}; do
	    if [ ${Type[$index]} == "ATAC" ] || [ ${Type[$index]} == "RNA" ]; then
		# count filtered reads
		if [ ! -f $Name.$Species2.filtered.counts.csv ]; then
		    echo "Count filtered reads" $Name.$Species2.rmdup.bam
		    if [ ! -f $Name.$Species2.rmdup.RG.bed ]; then
			samtools view -@ $cores $Name.$Species2.rmdup.bam | cut -f1 |sed 's/_/\t/g' | cut -f2 > $Name.$Species2.rmdup.RG.bed
		    fi
		    if [ ! -f $Name.$Species2.rmdup.RG.freq.bed ]; then
			cat $Name.$Species2.rmdup.RG.bed | sort --parallel=$cores -S 10G | uniq -c > $Name.$Species2.rmdup.RG.freq.bed
		    fi
		    Rscript $myPATH/sum_reads.R $dir/fastqs/ $Name.$Species2.rmdup.RG.freq.bed --save
		    mv $Name.$Species2.rmdup.RG.freq.bed.csv $Name.$Species2.filtered.counts.csv
		    rm $Name.$Species2.rmdup.RG.bed $Name.$Species2.rmdup.RG.freq.bed
		fi
		# remove barcode combinations with less then N reads
		if [ -f $Name.$Species2.rmdup.cutoff.bam.bai ]; then
		    echo "Skip removing low counts barcode combination"
		else
		    echo "Remove low counts barcode combination"
		    sed -e 's/,/\t/g' $Name.$Species2.filtered.counts.csv | \
			awk -v OFS=',' 'NR>=2 {if($5 >= '${ReadsPerBarcode[$index]}') print $1,$2,$3,$4} '  > $Name.$Species2.barcodes.txt
		    samtools view -@ $cores -R $Name.$Species2.barcodes.txt $Name.$Species2.rmdup.bam -o $Name.$Species2.rmdup.cutoff.bam
		    samtools index -@ $cores $Name.$Species2.rmdup.cutoff.bam
		fi
		if [ ${Type[$index]} == "ATAC" ]; then
		    if [ -f $Name.$Species2.rmdup.cutoff.namesort.bam ]; then
			echo "Skip sort bam on name"
		else
		    echo "Sort bam on name"
		    samtools sort -@ $cores -n -o $Name.$Species2.rmdup.cutoff.namesort.bam $Name.$Species2.rmdup.cutoff.bam
		    fi
		    # convert bam to bed
		    if [ -f $Name.$Species2.rmdup.cutoff.bed.gz ]; then
			echo "Skip converting rmdup.cutoff.bam to rmdup.cutoff.bed.gz"
		    else
			echo "Convert rmdup.cutoff.bam to rmdup.cutoff.bed.gz"
			bedtools bamtobed -i $Name.$Species2.rmdup.cutoff.namesort.bam -bedpe | \
			    sed 's/_/\t/g' | \
			    awk -v OFS="\t" '{if($10=="+"){print $1,$2+4,$6+4,$8}else if($10=="-"){print $1,$2-5,$6-5,$8}}' |\
			    gzip > $Name.$Species2.rmdup.cutoff.bed.gz
			rm $Name.$Species2.rmdup.cutoff.namesort.bam
		    fi
		fi
		
		if [ ${Type[$index]} == "RNA" ]; then
		    if [ -f $dir/fastqs/$Name.$Species2.UMIcounts.csv ]; then
			echo "Skip counting UMIs"
		    else
			echo "Count UMIs per cell"
			umi_tools count --extract-umi-method=read_id --per-gene --gene-tag=XT --per-cell \
				  -I $dir/fastqs/$Name.$Species2.rmdup.cutoff.bam --wide-format-cell-counts \
				  --skip-tags-regex=Unassigned \
				  -S $dir/fastqs/$Name.$Species2.UMIcounts.csv >>$dir/Run.log
		    fi
		fi
	    fi
	done



	# estimate lib size
	if [ -f $Name.counts.csv ]; then
	    echo "Found $Name.counts.csv, skip calculate lib size"
	else
	    echo "Estimate lib size"
            if [ ${Genomes[$index]} == "both" ]; then
		if [ -f $Name.hg19.unfiltered.counts.csv ] && [ ! -f $Name.hg19.filtered.counts.csv ]; then
		    cp $Name.hg19.unfiltered.counts.csv $Name.hg19.filtered.counts.csv
		    echo "Error: Could locate $Name.hg19.unfiltered.counts.csv"
		fi
		if [ -f $Name.mm10.unfiltered.counts.csv ] && [ ! -f $Name.mm10.filtered.counts.csv ]; then
		    cp $Name.mm10.unfiltered.counts.csv $Name.mm10.filtered.counts.csv
		    echo "Error: Could locate $Name.mm10.unfiltered.counts.csv"
		fi
		if [ -f $Name.hg19.unfiltered.counts.csv ] && [ -f $Name.mm10.unfiltered.counts.csv ]; then
                    echo "Calcuating library size for $Name"
		    Rscript $myPATH/lib_size_sc_V4_species_mixing.R ./ $Name ${ReadsPerBarcode[$index]} ${Type[$index]} --save
		fi
	    else
		if [ -f $Name.${Genomes[$index]}.unfiltered.counts.csv ] && [ ! -f $Name.${Genomes[$index]}.filtered.counts.csv ]; then
		    cp $Name.${Genomes[$index]}.unfiltered.counts.csv $Name.${Genomes[$index]}.filtered.counts.csv
		else
		    Rscript $myPATH/lib_size_sc_V4_single_species.R ./ $Name ${ReadsPerBarcode[$index]} ${Genomes[$index]} ${Type[$index]} --save
		fi
	    fi
	fi
	if [ ! -d "$dir/$Name/" ]; then 
	    mkdir $dir/$Name/ && mv $dir/fastqs/$Name.* $dir/$Name
	fi
    fi
    index=$((index + 1))
done



rm $dir/Useful/*filtered.counts.csv 
rm $dir/filesr1.xls $dir/filesr2.xls $dir/filesi1.xls $dir/filesi2.xls
rm -r $dir/temp/

gzip $dir/*/*.csv
gzip $dir/*/*.groups.tsv

# rm  $dir/*/*rmdup.bam* $dir/*/*wdup.bam* $dir/*/*grouped.bam* $dir/*/*wdup.all.bam*  $dir/*/*.fastq
# rm -r $dir/Indexed/* $dir/smallfastqs/* $dir/*/*Sub*fq.qz
touch $dir/smallfastqs/0001.$Run.1.I1.fastq.gz
touch $dir/Indexed/Sub.0001.$Run.1.R1.fastq.gz
touch $dir/temp4.fastq.gz

echo "The pipeline is completed!! Author: Sai Ma <sai@broadinstitute.org>"
exit

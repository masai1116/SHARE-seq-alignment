#!/bin/bash

# data dirs
dir1='/raid/USRdirs/Buenrostro_data/05_leukemia/02_mergeReplicates/00_Heme/'  #### merged and peakcalls
dir2='/raid/USRdirs/Buenrostro_data/12_schumanHeme/02_mergeReplicates/'
blacklist='/raid/USRdirs/Buenrostro_data/02_20130527_singleCell/23_mitoBlackList/JDB_blacklist.bed'
data1='/raid/USRdirs/Buenrostro_data/05_leukemia/00_HSCData/'
data2='/raid/USRdirs/Buenrostro_data/12_schumanHeme/00_data/'

# peakFile name
outPeak="PeakFile_"`date "+%Y%m%d"`".bed";echo $outPeak

# get top summits
echo "Getting top summits"
cat $dir1/peak*/*summits.bed | awk '{OFS="\t"}{print$1,$2-250,$3+250,$4,$5}' > allSummits.bed
cat $dir2/peak*/*summits.bed | awk '{OFS="\t"}{print$1,$2-250,$3+250,$4,$5}' >> allSummits.bed
pyGetTopSummits.py -i allSummits.bed -o allSummits.top.bed --top 99999999999
awk '{if(length($1)<6)print}' allSummits.top.bed | grep -v chrY | grep -v chrM | bedtools sort -i stdin | bedtools intersect -v -a stdin -b $blacklist > filteredPeakCalls.topSummits.bed
awk '{OFS="\t"}{print$1,$2,$3,$4":"$5}' filteredPeakCalls.topSummits.bed | sed 's/peakCalls\///g' | sed 's/peakCalls_singles\///g' > filteredPeakCalls.topSummits.withName.bed

# annotate peaks
echo "Annotating peaks"
annotatePeaks.pl <(cut -f1,2,3 filteredPeakCalls.topSummits.withName.bed) hg19 -size given -annStats peakStats.bed > annotateOutput.bed 2> annotate.log
awk 'NR>1' annotateOutput.bed | sed 's/\t\t/\tNA\t/g' | sed 's/ /,/g' | cut -f2-8,10-14,16,18 | bedtools sort -i stdin > $outPeak

# add read counts
echo "Getting peak counts data1"
outData1="Data_"`date "+%Y%m%d"`".hsc.bed";echo $outData1
cut -f2-4 annotateOutput.bed | head -1 >$outData1
cut -f1-3 $outPeak >>$outData1
for i in `ls $data1/*/*.flt.bam`
do
	echo $i
	pyBamBedCount.py -a $i -b $outData1 -o $outData1 --append --header --integer
done

# add read counts to Mega
data2='/raid/Tn5_chromatin/2015_ATACdata/151110_scATACheme_Mega/BM0828-MEGA-bulk/'
outPeak='PeakFile_20160523.bed'
echo "Getting peak counts data2"
outData2="Data_"`date "+%Y%m%d"`".mega.bed";echo $outData2
cut -f2-4 annotateOutput.bed | head -1 >$outData2
cut -f1-3 $outPeak >>$outData2
for i in `ls $data2/*/*.flt.bam`
do
        echo $i
        pyBamBedCount.py -a $i -b $outData2 -o $outData2 --append --header --integer
done

######### add new CD34 data
data4='/raid/Tn5_chromatin/2015_ATACdata/150731_scAR_mDCs_CD34/BM-CD34/'
outPeak='PeakFile_20160523.bed'
echo "Getting peak counts data4"
outData4="Data_"`date "+%Y%m%d"`".BM-CD34.bed";echo $outData4
cut -f2-4 annotateOutput.bed | head -1 >$outData4
cut -f1-3 $outPeak >>$outData4
for i in `ls $data4/*/*.flt.bam`
do
        echo $i
        pyBamBedCount.py -a $i -b $outData4 -o $outData4 --append --header --integer
done

######## MCP and UNK
data5='/raid/Tn5_chromatin/160201_MCP-Heme/BM-bulk/'
outPeak='PeakFile_20160523.bed'
echo "Getting peak counts data4"
outData5="Data_"`date "+%Y%m%d"`".MCP-UNK.bed";echo $outData5
cut -f2-4 annotateOutput.bed | head -1 >$outData5
cut -f1-3 $outPeak >>$outData5
for i in `ls $data5/*/*.flt.bam`
do
        echo $i
        pyBamBedCount.py -a $i -b $outData5 -o $outData5 --append --header --integer
done

######## mDC and pDC
data6='/raid/Tn5_chromatin/160422_TCR_heme/Bulk-DC/'
outPeak='PeakFile_20160523.bed'
echo "Getting peak counts data6"
outData6="Data_"`date "+%Y%m%d"`".DC.bed";echo $outData6
cut -f2-4 annotateOutput.bed | head -1 >$outData6
cut -f1-3 $outPeak >>$outData6
for i in `ls $data6/*/*.flt.bam`
do
        echo $i
        pyBamBedCount.py -a $i -b $outData6 -o $outData6 --append --header --integer
done

cd /media/ssd/sdc1/data/ljh/dnaHurt/

## narrowpeak 文件
FU1_peak=/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/peak/rep1/A5-FU-1_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz
FU2_peak=/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/peak/rep2/A5-FU-2_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz
FU3_peak=/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/peak/rep3/A5-FU-3_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz
PBS1_peak=/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/peak/rep1/PBS-1_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz
PBS2_peak=/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/peak/rep2/PBS-2_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz
PBS3_peak=/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/peak/rep3/PBS-3_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz
FUANA1_peak=/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_ANA/peak/rep1/A5-FU-Ana-1_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz
FUANA2_peak=/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_ANA/peak/rep2/A5-FU-Ana-2_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz
FUANA3_peak=/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_ANA/peak/rep3/A5-FU-Ana-3_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz

##
resDir=result/2025-10-31-样本间相关性/
mkdir -p $resDir

## 合并排序并merge
cat $FU1_peak \
    $FU2_peak \
    $FU3_peak \
    $PBS1_peak \
    $PBS2_peak \
    $PBS3_peak \
    $FUANA1_peak \
    $FUANA2_peak \
    $FUANA3_peak \
    | zcat \
    | cut -f1-3 \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin > $resDir/all.merge.narrowPeak


## 去掉不是chr开头的区域
grep "^chr" $resDir/all.merge.narrowPeak > $resDir/all.merge.narrowPeak.tmp 
mv $resDir/all.merge.narrowPeak.tmp $resDir/all.merge.narrowPeak

## sort
sort -k1,1 -k2,2n $resDir/all.merge.narrowPeak > $resDir/all.merge.sort.narrowPeak


## 排序
zcat $FU1_peak | sort -k1,1 -k2,2n > $resDir/FU1_peak
zcat $FU2_peak | sort -k1,1 -k2,2n > $resDir/FU2_peak
zcat $FU3_peak | sort -k1,1 -k2,2n > $resDir/FU3_peak
zcat $PBS1_peak | sort -k1,1 -k2,2n > $resDir/PBS1_peak
zcat $PBS2_peak | sort -k1,1 -k2,2n > $resDir/PBS2_peak
zcat $PBS3_peak | sort -k1,1 -k2,2n > $resDir/PBS3_peak
zcat $FUANA1_peak | sort -k1,1 -k2,2n > $resDir/FUANA1_peak
zcat $FUANA2_peak | sort -k1,1 -k2,2n > $resDir/FUANA2_peak
zcat $FUANA3_peak | sort -k1,1 -k2,2n > $resDir/FUANA3_peak

## bedtools map 将每个样本的narrowpeak文件与合并后的文件进行交集操作
## 并对落在同一个区域的信号值进行加和
## 并将结果保存到新的文件中
bedtools map -b $resDir/FU1_peak -a $resDir/all.merge.sort.narrowPeak -c 7 -o sum -null 0 > $resDir/FU1.map.narrowPeak
bedtools map -b $resDir/FU2_peak -a $resDir/all.merge.sort.narrowPeak -c 7 -o sum -null 0 > $resDir/FU2.map.narrowPeak
bedtools map -b $resDir/FU3_peak -a $resDir/all.merge.sort.narrowPeak -c 7 -o sum -null 0 > $resDir/FU3.map.narrowPeak
bedtools map -b $resDir/PBS1_peak -a $resDir/all.merge.sort.narrowPeak -c 7 -o sum -null 0 > $resDir/PBS1.map.narrowPeak
bedtools map -b $resDir/PBS2_peak -a $resDir/all.merge.sort.narrowPeak -c 7 -o sum -null 0 > $resDir/PBS2.map.narrowPeak
bedtools map -b $resDir/PBS3_peak -a $resDir/all.merge.sort.narrowPeak -c 7 -o sum -null 0 > $resDir/PBS3.map.narrowPeak
bedtools map -b $resDir/FUANA1_peak -a $resDir/all.merge.sort.narrowPeak -c 7 -o sum -null 0 > $resDir/FUANA1.map.narrowPeak
bedtools map -b $resDir/FUANA2_peak -a $resDir/all.merge.sort.narrowPeak -c 7 -o sum -null 0 > $resDir/FUANA2.map.narrowPeak
bedtools map -b $resDir/FUANA3_peak -a $resDir/all.merge.sort.narrowPeak -c 7 -o sum -null 0 > $resDir/FUANA3.map.narrowPeak




## check 
head -1 $resDir/FU1.map.narrowPeak
head -1 $resDir/FU2.map.narrowPeak
head -1 $resDir/FU3.map.narrowPeak
head -1 $resDir/PBS1.map.narrowPeak
head -1 $resDir/PBS2.map.narrowPeak
head -1 $resDir/PBS3.map.narrowPeak
head -1 $resDir/FUANA1.map.narrowPeak
head -1 $resDir/FUANA2.map.narrowPeak
head -1 $resDir/FUANA3.map.narrowPeak

## 剩下的拼接结果和计算相关性都在R语言中完成
## 2025-10-31-B-绘制样本间的相关性.R


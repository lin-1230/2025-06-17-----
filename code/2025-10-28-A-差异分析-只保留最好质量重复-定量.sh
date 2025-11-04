# @date 2025-10-28
# @author lin
# @description 因为这次D21 的ATAC 数据质量并不是特别好，
# 且这次数据没办法富集到细胞迁移的相关通路
# 怀疑有可能是数据质量的原因导致
# 所以决定只保留质量最好的数据进行分析看看

## 进行定量分析
## 根据质控结果分别选择最好质量的重复样本进行定量分析
## 5FU 重复样本1
## PBS 重复样本1
## 5FU_ANA 重复样本3
cd /media/ssd/sdc1/data/ljh/dnaHurt/

## bam文件
FU_bam='/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/align/rep1/A5-FU-1_R1.trim.srt.nodup.no_chrM_MT.bam'
PBS_bam='/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/align/rep1/PBS-1_R1.trim.srt.nodup.no_chrM_MT.bam'
FUANA_bam='/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_ANA/align/rep3/A5-FU-Ana-3_R1.trim.srt.nodup.no_chrM_MT.bam'

## narowpeak文件
FU_peak='/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/peak/rep1/A5-FU-1_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz'
PBS_peak='/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/peak/rep1/PBS-1_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz'
FUANA_peak='/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_ANA/peak/rep3/A5-FU-Ana-3_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz'

## merge 5FU_peak, PBS_peak, 5FUANA_peak
mkdir -p result/2025-10-28-chaYi_onlyBest/
merged_peak='result/2025-10-28-chaYi_onlyBest/merged_peak.narrowPeak'
zcat $FU_peak $PBS_peak $FUANA_peak | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i stdin  > $merged_peak
head -10 $merged_peak

## bedtools 定量
bedtools multicov -bams $PBS_bam $FU_bam $FUANA_bam -bed $merged_peak \
  > result/2025-10-28-chaYi_onlyBest/multicov_PBS_5FU_5FUANA.txt
head -10 result/2025-10-28-chaYi_onlyBest/multicov_PBS_5FU_5FUANA.txt


## homer 注释
annotatePeaks.pl $merged_peak  mm10 > result/2025-10-28-chaYi_onlyBest/annotatePeaks.homer.txt
head -5 result/2025-10-28-chaYi_onlyBest/annotatePeaks.homer.txt

## 标准化
## 公式: counts / (total_mapped_reads / 1e6)

## 差异分析

## GSEA



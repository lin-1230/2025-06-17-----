# @date 2025-10-11
# @description 对10-10找到的差异可及性区域进行注释

conda activate ATACseq
cd /media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/



## 输入文件
## narrow_DEseq2
## narrow_DEseq2
narrow_5FU_PBS_DEpeak_DEseq2='/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_DiffBind_res_narrowPeak_DEpeak1_DEseq2.bed'
narrow_5FUANA_5FU_DEpeak_DEseq2='/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_DiffBind_res_narrowPeak_DEpeak2_DEseq2.bed'
## narrow_edgeR
narrow_5FU_PBS_DEpeak_edgeR='/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_DiffBind_res_narrowPeak_DEpeak1_edgeR.bed'
narrow_5FUANA_5FU_DEpeak_edgeR='/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_DiffBind_res_narrowPeak_DEpeak2_edgeR.bed'

## 输出文件
# narrow_DEseq2
narrow_5FU_PBS_Ann_DEseq2='/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_5FU_PBS_DEpeak_Ann_DEseq2.txt'
narrow_5FUANA_5FU_Ann_DEseq2='/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_5FUANN_5FU_DEpeak_Ann_DEseq2.txt'
## narrow_edgeR
narrow_5FU_PBS_Ann_edgeR='/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_5FU_PBS_DEpeak_Ann_edgeR.txt'
narrow_5FUANA_5FU_Ann_edgeR='/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_5FUANN_5FU_DEpeak_Ann_edgeR.txt'


head -10 ${narrow_5FU_PBS_DEpeak_DEseq2}


## narrow_DEseq2
# 5FU PBS
cut -f1-3 ${narrow_5FU_PBS_DEpeak_DEseq2} | annotatePeaks.pl - mm10 > \
  ${narrow_5FU_PBS_Ann_DEseq2}
head -10 ${narrow_5FU_PBS_Ann_DEseq2}
# 5FU-ANA 5FU
cut -f1-3 ${narrow_5FUANA_5FU_DEpeak_DEseq2} | annotatePeaks.pl - mm10 > \
  ${narrow_5FUANA_5FU_Ann_DEseq2}

## narrow_edgeR
# 5FU PBS
cut -f1-3 ${narrow_5FU_PBS_DEpeak_edgeR} | annotatePeaks.pl - mm10 > \
  ${narrow_5FU_PBS_Ann_edgeR}
head -10 ${narrow_5FU_PBS_Ann_edgeR}
# 5FU-ANA 5FU
cut -f1-3 ${narrow_5FUANA_5FU_DEpeak_edgeR} | annotatePeaks.pl - mm10 > \
  ${narrow_5FUANA_5FU_Ann_edgeR}

########################
## IDR的结果
########################

## 输入文件
# DEseq2
idr_5FU_PBS_DEpeak_DEseq2=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_idrPeak_DEpeak1_DEseq2.bed
idr_5FUANA_5FU_DEpeak_DEseq2=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_idrPeak_DEpeak2_DEseq2.bed
# edgeR
idr_5FU_PBS_DEpeak_edgeR=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_idrPeak_DEpeak1_edgeR.bed
idr_5FUANA_5FU_DEpeak_edgeR=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_idrPeak_DEpeak2_edgeR.bed

## 输出文件
# DEseq2
idr_5FU_PBS_Ann_DEseq2=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_idr_5FU_PBS_DEpeak_Ann_DEseq2.txt
idr_5FUANA_5FU_Ann_DEseq2=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_idr_5FUANN_5FU_DEpeak_Ann_DEseq2.txt
# edgeR
idr_5FU_PBS_Ann_edgeR=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_idr_5FU_PBS_DEpeak_Ann_edgeR.txt
idr_5FUANA_5FU_Ann_edgeR=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_idr_5FUANN_5FU_DEpeak_Ann_edgeR.txt

## 5FU vs PBS 
# DEseq2
cut -f1-3 ${idr_5FU_PBS_DEpeak_DEseq2} | annotatePeaks.pl - mm10 > \
  ${idr_5FU_PBS_Ann_DEseq2}
# edgeR
cut -f1-3 ${idr_5FU_PBS_DEpeak_edgeR} | annotatePeaks.pl - mm10 > \
  ${idr_5FU_PBS_Ann_edgeR}


## 5FU-ANA vs 5FU
# DEseq2
cut -f1-3 ${idr_5FUANA_5FU_DEpeak_DEseq2} | annotatePeaks.pl - mm10 > \
  ${idr_5FUANA_5FU_Ann_DEseq2}
# edgeR
cut -f1-3 ${idr_5FUANA_5FU_DEpeak_edgeR} | annotatePeaks.pl - mm10 > \
  ${idr_5FUANA_5FU_Ann_edgeR}


########################
## 1.5M narrow的结果
########################

## 输入文件
# edgeR
narrow_5FU_PBS_DEpeak_edgeR_1p5M='/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/1.5M_DiffBind_res_narrowPeak_DEpeak1_edgeR.bed'
## 输出文件
narrow_5FU_PBS_Ann_edgeR_1p5M='/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/1.5M_DiffBind_res_narrowPeak_DEpeak1_Ann_edgeR.txt'

## 注释
cut -f1-3 ${narrow_5FU_PBS_DEpeak_edgeR_1p5M} | annotatePeaks.pl - mm10 > \
  ${narrow_5FU_PBS_Ann_edgeR_1p5M}
head -10 ${narrow_5FU_PBS_Ann_edgeR_1p5M}
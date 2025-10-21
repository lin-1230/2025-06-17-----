# @date 2025-10-11
# @description 对10-10找到的差异可及性区域进行注释

conda activate ATACseq
cd /media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/

## 5FU vs PBS
# head -10 /media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_narrowPeak_DEpeak1.bed
## 输入文件
narrow_5FU_PBS_DEpeak=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_narrowPeak_DEpeak1.bed
narrow_5FUANA_5FU_DEpeak=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_narrowPeak_DEpeak2.bed
## 输出文件
narrow_5FU_PBS_Ann=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_narrow_5FU_PBS_DEpeak_Ann.txt
narrow_5FUANA_5FU_Ann=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_narrow_5FUANN_5FU_DEpeak_Ann.txt

head -10 ${narrow_5FU_PBS_DEpeak}

## 5FU vs PBS 
# 取前三列
cut -f1-3 ${narrow_5FU_PBS_DEpeak} | annotatePeaks.pl - mm10 > \
  ${narrow_5FU_PBS_Ann}

## 5FU-ANA vs 5FU
cut -f1-3 ${narrow_5FUANA_5FU_DEpeak} | annotatePeaks.pl - mm10 > \
  ${narrow_5FUANA_5FU_Ann}

########################
## IDR的结果
########################

## 输入文件
idr_5FU_PBS_DEpeak=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_idrPeak_DEpeak1.bed
idr_5FUANA_5FU_DEpeak=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_idPeak_DEpeak2.bed
## 输出文件
idr_5FU_PBS_Ann=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_idr_5FU_PBS_DEpeak_Ann.txt
idr_5FUANA_5FU_Ann=/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_idr_5FUANN_5FU_DEpeak_Ann.txt

head -10 ${idr_5FU_PBS_DEpeak}

## 5FU vs PBS 
# 取前三列
cut -f1-3 ${idr_5FU_PBS_DEpeak} | annotatePeaks.pl - mm10 > \
  ${idr_5FU_PBS_Ann}

## 5FU-ANA vs 5FU
cut -f1-3 ${idr_5FUANA_5FU_DEpeak} | annotatePeaks.pl - mm10 > \
  ${idr_5FUANA_5FU_Ann}



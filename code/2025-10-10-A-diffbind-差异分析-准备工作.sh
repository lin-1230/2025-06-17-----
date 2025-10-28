# 1. 安装DiffBind 
# 检查哪个conda环境安装了diffbind包
# 因为ATACseq环境装不上包
for e in $(conda env list | awk 'NR>2 {print $1}'); \
do echo "== $e =="; conda list -n "$e" | grep -iE '^(bioconductor-)?diffbind' || echo "not installed"; \
done

# conda create -n diffbind_env r-base=4.2 bioconda::bioconductor-diffbind

# 2. 解压peak文件
## 5-FU 样本 IDR-peak文件
# -k 保留原始压缩文件
gunzip -k -c \
  /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz > \
  /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/peak/idr_reproducibility/idr.conservative_peak.narrowPeak

## PBS 样本 IDR-peak文件
gunzip -k -c \
  /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz > \
  /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/peak/idr_reproducibility/idr.conservative_peak.narrowPeak

## 5-FU-ANA IDR-peak文件
gunzip -k -c \
  /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_ANA/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz > \
  /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_ANA/peak/idr_reproducibility/idr.conservative_peak.narrowPeak


## 常规做法是将call的peak，不同重复进行merge成更大的peak文件
# PBS
PBS_peak1=/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/peak/rep1/PBS-1_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak
PBS_peak2=/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/peak/rep2/PBS-2_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak
PBS_peak3=/media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/peak/rep3/PBS-3_R1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak

wc -l $PBS_peak1
wc -l $PBS_peak2
wc -l $PBS_peak3



## 解压
gunzip -k -c \
  ${PBS_peak1}.gz > \
  ${PBS_peak1}
  
gunzip -k -c \
  ${PBS_peak2}.gz > \
  ${PBS_peak2}
  
gunzip -k -c \
  ${PBS_peak3}.gz > \
  ${PBS_peak3}

## 合并、排序和merge
cat ${PBS_peak1} ${PBS_peak2} ${PBS_peak3} | bedtools sort -i - | bedtools merge -i - > \
  /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/peak/PBS_merge.narrowPeak

head -20 /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/peak/PBS_merge.narrowPeak
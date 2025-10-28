# @date 2025-10-27
# @description 前面对5FU组和PBS组进行了差异分析
# 现在想看这些区域在ANA组是否会被挽回

conda activate ATACseq
cd /media/ssd/sdc1/data/ljh/dnaHurt/
mkdir -p result/2025-10-27-plotHeatMap/

## 生成矩阵
tail -n +2 DiffBind/2025_10_22_res/narrow_5FU_PBS_edgeR_gain_p_005_Fold_0.csv | cut -d ',' -f1-3 | tr ',' '\t' > result/2025-10-27-plotHeatMap/narrow_5FU_PBS_edgeR_gain_p_005_Fold_0.bed
head -10 result/2025-10-27-plotHeatMap/narrow_5FU_PBS_edgeR_gain_p_005_Fold_0.bed
## 去掉随机染色组区域
grep -v "random" result/2025-10-27-plotHeatMap/narrow_5FU_PBS_edgeR_gain_p_005_Fold_0.bed | \
grep -v "chrUn" | \
grep -v "hap" > result/2025-10-27-plotHeatMap/narrow_5FU_PBS_edgeR_gain_p_005_Fold_0.filtered.bed
## 去掉引号
tr -d '"' < result/2025-10-27-plotHeatMap/narrow_5FU_PBS_edgeR_gain_p_005_Fold_0.filtered.bed > result/2025-10-27-plotHeatMap/narrow_5FU_PBS_edgeR_gain_p_005_Fold_0.filtered.bed.noquote
head -10 result/2025-10-27-plotHeatMap/narrow_5FU_PBS_edgeR_gain_p_005_Fold_0.filtered.bed.noquote


computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -R result/2025-10-27-plotHeatMap/narrow_5FU_PBS_edgeR_gain_p_005_Fold_0.filtered.bed.noquote \
-S /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/signal/rep1/PBS-1_R1.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig\
    /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/signal/rep2/PBS-2_R1.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig \
    /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/signal/rep3/PBS-3_R1.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig \
    /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/signal/rep1/A5-FU-1_R1.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig \
    /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/signal/rep2/A5-FU-2_R1.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig \
    /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/signal/rep3/A5-FU-3_R1.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig \
    /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_ANA/signal/rep1/A5-FU-Ana-1_R1.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig \
    /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_ANA/signal/rep2/A5-FU-Ana-2_R1.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig \
    /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_ANA/signal/rep3/A5-FU-Ana-3_R1.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig \
    --skipZeros \
    -o result/2025-10-27-plotHeatMap/matrix_5FU_PBS_ANA_edgeR_gain_p_005_Fold_0.gz




plotHeatmap \
  -m result/2025-10-27-plotHeatMap/matrix_5FU_PBS_ANA_edgeR_gain_p_005_Fold_0.gz \
  -out result/2025-10-27-plotHeatMap/heatmap_5FU_PBS_ANA.pdf \
  --samplesLabel "PBS-1" "PBS-2" "PBS-3" "5-FU-1" "5-FU-2" "5-FU-3" "5-FU-ANA-1" "5-FU-ANA-2" "5-FU-ANA-3" \
  --colorMap RdBu \
  --sortRegions descend \
  --plotTitle ""



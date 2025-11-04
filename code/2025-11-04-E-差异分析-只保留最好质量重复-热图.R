# @date 2025-11-04
# @description 前面对5FU组和PBS组进行了差异分析（只保留了5-FU1和PBS1）
# 现在想看这些区域在ANA组是否会被挽回

conda activate ATACseq
cd /media/ssd/sdc1/data/ljh/dnaHurt/
mkdir -p result/2025-11-04-plotHeatMap/

## 生成矩阵
DEpeakDF='result/2025-10-28-chaYi_onlyBest/diff_analysis_5FUANA_5FU_gt5.csv'
head -5 $DEpeakDF

tail -n +2 $DEpeakDF | cut -d ',' -f2-4 | tr ',' '\t' > result/2025-11-04-plotHeatMap/narrow_5FU_1_PBS_1_gt5.bed
head -10 result/2025-11-04-plotHeatMap/narrow_5FU_1_PBS_1_gt5.bed

## 去掉随机染色组区域
grep -v "random" result/2025-11-04-plotHeatMap/narrow_5FU_1_PBS_1_gt5.bed | \
grep -v "chrUn" | \
grep -v "hap" > result/2025-11-04-plotHeatMap/narrow_5FU_1_PBS_1_gt5.filtered.bed
## 去掉引号
tr -d '"' < result/2025-11-04-plotHeatMap/narrow_5FU_1_PBS_1_gt5.filtered.bed > result/2025-11-04-plotHeatMap/narrow_5FU_1_PBS_1_gt5.filtered.bed.noquote
head -10 result/2025-11-04-plotHeatMap/narrow_5FU_1_PBS_1_gt5.filtered.bed.noquote


computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -R result/2025-11-04-plotHeatMap/narrow_5FU_1_PBS_1_gt5.filtered.bed.noquote \
-S /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/signal/rep1/PBS-1_R1.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig\
    /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/signal/rep1/A5-FU-1_R1.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig \
    /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_ANA/signal/rep3/A5-FU-Ana-3_R1.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig \
    --skipZeros \
    -o result/2025-11-04-plotHeatMap/matrix_5FU_1_PBS_1_ANA_3_gt5.gz




plotHeatmap \
  -m result/2025-11-04-plotHeatMap/matrix_5FU_1_PBS_1_ANA_3_gt5.gz \
  -out result/2025-11-04-plotHeatMap/heatmap_5FU_1_PBS_1_ANA_3_gt5.pdf \
  --samplesLabel "PBS-1" "5-FU-1" "5-FU-ANA-3" \
  --colorMap BuRd \
  --sortRegions descend \
  --plotTitle ""



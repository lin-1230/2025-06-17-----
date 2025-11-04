## 进入环境
# conda activate hicseq
# R

setwd('/media/ssd/sdc1/data/ljh/dnaHurt/')
resDir <- 'result/2025-11-02-chaYi_onlyBest/'
dir.create(resDir,showWarnings = F)

## 读取差异分析结果
library(data.table)
data <- fread('result/2025-10-28-chaYi_onlyBest/diff_analysis.csv')
head(data)

## 绘制log2FC_5FU_PBS的分布
library(ggplot2)
themeSet <- theme(
  axis.text = element_text(size = 12, color = "black"),
  axis.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14, face = "bold"),
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
)


## 绘制log2FC_5FU_PBS的分布
p <- ggplot(data, aes(x = log2FC_5FU_PBS)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(title = "Distribution of log2FC_5FU_PBS",
       x = "log2FC_5FU_PBS",
       y = "Frequency") +
  themeSet
ggsave(paste0(resDir,'log2FC_5FU_PBS_distribution.pdf'),p,width = 8,height = 6)


## 绘制log2FC_5FU_PBS的分布（仅保留promoter区域）
p <- ggplot(data_promoter, aes(x = log2FC_5FU_PBS)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(title = "Distribution of log2FC_5FU_PBS (promoter-TSS)",
       x = "log2FC_5FU_PBS",
       y = "Frequency") +
  themeSet
ggsave(paste0(resDir,'log2FC_5FU_PBS_distribution_promoter.pdf'),p,width = 8,height = 6)


## 取出其中promoter 区域的结果
data_promoter <- data[data$newAnnotation == "promoter-TSS",]
head(data_promoter)

## 定义细胞周期基因集
cell_cycle_geneSet <- c('Myc','Fbl','Cond1','Cond2',
                        'Dkc1','Nop10','Polr1a','Eif4a1',
                        'Mars','Nars')

## 提取和细胞周期基因集相关的结果
data_promoter_cell_cycle <- data_promoter[data_promoter$`Gene Name` %in% cell_cycle_geneSet,]
head(data_promoter_cell_cycle)


## 按照基因名称和FC值排序
# 先按 Gene Name 升序，再按 log2FoldChange 降序
data_promoter_cell_cycle <- data_promoter_cell_cycle[order(data_promoter_cell_cycle$`Gene Name`,
                                                             -data_promoter_cell_cycle$log2FC_5FU_PBS),]
head(data_promoter_cell_cycle)


## 保存结果
write.csv(data_promoter_cell_cycle, paste0(resDir,'DEpeak_promoter_cell_cycle.csv'),row.names = F)
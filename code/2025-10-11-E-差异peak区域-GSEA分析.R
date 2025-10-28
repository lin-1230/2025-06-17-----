## @date 2025-10-18
## @description 对ATAC数据做GSEA分析

# conda activate RNAseq
# R

## TODO
## 准备数据
# 1. 基因Rnak列表
# 2. 获取TCA和FAO通路
# 3. 富集分析
# 4. 可视化

dir.create('/media/ssd/sdc1/data/ljh/dnaHurt/GSEA/',recursive = T)
setwd('/media/ssd/sdc1/data/ljh/dnaHurt/')
source('code/config_seurat.R')

## 自定义读取并筛选差异基因的方程
# 定义读取差异基因的函数
# 加载库
library(msigdbr)
library(fgsea)
library(data.table)
library(Seurat)
library(qs)



resDir <- 'result/2025-10-18/'
dir.create(resDir,recursive = T)


## 制作rank向量
#  rank向量总共需要两列，一列是基因名，一列是基因的log2FC值

## 读取peak注释结果
narrow_5FU_PBS_Ann_dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/narrow_5FU_PBS_Ann_DEpeak.csv'
narrow_5FUANA_5FU_Ann_dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/narrow_5FUANA_5FU_Ann_DEpeak.csv'

narrow_5FU_PBS_Ann <- fread(narrow_5FU_PBS_Ann_dir)
narrow_5FUANA_5FU_Ann <- fread(narrow_5FUANA_5FU_Ann_dir)

colnames(narrow_5FU_PBS_Ann)
colnames(narrow_5FUANA_5FU_Ann)

head(narrow_5FU_PBS_Ann)
head(narrow_5FUANA_5FU_Ann)

## 对narrow_5FU_PBS_Ann按Gene Name去重，保留Fold最大行
FiveFU_PBS_max <- narrow_5FU_PBS_Ann[, .(Fold = max(Fold)), by = `Gene Name`]
rank_5FU_PBS_max <- FiveFU_PBS_max$Fold
names(rank_5FU_PBS_max) <- FiveFU_PBS_max$`Gene Name`
rank_5FU_PBS_max <- sort(rank_5FU_PBS_max,decreasing = TRUE)



## 对narrow_5FU_PBS_Ann按Gene Name分组，计算Fold均值
FiveFU_PBS_mean <- narrow_5FU_PBS_Ann[, .(Fold = mean(Fold)), by = `Gene Name`]
rank_5FU_PBS_mean <- FiveFU_PBS_mean$Fold
names(rank_5FU_PBS_mean) <- FiveFU_PBS_mean$`Gene Name`
rank_5FU_PBS_mean <- sort(rank_5FU_PBS_mean,decreasing = TRUE)


# 获取所有基因集(小鼠)
msig_data <- msigdbr(species = "Mus musculus") 


# 获取包含“MIGRATION”关键字的基因集
migration_genesets <- msig_data[grep("MIGRATION", msig_data$gs_name), ]
# 查看基因集
class(migration_genesets)
colnames(migration_genesets)
head(migration_genesets)

## 保存migration_genesets为csv
write.csv(migration_genesets,paste0(resDir,'migration_genesets.csv'),row.names = F)



migration_genesets_gs <- migration_genesets[,c('gs_name','gs_description')]
## 去重
migration_genesets_gs <- unique(migration_genesets_gs)
head(migration_genesets_gs)

## 保存migration_genesets_gs为csv
write.csv(migration_genesets_gs,paste0(resDir,'migration_genesets_gs.csv'),row.names = F)

gsea_migration_max <- fgsea(
  pathways = split(migration_genesets$gene_symbol, migration_genesets$gs_name),
  stats = rank_5FU_PBS_max, 
  nPermSimple = 10000  # 设置置换次数
)

head(gsea_migration_max)
dim(gsea_migration_max)
class(gsea_migration_max)

# gsea_migration_max2 <- gsea_migration_max[,c('pathway','pval','padj','NES','size')]
## 合并gsea_migration_max2和migration_genesets_gs
gsea_migration_max3 <- merge(gsea_migration_max, migration_genesets_gs, by.x = 'pathway', by.y = 'gs_name', all = FALSE)
head(gsea_migration_max3)

## 按照NES排序
gsea_migration_max3 <- gsea_migration_max3[order(gsea_migration_max3$NES,decreasing = TRUE),]
## 给显著的通路添加一个列
gsea_migration_max3$significant <- ifelse(gsea_migration_max3$pval < 0.05, 'significant', 'not significant')

library(dplyr)
library(purrr)
gsea_migration_max3_to_save <- gsea_migration_max3 %>% mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";"))

## 保存结果
write.csv(gsea_migration_max3_to_save,paste0(resDir,'fgsea_migration_max3.csv'),row.names = F)


## 可视化

## 全局ggplot2主题设置
axis_text_size <- 20  # 全局坐标轴字体大小
legend_text_size <- 18  # 全局图例字体大小
label_text_size <- 20  # 全局标签字体大小
title_text_size <- 20  # 全局标题字体大小

pdf_width <- 13 # PDF图形默认宽度
pdf_height <- 10 # PDF图形默认高度

themeSet <- theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),  # 设置面板背景
        plot.title = element_text(hjust = 0.5, size = title_text_size),  # 设置标题水平居中对齐并增加字体大小
        strip.text = element_text(size = title_text_size,  # 继承标题字号
        face = "bold",          # 加粗显示
        ),  # 设置strip文本样式
        axis.line = element_line(colour = "black", linewidth = 0.25),  # 设置坐标轴线样式
        axis.title = element_text(size = (axis_text_size + 2), face = "plain", color = "black"),  # 设置坐标轴标题样式
        axis.text = element_text(size = axis_text_size, face = "plain", color = "black"),  # 设置坐标轴标签样式
        legend.text = element_text(face = "plain", colour = "black", size = (legend_text_size)),  # 设置图例文本样式
        legend.title = element_text(face = "plain", colour = "black", size = (legend_text_size+2)),  # 设置图例标题样式
        panel.grid.major = element_blank(),  # 隐藏主要网格线
        panel.grid.minor = element_blank(),  # 隐藏次要网格线
        axis.text.x = element_text(angle = 0, hjust = 0.5))  # 居中对齐


# 使用ggplot2可视化结果
library(ggplot2)

p1 <- ggplot(gsea_migration_max3, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, color = pval)) +
  theme_minimal() +
  labs(title = "Migration Related Pathways Enrichment") +
  xlab("Normalized Enrichment Score (NES)") +
  ylab("Pathway")+ themeSet

ggsave(p1,filename = paste0(resDir,'fgsea_migration_max.pdf'),width = pdf_width*3,height = pdf_height*4)

colnames(gsea_migration_max3)
## 横坐标为NES，纵坐标为FDR值的可视化
p_fdr <- ggplot(gsea_migration_max3,
                aes(x = NES, y = padj)) +
  geom_point(aes(size = size, color = pval)) +
  scale_y_continuous(
                     breaks = c(0.05,0.1,0.5,1),
                     labels = c(  "0.05", "0.1", "0.5", "1")) +
  scale_color_gradient(low = "red", high = "blue", name = "FDR") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
  labs(title = "",
       x = "Normalized Enrichment Score (NES)",
       y = "FDR") +
  theme_minimal() +
  themeSet

ggsave(p_fdr,
       filename = paste0(resDir, "NES_vs_FDR_scatter.pdf"),
       width = pdf_width,
       height = pdf_height)


## 筛选所有GO的功能通路



## 绘制所有基因集的FDR-NES散点图
gsea_all_max <- fgsea(
  pathways = split(msig_data$gene_symbol, msig_data$gs_name),
  stats = rank_5FU_PBS_max, 
  nPermSimple = 10000  # 设置置换次数
)

p_fdr <- ggplot(gsea_all_max,
                aes(x = NES, y = padj)) +
  geom_point(aes(size = size, color = pval)) +
  scale_y_continuous(
                     breaks = c(0.05,0.1,0.5,1),
                     labels = c(  "0.05", "0.1", "0.5", "1")) +
  scale_color_gradient(low = "red", high = "blue", name = "FDR") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
  labs(title = "",
       x = "Normalized Enrichment Score (NES)",
       y = "FDR") +
  theme_minimal() +
  themeSet

ggsave(p_fdr,
       filename = paste0(resDir, "NES_vs_FDR_scatter_all.pdf"),
       width = pdf_width,
       height = pdf_height)


## 保存结果
gsea_all_max <- gsea_all_max[order(gsea_all_max$NES,decreasing = TRUE),]
head(gsea_all_max)
write.csv(gsea_all_max, paste0(resDir,'fgsea_all_max.csv'),row.names = F)
class(gsea_all_max)



# 同时筛选 GOBP、GOCC 和 GOMF 三类 GO 通路
gsea_all_max_go <- gsea_all_max[grep("GOBP|GOCC|GOMF", gsea_all_max$pathway), ]
head(gsea_all_max_go)
gsea_all_max_go_save <- gsea_all_max_go[order(gsea_all_max_go$NES,decreasing = TRUE),]
gsea_all_max_go_save <-  gsea_all_max_go_save[,-c('leadingEdge')]
## 保存结果
write.csv(gsea_all_max_go_save, paste0(resDir,'fgsea_all_max_go.csv'),row.names = F)


# HSC_quitancer <- msig_data[grep("GOBP_NEGATIVE_REGULATION_OF_HEMATOPOIETIC_STEM_CELL_PROLIFERATION", msig_data$gs_name), ]
gsea_HSC_quitancer <- gsea_all_max_go[gsea_all_max_go$pathway == 'GOBP_NEGATIVE_REGULATION_OF_HEMATOPOIETIC_STEM_CELL_PROLIFERATION', ]
gsea_HSC_quitancer


gsea_all_max_go[gsea_all_max_go$pathway %like% 'CELL_DIFFERENTIATION', ]

hist(rank_5FU_PBS_max)
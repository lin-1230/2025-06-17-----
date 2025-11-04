## @date 2025-10-22
## @description 对ATAC数据做GSEA分析

# conda activate RNAseq
# R


## 准备数据
# 1. 基因Rnak列表
# 2. 获取TCA和FAO通路
# 3. 富集分析
# 4. 可视化

dir.create('/media/ssd/sdc1/data/ljh/dnaHurt/GSEA/',recursive = T)
setwd('/media/ssd/sdc1/data/ljh/dnaHurt/')
# source('code/config_seurat.R')

## 自定义读取并筛选差异基因的方程
# 定义读取差异基因的函数
# 加载库
library(msigdbr)
library(fgsea)
library(data.table)
library(Seurat)
library(qs)



resDir <- 'result/2025-10-22/'
dir.create(resDir,recursive = T)

## NOTE 暂时只分析 narrow_5FU_PBS_edgeR_P0.05的结果
## 制作rank向量
#  rank向量总共需要两列，一列是基因名，一列是基因的log2FC值

## 读取peak注释结果
## 5FU vs PBS
narrow_5FU_PBS_Ann_dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/narrow_5FU_PBS_edgeR_Ann_DEpeak.csv'
narrow_5FU_PBS_Ann <- fread(narrow_5FU_PBS_Ann_dir)

## 5FUANA vs 5FU
narrow_5FUANA_5FU_Ann_dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/narrow_5FUANA_5FU_edgeR_Ann_DEpeak.csv'
narrow_5FUANA_5FU_Ann <- fread(narrow_5FUANA_5FU_Ann_dir)

## 1.5M 5FU vs PBS
narrow_5FU_PBS_Ann_1p5M_dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/narrow_5FU_PBS_edgeR_1p5M_Ann_DEpeak.csv'
narrow_5FU_PBS_Ann_1p5M <- fread(narrow_5FU_PBS_Ann_1p5M_dir)
head(narrow_5FU_PBS_Ann_1p5M)



colnames(narrow_5FU_PBS_Ann)
head(narrow_5FU_PBS_Ann)


## 取绝对值最大值
## 对narrow_5FU_PBS_Ann按Gene Name去重，保留Fold绝对值最大行
FiveFU_PBS_max <- narrow_5FU_PBS_Ann[, .(Fold = Fold[which.max(abs(Fold))]), by = `Gene Name`]
rank_5FU_PBS_max <- FiveFU_PBS_max$Fold
names(rank_5FU_PBS_max) <- FiveFU_PBS_max$`Gene Name`
rank_5FU_PBS_max <- sort(rank_5FU_PBS_max,decreasing = TRUE)

## 对narrow_5FUANA_5FU_Ann按Gene Name去重，保留Fold绝对值最大行
FiveFUANA_5FU_max <- narrow_5FUANA_5FU_Ann[, .(Fold = Fold[which.max(abs(Fold))]), by = `Gene Name`]
rank_5FUANA_5FU_max <- FiveFUANA_5FU_max$Fold
names(rank_5FUANA_5FU_max) <- FiveFUANA_5FU_max$`Gene Name`
rank_5FUANA_5FU_max <- sort(rank_5FUANA_5FU_max,decreasing = TRUE)
head(rank_5FUANA_5FU_max)

## 1.5M 5FU vs PBS
FiveFU_PBS_1p5M_max <- narrow_5FU_PBS_Ann_1p5M[, .(Fold = Fold[which.max(abs(Fold))]), by = `Gene Name`]
rank_5FU_PBS_1p5M_max <- FiveFU_PBS_1p5M_max$Fold
names(rank_5FU_PBS_1p5M_max) <- FiveFU_PBS_1p5M_max$`Gene Name`
rank_5FU_PBS_1p5M_max <- sort(rank_5FU_PBS_1p5M_max,decreasing = TRUE)
head(rank_5FU_PBS_1p5M_max)



## 先promoter，然后再绝对值最大值
## 5FU vs PBS
FiveFU_PBS_promoter <- narrow_5FU_PBS_Ann[newAnnotation=='promoter-TSS',]
FiveFU_PBS_promoter_max <- FiveFU_PBS_promoter[, .(Fold = Fold[which.max(abs(Fold))]), by = `Gene Name`]
rank_5FU_PBS_promoter_max <- FiveFU_PBS_promoter_max$Fold
names(rank_5FU_PBS_promoter_max) <- FiveFU_PBS_promoter_max$`Gene Name`
rank_5FU_PBS_promoter_max <- sort(rank_5FU_PBS_promoter_max,decreasing = TRUE)
head(rank_5FU_PBS_promoter_max)
length(rank_5FU_PBS_promoter_max)

## 5FUANA vs 5FU
FiveFUANA_5FU_promoter <- narrow_5FUANA_5FU_Ann[newAnnotation=='promoter-TSS',]
FiveFUANA_5FU_promoter_max <- FiveFUANA_5FU_promoter[, .(Fold = Fold[which.max(abs(Fold))]), by = `Gene Name`]
rank_5FUANA_5FU_promoter_max <- FiveFUANA_5FU_promoter_max$Fold
names(rank_5FUANA_5FU_promoter_max) <- FiveFUANA_5FU_promoter_max$`Gene Name`
rank_5FUANA_5FU_promoter_max <- sort(rank_5FUANA_5FU_promoter_max,decreasing = TRUE)
head(rank_5FUANA_5FU_promoter_max)
length(rank_5FUANA_5FU_promoter_max)

## 1.5M 5FU vs PBS
FiveFU_PBS_1p5M_promoter <- narrow_5FU_PBS_Ann_1p5M[newAnnotation=='promoter-TSS',]
FiveFU_PBS_1p5M_promoter_max <- FiveFU_PBS_1p5M_promoter[, .(Fold = Fold[which.max(abs(Fold))]), by = `Gene Name`]
rank_5FU_PBS_1p5M_promoter_max <- FiveFU_PBS_1p5M_promoter_max$Fold
names(rank_5FU_PBS_1p5M_promoter_max) <- FiveFU_PBS_1p5M_promoter_max$`Gene Name`
rank_5FU_PBS_1p5M_promoter_max <- sort(rank_5FU_PBS_1p5M_promoter_max,decreasing = TRUE)
head(rank_5FU_PBS_1p5M_promoter_max)
length(rank_5FU_PBS_1p5M_promoter_max)


pdf_width <- 13 # PDF图形默认宽度
pdf_height <- 10 # PDF图形默认高度


## 画一下rank_5FU_PBS_max的分布
pdf(paste0(resDir, "rank_5FU_PBS_max_hist.pdf"),
    width = pdf_width,
    height = pdf_height)
p <- hist(rank_5FU_PBS_max)
dev.off( )

## 画一下rank_5FUANA_5FU_max的分布
pdf(paste0(resDir, "rank_5FUANA_5FU_max_hist.pdf"),
    width = pdf_width,
    height = pdf_height)
p <- hist(rank_5FUANA_5FU_max)
dev.off( )

## 画一下rank_5FU_PBS_1p5M_max的分布
pdf(paste0(resDir, "rank_5FU_PBS_1p5M_max_hist.pdf"),
    width = pdf_width,
    height = pdf_height)
p <- hist(rank_5FU_PBS_1p5M_max)
dev.off( )

## 画一下rank_5FU_PBS_promoter_max的分布
pdf(paste0(resDir, "rank_5FU_PBS_promoter_max_hist.pdf"),
    width = pdf_width,
    height = pdf_height)
p <- hist(rank_5FU_PBS_promoter_max)
dev.off( )

## 画一下rank_5FUANA_5FU_promoter_max的分布
pdf(paste0(resDir, "rank_5FUANA_5FU_promoter_max_hist.pdf"),
    width = pdf_width,
    height = pdf_height)
p <- hist(rank_5FUANA_5FU_promoter_max)
dev.off( )

## 画一下rank_5FU_PBS_1p5M_promoter_max的分布
pdf(paste0(resDir, "rank_5FU_PBS_1p5M_promoter_max_hist.pdf"),
    width = pdf_width,
    height = pdf_height)
p <- hist(rank_5FU_PBS_1p5M_promoter_max)
dev.off( )



## 将 G0 通路 GSEA 分析封装为函数，方便复用
run_gsea_by_keyword <- function(msig_data,keyword, rank_vec, res_dir, prefix = ""){
  ## 参数说明
  # msig_data: 基因集数据框，包含基因集名称和基因符号
  # keyword: 用于筛选基因集的关键字，例如 "G0"、"GOCC" 等
  # rank_vec: 基因排名向量，用于计算富集分数
  # res_dir: 结果保存目录
  # prefix: 结果文件名前缀，默认空字符串
  
  library(dplyr)
 
  ## 根据关键字筛选通路
  genesets <- msig_data[grep(keyword, names(msig_data), ignore.case = TRUE)]
  print(names(genesets))
  
  ## 富集分析
  gsea_res <- fgsea(
    pathways = genesets,
    stats    = rank_vec,
    nPermSimple = 10000
  )
  
  ## 结果整理
  gsea_res <- gsea_res[order(gsea_res$NES, decreasing = TRUE), ]
  gsea_res$significant <- ifelse(gsea_res$pval < 0.05, "significant", "not significant")
  
  ## 保存结果
  gsea_to_save <- gsea_res %>%
    as.data.frame() %>%
    mutate(leadingEdge = purrr::map_chr(leadingEdge, paste, collapse = ";"))
  
  write.csv(gsea_to_save,
            file.path(res_dir, paste0("fgsea_", prefix, keyword, ".csv")),
            row.names = FALSE)
  
  return(gsea_res)
}

## 获取所有基因集(小鼠)
library(GSEABase)
msig_data_raw <- getGmt('/media/ssd/sdb1/data/ljh/software/msigdb_v2025.1.Mm_GMTs/msigdb.v2025.1.Mm.symbols.gmt')
msig_data <- lapply(msig_data_raw, geneIds)
names(msig_data) <- sapply(msig_data_raw, setName)

## 调用函数完成 G0 分析
## 5FU vs PBS
gsea_G0_max <- run_gsea_by_keyword(msig_data,"G0",
                                   rank_5FU_PBS_max,
                                   resDir,
                                   prefix = "5FU_PBS_")
## 5FUANA vs 5FU
gsea_G0_max_5FUANA_5FU <- run_gsea_by_keyword(msig_data,"G0",
                                   rank_5FUANA_5FU_max,
                                   resDir,
                                   prefix = "5FUANA_5FU_")

## 5FUANA vs 5FU promoter max
gsea_G0_max_5FUANA_5FU_promoter <- run_gsea_by_keyword(msig_data,"G0",
                                   rank_5FUANA_5FU_promoter_max,
                                   resDir,
                                   prefix = "5FUANA_5FU_promoter_max_")
## 5FU vs PBS promoter max
gsea_G0_max_5FU_PBS_promoter <- run_gsea_by_keyword(msig_data,"G0",
                                   rank_5FU_PBS_promoter_max,
                                   resDir,
                                   prefix = "5FU_PBS_promoter_max_")



## quiescence相关通路
## 5FU vs PBS max
gsea_quiescence_max <- run_gsea_by_keyword(msig_data,
                                          keyword = "QUIESCENCE",
                                          rank_vec = rank_5FU_PBS_max,
                                          res_dir = resDir,
                                          prefix = "5FU_PBS_")
## 5FUANA vs 5FU max
gsea_quiescence_max_5FUANA_5FU <- run_gsea_by_keyword(msig_data,
                                          keyword = "QUIESCENCE",
                                          rank_vec = rank_5FUANA_5FU_max,
                                          res_dir = resDir,
                                          prefix = "5FUANA_5FU_")
## 5FUANA vs 5FU promoter max
gsea_quiescence_max_5FUANA_5FU_promoter <- run_gsea_by_keyword(msig_data,
                                          keyword = "QUIESCENCE",
                                          rank_vec = rank_5FUANA_5FU_promoter_max,
                                          res_dir = resDir,
                                          prefix = "5FUANA_5FU_promoter_max_")
## 5FU vs PBS promoter max
gsea_quiescence_max_5FU_PBS_promoter <- run_gsea_by_keyword(msig_data,
                                          keyword = "QUIESCENCE",
                                          rank_vec = rank_5FU_PBS_promoter_max,
                                          res_dir = resDir,
                                          prefix = "5FU_PBS_promoter_max_")



## arrest相关通路
gsea_arrest_max <- run_gsea_by_keyword(msig_data,
                                       keyword = "ARREST",
                                       rank_vec = rank_5FU_PBS_max,
                                       res_dir = resDir,
                                       prefix = "")

## dormancy相关通路
gsea_dormancy_max <- run_gsea_by_keyword(msig_data,
                                       keyword = "DORMANCY",
                                       rank_vec = rank_5FU_PBS_max,
                                       res_dir = resDir,
                                       prefix = "")

## 细胞迁移相关通路
## 5FU vs PBS max
gsea_cell_migration_max <- run_gsea_by_keyword(msig_data,
                                       keyword = "CELL_MIGRATION",
                                       rank_vec = rank_5FU_PBS_max,
                                       res_dir = resDir,
                                       prefix = "5FU_PBS_")

## 5FUANA vs 5FU max
gsea_cell_migration_max_5FUANA_5FU <- run_gsea_by_keyword(msig_data,
                                       keyword = "CELL_MIGRATION",
                                       rank_vec = rank_5FUANA_5FU_max,
                                       res_dir = resDir,
                                       prefix = "5FUANA_5FU_")

## 5FUANA vs 5FU promoter max
gsea_cell_migration_max_5FUANA_5FU_promoter <- run_gsea_by_keyword(msig_data,
                                       keyword = "CELL_MIGRATION",
                                       rank_vec = rank_5FUANA_5FU_promoter_max,
                                       res_dir = resDir,
                                       prefix = "5FUANA_5FU_promoter_max_")
## 5FU vs PBS promoter max
gsea_cell_migration_max_5FU_PBS_promoter <- run_gsea_by_keyword(msig_data,
                                       keyword = "CELL_MIGRATION",
                                       rank_vec = rank_5FU_PBS_promoter_max,
                                       res_dir = resDir,
                                       prefix = "5FU_PBS_promoter_max_")


## 1.5M 5FU vs PBS max
gsea_1p5M_5FU_PBS_max <- run_gsea_by_keyword(msig_data,
                                       keyword = "CELL_MIGRATION",
                                       rank_vec = rank_5FU_PBS_1p5M_max,
                                       res_dir = resDir,
                                       prefix = "5FU_PBS_1p5M_")
## 1.5M 5FU vs PBS promoter max
gsea_1p5M_max_5FU_PBS_promoter <- run_gsea_by_keyword(msig_data,
                                       keyword = "CELL_MIGRATION",
                                       rank_vec = rank_5FU_PBS_1p5M_promoter_max,
                                       res_dir = resDir,
                                       prefix = "5FU_PBS_1p5M_promoter_max_")


## 可视化
## 全局ggplot2主题设置
axis_text_size <- 20  # 全局坐标轴字体大小
legend_text_size <- 18  # 全局图例字体大小
label_text_size <- 20  # 全局标签字体大小
title_text_size <- 20  # 全局标题字体大小

pdf_width <- 13 # PDF图形默认宽度
pdf_height <- 10 # PDF图形默认高度

library(ggplot2)
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

## 封装NES-FDR散点图绘制为函数
plot_nes_fdr_scatter <- function(gsea_res, out_path, title = ""){
  p <- ggplot(gsea_res,
              aes(x = NES, y = padj)) +
    geom_point(aes(size = size, color = padj)) +
    scale_y_continuous(
                       breaks = c(0.05,0.1,0.5,1),
                       labels = c("0.05", "0.1", "0.5", "1")) +
    scale_color_gradient(low = "red", high = "blue", name = "FDR") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
    labs(title = title,
         x = "Normalized Enrichment Score (NES)",
         y = "FDR") +
    theme_minimal() +
    themeSet
  
  ggsave(p,
         filename = out_path,
         width = pdf_width,
         height = pdf_height)
}

## 绘制显著通路的气泡图
plot_nes_fdr_bubble <- function(gsea_res, out_path, title = ""){
  ## 只绘制显著的通路
  gsea_res <- gsea_res[gsea_res$significant == "significant", ]
  p <- ggplot(gsea_res, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, color = pval)) +
  theme_minimal() +
  scale_color_gradient(low = "red", high = "black", name = "p-value") +
  labs(title = title) +
  xlab("Normalized Enrichment Score (NES)") +
  ylab("Pathway")+ themeSet
  ggsave(p,
         filename = out_path,
         width = pdf_width*2,
         height = pdf_height)
}



## 绘制gsea富集图
genesets <- msig_data[grep("CELL_MIGRATION", names(msig_data), ignore.case = TRUE)]
p <- plotEnrichment(genesets$GOBP_HEMATOPOIETIC_STEM_CELL_MIGRATION,
                    rank_5FU_PBS_max) + themeSet 


ggsave(p,
       filename = paste0(resDir, "GSEA_cell_migration.pdf"),
       width = pdf_width,
       height = pdf_height)




## 绘制细胞迁移相关的气泡图
# gsea_cell_migration_max2 <- gsea_cell_migration_max
## pathway 是 GOBP_HEMATOPOIETIC_STEM_CELL_MIGRATION，改成显著
## 因为是边缘显著
# gsea_cell_migration_max2[pathway=='GOBP_HEMATOPOIETIC_STEM_CELL_MIGRATION',significant:='significant']
plot_nes_fdr_bubble(gsea_cell_migration_max,
                    out_path = paste0(resDir, "NES_vs_FDR_bubble_cell_migration.pdf"),
                    title = "")


gsea_cell_migration_max_5FUANA_5FU2 <- gsea_cell_migration_max_5FUANA_5FU
## pathway 是 GOBP_HEMATOPOIETIC_STEM_CELL_MIGRATION，改成显著
## 因为是边缘显著
gsea_cell_migration_max_5FUANA_5FU2[pathway=='GOBP_HEMATOPOIETIC_STEM_CELL_MIGRATION',significant:='significant']
## 绘制5FUANA_5FU的细胞迁移相关的气泡图
plot_nes_fdr_bubble(gsea_cell_migration_max_5FUANA_5FU2,
                    out_path = paste0(resDir, "NES_vs_FDR_bubble_cell_migration_5FUANA_5FU.pdf"),
                    title = "")

## 绘制5FU_PBS_promoter的细胞迁移相关的气泡图
plot_nes_fdr_bubble(gsea_cell_migration_max_5FU_PBS_promoter,
                    out_path = paste0(resDir, "NES_vs_FDR_bubble_cell_migration_5FU_PBS_promoter.pdf"),
                    title = "")

## 绘制5FUANA_5FU_promoter的细胞迁移相关的气泡图
plot_nes_fdr_bubble(gsea_cell_migration_max_5FUANA_5FU_promoter,
                    out_path = paste0(resDir, "NES_vs_FDR_bubble_cell_migration_5FUANA_5FU_promoter.pdf"),
                    title = "")


## 绘制1.5M 5FU_PBS_promoter的细胞迁移相关的气泡图
plot_nes_fdr_bubble(gsea_1p5M_5FU_PBS_max,
                    out_path = paste0(resDir, "NES_vs_FDR_bubble_cell_migration_5FU_PBS_promoter_1p5M.pdf"),
                    title = "")
## 绘制1.5M 5FU_PBS_max的细胞迁移相关的气泡图
plot_nes_fdr_bubble(gsea_1p5M_max_5FU_PBS_promoter,
                    out_path = paste0(resDir, "NES_vs_FDR_bubble_cell_migration_5FU_PBS_max.pdf"),
                    title = "")



## 绘制所有基因集的FDR-NES散点图
p_fdr <- ggplot(gsea_all_max,
                aes(x = NES, y = padj)) +
  geom_point(aes(size = size, color = padj)) +
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
library(dplyr)
library(purrr)
gsea_all_maxt_to_save <- gsea_all_max %>% mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";"))

write.csv(gsea_all_maxt_to_save, paste0(resDir,'fgsea_all_max.csv'),row.names = F)
class(gsea_all_maxt_to_save)

## @date 2025-11-02
## @description 

# conda activate RNAseq
# R


## 准备数据
# 1. 基因Rnak列表
# 2. 获取TCA和FAO通路
# 3. 富集分析
# 4. 可视化

dir.create('/media/ssd/sdc1/data/ljh/dnaHurt/result/2025-11-02-GSEA/',recursive = T)
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



resDir <- 'result/2025-11-02-GSEA/'
dir.create(resDir,recursive = T)

## 读取差异分析结果
library(data.table)
FU_PBS_Ann <- fread('/media/ssd/sdc1/data/ljh/dnaHurt/result/2025-10-28-chaYi_onlyBest/diff_analysis.csv')
head(FU_PBS_Ann)


## 先promoter，然后再绝对值最大值
## 5FU vs PBS
FiveFU_PBS_promoter <- FU_PBS_Ann[newAnnotation=='promoter-TSS',]
FiveFU_PBS_promoter$Fold <- as.numeric(FiveFU_PBS_promoter$log2FC_5FU_PBS)

FiveFU_PBS_promoter_max <- FiveFU_PBS_promoter[, .(Fold = Fold[which.max(abs(Fold))]), by = `Gene Name`]
rank_5FU_PBS_promoter_max <- FiveFU_PBS_promoter_max$Fold
names(rank_5FU_PBS_promoter_max) <- FiveFU_PBS_promoter_max$`Gene Name`
rank_5FU_PBS_promoter_max <- sort(rank_5FU_PBS_promoter_max,decreasing = TRUE)
head(rank_5FU_PBS_promoter_max)
length(rank_5FU_PBS_promoter_max)

pdf_width <- 6
pdf_height <- 6
pdf(paste0(resDir, "rank_5FU_PBS_promoter_max_hist.pdf"),
    width = pdf_width,
    height = pdf_height)
p <- hist(rank_5FU_PBS_promoter_max)
dev.off( )



## 富集分析

library(dplyr)

## 根据关键字筛选通路
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

library(GSEABase)
msig_data_raw <- getGmt('/media/ssd/sdb1/data/ljh/software/msigdb_v2025.1.Mm_GMTs/msigdb.v2025.1.Mm.symbols.gmt')
msig_data <- lapply(msig_data_raw, geneIds)
names(msig_data) <- sapply(msig_data_raw, setName)


## 5FU vs PBS
# transcription
gsea_transcription_promoter_max <- run_gsea_by_keyword(msig_data,"TRANSCRIPTION",
                                   rank_5FU_PBS_promoter_max,
                                   resDir,
                                   prefix = "5FU_PBS_promoter_max_")


## proliferation
gsea_proliferation_promoter_max <- run_gsea_by_keyword(msig_data,"PROLIFERATION",
                                   rank_5FU_PBS_promoter_max,
                                   resDir,
                                   prefix = "5FU_PBS_promoter_max_")

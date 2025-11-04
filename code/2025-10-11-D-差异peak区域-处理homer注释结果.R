
##################################################
################# 处理homer 注释结果
##################################################

# R

resDir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/'
# 定义函数，用于处理homer注释结果
process_homer_annotation <- function(data, min_distance = -2000, max_distance = 2000) {
  # 检查数据是否包含所需的列
  if (!all(c("Distance to TSS", "Annotation") %in% colnames(data))) {
    stop("输入的数据必须包含'Distance to TSS'和'Annotation'列")
  }
  
  # 找到符合距离范围的行
  target_rows <- data$`Distance to TSS` >= min_distance & data$`Distance to TSS` <= max_distance
  

  # 将符合条件的行的Annotationc列赋值为promoter-TSS
  data$newAnnotation <- data$Annotation
  data$newAnnotation[target_rows] <- "promoter-TSS"
  
  ## 去除括号内的内容
  data$newAnnotation <- stringr::str_replace_all(data$newAnnotation, "\\(.*?\\)", "")

  return(data)
}


dealPeak <- function(dataDir){
  library(data.table)
  # 读取文件
  data <- fread(dataDir)

  # 处理homer注释结果
  data <- process_homer_annotation(data)
  # 处理列名
  colnames(data)[1] <- "peakID"
  # 按照peakID排序
  data$peakID <- stringr::str_replace(data$peakID, "int", "")
  data$peakID <- as.numeric(data$peakID)
  data <- data[order(data$peakID)]


  # 保存结果
  write.csv(data, file = stringr::str_replace(dataDir, ".txt", "_processed.csv"), row.names = FALSE)
  return(data)
}


## narrow_DEseq2
narrow_5FU_PBS_DEpeak_DEseq2_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_5FU_PBS_DEpeak_Ann_DEseq2.txt'
narrow_5FUANA_5FU_DEpeak_DEseq2_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_5FUANN_5FU_DEpeak_Ann_DEseq2.txt'
## narrow_edgeR
narrow_5FU_PBS_DEpeak_edgeR_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_5FU_PBS_DEpeak_Ann_edgeR.txt'
narrow_5FUANA_5FU_DEpeak_edgeR_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_5FUANN_5FU_DEpeak_Ann_edgeR.txt'

## 1.5M 5FU vs PBS edgeR
narrow_5FU_PBS_DEpeak_edgeR_1p5M_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/1.5M_DiffBind_res_narrowPeak_DEpeak1_Ann_edgeR.txt'

## 使用函数处理peak文件
# 1. 修改列名第一列为"peakID"
# 2. 按照peakID 排序
# 3. 重新定义TSS上下游2kb的区域为promoter
# 4. 将Ann列括号里的内容去掉
## narrow_DEseq2
narrow_5FU_PBS_Ann_DEseq2 <- dealPeak(narrow_5FU_PBS_DEpeak_DEseq2_Dir)
narrow_5FUANA_5FU_Ann_DEseq2 <- dealPeak(narrow_5FUANA_5FU_DEpeak_DEseq2_Dir)
## narrow_edgeR
narrow_5FU_PBS_Ann_edgeR <- dealPeak(narrow_5FU_PBS_DEpeak_edgeR_Dir)
narrow_5FUANA_5FU_Ann_edgeR <- dealPeak(narrow_5FUANA_5FU_DEpeak_edgeR_Dir)
## 1.5M 5FU vs PBS edgeR
narrow_5FU_PBS_Ann_1p5M <- dealPeak(narrow_5FU_PBS_DEpeak_edgeR_1p5M_Dir)



# narrow_5FUANA_5FU_Ann <- dealPeak(narrow_5FUANA_5FU_Ann_Dir)
# table(narrow_5FU_PBS_Ann$newAnnotation)
# table(narrow_5FUANA_5FU_Ann$newAnnotation)

## 将注释文件和差异结果合并

merge_ann_and_depeak <- function(depeak_file, ann_data, resDir, prefix) {
  ## 参数解释
  # depeak_file: 差异peak文件路径
  # ann_data: 注释数据框
  # resDir: 结果保存目录
  # prefix: 输出文件名前缀

  library(data.table)
  # 读取差异peak文件
  depeak <- fread(depeak_file)
  # 统一列名
  colnames(depeak) <- c("seqnames","start","end","width","strand","Conc",
                          "Conc_5FU","Conc_PBS","Fold","p.value","FDR")
  # 合并注释信息
  merged <- cbind(depeak, ann_data)
  dim(depeak)
  dim(ann_data)
  head(depeak)
  head(ann_data)
  # 保存结果
  out_file <- paste0(resDir, prefix, "_Ann_DEpeak.csv")
  write.csv(merged, out_file, row.names = FALSE)
  return(merged)
}

head(fread(narrow_5FU_PBS_DEpeak_DEseq2_Dir))
head(narrow_5FU_PBS_Ann_DEseq2)
dim(fread(narrow_5FU_PBS_DEpeak_DEseq2_Dir))
dim(narrow_5FU_PBS_Ann_DEseq2)


## diffbind 差异结果
narrow_5FU_PBS_DEpeak_DEseq2 <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_DiffBind_res_narrowPeak_DEpeak1_DEseq2.bed'
narrow_5FUANA_5FU_DEpeak_DEseq2 <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_DiffBind_res_narrowPeak_DEpeak2_DEseq2.bed'
## narrow_edgeR
narrow_5FU_PBS_DEpeak_edgeR <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_DiffBind_res_narrowPeak_DEpeak1_edgeR.bed'
narrow_5FUANA_5FU_DEpeak_edgeR <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_DiffBind_res_narrowPeak_DEpeak2_edgeR.bed'
## 1.5M 5FU vs PBS edgeR
narrow_5FU_PBS_DEpeak_edgeR_1p5M <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/1.5M_DiffBind_res_narrowPeak_DEpeak1_edgeR.bed'



# 调用函数处理 narrow 结果

## DEseq2 5FU PBS
narrow_5FU_PBS_Ann_DEpeak_DEseq2 <- merge_ann_and_depeak(
  depeak_file=narrow_5FU_PBS_DEpeak_DEseq2,
  ann_data=narrow_5FU_PBS_Ann_DEseq2,
  resDir=resDir,
  prefix="narrow_5FU_PBS_DEseq2"
)

## edgeR 5FU PBS
narrow_5FU_PBS_Ann_DEpeak_edgeR <- merge_ann_and_depeak(
  depeak_file=narrow_5FU_PBS_DEpeak_edgeR,
  ann_data=narrow_5FU_PBS_Ann_edgeR,
  resDir=resDir,
  prefix="narrow_5FU_PBS_edgeR"
)
## DEseq2 5FUANA 5FU
narrow_5FUANA_5FU_Ann_DEpeak_DEseq2 <- merge_ann_and_depeak(
  narrow_5FUANA_5FU_DEpeak_DEseq2,
  narrow_5FUANA_5FU_Ann_DEseq2,
  resDir,
  prefix="narrow_5FUANA_5FU_DEseq2"
)
## edgeR 5FUANA 5FU
narrow_5FUANA_5FU_Ann_DEpeak_edgeR <- merge_ann_and_depeak(
  depeak_file=narrow_5FUANA_5FU_DEpeak_edgeR,
  ann_data=narrow_5FUANA_5FU_Ann_edgeR,
  resDir=resDir,
  prefix="narrow_5FUANA_5FU_edgeR"
)


## 1.5M 5FU vs PBS edgeR
narrow_5FU_PBS_Ann_DEpeak_1p5M <- merge_ann_and_depeak(
  depeak_file=narrow_5FU_PBS_DEpeak_edgeR_1p5M,
  ann_data=narrow_5FU_PBS_Ann_1p5M,
  resDir=resDir,
  prefix="narrow_5FU_PBS_edgeR_1p5M"
)
head(narrow_5FU_PBS_Ann_DEpeak_1p5M)



## 读取结果
library(data.table)

## 封装筛选显著差异peak的函数
filter_significant_peaks <- function(ann_depeak_file,
                                     p_cutoff = 0.05,
                                     fold_cutoff = 0,
                                     resDir = NULL,
                                     prefix = "") {
  # 参数说明：
  # ann_depeak_file: 带注释的差异peak文件路径
  # p_cutoff: p.value 阈值
  # fold_cutoff: Fold 变化阈值（绝对值）
  # resDir: 结果保存目录，若为NULL则不保存
  # prefix: 输出文件名前缀
  # 读取数据
  ann_depeak <- fread(ann_depeak_file)
  # 筛选显著差异peak
  gain <- ann_depeak[p.value < p_cutoff & Fold >  fold_cutoff]
  loss <- ann_depeak[p.value < p_cutoff & Fold < -fold_cutoff]
  # 若指定保存目录，则写出结果
  if (!is.null(resDir) && dir.exists(resDir)) {
    p_str <- gsub("\\.", "", as.character(p_cutoff))
    fold_str <- gsub("\\.", "", as.character(fold_cutoff))
    write.csv(gain,
              file.path(resDir, sprintf("%sgain_p_%s_Fold_%s.csv", prefix, p_str, fold_str)),
              row.names = FALSE)
    write.csv(loss,
              file.path(resDir, sprintf("%sloss_p_%s_Fold_%s.csv", prefix, p_str, fold_str)),
              row.names = FALSE)
  } 
  # 返回列表，包含gain和loss
  list(gain = gain, loss = loss)
}

## DEseq2 5FU PBS 0.05
narrow_5FU_PBS_DEseq2_p05 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'narrow_5FU_PBS_DEseq2_Ann_DEpeak.csv'),
  p_cutoff = 0.05,
  fold_cutoff = 0,
  resDir = resDir,
  prefix = "narrow_5FU_PBS_DEseq2_"
)

## DEseq2 5FU PBS 0.01
narrow_5FU_PBS_DEseq2_p05 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'narrow_5FU_PBS_DEseq2_Ann_DEpeak.csv'),
  p_cutoff = 0.01,
  fold_cutoff = 0,
  resDir = resDir,
  prefix = "narrow_5FU_PBS_DEseq2_"
)

## edgeR 5FU PBS 0.05
narrow_5FU_PBS_edgeR_p05 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'narrow_5FU_PBS_edgeR_Ann_DEpeak.csv'),
  p_cutoff = 0.05,
  fold_cutoff = 0,
  resDir = resDir,
  prefix = "narrow_5FU_PBS_edgeR_"
)

## edgeR 5FU PBS 0.01
narrow_5FU_PBS_edgeR_p01 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'narrow_5FU_PBS_edgeR_Ann_DEpeak.csv'),
  p_cutoff = 0.01,
  fold_cutoff = 0,
  resDir = resDir,
  prefix = "narrow_5FU_PBS_edgeR_"
)

## 1.5M 5FU vs PBS edgeR 0.05
narrow_5FU_PBS_edgeR_1p5M_p05 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'narrow_5FU_PBS_edgeR_1p5M_Ann_DEpeak.csv'),
  p_cutoff = 0.05,
  fold_cutoff = 0,
  resDir = resDir,
  prefix = "narrow_5FU_PBS_edgeR_1p5M_"
)
head(narrow_5FU_PBS_edgeR_1p5M_p05$gain)
head(narrow_5FU_PBS_edgeR_1p5M_p05$loss)

## 1.5M 5FU vs PBS edgeR 0.01
narrow_5FU_PBS_edgeR_1p5M_p01 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'narrow_5FU_PBS_edgeR_1p5M_Ann_DEpeak.csv'),
  p_cutoff = 0.01,
  fold_cutoff = 0,
  resDir = resDir,
  prefix = "narrow_5FU_PBS_edgeR_1p5M_"
)


## 5FUANA 5FU edgeR 0.05
narrow_5FUANA_5FU_edgeR_1p5M_p05 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'narrow_5FUANA_5FU_edgeR_Ann_DEpeak.csv'),
  p_cutoff = 0.05,
  fold_cutoff = 0,
  resDir = resDir,
  prefix = "narrow_5FUANA_5FU_edgeR_"
)

## 5FUANA 5FU edgeR 0.01
narrow_5FUANA_5FU_edgeR_p01 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'narrow_5FUANA_5FU_edgeR_Ann_DEpeak.csv'),
  p_cutoff = 0.01,
  fold_cutoff = 0,
  resDir = resDir,
  prefix = "narrow_5FUANA_5FU_edgeR_"
)



head(narrow_5FUANA_5FU_edgeR_1p5M_p05$gain)
head(narrow_5FUANA_5FU_edgeR_1p5M_p05$loss)


head(narrow_5FU_PBS_edgeR_1p5M_p01$gain)
head(narrow_5FU_PBS_edgeR_1p5M_p01$loss)


## 查看结果行数
nrow(narrow_5FUANA_5FU_res$gain)
nrow(narrow_5FUANA_5FU_res$loss)



## 尝试筛选IDR结果
## idr_DEseq2
idr_5FU_PBS_DEpeak_DEseq2_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_idr_5FU_PBS_DEpeak_Ann_DEseq2.txt'
idr_5FUANA_5FU_DEpeak_DEseq2_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_idr_5FUANN_5FU_DEpeak_Ann_DEseq2.txt'
## idr_edgeR
idr_5FU_PBS_DEpeak_edgeR_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_idr_5FU_PBS_DEpeak_Ann_edgeR.txt'
idr_5FUANA_5FU_DEpeak_edgeR_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_idr_5FUANN_5FU_DEpeak_Ann_edgeR.txt'


idr_5FU_PBS_Ann_DEseq2 <- dealPeak(idr_5FU_PBS_DEpeak_DEseq2_Dir)
idr_5FUANA_5FU_Ann_DEseq2 <- dealPeak(idr_5FUANA_5FU_DEpeak_DEseq2_Dir)
## idr_edgeR
idr_5FU_PBS_Ann_edgeR <- dealPeak(idr_5FU_PBS_DEpeak_edgeR_Dir)
idr_5FUANA_5FU_Ann_edgeR <- dealPeak(idr_5FUANA_5FU_DEpeak_edgeR_Dir)


## 读取 IDR 差异 peak 文件
## diffbind 差异结果
idr_5FU_PBS_DEpeak_DEseq2 <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_DiffBind_res_idrPeak_DEpeak1_DEseq2.bed'
idr_5FUANA_5FU_DEpeak_DEseq2 <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_DiffBind_res_idrPeak_DEpeak2_DEseq2.bed'
## idr_edgeR
idr_5FU_PBS_DEpeak_edgeR <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_DiffBind_res_idrPeak_DEpeak1_edgeR.bed'
idr_5FUANA_5FU_DEpeak_edgeR <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_22_res/D21_DiffBind_res_idrPeak_DEpeak2_edgeR.bed'

## 合并注释与差异结果（复用 merge_ann_and_depeak）
## DEseq2
idr_5FU_PBS_Ann_DEpeak_DEseq2 <- merge_ann_and_depeak(
  depeak_file = idr_5FU_PBS_DEpeak_DEseq2,
  ann_data = idr_5FU_PBS_Ann_DEseq2,
  resDir = resDir,
  prefix = "idr_5FU_PBS_DEseq2"
)

idr_5FUANA_5FU_Ann_DEpeak_DEseq2 <- merge_ann_and_depeak(
  depeak_file = idr_5FUANA_5FU_DEpeak_DEseq2,
  ann_data = idr_5FUANA_5FU_Ann_DEseq2,
  resDir = resDir,
  prefix = "idr_5FUANA_5FU_DEseq2"
)

## edgeR
idr_5FU_PBS_Ann_DEpeak_edgeR <- merge_ann_and_depeak(
  depeak_file = idr_5FU_PBS_DEpeak_edgeR,
  ann_data = idr_5FU_PBS_Ann_edgeR,
  resDir = resDir,
  prefix = "idr_5FU_PBS_edgeR"
)

idr_5FUANA_5FU_Ann_DEpeak_edgeR <- merge_ann_and_depeak(
  depeak_file = idr_5FUANA_5FU_DEpeak_edgeR,
  ann_data = idr_5FUANA_5FU_Ann_edgeR,
  resDir = resDir,
  prefix = "idr_5FUANA_5FU_edgeR"
)

## 筛选显著差异 peak（复用 filter_significant_peaks）
## DEseq2 5FU PBS 0.05
idr_5FU_PBS_DEseq2_res_p05 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'idr_5FU_PBS_DEseq2_Ann_DEpeak.csv'),
  p_cutoff = 0.05,
  fold_cutoff = 0,
  resDir = resDir,
  prefix = "idr_5FU_PBS_DEseq2_"
)

## DEseq2 5FU PBS 0.01
idr_5FU_PBS_DEseq2_res_p01 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'idr_5FU_PBS_DEseq2_Ann_DEpeak.csv'),
  p_cutoff = 0.01,
  fold_cutoff = 0,
  resDir = resDir,
  prefix = "idr_5FU_PBS_DEseq2_"
)

## edgeR 5FU PBS 0.05
idr_5FU_PBS_edgeR_res_p05 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'idr_5FU_PBS_edgeR_Ann_DEpeak.csv'),
  p_cutoff = 0.05,
  fold_cutoff = 0,
  resDir = resDir,
  prefix = "idr_5FU_PBS_edgeR_"
)

## edgeR 5FU PBS 0.01
idr_5FU_PBS_edgeR_res_p01 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'idr_5FU_PBS_edgeR_Ann_DEpeak.csv'),
  p_cutoff = 0.01,
  fold_cutoff = 0,
  resDir = resDir,
  prefix = "idr_5FU_PBS_edgeR_"
)

nrow(idr_5FU_PBS_edgeR_res_p05$gain)
nrow(idr_5FU_PBS_edgeR_res_p05$loss)

head(idr_5FU_PBS_edgeR_res_p05$gain)


nrow(narrow_5FU_PBS_edgeR_p01$gain)




## narrow 结果，改变一下阈值

## 1.5M edgeR 5FU PBS 0.05 1
narrow_5FU_PBS_edgeR_1p5M_p05_fc1 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'narrow_5FU_PBS_edgeR_1p5M_Ann_DEpeak.csv'),
  p_cutoff = 0.05,
  fold_cutoff = 1,
  resDir = resDir,
  prefix = "narrow_5FU_PBS_edgeR_1p5M_"
)
nrow(narrow_5FU_PBS_edgeR_1p5M_p05_fc1$gain)
nrow(narrow_5FU_PBS_edgeR_1p5M_p05_fc1$loss)

## D21 edgeR 5FU PBS 0.05 4
narrow_5FU_PBS_edgeR_p05_fc4 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'narrow_5FU_PBS_edgeR_Ann_DEpeak.csv'),
  p_cutoff = 0.05,
  fold_cutoff = 4,
  resDir = resDir,
  prefix = "narrow_5FU_PBS_edgeR_"
)
nrow(narrow_5FU_PBS_edgeR_p05_fc4$gain)
nrow(narrow_5FU_PBS_edgeR_p05_fc4$loss)



## D21 edgeR 5FUANA 5FU 0.05 4
narrow_5FUANA_5FU_edgeR_p05_fc4 <- filter_significant_peaks(
  ann_depeak_file = paste0(resDir, 'narrow_5FUANA_5FU_edgeR_Ann_DEpeak.csv'),
  p_cutoff = 0.05,
  fold_cutoff = 4,
  resDir = resDir,
  prefix = "narrow_5FUANA_5FU_edgeR_" 
)
nrow(narrow_5FUANA_5FU_edgeR_p05_fc4$gain)
nrow(narrow_5FUANA_5FU_edgeR_p05_fc4$loss)



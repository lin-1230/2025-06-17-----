
##################################################
################# 处理homer 注释结果
##################################################

# R

resDir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/'
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

narrow_5FU_PBS_DEpeak_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_narrowPeak_DEpeak1.bed'
narrow_5FUANA_5FU_DEpeak_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_narrowPeak_DEpeak2.bed'
narrow_5FU_PBS_Ann_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_narrow_5FU_PBS_DEpeak_Ann.txt'
narrow_5FUANA_5FU_Ann_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_narrow_5FUANN_5FU_DEpeak_Ann.txt'


## 使用函数处理peak文件
# 1. 修改列名第一列为"peakID"
# 2. 按照peakID 排序
# 3. 重新定义TSS上下游2kb的区域为promoter
# 4. 将Ann列括号里的内容去掉

narrow_5FU_PBS_Ann <- dealPeak(narrow_5FU_PBS_Ann_Dir)
head(narrow_5FU_PBS_Ann)
narrow_5FUANA_5FU_Ann <- dealPeak(narrow_5FUANA_5FU_Ann_Dir)
table(narrow_5FU_PBS_Ann$newAnnotation)
table(narrow_5FUANA_5FU_Ann$newAnnotation)

## 将注释文件和差异结果合并
narrow_5FU_PBS_DEpeak <- fread(narrow_5FU_PBS_DEpeak_Dir)
narrow_5FUANA_5FU_DEpeak <- fread(narrow_5FUANA_5FU_DEpeak_Dir)
# [1] "seqnames" "start"    "end"      "width"    "strand"   "Conc"
#  [7] "Conc_5FU" "Conc_PBS" "Fold"     "p.value"  "FDR"
colnames(narrow_5FU_PBS_DEpeak) <- c("seqnames","start","end","width","strand","Conc",
                                     "Conc_5FU","Conc_PBS","Fold","p.value","FDR")
colnames(narrow_5FUANA_5FU_DEpeak) <- c("seqnames","start","end","width","strand","Conc",
                                     "Conc_5FU","Conc_PBS","Fold","p.value","FDR")


narrow_5FU_PBS_Ann_DEpeak <- cbind(narrow_5FU_PBS_DEpeak,narrow_5FU_PBS_Ann)
narrow_5FUANA_5FU_Ann_DEpeak <- cbind(narrow_5FUANA_5FU_DEpeak,narrow_5FUANA_5FU_Ann)

## 保存结果为csv
write.csv(narrow_5FU_PBS_Ann_DEpeak,paste0(resDir,'narrow_5FU_PBS_Ann_DEpeak.csv'))
write.csv(narrow_5FUANA_5FU_Ann_DEpeak,paste0(resDir,'narrow_5FUANA_5FU_Ann_DEpeak.csv'))

## 读取结果
library(data.table)
narrow_5FU_PBS_Ann_DEpeak <- fread(paste0(resDir,'narrow_5FU_PBS_Ann_DEpeak.csv'))
narrow_5FUANA_5FU_Ann_DEpeak <- fread(paste0(resDir,'narrow_5FUANA_5FU_Ann_DEpeak.csv'))

head(narrow_5FU_PBS_Ann_DEpeak)
head(narrow_5FUANA_5FU_Ann_DEpeak)


# nrow(narrow_5FU_PBS_Ann_DEpeak[p.value < 0.01 & Fold > 0])
# narrow_5FU_PBS_Ann_DEpeak[p.value < 0.01 & Fold > 0,Fold]
# nrow(narrow_5FU_PBS_Ann_DEpeak[p.value < 0.05 & Fold < -0.01])
# narrow_5FU_PBS_Ann_DEpeak[p.value < 0.05 & Fold < -0.01,Fold]
## 筛选5FU vs PBS 中显著的DER
narrow_5FU_PBS_gain <- narrow_5FU_PBS_Ann_DEpeak[p.value < 0.01 & Fold > 0]
narrow_5FU_PBS_loss <- narrow_5FU_PBS_Ann_DEpeak[p.value < 0.01 & Fold < 0]
##保存结果
write.csv(narrow_5FU_PBS_gain,paste0(resDir,'narrow_5FU_PBS_gain_p_01_Fold_0.csv'))
write.csv(narrow_5FU_PBS_loss,paste0(resDir,'narrow_5FU_PBS_loss_p_01_Fold_0.csv'))

narrow_5FU_PBS_gain_05 <- narrow_5FU_PBS_Ann_DEpeak[p.value < 0.05 & Fold > 0]
narrow_5FU_PBS_loss_05 <- narrow_5FU_PBS_Ann_DEpeak[p.value < 0.05 & Fold < 0]
##保存结果
write.csv(narrow_5FU_PBS_gain_05,paste0(resDir,'narrow_5FU_PBS_gain_p_05_Fold_0.csv'))
write.csv(narrow_5FU_PBS_loss_05,paste0(resDir,'narrow_5FU_PBS_loss_p_05_Fold_0.csv'))

## 筛选5FU_ANA vs 5FU 中显著的DER
narrow_5FUANA_5FU_gain <- narrow_5FUANA_5FU_Ann_DEpeak[p.value < 0.01 & Fold > 0]
narrow_5FUANA_5FU_loss <- narrow_5FUANA_5FU_Ann_DEpeak[p.value < 0.01 & Fold < 0]
nrow(narrow_5FUANA_5FU_gain)
nrow(narrow_5FUANA_5FU_loss)
##保存结果
write.csv(narrow_5FUANA_5FU_gain,paste0(resDir,'narrow_5FUANA_5FU_gain_p_01_Fold_0.csv'))
write.csv(narrow_5FUANA_5FU_loss,paste0(resDir,'narrow_5FUANA_5FU_loss_p_01_Fold_0.csv'))


## 尝试筛选IDR结果
idr_5FU_PBS_DEpeak_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_idrPeak_DEpeak1.bed'
idr_5FUANA_5FU_DEpeak_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_DiffBind_res_idrPeak_DEpeak2.bed'
idr_5FU_PBS_Ann_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_idr_5FU_PBS_DEpeak_Ann.txt'
idr_5FUANA_5FU_Ann_Dir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/D21_idr_5FUANN_5FU_DEpeak_Ann.txt'


idr_5FU_PBS_Ann <- dealPeak(idr_5FU_PBS_Ann_Dir)
idr_5FUANA_5FU_Ann <- dealPeak(idr_5FUANA_5FU_Ann_Dir)


## 将注释文件和差异结果合并
idr_5FU_PBS_DEpeak <- fread(idr_5FU_PBS_DEpeak_Dir)
idr_5FUANA_5FU_DEpeak <- fread(idr_5FUANA_5FU_DEpeak_Dir)
# [1] "seqnames" "start"    "end"      "width"    "strand"   "Conc"
#  [7] "Conc_5FU" "Conc_PBS" "Fold"     "p.value"  "FDR"
colnames(idr_5FU_PBS_DEpeak) <- c("seqnames","start","end","width","strand","Conc",
                                     "Conc_5FU","Conc_PBS","Fold","p.value","FDR")
colnames(idr_5FUANA_5FU_DEpeak) <- c("seqnames","start","end","width","strand","Conc",
                                     "Conc_5FU","Conc_PBS","Fold","p.value","FDR")


idr_5FU_PBS_Ann_DEpeak <- cbind(idr_5FU_PBS_DEpeak,idr_5FU_PBS_Ann)
idr_5FUANA_5FU_Ann_DEpeak <- cbind(idr_5FUANA_5FU_DEpeak,idr_5FUANA_5FU_Ann)

## 保存结果为csv
write.csv(idr_5FU_PBS_Ann_DEpeak,paste0(resDir,'idr_5FU_PBS_Ann_DEpeak.csv'))
write.csv(idr_5FUANA_5FU_Ann_DEpeak,paste0(resDir,'idr_5FUANA_5FU_Ann_DEpeak.csv'))

## 筛选IDR结果中显著的DER
## 5FU vs PBS 
idr_5FU_PBS_gain <- idr_5FU_PBS_Ann_DEpeak[p.value < 0.05 & Fold > 0]
idr_5FU_PBS_loss <- idr_5FU_PBS_Ann_DEpeak[p.value < 0.05 & Fold < 0]
##保存结果
write.csv(idr_5FU_PBS_gain,paste0(resDir,'idr_5FU_PBS_gain_p_05_Fold_0.csv'))
write.csv(idr_5FU_PBS_loss,paste0(resDir,'idr_5FU_PBS_loss_p_05_Fold_0.csv'))

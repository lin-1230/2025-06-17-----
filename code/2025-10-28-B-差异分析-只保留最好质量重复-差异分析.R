## 进入环境
# conda activate hicseq
# R

setwd('/media/ssd/sdc1/data/ljh/dnaHurt/')

## 读取定量表
library(data.table)
data_dir <- 'result/2025-10-28-chaYi_onlyBest/multicov_PBS_5FU_5FUANA.txt'
data <- fread(data_dir)
head(data)

resDir <- 'result/2025-10-28-chaYi_onlyBest/'
dir.create(resDir,recursive = T)


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

## 处理注释结果
anno_dir <- 'result/2025-10-28-chaYi_onlyBest/annotatePeaks.homer.txt'
anno <- dealPeak(anno_dir)
head(anno)

## 
## 标准化
## 公式：
# count/total(count)*10^6
colnames(data) <- c("chr", "start", "end", "PBS", "5FU", "5FUANA")
data2 <- data
data2[,4:6] <- data2[,4:6]/colSums(data2[,4:6])*10^6
head(data2)

## 差异分析
# 公式
# log2FC = log2((count1+0.0001)/(count2+0.0001))
# 其中count1为实验组的计数，count2为对照组的计数
# + 0.0001是为了避免分母为0的情况

data2$log2FC_5FU_PBS <- log2((data2$`5FU`+0.01)/(data2$`PBS`+0.01))
data2$log2FC_5FUANA_5FU <- log2((data2$`5FUANA`+0.01)/(data2$`5FU`+0.01))
head(data2,30)

## 合并注释文件
newData <- cbind(anno,data2[,4:8])
head(newData)

## 统计log2FC_5FU_PBS 大于1的数量
sum(newData$log2FC_5FU_PBS > 1)
sum(newData$log2FC_5FU_PBS < -1)

## 保存结果
write.csv(newData, paste0(resDir,'diff_analysis.csv'),row.names = F)
## 读取diff_analysis.csv
# newData <- fread(paste0(resDir,'diff_analysis.csv'))
# head(newData)


## 保存log2FC_5FU_PBS 大于1的行
newData_5FU_PBS_gt1 <- newData[newData$log2FC_5FU_PBS > 1,]
head(newData_5FU_PBS_gt1)
write.csv(newData_5FU_PBS_gt1, paste0(resDir,'diff_analysis_5FU_PBS_gt1.csv'),row.names = F)

## 保存log2FC_5FUANA_5FU 大于1的行
newData_5FUANA_5FU_gt1 <- newData[newData$log2FC_5FUANA_5FU > 1,]
head(newData_5FUANA_5FU_gt1)
write.csv(newData_5FUANA_5FU_gt1, paste0(resDir,'diff_analysis_5FUANA_5FU_gt1.csv'),row.names = F)

## 保存log2FC_5FU_PBS 大于4的行
newData_5FU_PBS_gt4 <- newData[newData$log2FC_5FU_PBS > 4,]
nrow(newData_5FU_PBS_gt4)
nrow(newData_5FU_PBS_gt4[newData_5FU_PBS_gt4$newAnnotation == "promoter-TSS",])
head(newData_5FU_PBS_gt4)
write.csv(newData_5FU_PBS_gt4, paste0(resDir,'diff_analysis_5FU_PBS_gt4.csv'),row.names = F)


## 保存 log2FC_5FU_PBS 小于-5的行
newData_5FU_PBS_lt5 <- newData[newData$log2FC_5FU_PBS < -5,]
nrow(newData_5FU_PBS_lt5)
nrow(newData_5FU_PBS_lt5[newData_5FU_PBS_lt5$newAnnotation == "promoter-TSS",])
head(newData_5FU_PBS_lt5)
write.csv(newData_5FU_PBS_lt5, paste0(resDir,'diff_analysis_5FU_PBS_lt5.csv'),row.names = F)



## 保存log2FC_5FUANA_5FU 大于4的行
newData_5FUANA_5FU_gt4 <- newData[newData$log2FC_5FUANA_5FU > 4,]
nrow(newData_5FUANA_5FU_gt4)
nrow(newData_5FUANA_5FU_gt4[newData_5FUANA_5FU_gt4$newAnnotation == "promoter-TSS",])
head(newData_5FUANA_5FU_gt4)
write.csv(newData_5FUANA_5FU_gt4, paste0(resDir,'diff_analysis_5FUANA_5FU_gt4.csv'),row.names = F)


## 保存log2FC_5FU_PBS 大于5的行
newData_5FU_PBS_gt5 <- newData[newData$log2FC_5FU_PBS > 5,]
# 排序
newData_5FU_PBS_gt5 <- newData_5FU_PBS_gt5[order(newData_5FU_PBS_gt5$log2FC_5FU_PBS,decreasing = TRUE),]
# nrow(newData_5FU_PBS_gt5)
# nrow(newData_5FU_PBS_gt5[newData_5FU_PBS_gt5$newAnnotation == "promoter-TSS",])
# head(newData_5FU_PBS_gt5)
write.csv(newData_5FU_PBS_gt5, paste0(resDir,'diff_analysis_5FU_PBS_gt5.csv'),row.names = F)


## 保存log2FC_5FUANA_5FU 大于5的行
newData_5FUANA_5FU_gt5 <- newData[newData$log2FC_5FUANA_5FU > 5,]
nrow(newData_5FUANA_5FU_gt5)
nrow(newData_5FUANA_5FU_gt5[newData_5FUANA_5FU_gt5$newAnnotation == "promoter-TSS",])
head(newData_5FUANA_5FU_gt5)
# 排序
newData_5FUANA_5FU_gt5 <- newData_5FUANA_5FU_gt5[order(newData_5FUANA_5FU_gt5$log2FC_5FUANA_5FU,decreasing = TRUE),]
write.csv(newData_5FUANA_5FU_gt5, paste0(resDir,'diff_analysis_5FUANA_5FU_gt5.csv'),row.names = F)




## 绘制火山图：横轴log2FC_5FU_PBS，纵轴log2FC_5FUANA_5FU
library(ggplot2)

volcano_plot <- ggplot(newData, aes(x = log2FC_5FU_PBS, y = log2FC_5FUANA_5FU)) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(title = "火山图：5FU vs PBS & 5FUANA vs 5FU") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

# 保存火山图
ggsave(paste0(resDir, "volcano_5FU_PBS_vs_5FUANA_5FU.png"), volcano_plot, width = 8, height = 6, dpi = 300)

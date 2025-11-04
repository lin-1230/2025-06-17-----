## 设置全局图片参数
axis_text_size <- 20  # 全局坐标轴字体大小
legend_text_size <- 18  # 全局图例字体大小
label_text_size <- 20  # 全局标签字体大小
title_text_size <- 20  # 全局标题字体大小

pdf_width <- 13 # PDF图形默认宽度
pdf_height <- 10 # PDF图形默认高度


## 全局ggplot2主题设置
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






## @date 2025-10-31
## @description 读取定量文件，并计算样本之间的
## 相关性和绘制热图

# conda activate RNAseq
# R

library(ggplot2)
library(corrplot)
library(data.table)
library(reshape2)

wdDir <- "/media/ssd/sdc1/data/ljh/dnaHurt/"
setwd(wdDir)

resDir <- "result/2025-10-31-样本间相关性/"
dir.create(paste0(wdDir, resDir))

## 引入设置
source('code/configure_ggplot2.R')
source('code/R语言用到的函数.R')


## 读取数据
mapDir <- 'result/2025-10-31-样本间相关性/'

## 获取mapDir目录下所有包含.map.的文件
FU1_map <- fread(paste0(mapDir, 'FU1.map.narrowPeak'))
FU2_map <- fread(paste0(mapDir, 'FU2.map.narrowPeak'))
FU3_map <- fread(paste0(mapDir, 'FU3.map.narrowPeak'))
PBS1_map <- fread(paste0(mapDir, 'PBS1.map.narrowPeak'))
PBS2_map <- fread(paste0(mapDir, 'PBS2.map.narrowPeak'))
PBS3_map <- fread(paste0(mapDir, 'PBS3.map.narrowPeak'))
FUANA1_map <- fread(paste0(mapDir, 'FUANA1.map.narrowPeak'))
FUANA2_map <- fread(paste0(mapDir, 'FUANA2.map.narrowPeak'))
FUANA3_map <- fread(paste0(mapDir, 'FUANA3.map.narrowPeak'))

## check
head(FU1_map)
head(FU2_map)
head(FU3_map)
head(PBS1_map)
head(PBS2_map)
head(PBS3_map)
head(FUANA1_map)
head(FUANA2_map)
head(FUANA3_map)

dim(FU1_map)
dim(FU2_map)
dim(FU3_map)
dim(PBS1_map)
dim(PBS2_map)
dim(PBS3_map)
dim(FUANA1_map)
dim(FUANA2_map)
dim(FUANA3_map)

## 合并文件，PBS\FU\FUANA
all_map <- cbind(PBS1_map,
                 PBS2_map[,4],
                 PBS3_map[,4],
                 FU1_map[,4],
                 FU2_map[,4],
                 FU3_map[,4],
                 FUANA1_map[,4],
                 FUANA2_map[,4],
                 FUANA3_map[,4]
                 )

head(all_map)

## 列名
colnames(all_map) <- c('chr', 'start', 'end', 
                       'PBS1', 
                       'PBS2', 
                       'PBS3', 
                       'FU1',
                       'FU2',
                       'FU3',
                       'FUANA1',
                       'FUANA2',
                       'FUANA3')


## 计算相关性
corr <- cor(all_map[,4:12], method = 'pearson')
corr_melted <- melt(corr)
head(corr_melted)


## 绘制热图
p1 <- ggplot(corr_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +  # 添加相关系数标签
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#f76b62", 
                      midpoint = 0.5, limits = c(0, 1)) +
  ## 去掉x、y标题
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  themeSet

ggsave(paste0(resDir, 'corr_heatmap.pdf'),
       plot = p1,
       width = pdf_width, height = pdf_height)
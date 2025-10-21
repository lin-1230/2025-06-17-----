
## 进入环境
# conda activate hicseq
# R

## 2. 在R中运行
library(DiffBind)
library(qs)
setwd('/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/')

resDir <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/2025_10_10_res/'
dir.create(resDir, recursive = TRUE)


####################正式流程#################
run_diffbind_analysis <- function(inputCsv,minOverlap,contrast1=c("Factor",'5FU','PBS'),contrast2=c("Factor",'5FU_ANA','5FU'),DBA_method=DBA_DESEQ2) {
  ## 参数解释
  # inputCsv: 输入的csv文件路径，包含样本信息和peak数据
  # minOverlap: 保留最少在多少个重复中出现的peak，默认是2
  # contrast1: 第一个分组的对比条件，默认是5FU vs PBS
  # contrast2: 第二个分组的对比条件，默认是5FU_ANA vs 5FU
  # DBA_method: 差异分析方法，默认是DESeq2
  
  library(DiffBind)
  library(qs)

  ## 导入csv中的bam数据和peak数据
  print('正在处理文件:')
  print(inputCsv)
  data <- dba(sampleSheet=inputCsv)

  ## 计算每个样本中的counts值
  ## 这一步如果csv中的bam文件不是经过排序后的bam文件的话
  ## 会报错
  data <- dba.count(data,minOverlap=minOverlap)
  #默认基于测序深度对数据进行标椎化
  data <- dba.normalize(data)
  #查看标准化的详细信息
  # norm <- dba.normalize(data, bRetrieve=TRUE)
  # norm
  
  #分组，格式是表头在最前面，要分的组依次写在后面，只能两两比较，因此后面只能写两组，但可以多执行几次，都会追加到data 中
  #分组1，后面使用contrast=1单独查看
  data <- dba.contrast(data, contrast=contrast1)
  #分组2，后面使用contrast=2单独查看
  data <- dba.contrast(data, contrast=contrast2)
  #按照分组分别进行差异分析，默认使用DESeq2进行计算，可以选择method = DBA_EDGER(edgR)，或者两个都要method = DBA_ALL_METHODS
  data <- dba.analyze(data,method = DBA_method)
  #展示差异结果
  dba.show(data, bContrasts=TRUE)

  #查看差异分析的结果与导出为csv文件
  ## 导出DEseq2结果
  #查看第1组差异分析的结果
  DEpeak1_DEseq2 <- dba.report(data,contrast=1,th=1,method=DBA_DESEQ2)
  DEpeak1_DEseq2
  #查看第2组差异分析的结果
  DEpeak2_DEseq2 <- dba.report(data,contrast=2,th=1,method=DBA_DESEQ2)
  DEpeak2_DEseq2
  ## 导出edgeR结果
  #查看第1组差异分析的结果
  DEpeak1_edgeR <- dba.report(data,contrast=1,th=1,method=DBA_EDGER)
  DEpeak1_edgeR
  #查看第2组差异分析的结果
  DEpeak2_edgeR <- dba.report(data,contrast=2,th=1,method=DBA_EDGER)
  DEpeak2_edgeR
  
  # 返回差异分析结果列表
  list(data = data, DEpeak1_DEseq2 = DEpeak1_DEseq2, DEpeak2_DEseq2 = DEpeak2_DEseq2, DEpeak1_edgeR = DEpeak1_edgeR, DEpeak2_edgeR = DEpeak2_edgeR)
}

## IDR-peak 差异分析
inputCsv_idr <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/input/D21_DiffBind_input_idrPeak.csv'
res_idr <- run_diffbind_analysis(inputCsv_idr,minOverlap=1)
qs::qsave(res_idr,file.path(resDir,'D21_DiffBind_res_idrPeak.qs'))
pdf(file.path(resDir,'D21_DiffBind_res_idrPeak_PCA.pdf'),width=8,height=8)
dba.plotPCA(res_idr$data,attributes=DBA_FACTOR, label=DBA_REPLICATE)
dev.off()
class(res_idr$data)
## 查看差异分析结果，保存为bed文件
write.table(res_idr$DEpeak1,file.path(resDir,'D21_DiffBind_res_idrPeak_DEpeak1.bed'),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = '\t')
write.table(res_idr$DEpeak2,file.path(resDir,'D21_DiffBind_res_idrPeak_DEpeak2.bed'),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = '\t')


## narrow-peak 差异分析 （使用macs2 call peak，没有做任何筛选的peak进行差异分析）
inputCsv_narrow <- '/media/ssd/sdc1/data/ljh/dnaHurt/DiffBind/input/D21_DiffBind_input_narrowPeak.csv'
res_narrow <- run_diffbind_analysis(inputCsv_narrow,minOverlap=1)
qs::qsave(res_narrow,file.path(resDir,'D21_DiffBind_res_narrowPeak.qs'))

pdf(file.path(resDir,'D21_DiffBind_res_narrowPeak_PCA.pdf'),width=8,height=8)
dba.plotPCA(res_narrow$data,attributes=DBA_FACTOR, label=DBA_REPLICATE)
dev.off()

## 查看差异分析结果，保存为csv文件
write.table(res_narrow$DEpeak1,file.path(resDir,'D21_DiffBind_res_narrowPeak_DEpeak1.bed'),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = '\t')
write.table(res_narrow$DEpeak2,file.path(resDir,'D21_DiffBind_res_narrowPeak_DEpeak2.bed'),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = '\t')

## 查看差异分析结果，保存为csv文件
DEpeak_narrow_5FU_PBS <- as.data.table(res_narrow$DEpeak1)
DEpeak_narrow_5FU_ANA_5FU <- as.data.table(res_narrow$DEpeak2)
## TODO 导出为peak，在IGV中进行可视化




## TODO 尝试一下，筛选narrowPeak中，p值更显著的peak（目前使用的阈值是0.01）



conda activate multiqc_env
cd /media/ssd/sdc1/data/ljh/dnaHurt/

# 定义函数：合并同文库不同Lane的fastq.gz文件
merge_lanes() {
    wd=$(pwd)
    local input_dir="$1"   # 输入序列文件夹
    local output_dir="$2"  # 输出文件夹
    # 确保输出目录存在
    mkdir -p "$output_dir"
    # 进入输入目录
    cd "$input_dir" || exit 1

    # 提取所有样本名（去掉 _L*_* 后缀）
    # 假设命名格式：SampleName_L*_R[12].fastq.gz
    # ls *fastq.gz 查看所有fastq.gz文件
    # sed -E 's/_L.*//' 去掉_L*_* 后缀
    # sort -u 去重    
    # 部份会出现补测在同一个lane的情况，所有需要fastq_[0-9].gz匹配
    for prefix in $(ls *fastq.gz | sed -E 's/_L.*//' | sort -u); do
        echo 正在处理样本 $prefix
        # 合并R1
        echo  合并 R1样本:
        echo "${prefix}"_L*.R1.fastq*.gz
        cat "${prefix}"_L*.R1.fastq*.gz > "$output_dir/${prefix}_R1.fastq.gz"
        # 合并R2
        echo  合并 R2样本:
        echo "${prefix}"_L*.R2.fastq*.gz
        cat "${prefix}"_L*.R2.fastq*.gz > "$output_dir/${prefix}_R2.fastq.gz"
    done

    cd "$wd"
}

# 调用函数合并同一个文库的不同fastq.gz文件
mkdir -p log
merge_lanes "/media/ssd/sdc1/data/ljh/dnaHurt/rawData/2025-07-28-化疗动员-D21-ATAC/" \
            "/media/ssd/sdc1/data/ljh/dnaHurt/rawData/2025-07-28-化疗动员-D21-ATAC_merged/" \
            > log/merge_lanes.log



#######################################
###### ATAC encode pipline 
#######################################

conda create -n caper -c conda-forge caper
conda activate caper

## 1. 下载参考基因组，因为 ENCODE ATAC pipline 需要的参考基因组文件
## 从github 网站上下载script文件夹中的sh脚本，加到服务器中
# https://github.com/ENCODE-DCC/atac-seq-pipeline/tree/master
## 也可以从已经有的参考基因组构建
#  文档：https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/build_genome_database.md

## 和平常自己用的格式有点不一样
## 直接使用官方的脚本下载比较方便
cd /media/ssd/sdb1/data/ljh/software/ref/
mkdir -p ENCODE_ATAC_pipline_ref_mm10/
cd ENCODE_ATAC_pipline_ref_mm10/

mkdir -p log
# bash scripts/download_genome_data.sh [GENOME] [DESTINATION_DIR]
conda activate caper

## 安装需要用的工具
conda install bioconda::bowtie2
conda install bioconda::samtools

nohup bash /media/ssd/sdc1/data/ljh/ENCODE_ATAC_pipline_code/atac-seq-pipeline-master/scripts/build_genome_data.sh \
        mm10 \
        ./ \
        >> log/build_genome_data.log 2>&1 &
tail -f log/build_genome_data.log


## 2. 设置输入json 文件
# 完整文档：https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/input.md
# 最少运行必须：https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/input_short.md


## 3. 使用conda 环境前的准备工作
## NOTE （失败，跳过）
# 参考 https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/install_conda.md
# 中的第五点
# 这一步会安装四个encode专用的环境
# 都是encd开头
# conda activate caper
# nohup bash /media/ssd/sdc1/data/ljh/ENCODE_ATAC_pipline_code/atac-seq-pipeline-master/scripts/uninstall_conda_env.sh \
#         >> log/uninstall_conda_env.log 2>&1 &
# tail -f log/uninstall_conda_env.log
# nohup bash /media/ssd/sdc1/data/ljh/ENCODE_ATAC_pipline_code/atac-seq-pipeline-master/scripts/install_conda_env.sh mamba \
#         >> log/install_conda_env.log 2>&1 &
# tail -10 log/install_conda_env.log
# cd /media/ssd/sdc1/data/ljh/dnaHurt/
# ## 5-FU组
# nohup caper run \
#       /media/ssd/sdc1/data/ljh/ENCODE_ATAC_pipline_code/atac-seq-pipeline-master/atac.wdl \
#       -i code/2025-10-01-Chemo_mobil_atac_A5_FU.json \
#       --conda --max-concurrent-tasks 1 \
#       > log/atac_A5_FU.log 2>&1 &
# tail -f log/atac_A5_FU.log
# less log/atac_A5_FU.log
# chmod 777 /media/ssd/sdb1/data/ljh/software/ref/ENCODE_ATAC_pipline_ref_mm10/mm10.tsv



## 4. 尝试用singularity 运行
# 下载sif镜像
# https://encode-pipeline-singularity-image.s3.us-west-2.amazonaws.com/atac-seq-pipeline_v2.2.3.sif

## 5. 设置json文件，设置不从服务器下载sif，使用本地sif文件
## 如果没有设置，则是默认代码从服务器上下载sif文件
## 但是可能受到网络因素影响，最终导致任务失败
## 所以建议还是手动下载好sif文件
## 然后在json文件中指定镜像地址

## 激活环境
conda activate singularity
## 定义caper 路径（因为caper是另一个环境中安装的）
caper=/home/ljh/miniconda3/envs/caper/bin/caper

## 定义编码方式（不定义的话，如果任务json的描述或者title中出现中文，就会导致质控结果出问题）
export SINGULARITYENV_LC_ALL=C.UTF-8


## 5-FU组 用singularity 运行
log_5_FU=log/atac_A5_FU_singularity_10_03.log
nohup $caper run \
     /media/ssd/sdc1/data/ljh/ENCODE_ATAC_pipline_code/atac-seq-pipeline-master/atac.wdl \
     -i code/2025-10-01-Chemo_mobil_atac_A5_FU.json \
     --singularity --max-concurrent-tasks 3 \
     > $log_5_FU 2>&1 &
tail -f $log_5_FU
less $log_5_FU

## PBS组 用singularity 运行
log_PBS=log/atac_PBS_singularity_10_03.log
nohup $caper run \
     /media/ssd/sdc1/data/ljh/ENCODE_ATAC_pipline_code/atac-seq-pipeline-master/atac.wdl \
     -i code/2025-10-01-Chemo_mobil_atac_PBS.json \
     --singularity --max-concurrent-tasks 3 \
     > $log_PBS 2>&1 &
tail -f $log_PBS
less $log_PBS



## 5-FU+Ana组 用singularity 运行
log_5_FU_ANS=log/atac_A5_FU_ANS_singularity_10_03.log
nohup $caper run \
     /media/ssd/sdc1/data/ljh/ENCODE_ATAC_pipline_code/atac-seq-pipeline-master/atac.wdl \
     -i code/2025-10-01-Chemo_mobil_atac_A5_FU_ANS.json \
     --singularity --max-concurrent-tasks 3 \
     > $log_5_FU_ANS 2>&1 &
tail -f $log_5_FU_ANS
less $log_5_FU_ANS



## 整理输出
conda install bioconda::croo

conda activate caper
cd /media/ssd/sdc1/data/ljh/dnaHurt/

# croo /media/ssd/sdc1/data/ljh/dnaHurt/atac/65df076a-7202-4fcf-944b-6139adcfe5ca/metadata.json \
#     --out-dir /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/65df076a-7202-4fcf-944b-6139adcfe5ca/

# croo /media/ssd/sdc1/data/ljh/dnaHurt/atac/c4d22ef5-eb0b-4af2-9bdb-370db112869a/metadata.json \
#     --out-dir /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/c4d22ef5-eb0b-4af2-9bdb-370db112869a/

# croo /media/ssd/sdc1/data/ljh/dnaHurt/atac/6d9d31c7-4030-4a30-8d92-0ae141ce89f7/metadata.json \
#     --out-dir /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/6d9d31c7-4030-4a30-8d92-0ae141ce89f7/

## D21 5FU-ANA
croo /media/ssd/sdc1/data/ljh/dnaHurt/atac/f5cc4ac5-2216-4b02-8958-7e502f226343/metadata.json \
    --out-dir /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_ANA/

## D21 5-FU (用下面的代替，重跑了)
# croo /media/ssd/sdc1/data/ljh/dnaHurt/atac/424e1a58-6eba-4ee2-9e7d-048402ad4d40/metadata.json \
#     --out-dir /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/

## D21 PBS
croo /media/ssd/sdc1/data/ljh/dnaHurt/atac/f2372544-8bc0-4fa8-85f1-157fc0d0f579/metadata.json \
    --out-dir /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_PBS/


## 整理qc文件
conda install conda-forge::qc2tsv

mkdir -p qc
qc2tsv /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/f5cc4ac5-2216-4b02-8958-7e502f226343/qc/qc.json \
    /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/424e1a58-6eba-4ee2-9e7d-048402ad4d40/qc/qc.json \
    /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/f2372544-8bc0-4fa8-85f1-157fc0d0f579/qc/qc.json \
    > qc/spreadsheet.tsv



## 查看 peak数量


#######################################
###### @date 2025-10-09
###### 因为5-FU-2 样本 之前数据不全，现在
###### 补测数据已经下来了，所以需要重新跑5-FU的ATAC上游
#######################################

## 1。 合并5-FU2 的样本
merge_lanes "/media/ssd/sdc1/data/ljh/dnaHurt/rawData/2025-07-28-化疗动员-D21-ATAC/" \
            "/media/ssd/sdc1/data/ljh/dnaHurt/rawData/2025-07-28-化疗动员-D21-ATAC_merged/" \
            > log/merge_lanes_2025_10_09.log
tail -f log/merge_lanes_2025_10_09.log

## 2. 重新跑5-FU的ATAC上游
## 激活环境
conda activate singularity
## 定义caper 路径（因为caper是另一个环境中安装的）
caper=/home/ljh/miniconda3/envs/caper/bin/caper


## 定义编码方式（不定义的话，如果任务json的描述或者title中出现中文，就会导致质控结果出问题）
export SINGULARITYENV_LC_ALL=C.UTF-8
## 5-FU组 用singularity 运行
log_5_FU=log/atac_A5_FU_singularity_10_09.log
nohup $caper run \
     /media/ssd/sdc1/data/ljh/ENCODE_ATAC_pipline_code/atac-seq-pipeline-master/atac.wdl \
     -i code/2025-10-01-Chemo_mobil_atac_A5_FU.json \
     --singularity --max-concurrent-tasks 3 \
     > $log_5_FU 2>&1 &
tail -f $log_5_FU
less $log_5_FU


## 整理输出
conda activate caper
cd /media/ssd/sdc1/data/ljh/dnaHurt/
## 5-FU
croo /media/ssd/sdc1/data/ljh/dnaHurt/atac/91f78c19-b7fa-4a66-9a84-e0fbeb82fced/metadata.json \
    --out-dir /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU/

## 整理qc报告





#######################################
###### @date 2025-10-09
###### 添加1.5M的ATAC上游分析
###### 因为之前的上游流程并不完全一致
###### 保证一致的流程
#######################################


## 激活环境
conda activate singularity
cd /media/ssd/sdc1/data/ljh/dnaHurt/

## 定义caper 路径（因为caper是另一个环境中安装的）
caper=/home/ljh/miniconda3/envs/caper/bin/caper


## 定义编码方式（不定义的话，如果任务json的描述或者title中出现中文，就会导致质控结果出问题）
export SINGULARITYENV_LC_ALL=C.UTF-8

## 5-FU组 用singularity 运行
log_5_FU=log/atac_A5_FU_1.5M_singularity_10_09.log
nohup $caper run \
     /media/ssd/sdc1/data/ljh/ENCODE_ATAC_pipline_code/atac-seq-pipeline-master/atac.wdl \
     -i code/2025-10-09-Chemo_mobil_1.5M_5_FU.json \
     --singularity --max-concurrent-tasks 3 \
     > $log_5_FU 2>&1 &
tail -f $log_5_FU
less $log_5_FU

## 整理输出
conda activate caper
cd /media/ssd/sdc1/data/ljh/dnaHurt/
## 5-FU
croo /media/ssd/sdc1/data/ljh/dnaHurt/atac/91f78c19-b7fa-4a66-9a84-e0fbeb82fced/metadata.json \
    --out-dir /media/ssd/sdc1/data/ljh/dnaHurt/atac_res/D21_5FU_1.5M/


## PBS组 用singularity 运行
log_PBS=log/atac_PBS_1.5M_singularity_10_09.log
nohup $caper run \
     /media/ssd/sdc1全（）/data/ljh/ENCODE_ATAC_pipline_code/atac-seq-pipeline-master/atac.wdl \
     -i code/2025-10-09-Chemo_mobil_1.5M_PBS.json \
     --singularity --max-concurrent-tasks 3 \
     > $log_PBS 2>&1 &
tail -f $log_PBS
less $log_PBS



# @date 2025-09-28
# @author lin
# @description 化疗动员课题 Drugseq 上游处理

conda activate RNAseq
cd /media/ssd/sdc1/data/ljh/dnaHurt

cd /media/ssd/sdc1/data/ljh/dnaHurt/rawData/2025-07-23-化疗动员-D21-Drugseq/


mkdir log
mkdir B1
mv /media/ssd/sdc1/data/ljh/dnaHurt/rawData/2025-07-23-化疗动员-D21-Drugseq/RNA-seq-1_L1_801.R1.fastq.gz B1/
mv /media/ssd/sdc1/data/ljh/dnaHurt/rawData/2025-07-23-化疗动员-D21-Drugseq/RNA-seq-1_L1_801.R2.fastq.gz B1/

## 开始上游
nohup python /media/ssd/sdc1/data/ljh/dnaHurt/code/DrugcellUMI_V2.py -f B1 -b /media/ssd/sdc1/data/ljh/dnaHurt/code/barcode24.txt -o B1_out >log/B1_log.log &
tail -f log/B1_log.log


## 开始上游（全部barcode）
nohup python /media/ssd/sdc1/data/ljh/dnaHurt/code/DrugcellUMI_V2.py -f B1 -b /media/ssd/sdc1/data/ljh/dnaHurt/code/barcode96.txt -o B1_out_barcode96 >log/B1_log_96.log &
tail -f log/B1_log.log


## test
nohup python /media/ssd/sdc1/data/ljh/dnaHurt/code/DrugcellUMI_V2.py -f test -b /media/ssd/sdc1/data/ljh/dnaHurt/code/barcode96.txt -o test_out_barcode96 >log/test_log_96.log &


tail -f log/test_log_96.log


zcat test/Iron_L4_805.R1.fastq.gz | head -1000000  | grep GTACCGAACTTA | wc -l


zcat test/Iron_L4_805.R1.fastq.gz | head -10000000  | grep GTACAGATCGCA | wc -l


## 匹配序列 



zcat B1/RNA-seq-1_L1_801.R1.fastq.gz | head -1000000  | grep GTACAAACATCG | wc -l

zcat B1/RNA-seq-1_L1_801.R1.fastq.gz | head -1000000  | grep GTACACATTGGC | wc -l

zcat B1/RNA-seq-1_L1_801.R1.fastq.gz | head -1000000  | grep GTACGCCAAGAC | wc -l


zcat B1/RNA-seq-1_L1_801.R1.fastq.gz | head -1000000  | grep GTACGACAGTGC | wc -l

zcat B1/RNA-seq-1_L1_801.R1.fastq.gz | head -1000000  | grep GTACCGCTGATC | wc -l





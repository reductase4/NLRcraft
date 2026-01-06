#!/usr/bin/env Rscript

# 读取命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <model_file> <input_file> <output_file>")
}

model_file <- args[1]
input_file <- args[2]
output_file<- args[3]

# load packages
library(tidyr)
library(dplyr)
library(stringr)
library(randomForest)

# load the random forest model
loaded_model <- readRDS(model_file)

# read data
aln_results <- read.csv(input_file,sep="\t",header=F)
df_cols <- c("query","target","fident","alnlen","qcov","tcov","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","prob","lddt","alntmscore","qtmscore","ttmscore")
colnames(aln_results) <- df_cols

# 创建NBS数据集（target以PF00931开头）
aln_results_NBS <- aln_results %>%
  filter(str_detect(target, "^PF00931"))

# 筛选NBS数据集 -----------------------------------------------------------
aln_results_filtered <- aln_results_NBS %>%
  # 按query分组处理
  group_by(query) %>%
  # 按优先级排序：alntmscore降序 -> tcov降序 -> lddt降序 
  arrange(desc(alntmscore), desc(tcov), desc(lddt)) %>%
  # 保留每个分组的第一行（即最优记录）
  slice(1) %>%
  # 解除分组
  ungroup()
# get 51815 obs (rows)

#predict data
aln_results_pred <- predict(loaded_model, newdata = aln_results_filtered)

aln_results_filtered$predicted <- aln_results_pred

#save data
write.table(aln_results_filtered, output_file,quote=F,row.names = F, sep="\t")

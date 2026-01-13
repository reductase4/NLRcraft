#!/usr/bin/env Rscript

# 读取命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript script.R <nlr_file> <NTD_rep_file> <cluster_file> <output_file>")
}

nlr_file <- args[1]
NTD_rep_file <- args[2]
cluster_file<- args[3]
output_file<- args[4]

# load packages
library(dplyr)

# read data
nlr <- read.csv(nlr_file,sep="\t",header=F)
nlr <- na.omit(nlr)
colnames(nlr) <- c("ID","domain_detection")

NTD_rep <- read.csv(NTD_rep_file, header=T)

cluster <- read.csv(cluster_file, sep="\t",header=F)
colnames(cluster) <- c("clusterID", "proteinID")

# sublass assignment
nlr_subclass <- nlr %>%
  left_join(cluster, by = c("ID" = "proteinID")) %>%
  left_join(NTD_rep %>% select(UniProtID, subclass),
            by = c("clusterID" = "UniProtID"))

# save data
write.table(nlr_subclass, output_file,quote=F,row.names = F, sep="\t")

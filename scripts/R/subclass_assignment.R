#!/usr/bin/env Rscript

# 读取命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript script.R <nlr_file> <NTD_rep_file> <cluster_file> <output_file>")
}

nlr_file <- args[1]
NTD_rep_file <- args[2]
cluster_file<- args[3]
plddt_file <- args[4]
output_file<- args[5]

# load packages
library(dplyr)

################## read data ##################
#nlr <- read.csv("step3_identification/results_all_rm_FPs.txt",sep="\t",header=F)
nlr <- read.csv(nlr_file,sep="\t",header=F)
nlr <- na.omit(nlr) # remove non-NLR
colnames(nlr) <- c("ID","domain_detection")
# ID domain_detection
# 1        fold_adnr1_model_0                N
# 2    fold_alpha_cmv_model_0                N
# 3 fold_alpha_wolf10_model_0                N

#NTD_rep <- read.csv("../NLRcraft/pre_nlr_subclass_assignment.csv", header=T)
NTD_rep <- read.csv(NTD_rep_file, header=T)
# UniProtID    all_members_community  subclass
# 1 A0A2T7F4Y2            Community1   CC-NLR
# 2 A0A2T8IHV2            Community1   CC-NLR
# 3 A0A2T7D9L2            Community1   CC-NLR

#cluster <- read.csv("step4_classification/cluster/NTD_cluster_with_max_pLDDT_ID.tsv", sep="\t",header=T)
cluster <- read.csv(cluster_file, sep="\t",header=T)
#head(cluster)
# clusterID     max_pLDDT  proteinID
# 1 A0A022PQE2 A0A022Q4E0 A0A022PQE2
# 2 A0A022PQE2 A0A022Q4E0 A0A022PSD9
# 3 A0A022PQE2 A0A022Q4E0 A0A022PU06

#plddt_scores <- read.csv("step4_classification/cluster/plddt_scores.tsv", sep="\t",header=F)
plddt_scores <- read.csv(plddt_file, sep="\t",header=F)
colnames(plddt_scores) <- c("proteinID","avg_pLDDT")
#head(plddt_scores)
# proteinID    avg_pLDDT
# 1 A0A022PN71     79.61
# 2 A0A022PPJ5     77.98
# 3 A0A022PPP8     83.00


################## sublass assignment ################## 
# merge cluster & plddt_scores
cluster_scored <- cluster %>%
  left_join(plddt_scores, by="proteinID")

# rank cluster members based on avg_pLDDT
cluster_ranked <- cluster_scored %>%
  arrange(clusterID, desc(avg_pLDDT))

# search subclass of cluster members from NTD_rep dataset, return "subclass" or "others"
assign_subclass <- function(cluster_id) {
  proteins <- cluster_ranked %>%
    filter(clusterID == cluster_id) %>%
    arrange(desc(avg_pLDDT)) %>%
    pull(proteinID)
  
  for (p in proteins) {
    hit <- NTD_rep %>% filter(UniProtID == p)
    if (nrow(hit) > 0) {
      return(hit$subclass[1])
    }
  }
  
  return("others")
}

# sublass assignment
nlr$subclass <- sapply(nlr$ID, function(x) {
  cluster_rows <- cluster %>% filter(proteinID == x)
  if (nrow(cluster_rows) == 0) return("others")
  cluster_id <- cluster_rows$clusterID[1] # get clusterID of each nlr
  assign_subclass(cluster_id) # search subclass of cluster members from NTD_rep dataset
})

nlr_subclass <- nlr

# save data
write.table(nlr_subclass, output_file,quote=F,row.names = F, sep="\t")

# save supplementary file
nlr_subclass_cluster <- left_join(nlr_subclass, cluster, by=c("ID"="proteinID"))

base <- sub("\\.[^.]+$", "", output_file) # get filename
output_file2 <- paste0(base, "_cluster.tsv")
write.table(nlr_subclass_cluster, paste0(output_file2),quote=F,row.names = F, sep="\t")

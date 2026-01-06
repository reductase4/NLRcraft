#!/bin/bash

# =====================
# check parameters
# =====================
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <PDB_DIR> <SCRIPTS_DIR>"
    exit 1
fi

PDB_DIR=$1
SCRIPTS_DIR=$2

current_path=$(pwd)
CLUSTER_DIR="$current_path/cluster"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
COMMAND_LOG="$current_path/Step03_cluster_command_log_$TIMESTAMP.txt"
DOWNLOAD_LIST="$current_path/Step03_cluster_download_paths.txt"

> $COMMAND_LOG
> $DOWNLOAD_LIST

DOMAINS=("N_terminal" "C_terminal" "NBS")

# =====================
# cluster & align parameters
# =====================
CLUSTER_CUTOFFS=(0.8 0.8 0.8) # coverage
SEARCH_CUTOFFS=(0.5 0.5 0.5) # coverage
CLUSTER_EVALUES=(0.001 0.001 0.001) # E value
SEARCH_EVALUES=(0.001 0.001 0.001) # E value

# =====================
# function: record + run
# =====================
log_and_run() {
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "[$TIMESTAMP] $1\n" | tee -a $COMMAND_LOG
    eval "$1"
}

# =====================
# Step 1: Foldseek Clustering
# =====================
echo -e "\n=== Step 1: Foldseek Clustering ===\n" | tee -a $COMMAND_LOG
for idx in "${!DOMAINS[@]}"; do
    DOMAIN=${DOMAINS[$idx]}
    DOMAIN_DIR="$CLUSTER_DIR/$DOMAIN"
    LOGFILE="$DOMAIN_DIR/log.txt"

    # 判断 nlr_cluster.tsv 是否存在且不为空
    if [ ! -f "$DOMAIN_DIR/nlr_cluster.tsv" ] || [ ! -s "$DOMAIN_DIR/nlr_cluster.tsv" ]; then
        log_and_run "mkdir -p $DOMAIN_DIR/foldseek/repDB"
        log_and_run "cd $DOMAIN_DIR"
        log_and_run "foldseek easy-cluster $PDB_DIR/$DOMAIN nlr tmp -c ${CLUSTER_CUTOFFS[$idx]} -e ${CLUSTER_EVALUES[$idx]} >> $LOGFILE 2>&1"
        log_and_run "sed -i 's/\.pdb//g' nlr_cluster.tsv"
    else
        echo "$DOMAIN: nlr_cluster.tsv already exists, skipping Step 1." | tee -a $COMMAND_LOG
    fi
done

# =====================
# Step 2: Max pLDDT Selection
# =====================
echo -e "\n=== Step 2: Max pLDDT Selection ===\n" | tee -a $COMMAND_LOG
for idx in "${!DOMAINS[@]}"; do
    DOMAIN=${DOMAINS[$idx]}
    DOMAIN_DIR="$CLUSTER_DIR/$DOMAIN"

    # 判断 nlr_cluster_with_max_pLDDT_ID.tsv 是否存在且不为空
    if [ ! -f "$DOMAIN_DIR/nlr_cluster_with_max_pLDDT_ID.tsv" ] || [ ! -s "$DOMAIN_DIR/nlr_cluster_with_max_pLDDT_ID.tsv" ]; then
        log_and_run "cd $DOMAIN_DIR"
        log_and_run "python $SCRIPTS_DIR/generate_max_pLDDT_ID.py -c nlr_cluster.tsv -s $PDB_DIR/$DOMAIN -o nlr_cluster_with_max_pLDDT_ID.tsv"
        log_and_run "echo \"$DOMAIN_DIR/nlr_cluster_with_max_pLDDT_ID.tsv\" >> $DOWNLOAD_LIST"
    else
        echo "$DOMAIN: nlr_cluster_with_max_pLDDT_ID.tsv already exists, skipping Step 2." | tee -a $COMMAND_LOG
        log_and_run "echo \"$DOMAIN_DIR/nlr_cluster_with_max_pLDDT_ID.tsv\" >> $DOWNLOAD_LIST"
    fi
done

# =====================
# Step 3: Copy Best PDBs
# =====================
echo -e "\n=== Step 3: Extract Best Representative Structures ===\n" | tee -a $COMMAND_LOG
for idx in "${!DOMAINS[@]}"; do
    DOMAIN=${DOMAINS[$idx]}
    DOMAIN_DIR="$CLUSTER_DIR/$DOMAIN/foldseek"
    LOGFILE="$DOMAIN_DIR/log.txt"

    # 判断 nlr_structs_max 文件夹是否存在且不为空
    if [ ! -d "$DOMAIN_DIR/nlr_structs_max" ] || [ ! "$(ls -A $DOMAIN_DIR/nlr_structs_max)" ]; then
        log_and_run "cd $DOMAIN_DIR"
        log_and_run "cut -f 2 $CLUSTER_DIR/$DOMAIN/nlr_cluster_with_max_pLDDT_ID.tsv | sort | uniq > best_id.txt"
        log_and_run "bash $SCRIPTS_DIR/cp_target_pdbs.sh $PDB_DIR/$DOMAIN best_id.txt nlr_structs_max >> $LOGFILE 2>&1"
        log_and_run "echo \"$DOMAIN_DIR/nlr_structs_max\" >> $DOWNLOAD_LIST"
    else
        echo "$DOMAIN: nlr_structs_max already exists, skipping Step 3." | tee -a $COMMAND_LOG
        log_and_run "echo \"$DOMAIN_DIR/nlr_structs_max\" >> $DOWNLOAD_LIST"
    fi
done

# =====================
# Step 4: Build Foldseek DB
# =====================
echo -e "\n=== Step 4: Build Foldseek Database ===\n" | tee -a $COMMAND_LOG
for idx in "${!DOMAINS[@]}"; do
    DOMAIN=${DOMAINS[$idx]}
    DOMAIN_DIR="$CLUSTER_DIR/$DOMAIN/foldseek/repDB"
    LOGFILE="$DOMAIN_DIR/log.txt"

    # 判断 rep_nlrDB 是否存在且不为空
    if [ ! -f "$DOMAIN_DIR/repDB/rep_nlrDB" ] || [ ! -s "$DOMAIN_DIR/repDB/rep_nlrDB" ]; then
        log_and_run "cd $DOMAIN_DIR"
        log_and_run "foldseek createdb ../nlr_structs_max/ rep_nlrDB >> $LOGFILE 2>&1"
        echo \"$DOMAIN_DIR/repDB/rep_nlrDB\" >> $COMMAND_LOG
    else
        echo "$DOMAIN: rep_nlrDB already exists, skipping Step 4." | tee -a $COMMAND_LOG
    fi
done

# =====================
# Step 5: Easy Search
# =====================
echo -e "\n=== Step 5: Foldseek Easy Search ===\n" | tee -a $COMMAND_LOG
for idx in "${!DOMAINS[@]}"; do
    DOMAIN=${DOMAINS[$idx]}
    SEARCH_CUTOFF=${SEARCH_CUTOFFS[$idx]}
    SEARCH_EVALUE=${SEARCH_EVALUES[$idx]}
    DOMAIN_DIR="$CLUSTER_DIR/$DOMAIN/foldseek"
    LOGFILE="$DOMAIN_DIR/log.txt"

    # 判断 aln_all.txt 是否存在且不为空
    if [ ! -f "$DOMAIN_DIR/aln_all.txt" ] || [ ! -s "$DOMAIN_DIR/aln_all.txt" ]; then
        log_and_run "cd $DOMAIN_DIR"
        log_and_run "foldseek easy-search nlr_structs_max/ repDB/rep_nlrDB aln_all.txt tmp -c $SEARCH_CUTOFF -e $SEARCH_EVALUE --format-output query,target,fident,alnlen,qcov,tcov,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,prob,lddt,alntmscore,qtmscore,ttmscore >> $LOGFILE 2>&1"
        log_and_run "echo \"$DOMAIN_DIR/aln_all.txt\" >> $DOWNLOAD_LIST"
    else
        echo "$DOMAIN: aln_all.txt already exists, skipping Step 5." | tee -a $COMMAND_LOG
        log_and_run "echo \"$DOMAIN_DIR/aln_all.txt\" >> $DOWNLOAD_LIST"
    fi
done

echo "=== Pipeline finished ==="
echo "All executed commands logged at: $COMMAND_LOG"
echo "Downloadable files listed at: $DOWNLOAD_LIST"

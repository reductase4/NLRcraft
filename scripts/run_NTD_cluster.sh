#!/bin/bash

set -euo pipefail

PDB_TARGET=$1
PDB_REF=$2

CLUSTER_CUTOFFS=0.8
CLUSTER_EVALUES=0.001
PDB_DIR="NTD_structs"
LOG="cluster.log"

mkdir -p "$PDB_DIR"

# copy pdbs
cp "$PDB_TARGET"/*.pdb "$PDB_DIR"/
cp "$PDB_REF"/*.pdb    "$PDB_DIR"/

echo "[INFO] Clustering NTD structures" >> "$LOG"
foldseek easy-cluster "$PDB_DIR" NTD tmp -c ${CLUSTER_CUTOFFS} -e ${CLUSTER_EVALUES} >> "$LOG" 2>&1

# remove .pdb suffix in cluster output
sed 's/\.pdb//g' NTD_cluster.tsv > NTD_cluster_new.tsv



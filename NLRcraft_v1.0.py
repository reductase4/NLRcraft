#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NLRcraft: A structural-homology-based pipeline for plant NLR detection and classification.

This script automates:
1. Domain database building by Foldseek
2. pLDDT filtering
3. Identification base on structural alignments
4. Classification base on structural clustering of N-terminal domains

This script serves as the main controller of the NLRcraft pipeline.

Usage example:

python NLRcraft.py \
  --structs protein_structures \
  --ids proteins_ids.txt \
  --plddt pLDDT_cutoff_value
"""

import os
from datetime import datetime
import subprocess
import argparse
import sys

def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def run(cmd, label):
    print(f"\n[{timestamp()}] {label}")
    print(f"[{timestamp()}]   $ {' '.join(cmd)}")
    
    res = subprocess.run(cmd)
    if res.returncode != 0:
        print(f"[{timestamp()}] ERROR: {label} failed.")
        sys.exit(1)
    print(f"[{timestamp()}] Completed: {label}")


def main():

    parser = argparse.ArgumentParser(description="NLRcraft: Identify and classify plant NLRs based on structures.")
    
    parser.add_argument("-i", "--structs", required=True)
    parser.add_argument("-d", "--ids", required=True)
    parser.add_argument("-p", "--plddt", default=60, type=int)

    args = parser.parse_args()

    structs = args.structs
    plddt = str(args.plddt)
    id_file = args.ids


    base_dir = os.path.dirname(os.path.abspath(__file__))
    script_dir = f"{base_dir}/scripts"
    work_dir = os.getcwd()

    domain_folds = f"{base_dir}/domains_pdb"

    # output files
    filtered_structs = "structs_filter"
    aln_results = "aln_all.txt"
    aln_results_filter = "results_all.txt"
    rf_predict_results = "rf_prediction.txt" # prediction by a random forest model
    aln_results_rm_FPs = "results_all_rm_FPs.txt" # identifiaction results after removing non-NLR homologs (false positives)
    nbs_position = "NBS_pos.tsv"
    

    # STEP 1: Foldseek DB building
    run(["mkdir", "step1_domain_db"], "mkdir step1_domain_db")

    run(["cd", "step1_domain_db"], "cd step1_domain_db")

    run(
        ["foldseek", "createdb", domain_folds, "nlrDB"],
        "STEP 1: Build Foldseek database"
    )

    # STEP 2 Filter predicted structures based on pLDDT values
    run(["cd", ".."], "")

    run(["mkdir", "step2_plddt_filter"], "mkdir step2_plddt_filter")

    run(["cd", "step2_plddt_filter"], "cd step2_plddt_filter")

    run(
        ["python", f"{script_dir}/run_plddt_filter.py", plddt, structs, filtered_structs],
        "STEP 2: Filter structures based on pLDDT values"
    )

    # STEP 3 Identification
    run(["cd", ".."], "")

    run(["mkdir", "step3_identification"], "mkdir step3_identification")

    run(["cd", "step3_identification"], "cd step3_identification")

    # STEP 3.1 Foldseek alignment
    run(
        [
            "foldseek", "easy-search", f"{work_dir}/step2_plddt_filter/{filtered_structs}", f"{work_dir}/step1_domain_db/nlrDB", aln_results, "tmp",
            "--format-output",
            "query,target,fident,alnlen,qcov,tcov,mismatch,gapopen,qstart,qend,"
            "tstart,tend,evalue,bits,prob,lddt,alntmscore,qtmscore,ttmscore"
        ],
        "STEP 3.1: Run Foldseek alignment"
    )

    # STEP 3.2 Alignment results extraction
    # python extract_align_results.py aln.txt ids.txt results_all.txt -p 1 -tc 0.8 -tm 0.5 -e 0.001
    # default prob_cutoff=1.0, tcov_cutoff=0.8, tm_cutoff=.5, evalue_cutoff=0.001
    run(
        ["python", f"{script_dir}/extract_align_results.py", aln_results, {id_file}, aln_results_filter, "-p", "0.99"],
        "STEP 3.2: Extract alignment results (filter by parameteres)"
    )

    # STEP 3.3 Remove non-NLR homologs
    run(
        ["Rscript", f"{script_dir}/R/rf_predict.R", f"{script_dir}/R/final_rf_model_undersampling.rds", aln_results_filter, rf_predict_results],
        "STEP 3.3.1: Random Forest prediction"
    )

    run(
        ["python", f"{script_dir}/rm_FPs.py", "-i", aln_results_filter, "-p", rf_predict_results, "-o", aln_results_rm_FPs],
        "STEP 3.3.2: Remove false positives"
    )

    # STEP 4 Classification
    run(["cd", ".."], "")

    run(["mkdir", "step4_classification"], "mkdir step4_classification")

    run(["cd", "step4_classification"], "cd step4_classification")

    # STEP 4.1 Get NBS position
    run(
        ["python", f"{script_dir}/extract_NBS_pos.py", "-a", f"{work_dir}/step3_identification/aln_filtered.txt", "-i", id_file, "-o", nbs_position],
        "STEP 4.1: Extract NBS location"
    )

    # STEP 4.2 Split domains based on NBS location
    os.makedirs("split_domains", exist_ok=True)

    run(
        ["python", f"{script_dir}/split_pdb_by_NBS.py", "-n", nbs_position, "-p", f"{work_dir}/step2_plddt_filter/{filtered_structs}", "-o", "split_domains"],
        "STEP 4.2: Split domains based on NBS location"
    )

'''
    # STEP 4.3 Structural clustering
    run(
        ["bash", f"{script_dir}/run_NLR_domains_cluster.sh", "split_domains", args.script_dir],
        "STEP 4.3: Structural clustering"
    )

    # STEP 4.4 Community detection & final classification
    run(
        ["Rscript", f"{args.script_dir}/R/network_community.R"],
        "STEP 4.4: Community detection & final classification"
    )
'''

    print("\n NLRcraft Pipeline Finished Successfully!")
    print("Output directory includes:")
    print(f" - NLR identification: {work_dir}/step3_identification/{aln_results_rm_FPs}")
    #print(f" - NLR classification: {work_dir}/step4_classification/{classification}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NLRcraft: A structure-homology-based pipeline for plant NLR identification
and classification.

Major steps:
1. Domain database building by Foldseek
2. pLDDT-based structure filtering
3. NLR identification based on structural alignments
4. NLR classification based on structural clustering of N-terminal domains

This script serves as the main controller of the NLRcraft pipeline.
"""

import os
import sys
import argparse
import subprocess
from datetime import datetime


# ================================================================
# Utility functions
# ================================================================

def timestamp():
    """Return current timestamp string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def run(cmd, label, cwd=None):
    """
    Execute a command-line program with logging and error handling.

    Parameters
    ----------
    cmd : list
        Command and arguments.
    label : str
        Description of the pipeline step.
    cwd : str, optional
        Working directory for the command.
    """
    print(f"\n[{timestamp()}] {label}")
    print(f"[{timestamp()}]   $ {' '.join(cmd)}")

    res = subprocess.run(cmd, cwd=cwd)
    if res.returncode != 0:
        print(f"[{timestamp()}] ERROR: {label} failed.")
        sys.exit(1)

    print(f"[{timestamp()}] Completed: {label}")


# ================================================================
# Step control utilities (resume / skip mechanism)
# ================================================================

def step_done(step_name):
    """
    Check whether a pipeline step has been successfully completed.

    A step is considered completed if a corresponding '.done' file exists.
    """
    return os.path.exists(f"{step_name}.done")


def mark_step_done(step_name):
    """
    Mark a pipeline step as completed by creating a '.done' file.
    """
    with open(f"{step_name}.done", "w") as f:
        f.write(f"Completed at {timestamp()}\n")


def should_skip(step_name, args):
    """
    Determine whether a pipeline step should be skipped.

    A step is skipped if:
    - It is explicitly listed in --skip, or
    - --resume is enabled and the step has already been completed.
    """
    if step_name in args.skip:
        print(f"[SKIP] {step_name} (explicitly skipped by user)")
        return True

    if args.resume and step_done(step_name):
        print(f"[SKIP] {step_name} (already completed)")
        return True

    return False


# ================================================================
# Main pipeline
# ================================================================

def main():

    parser = argparse.ArgumentParser(
        description="NLRcraft: Identify and classify plant NLRs based on protein structures."
    )

    parser.add_argument(
        "-i", "--structs", required=True,
        help="Directory containing predicted protein structures (PDB/mmCIF)."
    )
    parser.add_argument(
        "-d", "--ids", required=True,
        help="File containing target protein IDs."
    )
    parser.add_argument(
        "-p", "--plddt", type=int, default=60,
        help="pLDDT cutoff for structure filtering (default: 60)."
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume pipeline and skip completed steps."
    )
    parser.add_argument(
        "--skip", nargs="+", default=[],
        help="Skip specific steps (e.g. step1 step3.2 step4)."
    )

    args = parser.parse_args()

    structs = os.path.abspath(args.structs)
    id_file = os.path.abspath(args.ids)
    plddt = str(args.plddt)

    # directories
    base_dir = os.path.dirname(os.path.abspath(__file__))
    script_dir = os.path.join(base_dir, "scripts")
    work_dir = os.getcwd()

    domain_folds = os.path.join(base_dir, "domains_pdb")

    # --------------------------------------------------------------
    # Output files (explicit domain-level vs protein-level semantics)
    # --------------------------------------------------------------
    filtered_structs = "structs_filter"

    aln_domain = "aln_all.txt"               # raw domain-level alignments
    aln_filtered = "aln_filtered.txt"         # filtered domain-level alignments
    results_all = "results_all.txt"           # protein-level inference (N / T / TN / Na)
    rf_predict = "rf_prediction.txt"          # prediction by a randomForest model
    results_rm_FPs = "results_all_rm_FPs.txt" # inference results after FP removal

    nbs_position = "NBS_pos_alignment.txt"


    # ================================================================
    # STEP 1: Domain database building
    # ================================================================
    if not should_skip("step1", args):

        os.makedirs("step1_domain_db", exist_ok=True)

        run(
            ["foldseek", "createdb", domain_folds, "nlrDB"],
            "STEP 1: Build Foldseek domain database",
            cwd="step1_domain_db"
        )

        mark_step_done("step1")


    # ================================================================
    # STEP 2: pLDDT-based structure filtering
    # ================================================================
    if not should_skip("step2", args):

        os.makedirs("step2_plddt_filter", exist_ok=True)

        run(
            [
                "python", f"{script_dir}/run_plddt_filter.py",
                plddt, structs, filtered_structs
            ],
            "STEP 2: Filter structures based on pLDDT",
            cwd="step2_plddt_filter"
        )

        mark_step_done("step2")


    # ================================================================
    # STEP 3: NLR identification by structural alignment
    # ================================================================
    os.makedirs("step3_identification", exist_ok=True)

    # STEP 3.1 Domain-level structural alignment
    if not should_skip("step3.1", args):

        run(
            [
                "foldseek", "easy-search",
                f"{work_dir}/step2_plddt_filter/{filtered_structs}",
                f"{work_dir}/step1_domain_db/nlrDB",
                aln_domain,
                "tmp",
                "--format-output",
                "query,target,fident,alnlen,qcov,tcov,mismatch,gapopen,"
                "qstart,qend,tstart,tend,evalue,bits,prob,lddt,"
                "alntmscore,qtmscore,ttmscore"
            ],
            "STEP 3.1: Domain-level structural alignment (Foldseek)",
            cwd="step3_identification"
        )

        mark_step_done("step3.1")


    # STEP 3.2 Filter domain-level alignments and Protein-level inference
    if not should_skip("step3.2", args):

        # python extract_align_results.py aln.txt ids.txt results_all.txt -p 1 -tc 0.8 -tm 0.5 -e 0.001
        # default parameter: prob_cutoff=1.0, tcov_cutoff=0.8, tm_cutoff=.5, evalue_cutoff=0.001
        # Automatically generate "aln_filtered.txt"
        run(
            [
                "python", f"{script_dir}/extract_align_results.py",
                aln_domain, id_file, results_all,
                "-p", "0.99"
            ],
            "STEP 3.2: Filter domain-level alignments",
            cwd="step3_identification"
        )

        mark_step_done("step3.2")


    # STEP 3.3 NLR/non-NLR classification prediction and false positive removal
    if not should_skip("step3.3", args):

        run(
            [
                "Rscript", f"{script_dir}/R/rf_predict.R",
                f"{script_dir}/R/final_rf_model_undersampling.rds",
                aln_filtered, rf_predict
            ],
            "STEP 3.3.1: Random Forest prediction (NLR/non-NLR)",
            cwd="step3_identification"
        )
        # Automatically generate "NLR_ids.txt"
        run(
            [
                "python", f"{script_dir}/rm_FPs.py",
                "-i", results_all,
                "-p", rf_predict,
                "-o", results_rm_FPs
            ],
            "STEP 3.3.2: Remove false positive NLR candidates",
            cwd="step3_identification"
        )

        mark_step_done("step3.3")


    # ================================================================
    # STEP 4: NLR classification
    # ================================================================
    os.makedirs("step4_classification", exist_ok=True)

    # STEP 4.1 Extract NBS (NB-ARC) positions from alignment file
    if not should_skip("step4.1", args):

        run(
            [
                "python", f"{script_dir}/extract_NBS_pos.py",
                "-a", f"{work_dir}/step3_identification/{aln_filtered}",
                "-i", f"{work_dir}/step3_identification/NLR_ids.txt",
                "-o", nbs_position
            ],
            "STEP 4.1: Extract NBS (NB-ARC) positions from alignment file",
            cwd="step4_classification"
        )

        mark_step_done("step4.1")


    # STEP 4.2 Split domains based on NBS position
    # Automatically generate "NBS_pos_original_pdb.txt"
    if not should_skip("step4.2", args):

        os.makedirs("step4_classification/split_domains", exist_ok=True)

        run(
            [
                "python", f"{script_dir}/split_pdb_by_NBS.py",
                "-n", nbs_position,
                "-p", f"{work_dir}/step2_plddt_filter/{filtered_structs}",
                "-o", f"{work_dir}/step4_classification/split_domains",
            ],
            "STEP 4.2: Split N-terminal and NB-ARC domains",
            cwd="step4_classification"
        )

        mark_step_done("step4.2")


    """
    STEP 4.3: Structural clustering of N-terminal domains
    STEP 4.4: Community detection and final classification
    """

    print("\n[NLRcraft] Pipeline finished successfully.")
    print("Key outputs:")
    print(f" - NLR identification: step3_identification/{results_rm_FPs}")
    print(f" - NBS positions: step4_classification/{nbs_position}")


if __name__ == "__main__":
    main()

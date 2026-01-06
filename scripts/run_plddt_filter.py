#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: run_plddt_filter.py
# Description: Batch filter PDB or CIF files based on pLDDT cutoff

import os
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

cutoff = sys.argv[1]            # pLDDT cutoff
input_folder = sys.argv[2]
output_folder = sys.argv[3]

current_script_dir = os.path.dirname(os.path.abspath(__file__))
filter_script_path = os.path.join(current_script_dir, 'filter_low_plddt.py')

os.makedirs(output_folder, exist_ok=True)

def process_structure_file(filename):
    if not filename.endswith(('.pdb', '.cif')):
        print(f"Skipping non-structure file: {filename}", flush=True)
        return

    input_path = os.path.join(input_folder, filename)
    output_path = os.path.join(output_folder, filename)

    print(f"Processing {filename} ...", flush=True)
    command = ['python', filter_script_path, str(cutoff), input_path, output_path]

    try:
        subprocess.run(command, check=True)
        print(f"Finished processing: {filename}", flush=True)
    except subprocess.CalledProcessError as e:
        print(f"Error processing {filename}: {e}", flush=True)

def main():
    structure_files = [f for f in os.listdir(input_folder) if f.endswith(('.pdb', '.cif'))]

    for structure in structure_files:
        try:
            process_structure_file(structure)
        except Exception as exc:
            print(f"{structure} error: {exc}")

    print(f"\nAll files processed. Results saved in: {output_folder}")

if __name__ == "__main__":
    main()

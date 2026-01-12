#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: filter_low_plddt.py
# Description: Filter residues with low pLDDT and return a PDB file with only high pLDDT residues
# Usage: python filter_low_plddt.py <cutoff> <input_pdb> <output_pdb>
# Created: 2025-11-04 20:42:32
# Author: Adapted from Hyunbin Kim (khb7840@gmail.com) by Qian Jiang

from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
import sys
import os
import argparse

# Custom residue selector based on B-factor
class BFactorFilter(Select):
    def __init__(self, cutoff):
        self.cutoff = cutoff
    def accept_residue(self, residue):
        if "CA" in residue:
            return residue["CA"].get_bfactor() >= self.cutoff
        return False

def main():
    parser = argparse.ArgumentParser(description="Filter low-pLDDT residues from PDB or mmCIF files")
    parser.add_argument("cutoff", type=float, help="pLDDT cutoff value")
    parser.add_argument("input_file", help="Input structure file (.pdb or .cif)")
    parser.add_argument("output_file", help="Output filtered structure file (always .pdb)")
    args = parser.parse_args()

    # Auto-detect format based on file extension
    if args.input_file.endswith(".pdb"):
        structure_parser = PDBParser(QUIET=True)
    elif args.input_file.endswith(".cif"):
        structure_parser = MMCIFParser(QUIET=True)
    else:
        sys.exit("Unsupported file format. Only .pdb and .cif are accepted.")

    # Parse the structure
    structure = structure_parser.get_structure("input", args.input_file)

    # Write filtered structure to .pdb file
    #output_pdb = args.output_file if args.output_file.endswith(".pdb") else args.output_file + ".pdb" ==> error: ".cif.pdb"
    base = os.path.splitext(args.output_file)[0]
    output_pdb = base + ".pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=BFactorFilter(args.cutoff))

    print(f"Saved filtered structure to {output_pdb}")

if __name__ == "__main__":
    main()

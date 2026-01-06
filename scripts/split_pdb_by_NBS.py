#!/usr/bin/env python3
#v3: + --split-pdb
# v4: + -c/--cutdir for optional cutting of original PDBs

import os
import argparse

def parse_nbs_file(nbs_file):
    nbs_dict = {}
    with open(nbs_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            uniprot_id, qstart, qend = parts
            try:
                nbs_dict[uniprot_id] = (int(qstart), int(qend))
            except ValueError:
                continue
    return nbs_dict

def split_pdb_real_resi(pdb_path, qstart, qend):
    """
    1) 先按filtered文件中的真实residue编号生成虚拟编号，
    2) 将Foldseek的qstart/qend映射到真实pdb_start/pdb_end，
    3) 返回切割结果与真实编号。
    """
    n_term, nbs, c_term = [], [], []
    res_list = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    resnum = int(line[22:26].strip())
                    res_list.append((resnum, line))
                except ValueError:
                    continue

    # 虚拟编号 -> 真实编号映射
    virtual_to_pdb = []
    seen = set()
    for resnum, _ in res_list:
        if resnum not in seen:
            virtual_to_pdb.append(resnum)
            seen.add(resnum)

    if qstart <= len(virtual_to_pdb) and qend <= len(virtual_to_pdb):
        pdb_start = virtual_to_pdb[qstart - 1]
        pdb_end   = virtual_to_pdb[qend - 1]
    else:
        pdb_start, pdb_end = None, None

    # 仅在需要时才使用（真实分割）
    for resnum, line in res_list:
        if pdb_start is None or pdb_end is None:
            continue
        if resnum < pdb_start:
            n_term.append(line)
        elif pdb_start <= resnum <= pdb_end:
            nbs.append(line)
        elif resnum > pdb_end:
            c_term.append(line)

    return n_term, nbs, c_term, pdb_start, pdb_end

def save_to_file(lines, outpath):
    if lines:
        with open(outpath, 'w') as f:
            f.writelines(lines)

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def main(nbs_file, pdb_dir_filtered, outdir, do_split, orig_dir=None):
    """
    pdb_dir_filtered: 经过 pLDDT 过滤的 PDB 文件目录（用于坐标映射）
    orig_dir: 若提供，则切割原始 PDB 文件；否则切割 filtered 文件
    """
    if do_split:
        nbs_out   = os.path.join(outdir, "NBS")
        nterm_out = os.path.join(outdir, "N_terminal")
        cterm_out = os.path.join(outdir, "C_terminal")
        for p in [nbs_out, nterm_out, cterm_out]:
            ensure_dir(p)

    nbs_dict = parse_nbs_file(nbs_file)

    # 写新版NBS坐标表
    if do_split:
        v2_outfile = os.path.join(outdir, "NBS_pos_v2.txt")
    else:
        v2_outfile = os.path.join("NBS_pos_new.txt")

    with open(v2_outfile, "w") as out_f:
        out_f.write("UniProtID\tqstart\tqend\tpdb_start\tpdb_end\n")

        for uniprot_id, (qstart, qend) in nbs_dict.items():
            filtered_pdb = os.path.join(pdb_dir_filtered, f"{uniprot_id}.pdb")
            if not os.path.isfile(filtered_pdb):
                print(f"Warning: Missing filtered PDB for {uniprot_id}")
                continue

            # 先用 filtered 文件确定真实坐标
            _, _, _, pdb_start, pdb_end = split_pdb_real_resi(filtered_pdb, qstart, qend)

            if pdb_start is None or pdb_end is None:
                print(f"Warning: Invalid mapping for {uniprot_id}, skipped.")
                continue

            # 是否需要真正分割 PDB 文件
            if do_split:
                # 如果提供原始目录则从原始文件切割，否则从 filtered 文件切割
                pdb_source = os.path.join(orig_dir if orig_dir else pdb_dir_filtered,
                                          f"{uniprot_id}.pdb")
                if not os.path.isfile(pdb_source):
                    print(f"Warning: Missing source PDB for {uniprot_id}")
                else:
                    if orig_dir:
                        # 使用映射到真实残基编号的区间切割原始结构
                        n_term, nbs, c_term, _, _ = split_pdb_real_resi(pdb_source, pdb_start, pdb_end)
                    else:
                        # 直接用 Foldseek 的 qstart/qend 切割 filtered 结构
                        n_term, nbs, c_term, _, _ = split_pdb_real_resi(pdb_source, qstart, qend)

                    save_to_file(n_term, os.path.join(nterm_out, f"{uniprot_id}.pdb"))
                    save_to_file(nbs,   os.path.join(nbs_out,   f"{uniprot_id}.pdb"))
                    save_to_file(c_term,os.path.join(cterm_out, f"{uniprot_id}.pdb"))

            # 输出新的区间信息
            out_f.write(f"{uniprot_id}\t{qstart}\t{qend}\t{pdb_start}\t{pdb_end}\n")

    print(f"Finished. Output saved to {v2_outfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split PDB files using Foldseek-based NBS positions, "
                    "mapping virtual to real residue numbers."
    )
    parser.add_argument("-n", "--nbs", required=True,
                        help="Path to NBS_pos.txt (no header)")
    parser.add_argument("-p", "--pdbdir", required=True,
                        help="Directory with filtered UniProtID.pdb (for coordinate mapping)")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Output directory (will contain NBS_pos_v2.txt and optionally split PDBs)")
    parser.add_argument("--split-pdb", action="store_true",
                        help="If set, write split PDBs to subfolders")
    parser.add_argument("-c", "--cutdir",
                        help="Optional directory of original (unfiltered) PDB files to cut; "
                             "if omitted, cut filtered PDBs")

    args = parser.parse_args()
    main(args.nbs, args.pdbdir, args.outdir, args.split_pdb, args.cutdir)


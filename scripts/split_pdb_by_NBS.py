#!/usr/bin/env python3

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

    # 生成虚拟编号到真实residue编号映射（去重、保序）
    virtual_to_pdb = []
    seen_resnums = set()
    for resnum, _ in res_list:
        if resnum not in seen_resnums:
            virtual_to_pdb.append(resnum)
            seen_resnums.add(resnum)

    # 判断qstart和qend是否在有效范围内
    if qstart <= len(virtual_to_pdb) and qend <= len(virtual_to_pdb):
        pdb_start = virtual_to_pdb[qstart - 1]
        pdb_end = virtual_to_pdb[qend - 1]
    else:
        pdb_start, pdb_end = None, None

    # 按真实编号区间进行切割
    for resnum, line in res_list:
        if pdb_start is None or pdb_end is None:
            continue  # 跳过非法
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

def main(nbs_file, pdb_dir, outdir):
    # 输出子文件夹
    nbs_out = os.path.join(outdir, "NBS")
    nterm_out = os.path.join(outdir, "N_terminal")
    cterm_out = os.path.join(outdir, "C_terminal")
    for p in [nbs_out, nterm_out, cterm_out]:
            ensure_dir(p)

    # 读取 NBS_pos
    nbs_dict = parse_nbs_file(nbs_file)

    # 新版坐标输出
    nbs_pos_v2_lines = ["UniProtID\tqstart\tqend\tpdb_start\tpdb_end"]

    for uniprot_id, (qstart, qend) in nbs_dict.items():
        pdb_file = os.path.join(pdb_dir, f"{uniprot_id}.pdb")
        if not os.path.isfile(pdb_file):
            print(f"Warning: Missing PDB file for {uniprot_id}")
            continue

        n_term, nbs, c_term, pdb_start, pdb_end = split_pdb_real_resi(pdb_file, qstart, qend)

        # 保存切割pdb
        save_to_file(n_term, os.path.join(nterm_out, f"{uniprot_id}.pdb"))
        save_to_file(nbs, os.path.join(nbs_out, f"{uniprot_id}.pdb"))
        save_to_file(c_term, os.path.join(cterm_out, f"{uniprot_id}.pdb"))

        # 保存新版区间信息
        if pdb_start is not None and pdb_end is not None:
            nbs_pos_v2_lines.append(f"{uniprot_id}\t{qstart}\t{qend}\t{pdb_start}\t{pdb_end}")
        else:
            print(f"Warning: Invalid mapping for {uniprot_id}, skipped.")

    # 写新版NBS坐标表
    v2_outfile = os.path.join("NBS_pos_original_pdb.txt")
    with open(v2_outfile, "w") as out_f:
        out_f.write("\n".join(nbs_pos_v2_lines))

    print(f"Finished. Output saved to {v2_outfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split PDB files into N-terminal, NBS, and C-terminal regions based on real residue numbers (considering pLDDT filtering gaps).")
    parser.add_argument("-n", "--nbs", required=True, help="Path to original NBS_pos.txt (no header)")
    parser.add_argument("-p", "--pdbdir", required=True, help="Directory containing UniProtID.pdb files")
    parser.add_argument("-o", "--outdir", required=True, help="Base output directory (will create NBS/, N_terminal/, C_terminal/ and NBS_pos_v2.txt)")

    args = parser.parse_args()
    main(args.nbs, args.pdbdir, args.outdir)

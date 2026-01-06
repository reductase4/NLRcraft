#!/usr/bin/env python3

import argparse

# aln_file: query, target, fident, alnlen, qcov, tcov, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits, prob, lddt, alntmscore, qtmscore, ttmscore
def extract_best_hits(aln_file, id_file, output_file):
    best_hit_dict = {}

    # 读取目标ID（无.pdb后缀）
    with open(id_file, 'r') as f:
        target_ids = set(line.strip() for line in f if line.strip())

    with open(aln_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 19:
                continue  # 忽略非标准行

            query_raw = parts[0]
            query = query_raw.replace('.pdb', '')

            if query not in target_ids:
                continue

            target = parts[1]
            if not target.startswith("PF00931_"):
                continue

            try:
                bits = float(parts[13])
                qstart = parts[8]
                qend = parts[9]
            except ValueError:
                continue  # 非数字等异常情况

            if query not in best_hit_dict or bits > best_hit_dict[query]['bits']:
                best_hit_dict[query] = {'bits': bits, 'qstart': qstart, 'qend': qend}

    # 写结果
    with open(output_file, 'w') as out:
        out.write("query\tqstart\tqend\n")
        for query, hit in best_hit_dict.items():
            out.write(f"{query}\t{hit['qstart']}\t{hit['qend']}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract qstart/qend for top PF00931 hits from large aln file without headers.")
    parser.add_argument('-a', '--aln', required=True, help='Path to aln_filtered.txt')
    parser.add_argument('-i', '--ids', required=True, help='Path to nlr_id.txt')
    parser.add_argument('-o', '--output', required=True, help='Output file name')

    args = parser.parse_args()
    extract_best_hits(args.aln, args.ids, args.output)

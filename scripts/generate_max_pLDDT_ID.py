import os
import argparse
from Bio.PDB import PDBParser

def parse_pdb_plddt(pdb_file):
    """
    计算 PDB 文件中所有原子的 pLDDT 平均值。
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", pdb_file)
        plddt_values = [
            atom.bfactor for atom in structure.get_atoms() if atom.bfactor is not None
        ]
        if plddt_values:
            return sum(plddt_values) / len(plddt_values)
    except Exception as e:
        print(f"Error reading {pdb_file}: {e}")
    return None

def calculate_plddt(structs_dir, uniprot_ids):
    """
    根据 Uniprot IDs 计算对应 PDB 文件的平均 pLDDT 值。
    """
    plddt_scores = {}
    for uniprot_id in uniprot_ids:
        pdb_file = os.path.join(structs_dir, f"{uniprot_id}.pdb")
        if os.path.exists(pdb_file):
            avg_plddt = parse_pdb_plddt(pdb_file)
            #print(avg_plddt)
            if avg_plddt is not None:
                plddt_scores[uniprot_id] = avg_plddt
    return plddt_scores

def load_plddt_scores(plddt_file):
    """
    从文件中加载 pLDDT 值。
    """
    plddt_scores = {}
    with open(plddt_file, "r") as f:
        for line in f:
            uniprot_id, plddt = line.strip().split("\t")
            plddt_scores[uniprot_id] = float(plddt)
    return plddt_scores

def process_clusters(cluster_file, structs_dir, output_file, plddt_file=None):
    """
    处理聚类文件，替换 repID 为平均 pLDDT 值最大的蛋白。
    """
    # 读取聚类文件
    cluster_map = {}
    with open(cluster_file, "r") as f:
        for line in f:
            rep_id, uniprot_id = line.strip().split("\t")
            cluster_map.setdefault(rep_id, []).append(uniprot_id)

    # 收集所有 uniprotID
    all_uniprot_ids = {uniprot_id for ids in cluster_map.values() for uniprot_id in ids}

    # 计算或加载 pLDDT 值
    if plddt_file:
        print(f"Loading pLDDT values from {plddt_file}...")
        plddt_scores = load_plddt_scores(plddt_file)
    else:
        print("Calculating pLDDT values...")
        plddt_scores = calculate_plddt(structs_dir, all_uniprot_ids)

        # 输出所有 uniprotID 对应的平均 pLDDT 值
        with open("plddt_scores.tsv", "w") as f:
            for uniprot_id, plddt in sorted(plddt_scores.items()):
                f.write(f"{uniprot_id}\t{plddt:.2f}\n")

    # 替换 repID 为平均 pLDDT 值最大的蛋白
    print("Replacing repID with highest pLDDT...")
    with open(output_file, "w") as f:
        f.write(f"clusterID\tmax_pLDDT\tproteinID\n")
        for rep_id, uniprot_ids in cluster_map.items():
            # 计算该 cluster 中的最高 pLDDT 值
            best_id = max(uniprot_ids, key=lambda x: plddt_scores.get(x, -1))
            for uniprot_id in uniprot_ids:
                f.write(f"{rep_id}\t{best_id}\t{uniprot_id}\n")  # 写入文件

    print(f"Finished! Results saved to {output_file}.")

def main():
    parser = argparse.ArgumentParser(description="Process NLR clusters based on pLDDT values.")
    parser.add_argument("-c", "--cluster_file", required=True, help="Input cluster file (e.g., nlr_cluster.tsv).")
    parser.add_argument("-s", "--structs_dir", help="Directory containing PDB files. Not required if --plddt_scores is provided.")
    parser.add_argument("-p", "--plddt_scores", help="Precomputed pLDDT scores file. Skips pLDDT calculation if provided.")
    parser.add_argument("-o", "--output_file", required=True, help="Output file for updated clusters.")
    args = parser.parse_args()

    if not args.structs_dir and not args.plddt_scores:
        parser.error("Either --structs_dir or --plddt_scores must be provided.")

    process_clusters(args.cluster_file, args.structs_dir, args.output_file, args.plddt_scores)

if __name__ == "__main__":
    main()


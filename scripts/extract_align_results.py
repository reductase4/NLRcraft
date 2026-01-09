# v2: 1.'PF01582'='TIR', 'PF13676'='TIR2', 'PF00931'='NBD'; 2."prob = float(data[12])" changed to "prob = float(data[14])", for 'qcov' & 'tcov' added in align cmd; 3.results='' > results='NA'
# v3: + tcov_cutoff, tm_cutoff, evalue_cutoff
# v4: + prob_cutoff, use argparse & default values
# v5: nlr_cluster_file => id_file
import pandas as pd
import sys
import argparse

# Read the aln_ae file and extract the results
def extract_results_from_aln(file_path, prob_cutoff,tcov_cutoff, tm_cutoff, evalue_cutoff):
    results = {}
    with open('aln_filtered.txt', 'w') as filtered_file, open(file_path, 'r') as file:
        for line in file:
            data = line.strip().split()
            query_id = data[0].split('.')[0]
            target = data[1]
            tcov = float(data[5])
            evalue = float(data[12])
            prob = float(data[14])
            alnTM = float(data[16])
            
            # results[query_id] = '', 'T', 'N', 'TN'
            if query_id not in results.keys():
                results[query_id] = 'NA'
            # 'PF01582'='TIR', 'PF13676'='TIR2', 'PF00931'='NBD'
            if prob >= prob_cutoff and tcov >= tcov_cutoff and alnTM >= tm_cutoff and evalue < evalue_cutoff:
                filtered_file.write(line)
                if 'PF01582' in target or 'PF13676' in target:
                    if results[query_id] == 'NA':
                        results[query_id] = 'T'
                    elif results[query_id] == 'T':
                        results[query_id] = 'T'
                    else:
                        results[query_id] = 'TN'
                elif 'PF00931' in target:
                    if results[query_id] == 'NA':
                        results[query_id] = 'N'
                    elif results[query_id] == 'N':
                        results[query_id] = 'N'
                    else:
                        results[query_id] = 'TN'
    
    return results



# Read id_file and append a column based on TIR/NBD match
def update_nlr_domains(file_path, results, output_file):
    with open(output_file, 'w') as out_file:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                query_id = line.split('.')[0]
                nlr_domain = results.get(query_id, 'NA')
                out_file.write(query_id + '\t' + nlr_domain + '\n')

def main():
    parser = argparse.ArgumentParser(description="Extract results of structural alignments for NLR domains with custom cutoffs")
    
    # Positional arguments
    parser.add_argument("aln_file", help="Alignment file to process")
    parser.add_argument("id_file", help="ids of candidates")
    parser.add_argument("output_file", help="Output file to save updated results")
    
    # Optional arguments with defaults
    parser.add_argument("-p", "--prob_cutoff", type=float, default=1.0, help="Probability cutoff (default: 1.0)")
    parser.add_argument("-tc", "--tcov_cutoff", type=float, default=0.8, help="Target coverage cutoff (default: 0.8)")
    parser.add_argument("-tm", "--tm_cutoff", type=float, default=0.5, help="TM-score cutoff (default: 0.5)")
    parser.add_argument("-e", "--evalue_cutoff", type=float, default=0.001, help="E-value cutoff (default: 0.001)")
    
    args = parser.parse_args()

    aln_file = args.aln_file
    id_file = args.id_file
    output_file = args.output_file

    prob_cutoff = args.prob_cutoff
    tcov_cutoff = args.tcov_cutoff
    tm_cutoff = args.tm_cutoff
    evalue_cutoff = args.evalue_cutoff
    
    # Extract results from aln file
    results = extract_results_from_aln(aln_file, prob_cutoff, tcov_cutoff, tm_cutoff, evalue_cutoff)

    # Update NLR types based on extracted results and write to the output file
    update_nlr_domains(id_file, results, output_file)

if __name__ == '__main__':
    main()
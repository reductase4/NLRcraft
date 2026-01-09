# python ../scripts/rm_FPs.py -i ../03-foldseek_v2/results_cutoff_prob0.99_tcov0.8_tm0.5_filter/fragment_label/results_all_complete.txt
# -p ../06_remove_false_positives/results_rf/all_data_predicted.txt -o results_all_rm_FPs.txt
import argparse
import csv

def remove_false_positives(results_file, predicted_file, output_file):
    # read IDs of non-NLR homologs
    non_nlr_ids = set()
    with open(predicted_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['predicted'] == 'non-NLR':
                # remove '.pdb' & save IDs
                protein_id = row['query'].split('.')[0]
                non_nlr_ids.add(protein_id)

    # Replace results predicted as non-NLR homologs
    with open(results_file, 'r') as infile, open(output_file, 'w') as outfile, open("NLR_ids.txt", 'w') as outfile2:
        for line in infile:
            line = line.rstrip('\n')
            row = line.split('\t')
            
            # replace results
            if len(row) >= 2 and row[0] in non_nlr_ids:
                original = row[1]
                if original == 'N':
                    row[1] = 'Na'
                elif original == 'TN':
                    row[1] = 'T'
            
            # save
            new_line = '\t'.join(row) + '\n'
            outfile.write(new_line)

            if row[1] != 'Na':
                outfile2.write(row[0] + '\n')

def main():
    parser = argparse.ArgumentParser(description='Process protein prediction results')
    parser.add_argument('-i', "--resluts_file", required=True, help='Path to results_all.txt file')
    parser.add_argument('-p', "--predict_file", required=True, help='Path to all_data_predicted.txt file')
    parser.add_argument('-o', "--output_file", required=True, help='Path to output file')
    
    args = parser.parse_args()
    
    results = args.resluts_file
    predict  = args.predict_file
    output = args.output_file

    remove_false_positives(results, predict, output)

if __name__ == "__main__":
    main()
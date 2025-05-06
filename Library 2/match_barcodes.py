import numpy as np
import pandas as pd
import re
import os

def check_continuous_substrings(strings, max_consecutive=8):
    """
    Check for barcodes with long homopolymer repeats that may be sequencing errors.
    
    Parameters:
    -----------
    strings : list
        List of barcode sequences to check
    max_consecutive : int
        Maximum allowed consecutive identical bases
        
    Returns:
    --------
    numpy.ndarray
        Boolean array indicating sequences with homopolymer repeats
    """
    results = []
    for string in strings:
        consecutive_count = 1
        has_consecutive_substring = False
        for i in range(1, len(string)):
            if string[i] == string[i-1]:
                consecutive_count += 1
                if consecutive_count > max_consecutive:
                    has_consecutive_substring = True
                    break
            else:
                consecutive_count = 1
        results.append(has_consecutive_substring)
    return np.array(results)

def match_barcodes(barcode_table_path, mapped_bcs_path, effector_details_path, 
                  output_dir="./data", min_reads_per_bin=5, min_reads_per_rep=100):
    """
    Match barcodes between tables, identify effector combinations, and calculate enrichment scores.
    
    Parameters:
    -----------
    barcode_table_path : str
        Path to the barcode count table (BCTable.csv)
    mapped_bcs_path : str
        Path to the mapped barcodes file from Nanopore data
    effector_details_path : str
        Path to the effector details CSV file
    output_dir : str
        Directory to save output files
    min_reads_per_bin : int
        Minimum number of reads per bin to include in analysis
    min_reads_per_rep : int
        Minimum number of reads per replicate to include in analysis
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Reading in effector details and creating a list sorted by length
    efDeets = pd.read_csv(effector_details_path)
    efDeets.sort_values('Length (aa)', inplace=True)
    sortedEffs = list(efDeets['Protein'])
    efDeets = efDeets.set_index('Protein')
    efDeets['Mods'] = efDeets['Mods'].astype(str)
    
    # Reading in barcode to effector pairings
    mappedBCs = pd.read_csv(mapped_bcs_path, index_col=2, usecols=[1, 2, 3])
    
    # Reading in BC counts table
    countedBCs = pd.read_csv(barcode_table_path, index_col=0)
    
    # Joining counted BCs to mapped BCs, keeping only BCs that were mapped
    foundBCs = countedBCs.join(mappedBCs, how='inner')
    
    # Dropping BCs with homopolymer repeats
    foundBCs = foundBCs.loc[~check_continuous_substrings(foundBCs.index.to_list())]
    
    # Making EF1 and EF2 categorical variables sorted by length
    foundBCs['EF1'] = pd.Categorical(foundBCs['EF1'], categories=sortedEffs, ordered=True)
    foundBCs['EF2'] = pd.Categorical(foundBCs['EF2'], categories=sortedEffs, ordered=True)
    foundBCs.sort_values(['EF1', 'EF2', 'sum'], inplace=True)
    
    # Making indices EF1, EF2, then BC
    foundBCs = foundBCs.reset_index().set_index(['EF1', 'EF2', 'BC'])
    
    # Summing over EF1 and EF2 combinations
    EFTable = foundBCs.groupby(level=['EF1', 'EF2']).sum()
    
    # Save raw counts output
    EFTable.to_csv(f"{output_dir}/raw_counts.csv")
    
    # Dict for results
    results = {}
    
    # Parameters for enrichment calculation
    gene = ['81']  # Only CD81 data
    timepoint = ['5ddox', '12dmem']
    
    # Calculating enrichment scores for each experimental condition
    for g in gene:
        for t in timepoint:
            # Gets column names for specific gene and timepoint
            pattern = re.compile(g+r"rep\w+_"+t+"_")
            matches = [name for name in EFTable.columns.to_list() if pattern.match(name)]
            
            # Makes subtable for experimental condition
            sub = EFTable.loc[:, matches]
            
            # Renaming columns to be less wordy
            column_names = sub.columns.to_list()
            new_column_names = [col.split('_')[0][-4:] + '_' + col.split('_')[2] for col in column_names]
            sub.columns = new_column_names
            
            # Dropping pairs with < min_reads_per_bin in any bin
            sub = sub.loc[(sub >= min_reads_per_bin).all(axis=1), :]
            
            # Summing reads from each rep
            sub['sum1'] = sub.loc[:, ['rep1_HIGH', 'rep1_LOW']].sum(axis=1)
            sub['sum2'] = sub.loc[:, ['rep2_HIGH', 'rep2_LOW']].sum(axis=1)
            
            # Dropping pairs with < min_reads_per_rep then normalizing across each bin
            sub = sub.loc[(sub['sum1'] >= min_reads_per_rep) & (sub['sum2'] >= min_reads_per_rep), new_column_names]
            sub = sub.div(sub.sum(axis=0), axis=1)
            
            # Calculating HIGH/LOW ratio and geometric average
            sub['rep1_R'] = sub['rep1_HIGH'] / sub['rep1_LOW']
            sub['rep2_R'] = sub['rep2_HIGH'] / sub['rep2_LOW']
            sub['avg_R'] = (sub['rep1_R'] * sub['rep2_R']) ** (1/2)
            
            # Adding to results dict
            results[(g, t)] = sub
            
            # Saving normalized data with appropriate filename
            if t == '5ddox':
                day_label = 'Day0'
            else:
                day_label = 'Day12'
            sub.to_csv(f"{output_dir}/CD{g}_{day_label}_normalized.csv")
    
    print(f"Successfully matched barcodes and calculated enrichment scores.")
    print(f"Found {len(foundBCs.index.get_level_values(2).unique())} unique barcodes")
    print(f"Found {len(EFTable)} unique effector combinations")
    print(f"Output files saved to {output_dir}/")
    
    return EFTable, results

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Match barcodes between tables and calculate enrichment scores.')
    parser.add_argument('--barcode_table', required=True, help='Path to barcode count table (BCTable.csv)')
    parser.add_argument('--mapped_bcs', required=True, help='Path to mapped barcodes file from Nanopore data')
    parser.add_argument('--effector_details', required=True, help='Path to effector details CSV file')
    parser.add_argument('--output_dir', default='./data', help='Directory to save output files')
    parser.add_argument('--min_reads_per_bin', type=int, default=5, help='Minimum reads per bin')
    parser.add_argument('--min_reads_per_rep', type=int, default=100, help='Minimum reads per replicate')
    
    args = parser.parse_args()
    
    EFTable, results = match_barcodes(
        args.barcode_table,
        args.mapped_bcs,
        args.effector_details,
        args.output_dir,
        args.min_reads_per_bin,
        args.min_reads_per_rep
    )
import numpy as np
import pandas as pd
import subprocess
import os
import time
from pathlib import Path
from IPython.display import clear_output
import shutil
from collections import Counter
import tempfile

# Create a temporary directory for intermediate files
temp_dir = tempfile.TemporaryDirectory()
temp_path = temp_dir.name

# Functions for mapping

# defining a function to map EF1 sequences
def mapEF1(fastq_file, max_error_rate=0.25):
    # Define paths for temporary files
    cutadapt_output = os.path.join(temp_path, 'cutadaptOut.fastq')
    minimap_output = os.path.join(temp_path, 'minimapOut.paf')

    # finding and cutting out EF1 sequence
    subprocess.run(['cutadapt',
                    '-g', 'GAACAAGTACAGAACCCTCTGAGGGGTCTCCAGGA;max_error_rate=' + str(max_error_rate) +
                    '...GGGAGCGGCAGCGAGACACCAGGAACAAGCGAGTCAGCAACACCAGAGTCAGAGGGGTCTCCAGGA;max_error_rate=' + str(max_error_rate),
                    '-o', cutadapt_output,
                    '--trimmed-only',
                    '--revcomp',
                    '--minimum-length', '200',
                    '--maximum-length', '3500',
                    '--cores=0',
                    '--quiet',
                    fastq_file])

    # mapping EF1 to individual effector sequences
    f = open(minimap_output, 'w')
    subprocess.call(['minimap2',
                     'library/Library2_sequences.fasta',
                     cutadapt_output,
                     '--for-only',
                     '--secondary=no',
                     '-x', 'map-ont',
                     '-t', '31'],
                    stdout=f)

    # assigning EF1 identies to read names
    try:
        EF1_aln = pd.read_csv(minimap_output, sep='\t', header=None, index_col=0)
    except:
        return None
    EF1_aln['%ID'] = (EF1_aln[9] / EF1_aln[6]) * 100
    EF1 = EF1_aln.loc[EF1_aln['%ID'] > (100 - max_error_rate * 100)][5]
    EF1.rename('EF1', inplace=True)
    EF1.rename_axis('Read', inplace=True)

    # returning series
    return EF1


# defining a function to map EF2 sequences
def mapEF2(fastq_file, max_error_rate=0.25):
    # Define paths for temporary files
    cutadapt_output = os.path.join(temp_path, 'cutadaptOut.fastq')
    minimap_output = os.path.join(temp_path, 'minimapOut.paf')

    # finding and cutting out EF2 sequence
    subprocess.run(['cutadapt',
                    '-g', 'GGGAGCGGCAGCGAGACACCAGGAACAAGCGAGTCAGCAACACCAGAGTCAGAGGGGTCTCCAGGA;max_error_rate=' + str(max_error_rate) +
                    '...GGGAGCTGAGACCGAGTGATCTCAGGCTAATAACTCGTC;max_error_rate=' + str(max_error_rate),
                    '-o', cutadapt_output,
                    '--trimmed-only',
                    '--revcomp',
                    '--minimum-length', '200',
                    '--maximum-length', '3500',
                    '--cores=0',
                    '--quiet',
                    fastq_file])

    # mapping EF2 to individual effector sequences
    f = open(minimap_output, 'w')
    subprocess.call(['minimap2',
                     'library/Library2_sequences.fasta',
                     cutadapt_output,
                     '--for-only',
                     '--secondary=no',
                     '-x', 'map-ont',
                     '-t', '31'],
                    stdout=f)

    # assigning EF2 identies to read names
    try:
        EF2_aln = pd.read_csv(minimap_output, sep='\t', header=None, index_col=0)
    except:
        return None
    EF2_aln['%ID'] = (EF2_aln[9] / EF2_aln[6]) * 100
    EF2 = EF2_aln.loc[EF2_aln['%ID'] > (100 - max_error_rate * 100)][5]
    EF2.rename('EF2', inplace=True)
    EF2.rename_axis('Read', inplace=True)

    # returning series
    return EF2


# defining a function to find barcode sequences
def findBCs(fastq_file, max_error_rate=0.25):
    # Define path for temporary file
    cutadapt_output = os.path.join(temp_path, 'cutadaptOut.fastq')

    # finding barcodes
    subprocess.run(['cutadapt',
                    '-g', 'GGGAGCTGAGACCGAGTGATCTCAGGCTAATAACTCGTC;max_error_rate=' + str(max_error_rate) +
                    '...GCTTATGGCGAGTCTACTAATCCTC;max_error_rate=' + str(max_error_rate),
                    '-o', cutadapt_output,
                    '--trimmed-only',
                    '--revcomp',
                    '--minimum-length', '10',
                    '--maximum-length', '30',
                    '--cores=0',
                    '--quiet',
                    fastq_file, ])

    # assinging barcodes to read IDs
    try:
        BCs = pd.read_csv(cutadapt_output, sep='\\n@', header=None, engine='python')
    except:
        return None
    BCs = pd.DataFrame(BCs.values.reshape(-1, 4)[:, 0:2], columns=['Read', 'BC'])
    BCs['Read'] = BCs['Read'].str.split(pat='@| ').str[1]
    BCs.set_index('Read', inplace=True)

    # returning dataframe
    return BCs


# defining a function to map barcodes to effector combos by read ID
def mapBCs(EF1, EF2, BCs):
    # joining EF1 and EF2
    combs = EF1.to_frame().join(EF2, how='inner')

    # joining combinations to barcodes
    mappedBCs = combs.join(BCs, how='inner')

    # returning dataframe
    return mappedBCs


# Ensure mapped directory exists
os.makedirs('fastqs/mapped', exist_ok=True)

# Mapping BCs

# getting all fastq files in directory
files = list(Path('fastqs/nanopore').rglob('*.fastq.gz'))

# reading mappedBCs if it exists or creating a new dataframe to store barcodes
try:
    mappedBCs = pd.read_csv('mappedBCs.csv', index_col=0)
except:
    mappedBCs = pd.DataFrame()

# iterating through each file, mapping barcodes and appending to data frame
# moving files to mapped folder and saving mappedBCs along the way
for i, file in enumerate(files):
    print('Mapping file ' + str(i + 1) + '/' + str(len(files)))
    EF1 = mapEF1(file)
    EF2 = mapEF2(file)
    BCs = findBCs(file)
    if EF1 is None or EF2 is None or BCs is None:
        shutil.move(file, 'fastqs/mapped')
        clear_output()
        continue
    mappedBCs = pd.concat((mappedBCs, mapBCs(EF1, EF2, BCs)), axis=0)
    shutil.move(file, 'fastqs/mapped')
    #mappedBCs.to_csv('mappedBCs.csv')
    clear_output()

# removing any full duplicates since there were some duplication events with pod5 conversion or partial basecalling etc
mappedBCs = mappedBCs.reset_index().drop_duplicates().set_index('Read')
#mappedBCs.to_csv('mappedBCs.csv')

# --------------------------------
# Additional Filtering Section
# --------------------------------

print("Starting additional filtering steps...")

# Get the list of effectors and sort by length
try:
    effectors = pd.read_csv('library/Library2_members.csv')
    effectors.sort_values('Length (aa)', inplace=True)
    sortedEffs = list(effectors['Protein'])
except:
    print("Warning: Could not find effector details file. Using default sorting.")
    # If effector details file is not available, just use the existing effectors in the mappedBCs
    sortedEffs = sorted(list(set(mappedBCs['EF1'].unique()).union(set(mappedBCs['EF2'].unique()))))

# 1. Delete reads that mapped to more than one effector pair
print(f"Before removing duplicate read mappings: {len(mappedBCs)} entries")
readIDs = list(mappedBCs.index)
cnt = Counter(readIDs)
doubles = [k for k, v in cnt.items() if v > 1]
mappedBCs.drop(doubles, inplace=True)
print(f"After removing duplicate read mappings: {len(mappedBCs)} entries")

# 2. Sort by EF1 and EF2 lengths
mappedBCs['EF1'] = pd.Categorical(mappedBCs['EF1'], categories=sortedEffs, ordered=True)
mappedBCs['EF2'] = pd.Categorical(mappedBCs['EF2'], categories=sortedEffs, ordered=True)
mappedBCs.sort_values(['EF1', 'EF2', 'BC'], inplace=True)

# 3. Get rid of multiple observations of the same BC and EF combo
print(f"Before removing duplicate BC-EF combinations: {len(mappedBCs)} entries")
mappedBCs.drop_duplicates(inplace=True)
print(f"After removing duplicate BC-EF combinations: {len(mappedBCs)} entries")

# 4. Get rid of cases where the same BC maps to multiple effectors
print(f"Before removing barcodes mapping to multiple effectors: {len(mappedBCs)} entries")
mappedBCs.drop_duplicates(['BC'], keep=False, inplace=True)
print(f"After removing barcodes mapping to multiple effectors: {len(mappedBCs)} entries")

# Save the filtered mappedBCs
mappedBCs.to_csv('mappedBCs_filt.csv')
print(f"Filtering complete. {len(mappedBCs)} entries remain in the filtered dataset.")
print("Filtered barcode mapping saved to mappedBCs_filt.csv")

# Clean up temporary directory
temp_dir.cleanup()
print("Temporary files cleaned up. Process complete.")
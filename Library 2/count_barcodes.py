import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os
import glob
import multiprocessing as mp

# Merge reads for each sample

# getting directories in fastq folder
path = 'fastqs/unmerged'
directories = [path + '/' + d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

# merging function
def merge(d):
    # getting names of fastq files
    files = glob.glob(os.path.join(d, '*R*_001.fastq.gz'))
    
    # creating output path and file name
    out = 'fastqs/merged/' + '_'.join(files[0].split('/')[3].split('_')[0:3]) + '.fastq.gz'

    # merging the files
    subprocess.run(['fastp',
                    '-i', files[0],
                    '-I', files[1],
                    '-m', '-c', '-w', '16',
                    '--overlap_len_require', '20',
                    '--merged_out', out],
                   check=True)

    return None

# iterating through the directories and merging reads
if __name__ == '__main__':
    # creating process pool only 2 because each will take 16 threads
    pool = mp.Pool(2)
    
    # merging reads for each directory
    for d in directories:
        pool.apply_async(merge, (d,))
        
    # closing and joining
    pool.close()
    pool.join()

# Extracting barcode sequences

# getting merged fastq files
path = 'fastqs/merged'
files = glob.glob(os.path.join(path, '*.fastq.gz'))

# iterating through each merged fastq and cutting reads which is already parallelized
# trimming function
def trim(f):
    # creating output path and file name
    out = 'fastqs/trimmed/' + f.split('/')[-1]

    # trimming the primers
    subprocess.run(['cutadapt',
                    '-g', 'AGTGATCTCAGGCTAATAACTCGTC;max_error_rate=0' +
                    '...GCTTATGGCGAGTCTACTAATCCTC;max_error_rate=0',
                    '-o', out,
                    '--trimmed-only',
                    '--revcomp',
                    '--minimum-length', '10',
                    '--maximum-length', '30',
                    '--cores=32',
                    f])

    return None

# iterating through the directories and merging reads
if __name__ == '__main__':
    # iterating through the files
    for f in files:
        trim(f)

# Counting barcodes

# getting trimmed in fastq file names
path = 'fastqs/trimmed'
files = glob.glob(os.path.join(path, '*.fastq.gz'))

# merging function
def count(f):
    # creating output path and file name
    out = 'fastqs/counts/' + f.split('/')[2].split('.')[0] + '.tsv'

    # create a subprocess to execute the commands
    with open(out, 'w') as o:
        p1 = subprocess.Popen(['seqfu',
                               'derep',
                               f],
                              stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['seqfu',
                               'tabulate'],
                              stdin=p1.stdout,
                              stdout=o)
        p1.stdout.close()  # close the output of the first command to prevent deadlock
        p2.communicate()   # wait for the second command to finish

    return None

# iterating through the files
if __name__ == '__main__':
    # creating process pool
    pool = mp.Pool(8)
    
    # counting
    for f in files:
        pool.apply_async(count, (f,))
        
    # closing and joining
    pool.close()
    pool.join()

# Creating barcode table

# getting count file names
path = 'fastqs/counts'
files = np.sort(glob.glob(os.path.join(path, '*.tsv')))

# creating BCTable dataframe
BCTable = pd.DataFrame()

# iterating through each count file and adding the to table
for f in files:
    # reading in tsv
    counts = pd.read_csv(f, sep='\t', usecols=[0, 2], header=None, names=['name', 'BC'])

    # sample ID
    samp = f.split('/')[-1].split('.')[0]

    # finding the count info
    counts[samp] = counts['name'].str.split(pat=';', expand=True)[1]
    counts[samp] = pd.to_numeric(counts[samp].str.split(pat='=', expand=True)[1])
    counts = counts.drop(index=0)
    counts = counts.drop(columns='name')

    # setting the index
    counts = counts.set_index('BC')
    
    # adding to the table
    BCTable = BCTable.join(counts, how='outer')

# sorting table by overall sum of each BC
BCTable['sum'] = BCTable.iloc[:, :16].sum(axis=1)
BCTable = BCTable.sort_values(by='sum', ascending=False)

# adding zeros for nans
BCTable = BCTable.fillna(value=0)

# saving table
BCTable.to_csv('BCTable.csv')

import numpy as np
import pandas as pd
import subprocess
import os
import glob
import multiprocessing as mp
import tempfile
import shutil

# Create temporary directories for intermediate files
temp_base = tempfile.mkdtemp()
temp_merged = os.path.join(temp_base, 'merged')
temp_trimmed = os.path.join(temp_base, 'trimmed')
temp_counts = os.path.join(temp_base, 'counts')
temp_fastp = os.path.join(temp_base, 'fastp')  # For fastp output files

# Create the temporary subdirectories
os.makedirs(temp_merged, exist_ok=True)
os.makedirs(temp_trimmed, exist_ok=True)
os.makedirs(temp_counts, exist_ok=True)
os.makedirs(temp_fastp, exist_ok=True)  # Create fastp output directory

# Merge reads for each sample

# getting directories in fastq folder
path = 'fastqs/NGS'
if not os.path.exists(path):
    print(f"Warning: Directory {path} does not exist. Creating it.")
    os.makedirs(path, exist_ok=True)
    
directories = [path + '/' + d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

# merging function
def merge(d):
    # getting names of fastq files
    files = glob.glob(os.path.join(d, '*R*_001.fastq.gz'))
    
    # Skip if no matching files found
    if not files:
        print(f"No paired-end fastq files found in {d}")
        return None
    
    # creating output path and file name
    try:
        filename = '_'.join(files[0].split('/')[3].split('_')[0:3])
    except IndexError:
        # If path structure is different, use directory name as fallback
        filename = os.path.basename(d)
    
    out = os.path.join(temp_merged, f'{filename}.fastq.gz')
    
    # Define fastp output files
    html_report = os.path.join(temp_fastp, f'{filename}_fastp.html')
    json_report = os.path.join(temp_fastp, f'{filename}_fastp.json')

    # merging the files
    try:
        subprocess.run(['fastp',
                        '-i', files[0],
                        '-I', files[1],
                        '-m', '-c', '-w', '16',
                        '--overlap_len_require', '20',
                        '--merged_out', out,
                        '--html', html_report,  # Redirect HTML report
                        '--json', json_report,  # Redirect JSON report
                        '--report_title', f"Merge report for {filename}"],
                       check=True)
        print(f"Successfully merged files from {d}")
    except subprocess.CalledProcessError as e:
        print(f"Error merging files from {d}: {e}")
    except Exception as e:
        print(f"Unexpected error merging files from {d}: {e}")

    return None

# iterating through the directories and merging reads
if __name__ == '__main__':
    print("Starting paired-end read merging...")
    if not directories:
        print("No sample directories found. Skipping merge step.")
    else:
        # creating process pool only 2 because each will take 16 threads
        pool = mp.Pool(2)
        
        # merging reads for each directory
        for d in directories:
            pool.apply_async(merge, (d,))
            
        # closing and joining
        pool.close()
        pool.join()
        
        print("Merging complete.")

    # Extracting barcode sequences

    # getting merged fastq files
    merged_files = glob.glob(os.path.join(temp_merged, '*.fastq.gz'))
    
    if not merged_files:
        print("No merged files found. Looking for pre-merged files in fastqs/merged...")
        if os.path.exists('fastqs/merged'):
            original_files = glob.glob(os.path.join('fastqs/merged', '*.fastq.gz'))
            # Copy files to temp directory
            for f in original_files:
                shutil.copy(f, temp_merged)
            merged_files = glob.glob(os.path.join(temp_merged, '*.fastq.gz'))

    # trimming function
    def trim(f):
        # creating output path and file name
        out = os.path.join(temp_trimmed, os.path.basename(f))

        # trimming the primers
        try:
            subprocess.run(['cutadapt',
                            '-g', 'AGTGATCTCAGGCTAATAACTCGTC;max_error_rate=0' +
                            '...GCTTATGGCGAGTCTACTAATCCTC;max_error_rate=0',
                            '-o', out,
                            '--trimmed-only',
                            '--revcomp',
                            '--minimum-length', '10',
                            '--maximum-length', '30',
                            '--cores=32',
                            f],
                          check=True)
            print(f"Successfully trimmed {f}")
        except subprocess.CalledProcessError as e:
            print(f"Error trimming {f}: {e}")
        except Exception as e:
            print(f"Unexpected error trimming {f}: {e}")

        return None

    # iterating through the merged files for trimming
    print("Starting barcode extraction...")
    if not merged_files:
        print("No files found for barcode extraction. Skipping trim step.")
    else:
        for f in merged_files:
            trim(f)
        
        print("Barcode extraction complete.")

    # Counting barcodes

    # getting trimmed fastq file names
    trimmed_files = glob.glob(os.path.join(temp_trimmed, '*.fastq.gz'))

    # count function
    def count(f):
        # creating output path and file name
        filename = os.path.basename(f).split('.')[0]
        out = os.path.join(temp_counts, f'{filename}.tsv')

        # create a subprocess to execute the commands
        try:
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
            print(f"Successfully counted barcodes in {f}")
        except Exception as e:
            print(f"Error counting barcodes in {f}: {e}")

        return None

    # iterating through the files for counting
    print("Starting barcode counting...")
    if not trimmed_files:
        print("No trimmed files found. Skipping count step.")
    else:
        # creating process pool
        pool = mp.Pool(8)
        
        # counting
        for f in trimmed_files:
            pool.apply_async(count, (f,))
            
        # closing and joining
        pool.close()
        pool.join()
        
        print("Barcode counting complete.")

    # Creating barcode table

    # getting count file names - use regular list instead of numpy array
    count_files = glob.glob(os.path.join(temp_counts, '*.tsv'))
    count_files.sort()  # Sort in place

    # creating BCTable dataframe
    BCTable = pd.DataFrame()

    print("Creating barcode table...")
    if not count_files:
        print("No count files found. Cannot create barcode table.")
    else:
        # iterating through each count file and adding to the table
        for f in count_files:
            try:
                # reading in tsv
                counts = pd.read_csv(f, sep='\t', usecols=[0, 2], header=None, names=['name', 'BC'])

                # sample ID
                samp = os.path.basename(f).split('.')[0]

                # finding the count info
                counts[samp] = counts['name'].str.split(pat=';', expand=True)[1]
                counts[samp] = pd.to_numeric(counts[samp].str.split(pat='=', expand=True)[1])
                counts = counts.drop(index=0)
                counts = counts.drop(columns='name')

                # setting the index
                counts = counts.set_index('BC')
                
                # adding to the table
                BCTable = BCTable.join(counts, how='outer')
                
                print(f"Added counts from {f} to barcode table")
            except Exception as e:
                print(f"Error processing count file {f}: {e}")

        # Check if we have data in the table
        if BCTable.empty:
            print("Warning: No data found in any count files. BCTable is empty.")
        else:
            # Calculate the sum of all columns for sorting
            n_cols = min(16, BCTable.shape[1])  # Use actual number of columns if less than 16
            BCTable['sum'] = BCTable.iloc[:, :n_cols].sum(axis=1)
            BCTable = BCTable.sort_values(by='sum', ascending=False)

            # adding zeros for nans
            BCTable = BCTable.fillna(value=0)

            # saving table
            output_file = 'BCTable.csv'
            BCTable.to_csv(output_file)
            print(f"Barcode table created with {len(BCTable)} barcodes and saved to {output_file}")

    # Clean up temporary directory
    try:
        shutil.rmtree(temp_base)
        print("Temporary files cleaned up.")
    except Exception as e:
        print(f"Warning: Could not clean up temporary files: {e}")


import numpy as np
import pandas as pd
import subprocess
import os
from pathlib import Path
from IPython.display import clear_output
import shutil
import re
import pandas as pd
import numpy as np
import re
from pathlib import Path
from typing import Tuple, List

class BatchSequenceMapper:
    def __init__(self):
        # Previous initialization code remains the same
        self.cutadapt_path = str(Path('~/usr/bin/cutadapt').expanduser())
        self.minimap2_path = str(Path('./minimap2-2.28_x64-linux/minimap2').expanduser())
        self.ref_path_r1 = str(Path('./R1_Forward_P5_301cycles.fasta').expanduser())
        self.ref_path_r2 = str(Path('./R2_Reverse_P7_301cycles.fasta').expanduser())
        
        os.makedirs('mapped', exist_ok=True)
        os.makedirs('output_csv', exist_ok=True)
        os.makedirs('debug_logs', exist_ok=True)

    def determine_read_type(self, filename):
        """
        Determine if a file is R1 or R2 NGS read based on the pattern before '001.fastq.gz'
        
        Args:
            filename: Path object or string of the filename
            
        Returns:
            bool: True if R1, False if R2, None if undetermined
        """
        # Look specifically for _R1_001 or _R2_001 at the end of the filename
        name = str(filename)
        if '_R1_001.fastq.gz' in name:
            return True
        elif '_R2_001.fastq.gz' in name:
            return False
        return None

    def map_sequences(self, fastq_file, is_r1=True, max_error_rate=0.25):
        """
        Map sequences using minimap2
        """
        # Create unique PAF file names based on the input file
        base_name = Path(fastq_file).stem.replace('.fastq', '')
        read_type = "R1" if is_r1 else "R2"
        output_paf = f'minimapOut_{base_name}_{read_type}.paf'
        debug_log = Path('debug_logs') / f"{base_name}_{read_type}_debug.log"
        
        # Select appropriate reference
        ref_path = self.ref_path_r1 if is_r1 else self.ref_path_r2
        
        # Log the mapping operation
        with open(debug_log, 'w') as log:
            log.write(f"Processing file: {fastq_file}\n")
            log.write(f"Using reference: {ref_path}\n")
            log.write(f"Output PAF: {output_paf}\n")
            log.write(f"Read type: {read_type}\n")
        
        # Rest of the mapping code remains the same
        with open(output_paf, 'w') as f:
            subprocess.call([
                self.minimap2_path,
                ref_path,
                str(fastq_file),
                '--for-only',
                '--secondary=no',
                '-x', 'map-iclr',
                '-t', '31'
                ], stdout=f)
        
        try:
            aln = pd.read_csv(output_paf, sep='\t', header=None, index_col=0)
            aln['%ID'] = (aln[9] / aln[6]) * 100
            ef = aln.loc[aln['%ID'] > (100 - max_error_rate * 100)][5]
            ef.rename('EF1' if is_r1 else 'EF2', inplace=True)
            ef.rename_axis('Read', inplace=True)
            
            with open(debug_log, 'a') as log:
                log.write(f"\nFound {len(ef)} mapped sequences\n")
                if len(ef) > 0:
                    log.write(f"Example mapping: {ef.iloc[0]}\n")
            
            return ef
            
        except Exception as e:
            with open(debug_log, 'a') as log:
                log.write(f"\nError processing alignments: {str(e)}\n")
            print(f"Error processing alignments for {fastq_file.name}: {str(e)}")
            return None

    def process_all_files(self, input_dir='raw_fastq/'):
        """
        Process all FASTQ files in the input directory
        """
        fastq_files = list(Path(input_dir).rglob('*.fastq.gz'))
        
        for i, fastq_file in enumerate(fastq_files, 1):
            print(f'Processing file {i}/{len(fastq_files)}: {fastq_file.name}')
            
            # Use the new determine_read_type function
            is_r1 = self.determine_read_type(fastq_file.name)
            
            if is_r1 is None:
                print(f"Warning: Cannot determine if {fastq_file.name} is R1 or R2 NGS read")
                continue
                
            # Create output CSV name based on input filename
            output_csv = Path('output_csv') / f"{fastq_file.stem.replace('.fastq', '')}.csv"
            
            print(f"Processing as {'R1' if is_r1 else 'R2'} NGS read")
            ef = self.map_sequences(fastq_file, is_r1=is_r1)
            
            if ef is not None:
                ef.to_csv(output_csv)
                print(f"Saved results to {output_csv} with {len(ef)} sequences")
            
            shutil.move(fastq_file, Path('mapped') / fastq_file.name)
            print(f"Moved {fastq_file.name} to mapped directory")
            clear_output(wait=True)

def main():
    mapper = BatchSequenceMapper()
    print("Starting batch processing...")
    mapper.process_all_files()
    print("Processing complete!")
    
    output_files = list(Path('output_csv').glob('*.csv'))
    print(f"\nProcessed {len(output_files)} files:")
    for file in output_files:
        stats = os.path.getsize(file)
        print(f"- {file.name} ({stats} bytes)")

if __name__ == "__main__":
    main()

class PairedCSVProcessor:
    def __init__(self, input_dir: str = 'output_csv', output_dir: str = 'count_matrices'):
        """
        Initialize the processor
        
        Args:
            input_dir: Directory containing input CSV files
            output_dir: Directory to save count matrices
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
    def find_paired_files(self) -> List[Tuple[Path, Path]]:
        """
        Find pairs of R1 and R2 CSV files that match based on their prefix
        """
        # Get all CSV files
        csv_files = list(self.input_dir.glob('*.csv'))
        
        # Dictionary to store file pairs
        file_pairs = []
        
        # Group files by their common prefix
        for file in csv_files:
            # Extract the common part of the filename (everything before R1/R2)
            match = re.match(r'(.+?_L\d+_)R[12]', file.stem)
            if not match:
                continue
                
            common_prefix = match.group(1)
            r1_file = self.input_dir / f"{common_prefix}R1_001.csv"
            r2_file = self.input_dir / f"{common_prefix}R2_001.csv"
            
            if r1_file.exists() and r2_file.exists():
                if (r1_file, r2_file) not in file_pairs:
                    file_pairs.append((r1_file, r2_file))
        
        return file_pairs

    def extract_number_from_ef1(self, ef1_string: str, vector_type: int) -> int:
        """
        Extract number from EF1 string based on vector type
        """
        pattern = f"KanR_receiver_vector{vector_type}_-_(\d+)_"
        match = re.search(pattern, ef1_string)
        if match:
            return int(match.group(1))
        return None

    def create_count_matrix(self, file1: Path, file2: Path, matrix_size: int = 156) -> np.ndarray:
        """
        Process paired CSV files and create count matrix
        """
        # Read CSV files
        df1 = pd.read_csv(file1)
        df2 = pd.read_csv(file2)

        # Initialize count matrix
        count_matrix = np.zeros((matrix_size, matrix_size), dtype=int)

        # Merge dataframes on 'Read' column
        merged_df = pd.merge(df1, df2, on='Read')

        # Process each matched pair
        for _, row in merged_df.iterrows():
            # Extract numbers from EF1 and EF2 columns (note the change here)
            num1 = self.extract_number_from_ef1(row['EF1'], 1)  # From R1 file
            num2 = self.extract_number_from_ef1(row['EF2'], 2)  # From R2 file

            # Update count matrix if both numbers are valid
            if num1 is not None and num2 is not None:
                if 1 <= num1 <= matrix_size and 1 <= num2 <= matrix_size:
                    count_matrix[num1-1][num2-1] += 1

        return count_matrix

    def process_all_pairs(self):
        """
        Process all paired CSV files and generate count matrices
        """
        # Find all paired files
        file_pairs = self.find_paired_files()
        
        if not file_pairs:
            print("No paired CSV files found!")
            return
        
        print(f"Found {len(file_pairs)} pairs of CSV files to process")
        
        # Process each pair
        for file1, file2 in file_pairs:
            # Extract common prefix for output filename
            common_prefix = re.match(r'(.+?_L\d+_)', file1.stem).group(1)
            output_file = self.output_dir / f"{common_prefix}count_matrix.csv"
            
            print(f"\nProcessing pair:")
            print(f"R1: {file1.name}")
            print(f"R2: {file2.name}")
            
            # Create count matrix
            count_matrix = self.create_count_matrix(file1, file2)
            
            # Convert to DataFrame for easier saving
            df_matrix = pd.DataFrame(
                count_matrix,
                index=[i+1 for i in range(156)],
                columns=[i+1 for i in range(156)]
            )
            
            # Save matrix to CSV
            df_matrix.to_csv(output_file)
            
            # Print statistics
            print(f"Created count matrix: {output_file}")
            print(f"Matrix shape: {count_matrix.shape}")
            print(f"Total counts: {count_matrix.sum()}")
            print(f"Non-zero entries: {np.count_nonzero(count_matrix)}")
            print(f"Maximum count in any cell: {count_matrix.max()}")

def main():
    # Create processor instance
    processor = PairedCSVProcessor()
    
    # Process all pairs
    processor.process_all_pairs()

if __name__ == "__main__":
    main()

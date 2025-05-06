# Combinatorial domain screening platform reveals epigenetic effector interactions for transcriptional perturbation

Additional data and codes to accompany the preprint: [https://www.biorxiv.org/content/10.1101/2024.10.28.620683v2](https://www.biorxiv.org/content/10.1101/2024.10.28.620683v2)

# NGS Read Mapping and Analysis Pipeline (Library 1)

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-Apache%202.0-green.svg)](LICENSE)

A Python-based pipeline for processing and analyzing paired-end NGS (Next-Generation Sequencing) reads. The pipeline includes functionality for sequence mapping, processing paired CSV files, and generating count matrices.

## 📋 Dependencies

### Tested on
- Python 3.10 
- minimap2 v2.28

### Python Packages
```bash
numpy
pandas
pathlib
IPython
```

## 🔧 Installation

1. Clone this repository:
```bash
git clone https://github.com/hsulab-arc/COMBINE.git
cd COMBINE
```

2. Install Python dependencies:
```bash
pip install numpy pandas pathlib IPython
```

3. Install minimap2:

Follow the installation instructions at: https://github.com/lh3/minimap2
```bash
curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf -
./minimap2-2.28_x64-linux/minimap2
```
   - Note: Move the minimap2-2.28_x64-linux folder to the Library 1 directory to run the example dataset.

## 📁 Directory Structure

```
Library 1/
├── minimap2-2.28_x64-linux/
├── raw_fastq/  
├── Lib1_bivalent_effectors.png
├── Minimap_batch_Lib1.py
├── R1_Forward_P5_301cycles.fasta
└── R2_Reverse_P7_301cycles.fasta           
```

## 📝 File Naming Convention

The pipeline expects paired-end read files to follow this naming convention:
- R1 files: `*_R1_001.fastq.gz`
- R2 files: `*_R2_001.fastq.gz`

## 📊 Output Files

1. Mapping Results:
   - CSV files containing mapping information for each input FASTQ file
   - Debug logs for each mapping operation

2. Count Matrices:
   - CSV files containing count matrices for each paired set of reads
   - Matrix dimensions: 156 x 156
   - Each cell represents the count of specific bivalent domain combinations

## Description of additional files

![Description](./Library%201/Lib1_bivalent_effectors.png)



## 🚀 Quick Start with Example Dataset

### Configure File Paths
Before running the pipeline, adjust the following file paths in the code according to your system setup if necessary:

```python
# In BatchSequenceMapper class:
self.cutadapt_path = str(Path('~/.local/bin/cutadapt').expanduser())
self.minimap2_path = str(Path('./minimap2-2.28_x64-linux/minimap2').expanduser())
self.ref_path_r1 = str(Path('./R1_Forward_P5_301cycles.fasta').expanduser())
self.ref_path_r2 = str(Path('./R2_Reverse_P7_301cycles.fasta').expanduser())
```

### Example Dataset
Example dataset is provided in `Library 1/raw_fastq`:
```
Library 1/raw_fastq/
├── subset_D6_R1_CD81_High_S1_L001_R1_001.fastq.gz
├── subset_D6_R1_CD81_High_S1_L001_R2_001.fastq.gz
├── subset_D6_R1_CD81_Repressed_S2_L001_R1_001.fastq.gz
└── subset_D6_R1_CD81_Repressed_S2_L001_R2_001.fastq.gz
```

### Running the Example
1. Run `python Minimap_batch_Lib1.py`
2. Check results in `count_matrices` folder (process takes ~10 seconds)

### Citations
Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191

# Library 2 Analysis

## Overview

Library 2 consists of bivalent epigenetic editors with unique 20bp barcodes. The analysis workflow involves three main steps:
1. Extracting barcodes from Nanopore sequencing data to establish barcode-to-editor mappings
2. Counting barcodes from NGS data to quantify representation of each construct
3. Matching barcodes between tables to create a comprehensive mapping of barcode counts to effector combinations

This pipeline processes data from screens targeting CD81, calculating enrichment scores to identify effective epigenetic editor combinations.

## Installation and Setup

Before running the analysis pipeline, you'll need to set up your environment with the necessary packages and tools.

### Clone the Repository

First, clone the COMBINE repository to your local machine:

```bash
git clone https://github.com/hsulab-arc/COMBINE.git
cd COMBINE
```

### Install Dependencies

The analysis pipeline requires Python 3.9 and several packages. It is recommended to install these in a dedicated conda environment:

```bash
# First create base environment with Python and standard packages
conda create -n combine python=3.9
conda activate combine
conda install -c bioconda -c conda-forge cutadapt minimap2 fastp seqfu
conda install numpy pandas scipy
```

### Directory Structure

Make sure your directory structure follows this pattern:

```
COMBINE/
├── bc_extraction.py
├── count_barcodes.py
├── match_barcodes.py
├── effector_details.csv
└── fastqs/
    ├── nanopore/ # For nanopore barcode mapping fastqs
    └── unmerged/ # For raw NGS fastqs
```

## Usage

### Step 1: Barcode Extraction from Nanopore Data

The `bc_extraction.py` script maps Nanopore reads to identify both effector domains (EF1 and EF2) and their associated barcodes.

```bash
python bc_extraction.py
```

#### Process

1. **EF1 Mapping**: 
   - Identifies the EF1 sequence using flanking regions with cutadapt
   - Maps the extracted sequence to reference effectors using minimap2
   - Assigns effector identities to reads based on sequence identity

2. **EF2 Mapping**:
   - Similar process to identify and map EF2 sequences
   - Uses different flanking regions specific to EF2

3. **Barcode Extraction**:
   - Identifies the barcode sequence between EF2 and a downstream constant region
   - Extracts 10-30bp sequences representing the barcodes
   - Associates barcodes with their read IDs

4. **Barcode-to-Effector Mapping**:
   - Joins EF1, EF2, and barcode information using read IDs
   - Creates a comprehensive mapping table (mappedBCs.csv)
   - Removes duplicates to ensure data integrity

This step creates `mappedBCs.csv` containing barcode-to-effector mappings.

### Step 2: Barcode Counting from NGS Data

The `count_barcodes.py` script processes NGS data to quantify the representation of each barcode across different samples.

```bash
python count_barcodes.py
```

#### Process

1. **Read Merging**:
   - Merges paired-end reads using fastp
   - Requires at least 20bp overlap between R1 and R2
   - Uses multithreading for efficient processing

2. **Barcode Extraction**:
   - Trims reads to isolate barcode sequences using cutadapt
   - Extracts sequences between consistent flanking regions
   - Filters for sequences between 10-30bp in length

3. **Barcode Counting**:
   - Dereplicated sequences using seqfu
   - Tabulates unique barcode sequences and their counts
   - Processes each sample independently

4. **Barcode Table Creation**:
   - Combines counts from all samples into a single table
   - Sorts barcodes by total abundance across samples
   - Generates a comprehensive barcode count table (BCTable.csv)

This step produces `BCTable.csv` with barcode counts across samples.

### Step 3: Barcode Matching and Enrichment Calculation

The `match_barcodes.py` script integrates the data from Steps 1 and 2 and calculates enrichment scores for effector combinations.

```bash
python match_barcodes.py --barcode_table BCTable.csv --mapped_bcs mappedBCs.csv --effector_details effector_details.csv --output_dir ./data
```

#### Process

1. **Read and Filter Data**:
   - Combine barcode counts with mapped barcodes from Nanopore data
   - Filter out barcodes with homopolymer repeats (potential sequencing errors)
   - Organize by effector combinations

2. **Create Effector Combination Matrix**:
   - Index by EF1, EF2, and barcode
   - Group by effector combinations to calculate total counts
   - Output comprehensive raw counts table

3. **Calculate Enrichment Scores**:
   - Filter combinations based on minimum read counts
   - Normalize read counts across bins
   - Calculate HIGH/LOW ratio for each replicate
   - Compute geometric mean of replicates

#### Output Files

The script generates multiple output files:
- `raw_counts.csv`: Count of each barcode for each effector combination
- `CD81_Day0_normalized.csv`: Normalized counts and enrichment scores for CD81 at 5 days with doxycycline
- `CD81_Day12_normalized.csv`: Normalized counts and enrichment scores for CD81 at 12 days memory

## Requirements

- Python 3.6+
- Required packages: numpy, pandas, scipy, statsmodels, regex
- External tools:
  - cutadapt
  - minimap2
  - fastp
  - seqfu

## Notes

- The barcode extraction process is optimized for the specific structure of Library 2 constructs
- The pipeline can process large datasets efficiently using parallelization
- Error rate parameters can be adjusted for different levels of sequence stringency
- Homopolymer filtering removes barcodes that may contain sequencing errors
- Effector combinations are sorted by protein length for consistent organization

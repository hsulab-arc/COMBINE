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

**# BAM_mtrx**

## Overview
BAM_mtrx is a Python application designed to analyze and visualize quality distributions of sequencing reads from BAM files. It processes BAM files to compute various metrics and generates histograms for base and mapping quality distributions.

## Features
- **GUI for BAM file selection**: Users can select their input file using a simple file dialog.
- **Comprehensive metrics computation**: Total reads, mapped reads, unmapped reads, and more.
- **Quality distribution plots**: Generates histograms of base and mapping qualities.
- **Command-line interface**: Supports customization through command-line arguments.

## Bioinformatics Context and Logic
BAM_mtrx processes BAM files, which are binary versions of SAM (Sequence Alignment/Map) files, to assess the quality of DNA sequence alignments. Below is a detailed explanation of the program's logic, particularly focusing on CIGAR string processing and alignment data handling:

### CIGAR String Handling
- **Identifying Deletions**: The program includes a function to check if a CIGAR string represents a deletion in the alignment. If a deletion is detected during the read processing, the read is either skipped or specially handled depending on its impact on the alignment metrics.
- **Efficiency**: The deletion check is designed for efficiency, returning immediately when a deletion is found, which optimizes performance when processing large datasets.

### Alignment Quality Analysis
- **Quality Score Collection**: For reads that pass the initial CIGAR check (i.e., no problematic deletions), the program collects base quality scores and mapping quality scores. These scores are crucial for assessing the reliability of the alignment.
- **High-Quality Reads Identification**: Reads with base quality scores of 30 or higher are flagged as high-quality reads. Similarly, reads with a mapping quality score of 60 (perfect mapping) are separately counted to provide insights into the dataset's overall alignment quality.

### Tested Data
BAM_mtrx was extensively tested using mitochondrial DNA data that was pre-processed (mapped, indexed, and sorted). This type of data provides a robust framework for evaluating the program's effectiveness in handling complex genomic datasets.

## Requirements
- Python 3
- Libraries: pysam, matplotlib, tkinter

## Installation
Ensure Python 3 is installed along with the required libraries:
```bash
pip install pysam matplotlib
```

## Usage
To use BAM_mtrx, follow these instructions:

1. Run the script from the command line:
```bash
python BAM_mtrx.py
```

2. A file dialog will appear for you to select your BAM file.

3. The program will process the selected BAM file and generate two plots showing the base and mapping quality distributions, which will be saved as PNG files.

4. A text file containing detailed metrics about the processed reads will also be generated.

### Command-line Arguments
- `--title` (optional): Specify a base title for the graphs.

## Output
PNG files: Plots showing the distributions of base and mapping qualities.These plots are saved as PNG images with filenames indicating their content (e.g., your_bam_file_base_quality.png).

Text file: Detailed metrics of the processed reads, including total reads, mapped reads, and high-quality reads.It contains various metrics like total reads, mapped reads, unmapped reads, average alignment score, and the number of reads with high base quality (above 30) and perfect mapping quality (score of 60).

## Note:
The script filters out reads with missing base quality scores to avoid errors. This filtering is reported in the metrics file.

## Contributing
Contributions to BAM_mtrx are welcome. Please feel free to fork the repository, make your changes, and submit a pull request.

## Contact
For issues, questions, or contributions, please open an issue in the GitHub repository.

Version 1.0 - May 2024 

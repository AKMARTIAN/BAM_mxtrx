# BAMview

## Overview
BAMview is a Python application designed to analyze and visualize quality distributions of sequencing reads from BAM files. It processes BAM files to compute various metrics and generates histograms for base and mapping quality distributions.

## Features
- **GUI for BAM file selection**: Users can select their input file using a simple file dialog.
- **Comprehensive metrics computation**: Total reads, mapped reads, unmapped reads, and more.
- **Quality distribution plots**: Generates histograms of base and mapping qualities.
- **Command-line interface**: Supports customization through command-line arguments.

## Requirements
- Python 3
- Libraries: pysam, matplotlib, tkinter

## Installation
Ensure Python 3 is installed along with the required libraries:
```bash
pip install pysam matplotlib
```

## Usage
To use BAMview, follow these instructions:

1. Run the script from the command line:
```bash
python BAMview.py
```

2. A file dialog will appear for you to select your BAM file.

3. The program will process the selected BAM file and generate two plots showing the base and mapping quality distributions, which will be saved as PNG files.

4. A text file containing detailed metrics about the processed reads will also be generated.

### Command-line Arguments
- `--title` (optional): Specify a base title for the graphs.

## Output

PNG files: Plots showing the distributions of base and mapping qualities.
Text file: Detailed metrics of the processed reads, including total reads, mapped reads, and high-quality reads.

## Note:

The script filters out reads with missing base quality scores to avoid errors. This filtering is reported in the metrics file.

## Contributing
Contributions to BAMview are welcome. Please feel free to fork the repository, make your changes, and submit a pull request.

## Contact
For issues, questions, or contributions, please open an issue in the GitHub repository.

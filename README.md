# BAM Analyser

The BAM Analyser is a Python script designed to analyze base and mapping quality distributions from BAM files. It provides essential quality metrics and generates informative visual plots to help understand the data better.

## Dependencies

- **Python 3**: Tested with Python 3.10
- **pysam**: A tool for reading, manipulating, and writing BAM files. Available at [pysam](https://pypi.org/project/pysam/)
- **matplotlib**: A plotting library, usually included with most Python installations. If not, install it via pip:
  ```bash
  pip install matplotlib
  ```

## Installation

1. **Install required libraries**: Follow the guide in the Dependencies section to install Python 3 and pysam.
2. **(Optional) Create a virtual environment**: Using tools like `venv` or `conda` is recommended to isolate project dependencies.

## Usage

1. **Save the script**: Download and save the script as `BAM_analyser.py`.
2. **Make executable (optional)**: To run the script directly from the command line, prepend it with a shebang line (`#!/usr/bin/env python3`) and set the file permissions to executable:
   ```bash
   chmod +x BAM_analyser.py
   ```
3. **Run the script**:
   - **Basic execution**:
     ```bash
     ./BAM_analyser.py your_bam_file.bam
     ```
   - **With optional title for plots**:
     ```bash
     ./BAM_analyser.py your_bam_file.bam --title "My BAM Analysis"
     ```

## Output

- **Metrics Text File**: Outputs a file named `your_bam_file_metrics.txt`, providing metrics like total reads, mapped reads, unmapped reads, average alignment score, and counts of reads with high base quality (above 30) and perfect mapping quality (score of 60).
- **Quality Distribution Plots (Optional)**: If matplotlib is installed, the script generates plots visualizing the base and mapping quality distributions, saved as PNG images (e.g., `your_bam_file_base_quality.png`).

## Note

The script excludes reads lacking base quality scores to prevent errors, and these exclusions are documented in the metrics file.

---

Feel free to customize the wording or add additional sections as needed for your project!

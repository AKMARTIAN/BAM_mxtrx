BAM Analyser

This Python script analyzes base and mapping quality distributions from a BAM file. It provides basic quality metrics and generates informative plots.

Dependencies:

Python 3 (tested with 3.10)
pysam library (https://pypi.org/project/pysam/)
matplotlib library (usually included with most Python installations, but can be installed using pip install matplotlib if necessary)
Installation:

Install required libraries: Ensure you have Python 3 and pysam installed as mentioned in the Dependencies section.

(Optional) Create a virtual environment: It's recommended to create a virtual environment to isolate project dependencies. You can use tools like venv or conda for this purpose.

Usage:

Save the script: Save the script as BAM_analyser.py.

Make executable (optional): If you want to run the script directly from the command line, add a shebang line (e.g., #!/usr/bin/env python3) at the very beginning of the script and set the file permissions to executable using chmod +x BAM_analyser.py.

Run the script:

Basic execution:

Bash
./BAM_analyser.py your_bam_file.bam
Use code with caution.

 Replace your_bam_file.bam with the path to your actual BAM file.

With optional title (for plots):

Bash
./BAM_analyser.py your_bam_file.bam --title "My BAM Analysis"
Use code with caution.

 This creates plots with the provided title.

Output:

The script generates two main outputs:

Metrics Text File: A text file named your_bam_file_metrics.txt (where your_bam_file is the name of your BAM file) is created. It contains various metrics like total reads, mapped reads, unmapped reads, average alignment score, and the number of reads with high base quality (above 30) and perfect mapping quality (score of 60).

Quality Distribution Plots (Optional): If matplotlib is installed, the script also generates plots visualizing the base quality and mapping quality distributions of the reads in the BAM file. These plots are saved as PNG images with filenames indicating their content (e.g., your_bam_file_base_quality.png).

Note:

The script filters out reads with missing base quality scores to avoid errors. This filtering is reported in the metrics file.
Further Developments:

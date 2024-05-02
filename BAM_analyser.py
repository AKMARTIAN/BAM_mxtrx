#!/usr/bin/env python3
import pysam  # Library for working with SAM/BAM files
import matplotlib.pyplot as plt  # Library for creating plots
import tkinter as tk  # Library for the GUI file selection dialog
from tkinter import filedialog
import argparse  # Library for handling command-line arguments
import os  # For file path operations

# --- Function Definitions ---

def select_bam_file():
    """
    Creates a file selection dialog for the user to choose a BAM file.

    Returns:
        str: The filepath of the selected BAM file.
    """
    root = tk.Tk()
    root.withdraw()  # Hide the main tkinter window
    file_path = filedialog.askopenfilename(title="Select BAM File", filetypes=[("BAM files", "*.bam")])
    return file_path

def is_del(cigar):
    """
    Checks if a CIGAR string represents a deletion in the alignment.

    Args:
        cigar (str): The CIGAR string.

    Returns:
        bool: True if the CIGAR string indicates a deletion, False otherwise.
    """
    for operation, length in cigar:
        if operation == 2:  # 2 represents a deletion in the CIGAR format
            return True
    return False

def main():
    """
    Main function that performs analysis on the selected BAM file.
    """
    bam_file_path = select_bam_file()

    if not bam_file_path:  # Check if a file was selected
        print("No BAM file selected.")
        return

    parser = argparse.ArgumentParser(description='Plot base and mapping quality distributions from a BAM file.')
    parser.add_argument('--title', type=str,  # Optional title override
                        help='Base title for the graphs (optional)')
    args = parser.parse_args()

    base_filename = os.path.splitext(os.path.basename(bam_file_path))[0]

    # --- Metric Calculation ---
    total_reads = 0
    mapped_reads = 0
    unmapped_reads = 0
    skipped_reads = 0
    alignment_scores = []
    base_qualities = []
    mapping_qualities = []
    high_base_quality_reads = 0
    perfect_mapping_reads = 0

    with pysam.AlignmentFile(bam_file_path) as bam_file:
        for read in bam_file.fetch():
            total_reads += 1

            # Collect data for base qualities and filtering
            if not is_del(read.cigar) and not read.has_tag('NH'):
                for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    if query_pos is not None and read.query_qualities is not None:
                        base_qualities.append(read.query_qualities[query_pos])
                    else:
                        print(f"Read with CIGAR {read.cigar} has None query_pos or None query_qualities")
                        skipped_reads += 1
            else:
                skipped_reads += 1

            # Collect data for mapping qualities, alignment stats, and new metrics
            if not read.is_unmapped:
                mapped_reads += 1
                alignment_scores.append(read.mapping_quality)
                mapping_qualities.append(read.mapping_quality)

                if read.mapping_quality == 60:
                    perfect_mapping_reads += 1

            else:
                unmapped_reads += 1

            # Base Quality Metric (with fix)
            if read.query_qualities is not None:
                for q_score in read.query_qualities:
                    if q_score >= 30:
                        high_base_quality_reads += 1
                        break

    # --- Writing Metrics to Text File ---
    output_filename = base_filename + "_metrics.txt"
    with open(output_filename, "w") as metrics_file:
        metrics_file.write("Total Reads: {}\n".format(total_reads))
        metrics_file.write("Mapped Reads: {}\n".format(mapped_reads))
        metrics_file.write("Unmapped Reads: {}\n".format(unmapped_reads))
        metrics_file.write("Skipped Reads (CIGAR issues): {}\n".format(skipped_reads))
        metrics_file.write("Average Alignment Score: {}\n".format(sum(alignment_scores) / len(alignment_scores)))
        metrics_file.write("Reads with Base Quality >= 30: {}\n".format(high_base_quality_reads))
        metrics_file.write("Reads with Mapping Quality 60: {}\n".format(perfect_mapping_reads))
   # --- Create the Plots ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

    # Base quality plot
    ax1.hist(base_qualities, bins=40)
    ax1.set_xlabel("Base Quality Score")
    ax1.set_ylabel("Frequency")
    title_text = base_filename if not args.title else args.title
    ax1.set_title(title_text + " Base Quality Distribution")

    # Mapping quality plot
    ax2.hist(mapping_qualities, bins=60)
    ax2.set_xlabel("Mapping Quality Score")
    ax2.set_ylabel("Frequency")
    ax2.set_title(title_text + " Mapping Quality Distribution")

    plt.tight_layout()

    # Save the figure
    output_filename = base_filename + "_combined_plots.png"
    plt.savefig(output_filename)

    plt.show()

if __name__ == "__main__":
    main()

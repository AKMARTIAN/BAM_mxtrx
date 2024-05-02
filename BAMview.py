#!/usr/bin/env python3
import pysam  # Library for working with SAM/BAM files
import matplotlib.pyplot as plt  # Library for creating plots
import tkinter as tk  # Library for the GUI file selection dialog
from tkinter import filedialog
import argparse  # Library for handling command-line arguments
import os  # For file path operations

# --- Global tkinter root ---
root = tk.Tk()
root.withdraw()  # Hide the main tkinter window


# --- Function Definitions ---

def parse_arguments():
    """
    Parses command line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Analyze and plot base and mapping quality distributions from a BAM file.')
    parser.add_argument('--title', type=str, help='Base title for the graphs (optional)')
    return parser.parse_args()


def select_bam_file():
    """
    Creates a file selection dialog for the user to choose a BAM file.
    Returns:
        str: The filepath of the selected BAM file.
    """
    file_path = filedialog.askopenfilename(parent=root, title="Select BAM File", filetypes=[("BAM files", "*.bam")])
    return file_path


def is_del(cigar):
    """
    Checks if a CIGAR string represents a deletion in the alignment.
    Returns immediately if a deletion is found for efficiency.
    """
    return any(operation == 2 for operation, length in cigar)


def process_bam_file(bam_file_path):
    """
    Process the BAM file to compute metrics and collect quality scores.
    """
    metrics = {
        'total_reads': 0,
        'mapped_reads': 0,
        'unmapped_reads': 0,
        'skipped_reads': 0,
        'high_base_quality_reads': 0,
        'perfect_mapping_reads': 0,
        'alignment_scores': [],
        'base_qualities': [],
        'mapping_qualities': []
    }

    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for read in bam_file:
            metrics['total_reads'] += 1

            if not is_del(read.cigar) and not read.is_unmapped and not read.has_tag('NH'):
                for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    if query_pos is not None and read.query_qualities is not None:
                        metrics['base_qualities'].append(read.query_qualities[query_pos])
            else:
                metrics['skipped_reads'] += 1

            if not read.is_unmapped:
                metrics['mapped_reads'] += 1
                metrics['alignment_scores'].append(read.mapping_quality)
                metrics['mapping_qualities'].append(read.mapping_quality)

                if read.mapping_quality == 60:
                    metrics['perfect_mapping_reads'] += 1
            else:
                metrics['unmapped_reads'] += 1

            if read.query_qualities:
                if any(q_score >= 30 for q_score in read.query_qualities):
                    metrics['high_base_quality_reads'] += 1

    return metrics


def plot_metrics(metrics, base_filename, title=None):
    """
    Plot the base and mapping quality distributions.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
    title_text = title if title else base_filename

    ax1.hist(metrics['base_qualities'], bins=40)
    ax1.set_xlabel("Base Quality Score")
    ax1.set_ylabel("Frequency")
    ax1.set_title(title_text + " Base Quality Distribution")

    ax2.hist(metrics['mapping_qualities'], bins=60)
    ax2.set_xlabel("Mapping Quality Score")
    ax2.set_ylabel("Frequency")
    ax2.set_title(title_text + " Mapping Quality Distribution")

    plt.tight_layout()
    plt.savefig(base_filename + "_combined_plots.png")
    plt.show()


def main():
    args = parse_arguments()
    bam_file_path = select_bam_file()
    if not bam_file_path:
        print("No BAM file selected.")
        return

    metrics = process_bam_file(bam_file_path)
    base_filename = os.path.splitext(os.path.basename(bam_file_path))[0]

    plot_metrics(metrics, base_filename, args.title)

    output_filename = base_filename + "_metrics.txt"
    with open(output_filename, "w") as metrics_file:
        metrics_file.write("Total Reads: {}\n".format(metrics['total_reads']))
        metrics_file.write("Mapped Reads: {}\n".format(metrics['mapped_reads']))
        metrics_file.write("Unmapped Reads: {}\n".format(metrics['unmapped_reads']))
        metrics_file.write("Skipped Reads (CIGAR issues): {}\n".format(metrics['skipped_reads']))
        metrics_file.write(
            "Average Alignment Score: {}\n".format(sum(metrics['alignment_scores']) / len(metrics['alignment_scores'])))
        metrics_file.write("Reads with Base Quality >= 30: {}\n".format(metrics['high_base_quality_reads']))
        metrics_file.write("Reads with Perfect Mapping Quality: {}\n".format(metrics['perfect_mapping_reads']))


if __name__ == "__main__":
    main()

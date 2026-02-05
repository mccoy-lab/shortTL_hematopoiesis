#!/usr/bin/env python
"""
aggregate_bam_stats.py

This script aggregates summary statistics from multiple samtools stats files into a single DataFrame.
It parses the summary data from each provided samtools stats file and merges them based on the metric names.

Functions:
    parse_samtools_stats(file_path)
        Parses the samtools stats file and extracts summary data.

        Main function to aggregate BAM stats from multiple samtools stats files.

Usage:
    aggregate_bam_stats.py <samtools_stats_file1> [<samtools_stats_file2> ...]

Arguments:
    <samtools_stats_file1> [<samtools_stats_file2> ...] : Paths to one or more samtools stats files to be aggregated.

Example:
    python aggregate_bam_stats.py stats1.txt stats2.txt stats3.txt
"""

import sys
import os
import pandas as pd

# Function to parse the samtools stats file and extract summary data
def parse_samtools_stats(file_path):
    summary_data = []
    with open(file_path, 'r') as file:
        in_summary = False
        for line in file:
            # Check if the line starts with "SN" indicating summary data
            if line.startswith("SN"):
                in_summary = True
                parts = line.strip().split('\t')
                # Append the metric and its value to the summary_data list
                summary_data.append((parts[1].strip(':'), parts[2]))
            # If we are in the summary section and the line does not start with "SN", break the loop
            elif in_summary and not line.startswith("SN"):
                break
    return summary_data

# Main function to aggregate BAM stats from multiple samtools stats files
def main():
    # Check if at least one file path is provided as an argument
    if len(sys.argv) < 2:
        print("Usage: aggregate_bam_stats.py <samtools_stats_file1> [<samtools_stats_file2> ...]")
        sys.exit(1)

    all_data = []
    file_names = []
    # Parse each provided samtools stats file and collect the summary data
    for file_path in sys.argv[1:]:
        summary_data = parse_samtools_stats(file_path)
        all_data.append(summary_data)
        file_names.append(os.path.basename(file_path))

    # Create a DataFrame from the first file's summary data
    df = pd.DataFrame(all_data[0], columns=["Metric", "Value"])
    # Merge the summary data from subsequent files into the DataFrame
    for i, data in enumerate(all_data[1:], start=1):
        df = df.merge(pd.DataFrame(data, columns=["Metric", f"Value_{i+1}"]), on="Metric")

    # Transpose the DataFrame and set the first column as headers
    df_transposed = df.set_index("Metric").T.reset_index()
    df_transposed.columns = ["File"] + df["Metric"].tolist()

    # Replace spaces with underscores in column headers
    df_transposed.columns = [col.replace(" ", "_") for col in df_transposed.columns]

    # Update the "File" column with the filenames
    df_transposed["File"] = file_names

    # Print the transposed DataFrame
    print(df_transposed.to_string(index=False))

# Entry point of the script
if __name__ == "__main__":
    main()
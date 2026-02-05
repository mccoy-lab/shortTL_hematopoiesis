#!/usr/bin/env python
"""
Usage:
    aggregate_bam_coverage_stats.py <file1> <file2> ... <fileN>

Description:
    This script aggregates coverage statistics from multiple samtools stats files.
    It parses each file to extract coverage data and combines the data into a single
    pandas DataFrame, which is then printed as TSV to standard output. Note, only coverage 
    values from [1:1000] are presented. The final field (>1000) is not included in the output.

Arguments:
    <file1> <file2> ... <fileN>  Paths to the samtools stats files to be aggregated.
"""
import sys
import os
import pandas as pd

def parse_samtools_stats(file_path):
    """
    Parses a samtools stats file to extract coverage data.

    Args:
        file_path (str): Path to the samtools stats file.

    Returns:
        list: A list of lists, where each sublist contains coverage data from a line starting with 'COV'.
    """
    cov_data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('COV'):
                parts = line.strip().split('\t')
                cov_data.append(parts)
    return cov_data

def aggregate_coverage_stats(file_paths):
    """
    Aggregates coverage statistics from multiple samtools stats files.

    Args:
        file_paths (list): List of paths to samtools stats files.

    Returns:
        DataFrame: A pandas DataFrame containing aggregated coverage data.
    """
    all_cov_data = []
    for file_path in file_paths:
        cov_data = parse_samtools_stats(file_path)
        for row in cov_data[:-1]:
            all_cov_data.append([os.path.basename(file_path)] + row)
    
    columns = ['File', 'Type', 'Coverage', 'Cov', 'Count']
    df = pd.DataFrame(all_cov_data, columns=columns)
    df = df[['File', 'Type', 'Cov', 'Count']]
    return df

if __name__ == "__main__":
    # Check if at least one file path is provided as an argument
    if len(sys.argv) < 2:
        print("Usage: python aggregate_bam_coverage_stats.py <file1> <file2> ... <fileN>")
        sys.exit(1)
    
    # Get the list of file paths from command line arguments
    file_paths = sys.argv[1:]
    
    # Aggregate coverage statistics from the provided files
    df = aggregate_coverage_stats(file_paths)
    
    # Print the aggregated data as TSV to standard output
    print(df.to_csv(sep='\t', index=False))
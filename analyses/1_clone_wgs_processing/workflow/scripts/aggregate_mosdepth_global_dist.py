#!/usr/bin/env python

"""
This script aggregates mosdepth global distribution data from multiple files.
Usage:
    python aggregate_mosdepth_summaries.py <file1> <file2> ... <fileN>
Functions:
    read_mosdepth_files(file_paths):
        Reads and processes mosdepth output files.
        Parameters:
            file_paths (list of str): List of file paths to mosdepth output files.
        Returns:
            data_frames (list of pandas.DataFrame): List of processed data frames.
    main():
        Main function to execute the script.
        - Reads file paths from command line arguments.
        - Aggregates data from the provided files.
        - Prints the combined data frame as a tab-separated values (TSV) string.
"""

import sys
import pandas as pd

def read_mosdepth_files(file_paths):
    data_frames = []
    for file_path in file_paths:
        df = pd.read_csv(file_path, sep='\t', header=None, names=['chrom', 'depth', 'proportion'])
        df = df[df['chrom'].isin(['chr' + str(i) for i in range(1, 23)])]
        df = df[df['depth'] <= 100]
        df['sample'] = file_path.split('/')[-1].split('.')[0]
        df = df[['sample', 'chrom', 'depth', 'proportion']]
        df = df.sort_values(by=['chrom', 'depth'], ascending=True)

        data_frames.append(df)
    return data_frames

def main():
    if len(sys.argv) < 2:
        print("Usage: python aggregate_mosdepth_summaries.py <file1> <file2> ... <fileN>")
        sys.exit(1)
    
    file_paths = sys.argv[1:]
    data_frames = read_mosdepth_files(file_paths)
    combined_df = pd.concat(data_frames)
    
    print(combined_df.to_csv(sep='\t', index=False))

if __name__ == "__main__":
    main()
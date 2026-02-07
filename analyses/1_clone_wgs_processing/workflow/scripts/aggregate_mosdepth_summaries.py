#!/usr/bin/env python
"""
This script aggregates mosdepth summary files for chromosomes 1 through 22.
Usage:
    python aggregate_mosdepth_summaries.py <file1> <file2> ... <fileN>
Arguments:
    <file1> <file2> ... <fileN> : List of mosdepth summary files to aggregate.
The script reads each input file, filters for chromosomes 1 through 22, adds a sample identifier based on the filename,
and concatenates the data into a single DataFrame. The resulting DataFrame is then output to stdout in tab-separated format.
Functions:
    aggregate_mosdepth_summaries(files):
        Aggregates mosdepth summary files into a single DataFrame and writes the result to stdout.
        Parameters:
            files (list of str): List of file paths to mosdepth summary files.
        Returns:
            None
"""

import sys
import pandas as pd

def aggregate_mosdepth_summaries(files):
    dataframes = []
    for file in files:
        df = pd.read_csv(file, sep='\t')
        df = df[df['chrom'].isin(['chr' + str(i) for i in range(1, 23)])]
        df['sample'] = file.split('/')[-1].split('.')[0]
        df = df[['sample','chrom','length','bases','mean','min','max']]
        dataframes.append(df)
    result = pd.concat(dataframes)
    result.to_csv(sys.stdout, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python aggregate_mosdepth_summaries.py <file1> <file2> ... <fileN>")
        sys.exit(1)
    
    files = sys.argv[1:]
    aggregate_mosdepth_summaries(files)
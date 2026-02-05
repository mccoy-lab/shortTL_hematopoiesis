#!/usr/bin/env python
"""
Filter TSV by values based on specified criteria.
This script reads a TSV file, filters it based on a numeric threshold applied to a specified field,
and outputs values from another field that meet the filter criteria. Here, we use this to remove samples
that have low germline sensitivity estimates.

Usage:
    python filter_tsv_by_values.py --tsv_path <path> --output_field <field> \
                                     --filter_field <field> --filter_threshold <value> \
                                     --output_path <path> --mode <mode> \
                                     [--remove_consolidated_germline]
"""

import argparse
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter TSV by values")
    parser.add_argument("--tsv_path", required=True, help="Path to the input TSV file")
    parser.add_argument("--output_field", required=True, help="Field to output")
    parser.add_argument("--filter_field", required=True, help="Field to filter by")
    parser.add_argument("--filter_threshold", type=float, required=True, help="Threshold value for filtering")
    parser.add_argument("--output_path", required=True, help="Path to the output text file")
    parser.add_argument("--mode", choices=["greater", "less", "greater_equal", "less_equal"], required=True, help="Mode of filtering: 'greater', 'less', 'greater_equal', or 'less_equal'")
    parser.add_argument("--remove_consolidated_germline", action="store_true", help="Remove 'consolidated_germline' if present in the output field")

    args = parser.parse_args()
    
    df = pd.read_csv(args.tsv_path, sep='\t')
    
    if args.mode == "greater":
        filtered_df = df[df[args.filter_field] > args.filter_threshold]
    elif args.mode == "less":
        filtered_df = df[df[args.filter_field] < args.filter_threshold]
    elif args.mode == "greater_equal":
        filtered_df = df[df[args.filter_field] >= args.filter_threshold]
    elif args.mode == "less_equal":
        filtered_df = df[df[args.filter_field] <= args.filter_threshold]
    
    output_series = filtered_df[args.output_field]
    
    if args.remove_consolidated_germline:
        output_series = output_series[output_series != "consolidated_germline"]
    
    output_series.to_csv(args.output_path, sep='\n', index=False, header=False)
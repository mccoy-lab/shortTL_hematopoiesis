#!/usr/bin/env python
"""
This script reads specified JSON files, extracts the "before filtering" summary section,
converts them into pandas DataFrames, concatenates them, and writes the output as a TSV to standard output.

Functions:
    read_json_files(file_list):
        Reads the given JSON files, extracts the "before filtering" summary section,
        converts them into pandas DataFrames, and returns a list of these DataFrames.
        
Example:
    python read_json_stats.py /path/to/file1.json /path/to/file2.json

Notes:
    - The script expects the provided paths to be JSON files.
    - Each JSON file is read and the "before filtering" summary section is normalized into a pandas DataFrame.
    - The DataFrames are concatenated and written as a TSV to standard output.
"""

import sys  # Import the sys module for accessing command-line arguments
import pandas as pd  # Import the pandas library for data manipulation and analysis
import json  # Import the json module for parsing JSON files

def read_json_files(file_list):
    dataframes = []  # Initialize an empty list to store DataFrames
    for filepath in file_list:  # Iterate over all provided file paths
        if filepath.endswith('.json'):  # Check if the file has a .json extension
            with open(filepath, 'r') as file:  # Open the JSON file for reading
                data = json.load(file)  # Load the JSON data from the file
                if "summary" in data and "before_filtering" in data["summary"]:  # Check if "summary" and "before filtering" sections exist
                    df = pd.json_normalize(data["summary"]["before_filtering"])  # Normalize the "before filtering" section into a pandas DataFrame
                    df['filename'] = filepath  # Add a column with the filename
                    df = df[['filename'] + [col for col in df.columns if col != 'filename']]  # Reorder columns to put 'filename' at the front
                    dataframes.append(df)  # Append the DataFrame to the list
    return dataframes  # Return the list of DataFrames

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) < 2:  # Ensure there is at least one command-line argument
        print("Usage: python read_json_stats.py <file1> <file2> ... <fileN>")  # Print usage information
        sys.exit(1)  # Exit the script with an error code

    # Get the list of files from the command line arguments
    file_list = sys.argv[1:]  # Assign the command-line arguments to 'file_list'
    
    # Read JSON files and concatenate DataFrames
    dataframes = read_json_files(file_list)  # Call the function to read JSON files into DataFrames
    if dataframes:  # Check if there are any DataFrames to concatenate
        concatenated_df = pd.concat(dataframes, ignore_index=True)  # Concatenate all DataFrames
        concatenated_df.to_csv(sys.stdout, sep='\t', index=False)  # Write the concatenated DataFrame as TSV to standard output
    else:
        print("No valid JSON files found or no 'before filtering' section found.")  # Print a message if no valid JSON files are found

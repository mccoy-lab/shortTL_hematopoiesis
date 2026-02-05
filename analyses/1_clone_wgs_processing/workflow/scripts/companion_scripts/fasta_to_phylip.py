#!/usr/bin/env python

"""
Convert FASTA sequence alignment to PHYLIP format.
"""

from Bio import AlignIO
import sys

def fasta_to_phylip(input_file, output_file):
    # Read the FASTA file
    alignment = AlignIO.read(input_file, "fasta")
    
    # Write the alignment to a PHYLIP file
    AlignIO.write(alignment, output_file, "phylip-relaxed")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fasta_to_phylip.py <input_fasta> <output_phylip>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    fasta_to_phylip(input_file, output_file)
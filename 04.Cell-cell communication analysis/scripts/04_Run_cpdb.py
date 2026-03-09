#!/usr/bin/env python3

import sys
import os
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

if len(sys.argv) < 4:
    print("Error: Missing arguments.")
    print("Usage: python 04_Run_cpdb.py <meta_file> <counts_file> <database.zip> <output_dir>")
    sys.exit(1)

meta_file_path = sys.argv[1]
counts_file_path = sys.argv[2]
cpdb_file_path = sys.argv[3]
output_path = sys.argv[4]

# Ensure output directory exists
os.makedirs(output_path, exist_ok=True)

# Call CellPhoneDB statistical analysis method
cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,
    meta_file_path = meta_file_path,
    counts_file_path = counts_file_path,
    counts_data = 'hgnc_symbol',
    separator = '|', 
    score_interactions = True,
    output_path = output_path,
    threads = 8,
    threshold = 0.1
)

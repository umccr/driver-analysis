"""
Reads elements_results.txt files in a root dircetory. 
Reads COSMIC gene list.
Sets CGC column to True, if GE in elements_results.txt exists in COSMIC gene list.

Inputs:
    cosmic_df: path to file with elements results
    root_directory: p-value column name in dataframe to calculate KS

Returns:
    An updated elements_results.txt file written adjacent to original input.
"""

import pandas as pd
import os

# Read the cosmic file into a dataframe
cosmic_df = pd.read_csv('cosmic_grch38_v96.txt', sep='\t', header=None, names=['SYMBOL'])

# Specify the root directory containing subdirectories with "elements_results.txt" files
root_directory = "/Users/kanwals/UMCCR/research/projects/PAAD_atlas/driver_analysis/oncodriveclustl/model_adjustment/results/"

for root, dirs, files in os.walk(root_directory):
    if "elements_results.txt" in files:
        file_path = os.path.join(root, "elements_results.txt")

        # Read the elements_results.txt file into a dataframe
        elements_df = pd.read_csv(file_path, sep='\t')

        # Merge the dataframes based on the SYMBOL column
        merged_df = pd.merge(elements_df, cosmic_df, on='SYMBOL', how='left')

        # Update the CGC column to True where SYMBOL is present in cosmic_df
        merged_df.loc[merged_df['SYMBOL'].isin(cosmic_df['SYMBOL']), 'CGC'] = 'True'

        # Write the updated dataframe to a new file
        merged_df.to_csv(f'{file_path}_updated', sep='\t', index=False)


print("Results updated")